#!/usr/bin/env python3
"""Validate a revocable-credential manifest and generate C/Python bindings.

The JSON manifest is the only hand-edited source for parameters shared by the
accumulator, network protocol, and Labrador MemVer backend.  A reproducible
candidate is deliberately distinct from a security-reviewed frozen release.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import pprint
import sys
from typing import Any


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = ROOT / "params" / "revac128-candidate-v1.json"
C_HEADER = ROOT / "code/entities/common/core/include/revac_profile_generated.h"
PY_MODULE = ROOT / "code/entities/zk_lazer/revac_profile_generated.py"


class ProfileError(ValueError):
    pass


def require_int(obj: dict[str, Any], key: str, minimum: int = 0) -> int:
    value = obj.get(key)
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ProfileError(f"{key} must be an integer >= {minimum}")
    return value


def canonical_parameters(manifest: dict[str, Any]) -> bytes:
    bound = {
        "schema_version": manifest["schema_version"],
        "profile_name": manifest["profile_name"],
        "security_target": {
            "computational_bits": manifest["security"]["target_computational_bits"],
            "statistical_bits": manifest["security"]["target_statistical_bits"],
            "signature_query_bound": manifest["security"]["signature_query_bound"],
        },
        "source_bindings": manifest["source_bindings"],
        "parameters": manifest["parameters"],
    }
    return json.dumps(bound, sort_keys=True, separators=(",", ":")).encode()


def profile_id(manifest: dict[str, Any]) -> bytes:
    return hashlib.sha3_256(canonical_parameters(manifest)).digest()


def validate(manifest: dict[str, Any], require_release: bool = False) -> list[str]:
    warnings: list[str] = []
    if require_int(manifest, "schema_version", 1) != 1:
        raise ProfileError("unsupported schema version")
    name = manifest.get("profile_name")
    if not isinstance(name, str) or not name or not name.isascii():
        raise ProfileError("profile_name must be non-empty ASCII")
    status = manifest.get("release_status")
    if status not in {"candidate", "frozen"}:
        raise ProfileError("release_status must be candidate or frozen")

    security = manifest.get("security")
    params = manifest.get("parameters")
    if not isinstance(security, dict) or not isinstance(params, dict):
        raise ProfileError("security and parameters objects are required")
    require_int(security, "target_computational_bits", 1)
    require_int(security, "target_statistical_bits", 1)
    require_int(security, "signature_query_bound", 1)
    blockers = security.get("release_blockers")
    if not isinstance(blockers, list) or any(
        not isinstance(item, str) or not item for item in blockers
    ):
        raise ProfileError("release_blockers must be a list of strings")
    complete = security.get("derivation_complete")
    if not isinstance(complete, bool):
        raise ProfileError("derivation_complete must be boolean")
    if status == "frozen" and (not complete or blockers):
        raise ProfileError("frozen profile requires completed derivation and no blockers")
    if require_release and status != "frozen":
        raise ProfileError("candidate profile is not a security-reviewed release")

    bindings = manifest.get("source_bindings")
    if not isinstance(bindings, dict):
        raise ProfileError("source_bindings object is required")
    params_h_rel = bindings.get("params_h_path")
    expected_params_h_hash = bindings.get("params_h_sha3_256")
    if not isinstance(params_h_rel, str) or not params_h_rel:
        raise ProfileError("params_h_path must be a non-empty string")
    if (not isinstance(expected_params_h_hash, str) or
            len(expected_params_h_hash) != 64 or
            any(c not in "0123456789abcdef" for c in expected_params_h_hash)):
        raise ProfileError("params_h_sha3_256 must be lowercase SHA3-256 hex")
    params_h_path = (ROOT / params_h_rel).resolve()
    try:
        params_h_path.relative_to(ROOT.resolve())
    except ValueError as exc:
        raise ProfileError("params_h_path escapes the repository") from exc
    actual_params_h_hash = hashlib.sha3_256(params_h_path.read_bytes()).hexdigest()
    if actual_params_h_hash != expected_params_h_hash:
        raise ProfileError(
            "params.h changed without a new parameter profile: "
            f"expected {expected_params_h_hash}, got {actual_params_h_hash}"
        )

    ring = params["ring"]
    n = require_int(ring, "n", 1)
    d = require_int(ring, "d", 1)
    q = require_int(ring, "q", 3)
    if n & (n - 1):
        raise ProfileError("ring degree must be a power of two")
    if require_int(ring, "q_bitlen", 1) != q.bit_length():
        raise ProfileError("ring q_bitlen does not match q")

    signature = params["signature"]
    q_low = require_int(signature, "q_low", 2)
    q_high = require_int(signature, "q_high", 2)
    if q_low * q_high != q:
        raise ProfileError("signature q_low * q_high must equal q")
    if require_int(signature, "base_low", 2) ** require_int(
        signature, "dim_low", 1
    ) != q_low:
        raise ProfileError("low signature gadget does not reconstruct q_low")
    if require_int(signature, "base_high", 2) ** require_int(
        signature, "dim_high", 1
    ) != q_high:
        raise ProfileError("high signature gadget does not reconstruct q_high")

    acc = params["accumulator"]
    handle_bits = require_int(acc, "handle_bits", 1)
    base = require_int(acc, "gadget_base", 2)
    gadget_dim = require_int(acc, "gadget_dim", 1)
    signed = acc.get("gadget_signed")
    if not isinstance(signed, bool):
        raise ProfileError("gadget_signed must be boolean")
    random_update_signs = acc.get("authenticated_random_update_signs")
    if not isinstance(random_update_signs, bool):
        raise ProfileError("authenticated_random_update_signs must be boolean")
    randomized_witness_refresh = acc.get("randomized_witness_refresh")
    if not isinstance(randomized_witness_refresh, bool):
        raise ProfileError("randomized_witness_refresh must be boolean")
    if signed and base % 2 == 0:
        raise ProfileError("balanced signed gadget requires odd base")
    digit_bound = (base - 1) // 2 if signed else base - 1
    coverage = digit_bound * (base**gadget_dim - 1) // (base - 1)
    if coverage < q // 2:
        raise ProfileError("accumulator gadget does not cover centered q")
    trunc_k = require_int(acc, "truncated_sampler_k", 1)
    trunc_l = require_int(acc, "truncated_sampler_l", 0)
    require_int(acc, "truncated_sampler_base", 2)
    if trunc_l >= trunc_k:
        raise ProfileError("truncated sampler requires l < k")
    require_int(acc, "max_dictionary_entries", 1)
    max_storage = require_int(acc, "max_deletions_storage", 1)

    system = params["system"]
    public_delta = require_int(system, "public_delta_layers", 1)
    private_delta = require_int(system, "private_delta_layers", 1)
    bml_capacity = require_int(system, "bml_capacity", 1)
    if not public_delta <= private_delta <= handle_bits:
        raise ProfileError("DELTA layers must satisfy public <= private <= handle bits")
    if bml_capacity > max_storage:
        raise ProfileError("BML capacity exceeds deletion storage")
    for key in (
        "seed_bytes", "show_nonce_bytes", "sync_nonce_bytes", "parameter_hash_bytes"
    ):
        require_int(system, key, 16)

    lab = params["labrador_memver"]
    proof_mod = require_int(lab, "proof_modulus", 3)
    if proof_mod.bit_length() != require_int(lab, "proof_modulus_bits", 2):
        raise ProfileError("Labrador proof modulus bit length mismatch")
    for key in (
        "linear_decomposition_base",
        "linear_decomposition_chunks",
        "small_chunk_norm_sq",
        "tag_norm_sq",
        "v2_norm_sq",
        "right_fresh_group_norm_sq",
        "right_group_norm_sq_per_deletion",
        "right_group_norm_sq_max",
    ):
        require_int(lab, key, 1)
    if lab.get("prover_zero_knowledge") is not True:
        raise ProfileError("accumulator prover must use zk=True")
    if lab.get("right_layout") != "grouped":
        raise ProfileError("this implementation requires grouped right layout")

    issue = params["issue_proof_compiled"]
    show = params["show_proof_compiled"]
    if require_int(issue, "q", 3).bit_length() != require_int(
        issue, "q_bitlen", 1
    ):
        raise ProfileError("Issue proof modulus bit length mismatch")
    if require_int(show, "q", 3).bit_length() != require_int(
        show, "q_bitlen", 1
    ):
        raise ProfileError("Show proof modulus bit length mismatch")
    expected_link = 2 * d * (n // require_int(show, "n", 1))
    if show["same_handle_link_randomness_polynomials"] != expected_link:
        raise ProfileError("same-handle polynomial count mismatch")

    link = params["same_handle_commitment"]
    if require_int(link, "opening_vectors", 1) != 2:
        raise ProfileError("current same-handle relation requires two opening vectors")
    if require_int(link, "opening_dimension", 1) != d:
        raise ProfileError("same-handle opening dimension must equal ring module rank")
    if link.get("opening_distribution") != "binary":
        raise ProfileError("current same-handle opening distribution must be binary")
    if require_int(link, "opening_poly_norm_sq", 1) != n:
        raise ProfileError("binary same-handle polynomial norm bound must equal n")
    for key in (
        "matrix_a_domain_separator",
        "matrix_d_domain_separator",
        "randomness_domain_separator",
    ):
        require_int(link, key, 0)
    versions = (
        "context_wire_version",
        "sync_wire_version",
        "show_packet_encoding_version",
        "labrador_raw_encoding_version",
        "labrador_compact_encoding_version",
        "labrador_worker_protocol_version",
    )
    for key in versions:
        if require_int(system, key, 1) > 255:
            raise ProfileError(f"{key} must fit one byte")
    if system["show_packet_encoding_version"] > 9:
        raise ProfileError("show packet encoding version must be one decimal digit")
    if system["labrador_worker_protocol_version"] > 99:
        raise ProfileError("worker protocol version must be at most two decimal digits")

    warnings.extend(str(item) for item in blockers)
    if security["target_statistical_bits"] > 40:
        warnings.append(
            "source parameter scripts still use LOG2_EPS=-40; "
            "the composed statistical-distance budget needs rederivation"
        )
    if show["q_bitlen"] <= 60:
        warnings.append(
            "compiled Show modulus fits the single-limb backend, while the "
            "current ordinary-rejection parameter model selects a larger modulus"
        )
    if acc["witness_beta_sq_policy"] >= (1 << 63) - 1:
        warnings.append("accumulator witness beta remains the permissive prototype bound")
    return warnings


def c_bytes(data: bytes) -> str:
    return "{ " + ", ".join(f"0x{x:02x}" for x in data) + " }"


def generated_header(
    manifest: dict[str, Any], digest: bytes, warnings: list[str]
) -> str:
    p = manifest["parameters"]
    r, a = p["ring"], p["accumulator"]
    lab, link, system = (
        p["labrador_memver"], p["same_handle_commitment"], p["system"]
    )
    sec = manifest["security"]
    release_ready = int(
        manifest["release_status"] == "frozen"
        and sec["derivation_complete"]
        and not sec["release_blockers"]
        and not warnings
    )
    coverage = ((a["gadget_base"] - 1) // 2) * (
        a["gadget_base"] ** a["gadget_dim"] - 1
    ) // (a["gadget_base"] - 1)
    return f"""/* Generated by scripts/params/generate_revac_profile.py.
 * Do not edit this file manually.
 */
#ifndef REVAC_PROFILE_GENERATED_H
#define REVAC_PROFILE_GENERATED_H

#define REVAC_PARAM_SCHEMA_VERSION {manifest["schema_version"]}
#define REVAC_PARAM_PROFILE_NAME "{manifest["profile_name"]}"
#define REVAC_PARAM_PROFILE_ID_HEX "{digest.hex()}"
#define REVAC_PARAM_PROFILE_ID_BYTES {len(digest)}
#define REVAC_PARAM_PROFILE_ID_BYTES_INIT {c_bytes(digest)}
#define REVAC_PARAM_RELEASE_READY {release_ready}
#define REVAC_PARAM_TARGET_COMPUTATIONAL_BITS {sec["target_computational_bits"]}
#define REVAC_PARAM_TARGET_STATISTICAL_BITS {sec["target_statistical_bits"]}
#define REVAC_PARAM_SIGNATURE_QUERY_BOUND UINT64_C({sec["signature_query_bound"]})

#define REVAC_PARAM_RING_N {r["n"]}
#define REVAC_PARAM_RING_D {r["d"]}
#define REVAC_PARAM_RING_Q {r["q"]}L
#define REVAC_PARAM_RING_Q_BITLEN {r["q_bitlen"]}

#define REVAC_PARAM_ACC_HANDLE_BITS {a["handle_bits"]}
#define REVAC_PARAM_ACC_GADGET_BASE {a["gadget_base"]}
#define REVAC_PARAM_ACC_GADGET_DIM {a["gadget_dim"]}
#define REVAC_PARAM_ACC_GADGET_SIGNED {int(a["gadget_signed"])}
#define REVAC_PARAM_ACC_RANDOM_UPDATE_SIGNS {int(a["authenticated_random_update_signs"])}
#define REVAC_PARAM_ACC_RANDOMIZED_WITNESS_REFRESH {int(a["randomized_witness_refresh"])}
#define REVAC_PARAM_ACC_GADGET_MAX_CENTERED {coverage}L
#define REVAC_PARAM_ACC_MAX_DICT_ENTRIES {a["max_dictionary_entries"]}
#define REVAC_PARAM_ACC_MAX_DELETIONS_STORAGE {a["max_deletions_storage"]}
#define REVAC_PARAM_ACC_TRUNC_BASE {a["truncated_sampler_base"]}
#define REVAC_PARAM_ACC_TRUNC_K {a["truncated_sampler_k"]}
#define REVAC_PARAM_ACC_TRUNC_L {a["truncated_sampler_l"]}
#define REVAC_PARAM_ACC_TRUNC_BASE_POW_L {a["truncated_sampler_base"] ** a["truncated_sampler_l"]}
#define REVAC_PARAM_ACC_WIT_BETA_SQ UINT64_C({a["witness_beta_sq_policy"]})

#define REVAC_PARAM_PUBLIC_DELTA_LAYERS {system["public_delta_layers"]}
#define REVAC_PARAM_PRIVATE_DELTA_LAYERS {system["private_delta_layers"]}
#define REVAC_PARAM_BML_CAPACITY {system["bml_capacity"]}
#define REVAC_PARAM_SEED_BYTES {system["seed_bytes"]}
#define REVAC_PARAM_SHOW_NONCE_BYTES {system["show_nonce_bytes"]}
#define REVAC_PARAM_SYNC_NONCE_BYTES {system["sync_nonce_bytes"]}
#define REVAC_PARAM_HASH_BYTES {system["parameter_hash_bytes"]}
#define REVAC_PARAM_CONTEXT_WIRE_VERSION {system["context_wire_version"]}
#define REVAC_PARAM_SYNC_WIRE_VERSION {system["sync_wire_version"]}
#define REVAC_PARAM_SHOW_PACKET_ENCODING_VERSION {system["show_packet_encoding_version"]}
#define REVAC_PARAM_LABRADOR_RAW_ENCODING_VERSION {system["labrador_raw_encoding_version"]}
#define REVAC_PARAM_LABRADOR_COMPACT_ENCODING_VERSION {system["labrador_compact_encoding_version"]}
#define REVAC_PARAM_LABRADOR_WORKER_PROTOCOL_VERSION {system["labrador_worker_protocol_version"]}

#define REVAC_PARAM_LINK_OPENING_VECTORS {link["opening_vectors"]}
#define REVAC_PARAM_LINK_OPENING_DIMENSION {link["opening_dimension"]}
#define REVAC_PARAM_LINK_OPENING_POLY_NORM_SQ {link["opening_poly_norm_sq"]}ULL

#define REVAC_PARAM_LABRADOR_PROOF_MODULUS {lab["proof_modulus"]}ULL
#define REVAC_PARAM_LABRADOR_PROOF_MODULUS_BITS {lab["proof_modulus_bits"]}
#define REVAC_PARAM_LABRADOR_LINEAR_DECOMP_BASE {lab["linear_decomposition_base"]}
#define REVAC_PARAM_LABRADOR_LINEAR_DECOMP_CHUNKS {lab["linear_decomposition_chunks"]}
#define REVAC_PARAM_LABRADOR_SMALL_CHUNK_NORM_SQ {lab["small_chunk_norm_sq"]}ULL
#define REVAC_PARAM_LABRADOR_TAG_NORM_SQ {lab["tag_norm_sq"]}ULL
#define REVAC_PARAM_LABRADOR_V2_NORM_SQ {lab["v2_norm_sq"]}ULL
#define REVAC_PARAM_LABRADOR_RIGHT_FRESH_GROUP_NORM_SQ {lab["right_fresh_group_norm_sq"]}ULL
#define REVAC_PARAM_LABRADOR_RIGHT_GROUP_NORM_SQ_PER_DELETION {lab["right_group_norm_sq_per_deletion"]}ULL
#define REVAC_PARAM_LABRADOR_RIGHT_GROUP_NORM_SQ_MAX {lab["right_group_norm_sq_max"]}ULL

#if PARAM_N_ISS != {p["issue_proof_compiled"]["n"]} || \
    PARAM_D_ISS != {p["issue_proof_compiled"]["module_rank"]} || \
    PARAM_Q_ISS != {p["issue_proof_compiled"]["q"]}L || \
    PARAM_Q_ISS_BITLEN != {p["issue_proof_compiled"]["q_bitlen"]} || \
    PARAM_M1_ISS != {p["issue_proof_compiled"]["witness_polynomials"]} || \
    PARAM_M2_ISS != {p["issue_proof_compiled"]["randomness_polynomials"]} || \
    PARAM_L_ISS != {p["issue_proof_compiled"]["soundness_repetitions"]} || \
    PARAM_RHO_ISS != {p["issue_proof_compiled"]["challenge_rho"]} || \
    PARAM_ETA_ISS != {p["issue_proof_compiled"]["challenge_eta"]}
#error "params.h Issue proof does not match the selected profile"
#endif
#if PARAM_N_SHOW != {p["show_proof_compiled"]["n"]} || \
    PARAM_D_SHOW != {p["show_proof_compiled"]["module_rank"]} || \
    PARAM_Q_SHOW != {p["show_proof_compiled"]["q"]}L || \
    PARAM_Q_SHOW_BITLEN != {p["show_proof_compiled"]["q_bitlen"]} || \
    PARAM_M1_SHOW != {p["show_proof_compiled"]["witness_polynomials"]} || \
    PARAM_M2_SHOW != {p["show_proof_compiled"]["randomness_polynomials"]} || \
    PARAM_L_SHOW != {p["show_proof_compiled"]["soundness_repetitions"]} || \
    PARAM_RHO_SHOW != {p["show_proof_compiled"]["challenge_rho"]} || \
    PARAM_ETA_SHOW != {p["show_proof_compiled"]["challenge_eta"]} || \
    PARAM_GAMMA_SHOW != {p["show_proof_compiled"]["compression_gamma"]}L
#error "params.h Show proof does not match the selected profile"
#endif

#if PARAM_N != REVAC_PARAM_RING_N || PARAM_D != REVAC_PARAM_RING_D || \\
    PARAM_Q != REVAC_PARAM_RING_Q || PARAM_Q_BITLEN != REVAC_PARAM_RING_Q_BITLEN
#error "params.h does not match the selected revocable-credential profile"
#endif
#if (PARAM_Q / 2) > REVAC_PARAM_ACC_GADGET_MAX_CENTERED
#error "Accumulator gadget does not cover every centered coefficient"
#endif
#if DOMAIN_SEPARATOR_REVAC_LINK_A != {link["matrix_a_domain_separator"]} || \
    DOMAIN_SEPARATOR_REVAC_LINK_D != {link["matrix_d_domain_separator"]} || \
    DOMAIN_SEPARATOR_REVAC_LINK_R != {link["randomness_domain_separator"]}
#error "same-handle commitment domains do not match the selected profile"
#endif
#if SEED_BYTES != REVAC_PARAM_SEED_BYTES || \
    CRS_SEED_BYTES != REVAC_PARAM_SEED_BYTES
#error "seed lengths do not match the selected profile"
#endif

#endif
"""


def generated_python(
    manifest: dict[str, Any], digest: bytes, warnings: list[str]
) -> str:
    # repr/pformat emits Python booleans (True/False); JSON output would leave
    # lowercase true/false tokens and make the generated module invalid.
    payload = pprint.pformat(
        manifest["parameters"], sort_dicts=True, width=100
    )
    ready = (
        manifest["release_status"] == "frozen"
        and manifest["security"]["derivation_complete"]
        and not manifest["security"]["release_blockers"]
        and not warnings
    )
    return f'''"""Generated revocable-credential profile. Do not edit."""

from __future__ import annotations

PROFILE_NAME = {manifest["profile_name"]!r}
PROFILE_ID = bytes.fromhex({digest.hex()!r})
PROFILE_ID_HEX = PROFILE_ID.hex()
RELEASE_READY = {ready}
TARGET_COMPUTATIONAL_BITS = {manifest["security"]["target_computational_bits"]}
TARGET_STATISTICAL_BITS = {manifest["security"]["target_statistical_bits"]}
PARAMETERS = {payload}

RING = PARAMETERS["ring"]
ACCUMULATOR = PARAMETERS["accumulator"]
LABRADOR_MEMVER = PARAMETERS["labrador_memver"]
SYSTEM = PARAMETERS["system"]
SAME_HANDLE_COMMITMENT = PARAMETERS["same_handle_commitment"]
N, D, Q, Q_BITLEN = RING["n"], RING["d"], RING["q"], RING["q_bitlen"]
HANDLE_BITS = ACCUMULATOR["handle_bits"]
GADGET_BASE = ACCUMULATOR["gadget_base"]
GADGET_DIM = ACCUMULATOR["gadget_dim"]
GADGET_SIGNED = ACCUMULATOR["gadget_signed"]
RANDOM_UPDATE_SIGNS = ACCUMULATOR["authenticated_random_update_signs"]
RANDOMIZED_WITNESS_REFRESH = ACCUMULATOR["randomized_witness_refresh"]
TRUNC_K = ACCUMULATOR["truncated_sampler_k"]
TRUNC_L = ACCUMULATOR["truncated_sampler_l"]
TRUNC_K_L = TRUNC_K - TRUNC_L
PUBLIC_DELTA_LAYERS = SYSTEM["public_delta_layers"]
PRIVATE_DELTA_LAYERS = SYSTEM["private_delta_layers"]
BML_CAPACITY = SYSTEM["bml_capacity"]
SEED_BYTES = SYSTEM["seed_bytes"]
HASH_BYTES = SYSTEM["parameter_hash_bytes"]
PROOF_MODULUS = LABRADOR_MEMVER["proof_modulus"]
PROOF_MODULUS_BITS = LABRADOR_MEMVER["proof_modulus_bits"]
LINEAR_DECOMP_BASE = LABRADOR_MEMVER["linear_decomposition_base"]
LINEAR_DECOMP_CHUNKS = LABRADOR_MEMVER["linear_decomposition_chunks"]
SMALL_CHUNK_NORM_SQ = LABRADOR_MEMVER["small_chunk_norm_sq"]
TAG_NORM_SQ = LABRADOR_MEMVER["tag_norm_sq"]
V2_NORM_SQ = LABRADOR_MEMVER["v2_norm_sq"]
RIGHT_FRESH_GROUP_NORM_SQ = LABRADOR_MEMVER["right_fresh_group_norm_sq"]
RIGHT_GROUP_NORM_SQ_PER_DELETION = LABRADOR_MEMVER["right_group_norm_sq_per_deletion"]
RIGHT_GROUP_NORM_SQ_MAX = LABRADOR_MEMVER["right_group_norm_sq_max"]
LINK_OPENING_POLY_NORM_SQ = SAME_HANDLE_COMMITMENT["opening_poly_norm_sq"]
CONTEXT_WIRE_VERSION = SYSTEM["context_wire_version"]
SYNC_WIRE_VERSION = SYSTEM["sync_wire_version"]
SHOW_PACKET_ENCODING_VERSION = SYSTEM["show_packet_encoding_version"]
LABRADOR_RAW_ENCODING_VERSION = SYSTEM["labrador_raw_encoding_version"]
LABRADOR_COMPACT_ENCODING_VERSION = SYSTEM["labrador_compact_encoding_version"]
LABRADOR_WORKER_PROTOCOL_VERSION = SYSTEM["labrador_worker_protocol_version"]
'''


def generated_report(
    manifest: dict[str, Any], digest: bytes, warnings: list[str]
) -> str:
    p, a = manifest["parameters"], manifest["parameters"]["accumulator"]
    r, issue, show = p["ring"], p["issue_proof_compiled"], p["show_proof_compiled"]
    lab, system = p["labrador_memver"], p["system"]
    coverage = ((a["gadget_base"] - 1) // 2) * (
        a["gadget_base"] ** a["gadget_dim"] - 1
    ) // (a["gadget_base"] - 1)
    q_poly_bits = r["n"] * r["q_bitlen"]
    acc_bytes = (r["d"] * q_poly_bits + 7) // 8
    compact_polys = (2 + a["truncated_sampler_k"] - a["truncated_sampler_l"]) * r["d"]
    full_polys = compact_polys + a["handle_bits"] * a["gadget_dim"] * r["d"]
    compact_witness_bytes = (compact_polys * q_poly_bits + r["n"] + 7) // 8
    full_witness_bytes = (full_polys * q_poly_bits + r["n"] + 7) // 8
    ba_coeffs = a["handle_bits"] * a["gadget_dim"] * r["d"] * r["d"] * r["n"]
    ba_bytes = (ba_coeffs * r["q_bitlen"] + 7) // 8
    lines = [
        "# Generated parameter-profile audit",
        "",
        f"- Profile: {manifest['profile_name']}",
        f"- Status: {manifest['release_status']}",
        f"- Profile ID (SHA3-256): {digest.hex()}",
        f"- params.h SHA3-256: {manifest['source_bindings']['params_h_sha3_256']}",
        f"- Computational target: {manifest['security']['target_computational_bits']} bits",
        f"- Statistical target: {manifest['security']['target_statistical_bits']} bits",
        f"- Balanced gadget centered coverage: {coverage}",
        f"- Centered modulus requirement: {p['ring']['q'] // 2}",
        "",
        "## Active implementation parameters",
        "",
        "| Component | Frozen candidate values |",
        "|---|---|",
        f"| Ring | n={r['n']}, d={r['d']}, q={r['q']} ({r['q_bitlen']} bits) |",
        f"| Accumulator gadget | balanced base {a['gadget_base']}, dimension {a['gadget_dim']}, {a['handle_bits']}-bit handles |",
        f"| Update signs | {'authenticated random (research profile)' if a['authenticated_random_update_signs'] else 'fixed +1 (PKC theorem-compatible)'} |",
        f"| Truncated sampler | b={a['truncated_sampler_base']}, k={a['truncated_sampler_k']}, l={a['truncated_sampler_l']} |",
        f"| Issue proof | n={issue['n']}, rank={issue['module_rank']}, q={issue['q']} ({issue['q_bitlen']} bits), m1={issue['witness_polynomials']}, m2={issue['randomness_polynomials']}, l={issue['soundness_repetitions']} |",
        f"| Show proof | n={show['n']}, rank={show['module_rank']}, q={show['q']} ({show['q_bitlen']} bits), m1={show['witness_polynomials']}, m2={show['randomness_polynomials']}, l={show['soundness_repetitions']} |",
        f"| LaBRADOR | q_p={lab['proof_modulus']}, grouped right, zk=True |",
        f"| Synchronization | public DELTA{system['public_delta_layers']}, private DELTA{system['private_delta_layers']}, BML={system['bml_capacity']} |",
        "",
        "## Deterministic structural sizes",
        "",
        f"- Accumulator value: {acc_bytes} bytes ({acc_bytes / 1024:.2f} KiB).",
        f"- Compact accumulator witness (without right): {compact_witness_bytes} bytes ({compact_witness_bytes / 1024:.2f} KiB).",
        f"- Full accumulator witness: {full_witness_bytes} bytes ({full_witness_bytes / 1024:.2f} KiB).",
        f"- Accumulator Ba matrix: {ba_bytes} bytes ({ba_bytes / 1024:.2f} KiB) in q-bit packed form.",
        "- These are algebraic packed sizes; codec framing, signatures, compression and allocator memory are excluded.",
        "",
        "## Release blockers and warnings",
        "",
    ]
    lines.extend(f"- {warning}" for warning in warnings)
    if not warnings:
        lines.append("- None.")
    lines.extend(
        [
            "",
            "## Rule",
            "",
            "This records a reproducible candidate, not a 128-bit security "
            "claim, until the generator succeeds with --require-release.",
            "",
        ]
    )
    return "\n".join(lines)


def write_or_check(path: pathlib.Path, content: str, check: bool) -> None:
    content = content.rstrip() + "\n"
    if check:
        try:
            current = path.read_text(encoding="utf-8")
        except FileNotFoundError as exc:
            raise ProfileError(f"missing generated file: {path}") from exc
        if current != content:
            raise ProfileError(f"stale generated file: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--require-release", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    report = args.manifest.with_suffix(".generated.md")
    warnings = validate(manifest, require_release=args.require_release)
    digest = profile_id(manifest)
    write_or_check(C_HEADER, generated_header(manifest, digest, warnings), args.check)
    write_or_check(PY_MODULE, generated_python(manifest, digest, warnings), args.check)
    write_or_check(report, generated_report(manifest, digest, warnings), args.check)
    if not args.quiet:
        print(f"profile={manifest['profile_name']}")
        print(f"profile_id={digest.hex()}")
        print(f"release_status={manifest['release_status']}")
        print(f"release_ready={'yes' if not warnings and manifest['release_status'] == 'frozen' else 'no'}")
        for warning in warnings:
            print(f"warning: {warning}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, json.JSONDecodeError, KeyError, ProfileError) as exc:
        print(f"parameter profile error: {exc}", file=sys.stderr)
        raise SystemExit(1)
