#!/usr/bin/env python3
"""Import real AccGenpp MemVer coefficients and build the Labrador statement.

Input is produced by the C target `revac_acc_memver_export`.

The script first checks the exact integer lifting of the C relation

    left + Atrap*v12
    - sum_i trunc_B_i*v2_i
    + sum_i tag*b^e_i*v2_i
    + sum_{h,j} Ba[h][j]*right[h][j]
    - sum_{h,j} x_h*G_j*right[h][j]
    + q*lift
    = acc

over Z[X]/(X^256+1).  With `--mode labrador-simple`, it also appends the real
witnesses, public coefficients, quadratic terms, and lifting witnesses into a
Labrador `proof_statement` and runs `simple_verify`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import sys
import time
from io import BytesIO
from typing import Any

from revac_profile_generated import (
    D as PROFILE_D,
    GADGET_BASE as PROFILE_GADGET_BASE,
    GADGET_DIM as PROFILE_GADGET_DIM,
    GADGET_SIGNED as PROFILE_GADGET_SIGNED,
    HANDLE_BITS as PROFILE_HANDLE_BITS,
    HASH_BYTES,
    LINEAR_DECOMP_BASE,
    LINEAR_DECOMP_CHUNKS,
    LINK_OPENING_POLY_NORM_SQ,
    LABRADOR_COMPACT_ENCODING_VERSION,
    LABRADOR_RAW_ENCODING_VERSION,
    N as PROFILE_N,
    PROFILE_ID_HEX,
    PROOF_MODULUS as PROOF_MOD_38,
    Q as PROFILE_Q,
    RIGHT_FRESH_GROUP_NORM_SQ,
    RIGHT_GROUP_NORM_SQ_MAX,
    RIGHT_GROUP_NORM_SQ_PER_DELETION,
    SMALL_CHUNK_NORM_SQ,
    TAG_NORM_SQ,
    TRUNC_K_L as PROFILE_TRUNC_K_L,
    V2_NORM_SQ,
)


LAZER_PYTHON = pathlib.Path("/home/wkk/pq_ac/lazer/python")
sys.path.append(str(LAZER_PYTHON))


# MemWitUp adds EvalFX terms to the initially ternary right-side witness.
# Labrador groups the 32 handle-indexed polynomials for each (gadget, row)
# pair into one exact-norm vector.  The public deletion counter selects a fixed
# bound, so the norm is neither selected by nor leaked from an individual
# prover.  The cap is the largest policy supported reliably by the current
# 38-bit Labrador parameter generator.  This remains a prototype policy and
# must be replaced by the concrete accumulator beta from the final
# parameter/security analysis.
PROOF_MAGIC = (b"REVAC-LABRADOR38-MEMVER" +
               bytes([LABRADOR_RAW_ENCODING_VERSION]))
PROOF_MAGIC_COMPACT = (b"REVAC-LABRADOR38-MEMVER" +
                       bytes([LABRADOR_COMPACT_ENCODING_VERSION]))
PROOF_PROFILE_ID_BYTES = HASH_BYTES


def print_timing(label: str, started: float) -> None:
    """输出便于与 C 端阶段日志对应的细粒度耗时。"""
    print(f"[Labrador 计时] {label:<36} {time.perf_counter() - started:10.6f} 秒",
          flush=True)


def zero(n: int) -> list[int]:
    return [0] * n


def add(a: list[int], b: list[int]) -> list[int]:
    return [x + y for x, y in zip(a, b)]


def sub(a: list[int], b: list[int]) -> list[int]:
    return [x - y for x, y in zip(a, b)]


def scalar_mul(c: int, a: list[int]) -> list[int]:
    return [c * x for x in a]


def negacyclic_mul(a: list[int], b: list[int]) -> list[int]:
    n = len(a)
    out = [0] * n
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0:
                continue
            k = i + j
            if k >= n:
                out[k - n] -= ai * bj
            else:
                out[k] += ai * bj
    return out


def mat_vec_row(mat: list[list[list[int]]], row: int, vec: list[list[int]]) -> list[int]:
    acc = zero(len(vec[0]))
    for col, v in enumerate(vec):
        acc = add(acc, negacyclic_mul(mat[row][col], v))
    return acc


def const_poly(n: int, c: int) -> list[int]:
    out = [0] * n
    out[0] = c
    return out


def accumulator_gadget(params: dict[str, Any]) -> tuple[int, int, bool]:
    """Return and validate the accumulator gadget description.

    The defaults retain import compatibility with JSON emitted by the earlier
    binary-gadget exporter, which only included gadget_dim.
    """
    base = int(params.get("gadget_base", 2))
    dimension = int(params["gadget_dim"])
    signed = bool(params.get("gadget_signed", False))
    if base < 2:
        raise ValueError("accumulator gadget base must be at least two")
    if dimension < 1:
        raise ValueError("accumulator gadget dimension must be positive")
    if signed and base % 2 == 0:
        raise ValueError("balanced signed decomposition requires an odd gadget base")
    return base, dimension, signed


def right_group_norm_bound(params: dict[str, Any]) -> int:
    deletions = int(params.get("deletion_count", 0))
    if deletions < 0:
        raise ValueError("deletion_count cannot be negative")
    if deletions == 0:
        return RIGHT_FRESH_GROUP_NORM_SQ
    return min(
        RIGHT_GROUP_NORM_SQ_MAX,
        deletions * RIGHT_GROUP_NORM_SQ_PER_DELETION,
    )


def right_group_width(params: dict[str, Any]) -> int:
    """Return the experimental handle-axis width of one right witness group.

    Width 32 is the production layout.  Smaller divisors are an isolated
    parameter experiment: they leave the polynomial witness and relation
    unchanged, but split each Labrador norm vector into more groups.
    """
    handle_bits = int(params["handle_bits"])
    raw = os.environ.get("REVAC_RIGHT_GROUP_WIDTH", str(handle_bits))
    try:
        width = int(raw)
    except ValueError as exc:
        raise ValueError("REVAC_RIGHT_GROUP_WIDTH must be an integer") from exc
    if width < 1 or width > handle_bits or handle_bits % width != 0:
        raise ValueError(
            "REVAC_RIGHT_GROUP_WIDTH must be a positive divisor of handle_bits")
    return width


def right_group_ranges(params: dict[str, Any]) -> list[list[int]]:
    """Return the handle-axis partition used by grouped-right witnesses.

    ``REVAC_RIGHT_GROUP_HOT_PREFIX=k`` isolates the first k low handle bits
    and combines all colder suffix bits.  EvalFX activates bit h only after h
    matching low bits, so this non-uniform partition avoids combining the
    always-active bit zero with the geometrically colder suffix.
    """
    handle_bits = int(params["handle_bits"])
    raw_hot = os.environ.get("REVAC_RIGHT_GROUP_HOT_PREFIX", "0")
    try:
        hot = int(raw_hot)
    except ValueError as exc:
        raise ValueError(
            "REVAC_RIGHT_GROUP_HOT_PREFIX must be an integer") from exc
    if hot < 0 or hot >= handle_bits:
        raise ValueError(
            "REVAC_RIGHT_GROUP_HOT_PREFIX must be in [0, handle_bits)")
    if hot:
        return ([[bit] for bit in range(hot)] +
                [list(range(hot, handle_bits))])
    width = right_group_width(params)
    if os.environ.get("REVAC_RIGHT_GROUP_BALANCED", "0") != "0":
        group_count = handle_bits // width
        groups: list[list[int]] = [[] for _ in range(group_count)]
        weights = [0] * group_count
        for bit in range(handle_bits):
            candidates = [index for index, group in enumerate(groups)
                          if len(group) < width]
            target = min(candidates, key=lambda index: (weights[index], index))
            groups[target].append(bit)
            weights[target] += 1 << (handle_bits - 1 - bit)
        return groups
    return [list(range(first, first + width))
            for first in range(0, handle_bits, width)]


def memver_labrador_profile_id(data: dict[str, Any]) -> bytes:
    """Return a canonical ID for locally reproducible public proof params.

    Native parameter serialization contains process-local pointer fields, so
    hashing those bytes would not be stable across a UAV and GCS.  This ID is
    instead derived from every public dimension/norm policy that feeds the
    deterministic parameter generator plus an explicit compiler-profile
    version.  A verifier still generates the expected parameters itself; the
    ID is version/profile binding, not a substitute for verification.
    """
    params = data["params"]
    profile = {
        "profile": "revac-memver-labrador38-grouped-zk-native-sparse-v4",
        "revocable_parameter_profile_id": PROFILE_ID_HEX,
        "n": int(params["n"]),
        "d": int(params["d"]),
        "q": int(params["q"]),
        "handle_bits": int(params["handle_bits"]),
        "gadget_base": int(params["gadget_base"]),
        "gadget_dim": int(params["gadget_dim"]),
        "gadget_signed": bool(params.get("gadget_signed", False)),
        "trunc_k_l": int(params["trunc_k_l"]),
        "deletion_count": int(params.get("deletion_count", 0)),
        "refresh_bits": int(params.get("refresh_bits", 0)),
        "proof_mod": PROOF_MOD_38,
        "linear_decomp_base": LINEAR_DECOMP_BASE,
        "linear_decomp_chunks": LINEAR_DECOMP_CHUNKS,
        "small_chunk_norm_sq": SMALL_CHUNK_NORM_SQ,
        "tag_norm_sq": TAG_NORM_SQ,
        "v2_norm_sq": V2_NORM_SQ,
        "right_group_norm_sq": right_group_norm_bound(params),
        "lnp_z1_sigma_rounding": os.environ.get(
            "REVAC_LNP_Z1_SIGMA", "nearest"),
    }
    ranges = right_group_ranges(params)
    if ranges != [list(range(int(params["handle_bits"])))]:
        # Experimental non-default norm partitions must obtain independent
        # proof parameters.  Omitting this field for the production 32-bit
        # group preserves the established profile and wire compatibility.
        profile["right_group_ranges"] = ranges
    canonical = json.dumps(profile, sort_keys=True,
                           separators=(",", ":")).encode("ascii")
    return hashlib.sha256(canonical).digest()


class MemverStaticTemplate:
    """Reusable public part of the grouped MemVer statement.

    A Labrador statement cannot be copied safely through the current public
    CFFI API: its right-hand side and all witness vectors are request-specific.
    This template therefore caches the immutable CFFI polynomials, vectors,
    selector relations, and decomposition coefficients.  Each request still
    creates a fresh proof_statement/witness and fresh ZK randomness.
    """

    def __init__(self, data: dict[str, Any]):
        from lazer import polyring_t, poly_t, polyvec_t

        params = data["params"]
        self.right_group_ranges = right_group_ranges(params)
        if (params["n"], params["d"], params["q"],
                params["handle_bits"], params["gadget_base"],
                params["gadget_dim"], bool(params.get("gadget_signed", False)),
                params["trunc_k_l"]) != (
                    PROFILE_N, PROFILE_D, PROFILE_Q, PROFILE_HANDLE_BITS,
                    PROFILE_GADGET_BASE, PROFILE_GADGET_DIM,
                    PROFILE_GADGET_SIGNED, PROFILE_TRUNC_K_L):
            raise ValueError("MemVer static template parameters do not match this build")
        self.ring = polyring_t(params["n"], PROOF_MOD_38)
        self.one = poly_t(self.ring, const_poly(params["n"], 1))
        self.minus_one = poly_t(self.ring, const_poly(params["n"], -1))
        self.zero = poly_t(self.ring, zero(params["n"]))
        self.bit_poly = (self.zero, self.one)
        self.selectors = []
        for group_range in self.right_group_ranges:
            group_selectors = []
            for local_h in range(len(group_range)):
                selector = [self.zero] * len(group_range)
                selector[local_h] = self.one
                group_selectors.append(polyvec_t(
                    self.ring, len(group_range), selector))
            self.selectors.append(group_selectors)

        self.scale_poly: dict[int, Any] = {}
        self.q_scale_poly: dict[int, Any] = {}
        scale = 1
        for _ in range(LINEAR_DECOMP_CHUNKS):
            self.scale_poly[scale] = poly_t(self.ring, const_poly(params["n"], scale))
            self.q_scale_poly[scale] = poly_t(self.ring, const_poly(params["n"], params["q"] * scale))
            scale *= LINEAR_DECOMP_BASE

        self.atrap_scaled = []
        for scale in self.scale_poly:
            self.atrap_scaled.append([
                [poly_t(self.ring, scalar_mul(scale, data["Atrap"][row][col]))
                 for col in range(params["d"])]
                for row in range(params["d"])
            ])
        self.trunc_neg = [
            [[poly_t(self.ring, scalar_mul(-1, data["trunc_B"][i][row][col]))
              for col in range(params["d"])]
             for row in range(params["d"])]
            for i in range(params["trunc_k_l"])
        ]
        self.ba_grouped = [
            [[[
                polyvec_t(
                    self.ring, len(group_range),
                    [poly_t(self.ring, data["Ba"][h][j][row][col])
                     for h in group_range],
                )
                for col in range(params["d"])]
              for row in range(params["d"])]
             for group_range in self.right_group_ranges]
            for j in range(params["gadget_dim"])
        ]
        self.link_a = [
            [poly_t(self.ring, data["link_A"][row][col])
             for col in range(params["d"])]
            for row in range(params["d"])
        ]
        self.link_d_bits = [
            [poly_t(self.ring,
                    negacyclic_monomial_mul(data["link_D"][row][0], bit))
             for bit in range(params["handle_bits"])]
            for row in range(params["d"])
        ]
        refresh_bits = int(params.get("refresh_bits", 0))
        self.refresh_coeff = []
        self.refresh_selectors = []
        refresh_a = data.get("refresh_A")
        native_sparse_enabled = os.environ.get(
            "REVAC_NATIVE_SPARSE_REFRESH", "1") != "0"
        if (refresh_bits and
                (not hasattr(refresh_a, "coeffs") or
                 not native_sparse_enabled)):
            for row in range(params["d"]):
                row_coeff = []
                row_selectors = []
                for gadget in range(params["gadget_dim"]):
                    for col in range(params["d"]):
                        row_coeff.append(polyvec_t(
                            self.ring,
                            refresh_bits,
                            [poly_t(self.ring,
                                    refresh_matrix_poly(
                                        refresh_a, i + 1, gadget, row, col))
                             for i in range(refresh_bits)],
                        ))
                        selector = [self.zero] * (
                            params["gadget_dim"] * params["d"])
                        selector[gadget * params["d"] + col] = self.one
                        row_selectors.append(polyvec_t(
                            self.ring, len(selector), selector))
                self.refresh_coeff.append(row_coeff)
                self.refresh_selectors.append(row_selectors)
        # Filled on the first request after deterministic witness indices have
        # been assigned.  It is the compiler output for the fixed MemVer
        # topology: static linear/quadratic operands plus witness-index wiring.
        # Dynamic accumulator/commitment RHS values are deliberately excluded.
        self.compiled_relation: dict[str, Any] | None = None
        # Raw CFFI form of the same topology.  This avoids rebuilding large
        # temporary polyvec_t objects and converting them back to int64 arrays
        # on every prove/verify request.
        self.compiled_raw_relation: dict[str, Any] | None = None
        # Stable native statements keyed by the complete witness/norm layout.
        # Each entry owns its statement and witness storage for the lifetime of
        # the persistent local worker.
        self.live_statements: dict[tuple[Any, ...], tuple[Any, bytes, int]] = {}
        # Exact dynamic RHS values currently bound to each live statement.
        # Prove is immediately followed by verify for the same show packet, so
        # the verifier can normally reuse the already-bound native statement.
        # Store full integer tuples (not a hash) so a collision can never turn
        # this performance cache into a verification bypass.
        self.live_dynamic_rhs: dict[tuple[Any, ...], tuple[Any, ...]] = {}
        # Exact ordinary-witness snapshots associated with each native
        # statement.  Consecutive shows for one credential keep almost all
        # accumulator witness vectors unchanged; retaining exact tuples and
        # their already-computed norms lets the worker update only changed
        # vectors.  Native-generated refresh/lift slots are never skipped.
        self.live_witness_keys: dict[tuple[Any, ...], list[Any]] = {}
        self.live_witness_norms: dict[tuple[Any, ...], list[int | None]] = {}
        self.sealed = False

    def compatible(self, data: dict[str, Any]) -> bool:
        params = data["params"]
        return (params["n"], params["d"], params["q"],
                params["handle_bits"], params["gadget_base"],
                params["gadget_dim"], bool(params.get("gadget_signed", False)),
                params["trunc_k_l"], params.get("refresh_bits", 0),
                tuple(tuple(group)
                      for group in right_group_ranges(params))) == (
                    PROFILE_N, PROFILE_D, PROFILE_Q, PROFILE_HANDLE_BITS,
                    PROFILE_GADGET_BASE, PROFILE_GADGET_DIM,
                    PROFILE_GADGET_SIGNED, PROFILE_TRUNC_K_L,
                    int(data["params"].get("refresh_bits", 0)),
                    tuple(tuple(group)
                          for group in self.right_group_ranges))

    def seal_native_statement(self, layout: tuple[int, ...]) -> None:
        """Release compiler-only Python/CFFI objects after native freezing.

        The resident native statement owns copies of every immutable
        coefficient.  Retaining the Python polynomial graph, the high-level
        relation plan and the static raw arrays after that point only keeps
        duplicate representations alive.
        """
        if self.sealed:
            return
        self.compiled_relation = {"layout": layout}
        if self.compiled_raw_relation is not None:
            self.compiled_raw_relation["bit_links"] = []
            self.compiled_raw_relation["refresh_links"] = []
        self.selectors = []
        self.scale_poly = {}
        self.q_scale_poly = {}
        self.atrap_scaled = []
        self.trunc_neg = []
        self.ba_grouped = []
        self.link_a = []
        self.link_d_bits = []
        self.refresh_coeff = []
        self.refresh_selectors = []
        self.bit_poly = ()
        self.one = None
        self.minus_one = None
        self.zero = None
        self.sealed = True


def refresh_matrix_poly(bank: Any, matrix: int, gadget: int,
                        row: int, col: int) -> list[int]:
    """Read one refresh-bank polynomial from compact or legacy storage."""
    if hasattr(bank, "poly"):
        return bank.poly(matrix, gadget, row, col)
    return bank[matrix][gadget][row][col]


class _NativeWitness:
    """Placeholder for a vector filled directly by the native bridge."""

    __slots__ = ()


def _dynamic_rhs_key(data: dict[str, Any]) -> tuple[Any, ...]:
    """Return an exact, immutable key for all request-specific RHS rows."""
    return (
        tuple(tuple(int(value) for value in row)
              for row in data["link_commitment"]),
        tuple(tuple(int(value) for value in row)
              for row in data["acc"]),
    )


def _dynamic_bind_cache_enabled() -> bool:
    value = os.environ.get("REVAC_LABRADOR_DYNAMIC_BIND_CACHE", "1")
    return value != "0"


def _compile_raw_linear(ps: Any, stat_list: list[Any],
                        wit_list: list[int], right_pol: Any) -> dict[str, Any]:
    """Flatten one immutable linear relation once for the native C API."""
    from lazer import ffi, poly_t, polyvec_t
    from labrados import poly_to_ar, polyvec_to_ar

    lengths: list[int] = []
    stat_size = 0
    for item in stat_list:
        if type(item) is poly_t:
            lengths.append(1)
            stat_size += 1
        elif type(item) is polyvec_t:
            lengths.append(item.dim)
            stat_size += item.dim
        else:
            raise TypeError("unsupported Labrador linear coefficient")
    ring = right_pol.ring
    stat_vec = polyvec_t(ring, stat_size, stat_list)
    return {
        "nvec": len(stat_list),
        "idx": ffi.new("size_t[]", [ps.map[idx] for idx in wit_list]),
        "lengths": ffi.new("size_t[]", lengths),
        "degree": ring.deg // 256,
        "phi": polyvec_to_ar(stat_vec),
        "rhs": poly_to_ar(right_pol),
    }


def _append_raw_linear(ps: Any, raw: dict[str, Any],
                       rhs_coeffs: list[int] | None = None) -> None:
    """Append a pre-flattened relation, optionally replacing its RHS."""
    from lazer import ffi, lib

    rhs = (raw["rhs"] if rhs_coeffs is None else
           ffi.new("int64_t[]", rhs_coeffs))
    output = lib.labrador38_py_append_constraint(
        ps.statement, raw["nvec"], raw["idx"], raw["lengths"],
        raw["degree"], raw["phi"], rhs, 1)
    if output != 0:
        raise ValueError(f"native cached linear constraint failed: {output}")
    ps.cur_constraint_num += 1


def _compile_raw_quadratic(ps: Any, relation: tuple[Any, ...],
                           right_pol: Any) -> dict[str, Any]:
    """Flatten immutable linear/quadratic operands once for the native API."""
    from lazer import ffi, poly_t, polyvec_t
    from labrados import poly_to_ar, polyvec_to_ar

    stat_lin, wit_lin, stat_quad, wit_quad1, wit_quad2 = relation
    lengths: list[int] = []
    stat_size = 0
    for item in stat_lin:
        if type(item) is poly_t:
            lengths.append(1)
            stat_size += 1
        elif type(item) is polyvec_t:
            lengths.append(item.dim)
            stat_size += item.dim
        else:
            raise TypeError("unsupported Labrador quadratic linear coefficient")
    ring = right_pol.ring
    stat_vec = polyvec_t(ring, stat_size, stat_lin)
    return {
        "nlin": len(stat_lin),
        "nprod": len(stat_quad),
        "wit_lin": ffi.new("size_t[]", [ps.map[idx] for idx in wit_lin]),
        "wit_quad1": ffi.new("size_t[]", [ps.map[idx] for idx in wit_quad1]),
        "wit_quad2": ffi.new("size_t[]", [ps.map[idx] for idx in wit_quad2]),
        "lengths": ffi.new("size_t[]", lengths),
        "degree": ring.deg // 256,
        "quad": ffi.new("int64_t[]", stat_quad),
        "phi": polyvec_to_ar(stat_vec),
        "rhs": poly_to_ar(right_pol),
    }


def _append_raw_quadratic(ps: Any, raw: dict[str, Any],
                          rhs_coeffs: list[int]) -> None:
    from lazer import ffi, lib

    rhs = ffi.new("int64_t[]", rhs_coeffs)
    output = lib.labrador38_py_append_quadratic(
        ps.statement, raw["nlin"], raw["nprod"], raw["wit_lin"],
        raw["wit_quad1"], raw["wit_quad2"], raw["lengths"],
        raw["degree"], raw["quad"], raw["phi"], rhs)
    if output != 0:
        raise ValueError(f"native cached quadratic constraint failed: {output}")
    ps.cur_constraint_num += 1


def signed_decompose_poly(poly: list[int],
                          base: int = LINEAR_DECOMP_BASE,
                          chunks: int = LINEAR_DECOMP_CHUNKS) -> list[list[int]]:
    """Return small signed base chunks whose weighted sum is poly."""
    out: list[list[int]] = []
    vals = list(poly)
    for _ in range(chunks):
        digit_poly: list[int] = []
        next_vals: list[int] = []
        for v in vals:
            digit = ((v + base // 2) % base) - base // 2
            digit_poly.append(digit)
            next_vals.append((v - digit) // base)
        out.append(digit_poly)
        vals = next_vals
    if any(vals):
        raise ValueError("linear witness coefficient is too large for fixed decomposition")
    return out


def negacyclic_monomial_mul(poly: list[int], degree: int) -> list[int]:
    """Return X^degree * poly in Z[X]/(X^n + 1)."""
    n = len(poly)
    out = [0] * n
    for idx, value in enumerate(poly):
        target = idx + degree
        if target < n:
            out[target] += value
        else:
            out[target - n] -= value
    return out


def l2sq(poly: list[int]) -> int:
    return sum(c * c for c in poly)


def right_witness_norms(data: dict[str, Any]) -> tuple[int, int, int]:
    """Return max polynomial/group norms and the aggregate right-side norm."""
    right = data["witness"]["right"]
    group_ranges = right_group_ranges(data["params"])
    norms = [
        l2sq(poly)
        for handle_row in right
        for gadget_vec in handle_row
        for poly in gadget_vec
    ]
    gadget_dim = len(right[0])
    rows = len(right[0][0])
    group_norms = [
        sum(l2sq(right[h][j][row]) for h in group_range)
        for j in range(gadget_dim)
        for row in range(rows)
        for group_range in group_ranges
    ]
    return max(norms, default=0), max(group_norms, default=0), sum(norms)


def append_poly(ps: Any, poly_t: Any, ring: Any, polys: list[list[int]], norms: list[int], p: list[int], norm: int) -> int:
    polys.append(p)
    norms.append(norm)
    return len(polys) - 1


def compute_lift(data: dict[str, Any]) -> tuple[list[list[int]], list[list[int]]]:
    params = data["params"]
    n = params["n"]
    q = params["q"]
    handle_bits = params["handle_bits"]
    gadget_base, gadget_dim, _ = accumulator_gadget(params)
    trunc_k_l = params["trunc_k_l"]
    w = data["witness"]
    acc = data["acc"]
    x_bits = data["x_bits"]
    tag = w["tag"]

    expr_rows: list[list[int]] = []
    lifts: list[list[int]] = []
    for row in range(params["d"]):
        bexpi = 196
        expr = list(w["left"][row])
        expr = add(expr, mat_vec_row(data["Atrap"], row, w["v12"]))
        for i in range(trunc_k_l):
            expr = sub(expr, mat_vec_row(data["trunc_B"][i], row, w["v2"][i]))
            expr = add(expr, scalar_mul(bexpi, negacyclic_mul(tag, w["v2"][i][row])))
            bexpi *= 14
        for h in range(handle_bits):
            for j in range(gadget_dim):
                expr = add(expr, mat_vec_row(data["Ba"][h][j], row, w["right"][h][j]))
                expr = sub(
                    expr,
                    scalar_mul((gadget_base**j) * x_bits[h], w["right"][h][j][row]),
                )
        diff = sub(acc[row], expr)
        if "memver_lhs" in data:
            lhs_delta = sub(data["memver_lhs"][row], expr)
            lhs_bad = [c for c in lhs_delta if c % q != 0]
            if lhs_bad:
                raise ValueError(
                    f"row {row}: Python lhs disagrees with exported C lhs; "
                    f"first non-q-multiple delta {lhs_bad[0]}"
                )
        bad = [c for c in diff if c % q != 0]
        if bad:
            raise ValueError(f"row {row}: lifting is not divisible by q; first bad coefficient {bad[0]}")
        lifts.append([c // q for c in diff])
        expr_rows.append(expr)
    return expr_rows, lifts


def compute_link_lift(data: dict[str, Any]) -> list[list[int]]:
    """Exact lift of r0 + A_link*r1 + D_link*x = C_link over R_q."""
    params = data["params"]
    q = int(params["q"])
    d = int(params["d"])
    n = int(params["n"])
    w = data["witness"]
    r0, r1 = w["link_r"]
    x_poly = [0] * n
    for bit, value in enumerate(data["x_bits"]):
        x_poly[bit] = int(value)
    x_vec = [x_poly] + [[0] * n for _ in range(d - 1)]
    lifts: list[list[int]] = []
    for row in range(d):
        expr = list(r0[row])
        expr = add(expr, mat_vec_row(data["link_A"], row, r1))
        expr = add(expr, mat_vec_row(data["link_D"], row, x_vec))
        diff = sub(data["link_commitment"][row], expr)
        bad = [value for value in diff if value % q != 0]
        if bad:
            raise ValueError(
                "same-handle commitment lifting is not divisible by q; "
                f"first bad coefficient {bad[0]}")
        lifts.append([value // q for value in diff])
    return lifts


def load_or_compute_lift(data: dict[str, Any]) -> tuple[list[list[int]], str]:
    """Prefer the exact C lift; retain Python as a slow reference backend."""
    raw = data.get("lift")
    if raw is None:
        _, lifts = compute_lift(data)
        return lifts, "Python 参考实现"

    d = int(data["params"]["d"])
    n = int(data["params"]["n"])
    if not isinstance(raw, list) or len(raw) != d:
        raise ValueError("C exact-lift 的行数错误")
    lifts: list[list[int]] = []
    for row, poly in enumerate(raw):
        if not isinstance(poly, list) or len(poly) != n:
            raise ValueError(f"C exact-lift 第 {row} 行的系数数目错误")
        if any(not isinstance(value, int) or isinstance(value, bool)
               for value in poly):
            raise ValueError(f"C exact-lift 第 {row} 行包含非整数系数")
        lifts.append(list(poly))

    if os.getenv("REVAC_LIFT_PYTHON_CROSSCHECK", "0") == "1":
        _, reference = compute_lift(data)
        if lifts != reference:
            raise ValueError("C exact-lift 与 Python 参考实现不一致")
        return lifts, "C exact-lift（已与 Python 参考实现一致性检查）"
    return lifts, "C exact-lift"


def _write_u64(out: BytesIO, value: int) -> None:
    out.write(int(value).to_bytes(8, "little", signed=False))


def _read_u64(inp: BytesIO) -> int:
    raw = inp.read(8)
    if len(raw) != 8:
        raise ValueError("truncated proof bytes")
    return int.from_bytes(raw, "little", signed=False)


def _write_blob(out: BytesIO, blob: bytes) -> None:
    _write_u64(out, len(blob))
    out.write(blob)


def _read_blob(inp: BytesIO) -> bytes:
    size = _read_u64(inp)
    raw = inp.read(size)
    if len(raw) != size:
        raise ValueError("truncated proof bytes")
    return raw


def _lab_polz_len(pp: Any) -> int:
    if pp.compressed:
        return 4 * pp.kappa[2]
    if not pp.tail:
        return 2 * pp.kappa[1] + 3
    return pp.len[1] + pp.len[2] + pp.len[8] + pp.len[4]


def _set_lab_m_pointers(lab_proof: Any, pp: Any, arr: Any) -> None:
    lab_proof.m[0] = arr
    if pp.compressed:
        lab_proof.m[1] = arr + pp.kappa[2]
        lab_proof.m[2] = lab_proof.m[1] + pp.kappa[2]
        lab_proof.m[3] = lab_proof.m[2] + pp.kappa[2]
    elif not pp.tail:
        lab_proof.m[1] = None
        lab_proof.m[2] = arr + pp.kappa[1]
        lab_proof.m[3] = lab_proof.m[2] + 3
    else:
        lab_proof.m[1] = None
        lab_proof.m[2] = arr + pp.len[1] + pp.len[2]
        lab_proof.m[3] = lab_proof.m[2] + pp.len[8]


def _lnp_polz_len(pp: Any) -> int:
    return 4 * pp.kappa_linfmsis + 4


def _set_lnp_m_pointers(lnp_proof: Any, pp: Any, arr: Any) -> None:
    lnp_proof.m[0] = arr
    lnp_proof.m[1] = arr + pp.kappa_linfmsis
    lnp_proof.m[2] = arr + 2 * pp.kappa_linfmsis
    lnp_proof.m[3] = arr + 2 * pp.kappa_linfmsis + 4
    lnp_proof.m[4] = arr + 3 * pp.kappa_linfmsis + 4


def serialize_labrador38_params(params: Any, require_zk: bool = True) -> bytes:
    """Serialize reusable public proof parameters only, never a witness/proof."""
    from lazer import ffi, lib

    if require_zk and lib.labrador38_py_params_is_zk(params) != 1:
        raise ValueError("拒绝缓存非零知识 Labrador 参数")
    params_nbytes = lib.labrador38_py_params_nbytes(params)
    out = ffi.new("uint8_t[]", params_nbytes)
    if lib.labrador38_py_params_serialize(out, params_nbytes, params) != 0:
        raise ValueError("Labrador 参数序列化失败")
    return bytes(ffi.buffer(out, params_nbytes))


class LiveLabrador38Params:
    """Own one deserialized/generated public parameter object in a worker.

    The old cache retained an immutable byte string and deserialized the same
    public object for every proof.  A persistent worker is single-requested, so
    it can safely retain the native read-only parameter tree itself.  Proofs and
    all ZK randomness remain fresh per request.
    """

    def __init__(self) -> None:
        from lazer import ffi

        self.params = ffi.new("labrador38_dch_pack_params")
        self.refs: list[Any] = []
        self._serialized: bytes | None = None
        self.closed = False

    @classmethod
    def from_bytes(cls, params_bytes: bytes,
                   require_zk: bool = True,
                   reset_comkey: bool = False) -> "LiveLabrador38Params":
        from lazer import ffi, lib

        obj = cls()
        if reset_comkey:
            lib.labrador38_py_comkey_reset()
        raw = ffi.new("uint8_t[]", params_bytes)
        obj.refs.append(raw)
        if lib.labrador38_py_params_deserialize(
                obj.params, raw, len(params_bytes)) != 0:
            obj.closed = True  # Deserializer owns cleanup on its error path.
            raise ValueError("缓存的 Labrador 参数反序列化失败")
        if require_zk and lib.labrador38_py_params_is_zk(obj.params) != 1:
            obj.close()
            raise ValueError("缓存参数不是 zk=True Labrador 参数")
        obj._serialized = bytes(params_bytes)
        return obj

    def checked_params(self, require_zk: bool = True) -> Any:
        from lazer import lib

        if self.closed:
            raise ValueError("尝试复用已经释放的 Labrador 参数")
        if require_zk and lib.labrador38_py_params_is_zk(self.params) != 1:
            raise ValueError("常驻参数不是 zk=True Labrador 参数")
        return self.params

    def serialized_bytes(self, require_zk: bool = True) -> bytes:
        self.checked_params(require_zk=require_zk)
        if self._serialized is None:
            self._serialized = serialize_labrador38_params(
                self.params, require_zk=require_zk)
        return self._serialized

    def close(self) -> None:
        if self.closed:
            return
        from lazer import lib

        lib.labrador38_py_free_params(self.params)
        self.refs.clear()
        self._serialized = None
        self.closed = True

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            # Interpreter shutdown may unload the CFFI module first.
            pass


def deserialize_labrador38_params(params_bytes: bytes,
                                  require_zk: bool = True,
                                  reset_comkey: bool = True) -> tuple[Any, list[Any]]:
    """Create a fresh CFFI parameter object from a cached immutable blob."""
    from lazer import ffi, lib

    if reset_comkey:
        lib.labrador38_py_comkey_reset()
    params = ffi.new("labrador38_dch_pack_params")
    raw = ffi.new("uint8_t[]", params_bytes)
    if lib.labrador38_py_params_deserialize(params, raw, len(params_bytes)) != 0:
        raise ValueError("缓存的 Labrador 参数反序列化失败")
    if require_zk and lib.labrador38_py_params_is_zk(params) != 1:
        lib.labrador38_py_free_params(params)
        raise ValueError("缓存参数不是 zk=True Labrador 参数")
    return params, [raw]


def serialize_labrador38_proof(params: Any, proof: Any,
                               require_zk: bool = True,
                               cached_params_bytes: bytes | None = None,
                               compact_profile_id: bytes | None = None) -> bytes:
    from lazer import ffi, lib

    if require_zk and lib.labrador38_py_params_is_zk(params) != 1:
        raise ValueError("拒绝序列化非零知识 Labrador proof")

    if compact_profile_id is not None:
        if len(compact_profile_id) != PROOF_PROFILE_ID_BYTES:
            raise ValueError("Labrador compact profile ID 长度错误")
        params_bytes = b""
    elif cached_params_bytes is None:
        params_nbytes = lib.labrador38_py_params_nbytes(params)
        params_out = ffi.new("uint8_t[]", params_nbytes)
        err = lib.labrador38_py_params_serialize(
            params_out, params_nbytes, params)
        if err != 0:
            raise ValueError(f"Labrador params serialization failed: {err}")
        params_bytes = bytes(ffi.buffer(params_out, params_nbytes))
    else:
        params_bytes = cached_params_bytes

    proof_nbytes = lib.labrador38_py_proof_nbytes(params, proof)
    proof_out = ffi.new("uint8_t[]", proof_nbytes)
    err = lib.labrador38_py_proof_serialize(proof_out, proof_nbytes, params, proof)
    if err != 0:
        raise ValueError(f"Labrador proof serialization failed: {err}")

    out = BytesIO()
    if compact_profile_id is not None:
        out.write(PROOF_MAGIC_COMPACT)
        out.write(compact_profile_id)
    else:
        out.write(PROOF_MAGIC)
        _write_blob(out, params_bytes)
    _write_blob(out, bytes(ffi.buffer(proof_out, proof_nbytes)))
    return out.getvalue()


def deserialize_labrador38_proof(
        proof_bytes: bytes, require_zk: bool = True,
        reset_comkey: bool = True,
        cached_params_obj: LiveLabrador38Params | None = None,
        expected_profile_id: bytes | None = None,
) -> tuple[Any, Any, list[Any]]:
    from lazer import ffi, lib

    inp = BytesIO(proof_bytes)
    magic = inp.read(len(PROOF_MAGIC))
    compact = magic == PROOF_MAGIC_COMPACT
    if magic == PROOF_MAGIC:
        params_bytes = _read_blob(inp)
        profile_id = None
    elif compact:
        profile_id = inp.read(PROOF_PROFILE_ID_BYTES)
        if len(profile_id) != PROOF_PROFILE_ID_BYTES:
            raise ValueError("truncated Labrador compact profile ID")
        params_bytes = b""
    else:
        raise ValueError("bad revocable Labrador proof envelope magic")
    inner_proof_bytes = _read_blob(inp)
    if inp.read(1):
        raise ValueError("trailing data after revocable Labrador proof envelope")

    refs: list[Any] = []
    if reset_comkey:
        lib.labrador38_py_comkey_reset()

    if cached_params_obj is not None:
        if compact:
            if expected_profile_id is None or profile_id != expected_profile_id:
                raise ValueError("compact proof profile ID 与 GCS 公开参数策略不一致")
        elif params_bytes != cached_params_obj.serialized_bytes(
                require_zk=require_zk):
                raise ValueError("proof 携带的 Labrador 参数与 GCS 常驻缓存不一致")
        params = cached_params_obj.checked_params(require_zk=require_zk)
    else:
        if compact:
            raise ValueError("compact proof 需要 GCS 本地生成的公开参数")
        params = ffi.new("labrador38_dch_pack_params")
        raw_params = ffi.new("uint8_t[]", params_bytes)
        refs.append(raw_params)
        err = lib.labrador38_py_params_deserialize(
            params, raw_params, len(params_bytes))
        if err != 0:
            raise ValueError(f"Labrador params deserialization failed: {err}")
        if require_zk and lib.labrador38_py_params_is_zk(params) != 1:
            lib.labrador38_py_free_params(params)
            raise ValueError("GCS 拒绝：收到的 Labrador MemVer proof 不是零知识证明")

    proof = ffi.new("labrador38_dch_pack_proof")
    # The deserializer treats its input as read-only.  Borrow the immutable
    # Python buffer instead of allocating and copying the whole (roughly
    # 165-KiB) proof at every GCS verification.  Keep both objects alive until
    # native deserialization has completed.
    raw_owner = inner_proof_bytes
    raw = ffi.from_buffer("const uint8_t[]", raw_owner)
    refs.extend((raw_owner, raw))
    err = lib.labrador38_py_proof_deserialize(proof, params, raw, len(inner_proof_bytes))
    if err != 0:
        raise ValueError(f"Labrador proof deserialization failed: {err}")
    return params, proof, refs


def verify_labrador_live_statement(
        data: dict[str, Any], proof_bytes: bytes,
        static_template: MemverStaticTemplate,
        cached_params_obj: LiveLabrador38Params,
        reuse_comkey: bool = True) -> bool | None:
    """Verify against an already-compiled native MemVer statement.

    The generic compiler must run once for each statement/parameter profile to
    establish witness dimensions, norm bounds and the immutable constraint
    skeleton.  Once that profile is active, a verifier has no witness to
    decompose: only four link-commitment RHS rows and four accumulator RHS rows
    vary per request.  Rebind those rows directly and enter Pack Verify.

    ``None`` means that the reusable native state is unavailable and the
    caller must use ``build_labrador``.  Cryptographic rejection is returned as
    ``False`` and must never trigger a weaker fallback verifier.
    """
    if (not static_template.compatible(data) or
            static_template.compiled_raw_relation is None or
            len(static_template.live_statements) != 1):
        return None

    from lazer import ffi, lib
    from labrados import pack_verify

    statement_key, live_entry = next(
        iter(static_template.live_statements.items()))
    ps, static_hash, dynamic_first = live_entry
    raw_relation = static_template.compiled_raw_relation
    link_rows = raw_relation.get("link_rows", [])
    accumulator_rows = raw_relation.get("accumulator_rows", [])
    if len(link_rows) != len(data["link_commitment"]):
        return None
    if len(accumulator_rows) != len(data["acc"]):
        return None

    dynamic_key = _dynamic_rhs_key(data)
    bind_hit = (_dynamic_bind_cache_enabled() and
                static_template.live_dynamic_rhs.get(statement_key) ==
                dynamic_key)
    if not bind_hit:
        lib.labrador38_py_statement_hash_reset(
            ps.statement, ffi.new("uint8_t[16]", static_hash))
        constraint_index = dynamic_first
        for row, raw in enumerate(link_rows):
            rhs = ffi.new("int64_t[]", data["link_commitment"][row])
            output = lib.labrador38_py_statement_rebind_constraint_raw(
                ps.statement, constraint_index, raw["nvec"], raw["lengths"],
                raw["degree"], raw["phi"], rhs, 0, ffi.NULL)
            if output != 0:
                raise ValueError(
                    f"常驻 MemVer link statement 快速绑定失败：{output}")
            constraint_index += 1
        for row, raw in enumerate(accumulator_rows):
            rhs = ffi.new("int64_t[]", data["acc"][row])
            output = lib.labrador38_py_statement_rebind_constraint_raw(
                ps.statement, constraint_index, raw["nlin"], raw["lengths"],
                raw["degree"], raw["phi"], rhs, raw["nprod"], raw["quad"])
            if output != 0:
                raise ValueError(
                    f"常驻 MemVer accumulator statement 快速绑定失败：{output}")
            constraint_index += 1
        static_template.live_dynamic_rhs[statement_key] = dynamic_key
    if os.environ.get("REVAC_LABRADOR_PROFILE", "0") != "0":
        print("[Labrador worker] 动态 RHS 绑定缓存：" +
              ("命中" if bind_hit else "未命中"))

    expected_profile_id = memver_labrador_profile_id(data)
    if (labrador38_proof_is_compact(proof_bytes) and
            labrador38_proof_profile_id(proof_bytes) != expected_profile_id):
        return False
    params, proof, refs = deserialize_labrador38_proof(
        proof_bytes, reset_comkey=not reuse_comkey,
        cached_params_obj=cached_params_obj,
        expected_profile_id=expected_profile_id)
    try:
        return bool(pack_verify(
            (params, proof), ps.output_statement(), "38"))
    finally:
        lib.labrador38_py_free_proof(proof)
        _ = refs


def labrador38_proof_params_bytes(proof_bytes: bytes) -> bytes:
    """Extract and validate the immutable public-parameter envelope field."""
    inp = BytesIO(proof_bytes)
    magic = inp.read(len(PROOF_MAGIC))
    if magic == PROOF_MAGIC_COMPACT:
        raise ValueError("compact Labrador proof 不携带公共参数字节")
    if magic != PROOF_MAGIC:
        raise ValueError("bad revocable Labrador proof envelope magic")
    params_bytes = _read_blob(inp)
    _ = _read_blob(inp)
    if inp.read(1):
        raise ValueError("trailing data after revocable Labrador proof envelope")
    return params_bytes


def labrador38_proof_is_compact(proof_bytes: bytes) -> bool:
    return proof_bytes.startswith(PROOF_MAGIC_COMPACT)


def labrador38_proof_profile_id(proof_bytes: bytes) -> bytes | None:
    inp = BytesIO(proof_bytes)
    magic = inp.read(len(PROOF_MAGIC))
    if magic == PROOF_MAGIC:
        return None
    if magic != PROOF_MAGIC_COMPACT:
        raise ValueError("bad revocable Labrador proof envelope magic")
    profile_id = inp.read(PROOF_PROFILE_ID_BYTES)
    if len(profile_id) != PROOF_PROFILE_ID_BYTES:
        raise ValueError("truncated Labrador compact profile ID")
    _ = _read_blob(inp)
    if inp.read(1):
        raise ValueError("trailing data after revocable Labrador proof envelope")
    return profile_id


def build_labrador(data: dict[str, Any],
                   lifts: list[list[int]] | None,
                   proof_out: pathlib.Path | None = None,
                   proof_in: pathlib.Path | None = None,
                   debug_hash: bool = False,
                   self_verify: bool = True,
                   right_layout: str = "grouped",
                   prover_zk: bool = True,
                   proof_bytes_out: list[bytes] | None = None,
                   proof_bytes_in: bytes | None = None,
                   static_template: MemverStaticTemplate | None = None,
                   cached_params_bytes: bytes | None = None,
                   params_bytes_out: list[bytes] | None = None,
                   cached_params_obj: LiveLabrador38Params | None = None,
                   params_obj_out: list[LiveLabrador38Params] | None = None,
                   reuse_comkey: bool = False,
                   compact_params: bool = False) -> bool:
    total_started = time.perf_counter()
    step_started = time.perf_counter()
    from lazer import ffi, lib, polyring_t, poly_t, polyvec_t
    from labrados import pack_verify, proof_statement

    if not reuse_comkey:
        lib.labrador38_py_comkey_reset()
    print_timing("导入 LaZer/Labrador 与%s comkey" %
                 ("复用" if reuse_comkey else "重置"), step_started)

    step_started = time.perf_counter()
    params = data["params"]
    n = params["n"]
    q = params["q"]
    handle_bits = params["handle_bits"]
    gadget_base, gadget_dim, _ = accumulator_gadget(params)
    right_norm_bound = right_group_norm_bound(params)
    grouped_ranges = right_group_ranges(params)
    if right_layout not in {"scalar", "grouped"}:
        raise ValueError("right_layout must be 'scalar' or 'grouped'")
    grouped_right = right_layout == "grouped"
    trunc_k_l = params["trunc_k_l"]
    if static_template is not None:
        if not static_template.compatible(data):
            raise ValueError("MemVer static template 与当前公开参数不兼容")
        ring = static_template.ring
        one = static_template.one
    else:
        ring = polyring_t(n, PROOF_MOD_38)
        one = poly_t(ring, const_poly(n, 1))

    have_witness = lifts is not None and "witness" in data
    witness_vectors: list[list[list[int]] | _NativeWitness | None] = []
    witness_dimensions: list[int] = []
    norms: list[int] = []
    w = data.get("witness", {})

    def append_vector(polys: list[list[int]] | _NativeWitness | None,
                      dimension: int,
                      norm: int) -> int:
        if (polys is not None and not isinstance(polys, _NativeWitness)
                and len(polys) != dimension):
            raise ValueError("witness vector has the wrong dimension")
        witness_vectors.append(polys)
        witness_dimensions.append(dimension)
        norms.append(norm)
        return len(witness_vectors) - 1

    def append_slot(poly: list[int] | None, norm: int) -> int:
        return append_vector([poly] if poly is not None else None, 1, norm)

    x_bits = data.get("x_bits", [0] * handle_bits)
    if len(x_bits) != handle_bits or any(bit not in (0, 1) for bit in x_bits):
        raise ValueError("撤销句柄的二进制分解不是规范的 0/1 向量")

    # These are witness values, even though their expected values are public.
    # Keep them as coefficient lists here.  ``poly_t`` objects are iterable in
    # a way that is not suitable for the Python-side 0/1 validation below;
    # passing them here caused the new worker and the JSON fallback to reject
    # every valid request as a non-binary witness.
    x_idx = [
        append_slot(const_poly(n, b) if have_witness else None, 0)
        for b in x_bits
    ]
    x_vector_idx: list[int] = []
    if grouped_right:
        for group_range in grouped_ranges:
            x_vector_idx.append(append_vector(
                [const_poly(n, x_bits[h]) for h in group_range]
                if have_witness else None,
                len(group_range),
                0,
            ))
    tag_idx = append_slot(w["tag"] if have_witness else None, TAG_NORM_SQ)

    def append_decomposed(poly: list[int]) -> list[tuple[int, int]]:
        terms: list[tuple[int, int]] = []
        scale = 1
        for chunk in signed_decompose_poly(poly):
            idx = append_slot(chunk if have_witness else None, SMALL_CHUNK_NORM_SQ)
            terms.append((scale, idx))
            scale *= LINEAR_DECOMP_BASE
        return terms

    def append_native_decomposed() -> list[tuple[int, int]]:
        terms: list[tuple[int, int]] = []
        scale = 1
        for _ in range(LINEAR_DECOMP_CHUNKS):
            terms.append((scale, append_vector(
                _NativeWitness(), 1, SMALL_CHUNK_NORM_SQ)))
            scale *= LINEAR_DECOMP_BASE
        return terms

    zero_poly = zero(n)
    left_idx = [append_decomposed(w["left"][i] if have_witness else zero_poly) for i in range(params["d"])]
    v12_idx = [append_decomposed(w["v12"][i] if have_witness else zero_poly) for i in range(params["d"])]
    v2_idx: list[list[int]] = []
    for i in range(trunc_k_l):
        v2_idx.append([
            append_slot(w["v2"][i][row] if have_witness else None, V2_NORM_SQ)
            for row in range(params["d"])
        ])
    right_group_idx: list[list[list[int]]] = []
    right_scalar_idx: list[list[list[int]]] = []
    if grouped_right:
        for j in range(gadget_dim):
            digit_groups = []
            for group_range in grouped_ranges:
                digit_groups.append([
                    append_vector(
                        [w["right"][h][j][row] for h in group_range]
                        if have_witness else None,
                        len(group_range),
                        right_norm_bound,
                    )
                    for row in range(params["d"])
                ])
            right_group_idx.append(digit_groups)
    else:
        scalar_bound = 256 if int(params.get("deletion_count", 0)) == 0 else right_norm_bound
        for h in range(handle_bits):
            h_idx = []
            for j in range(gadget_dim):
                h_idx.append([
                    append_slot(
                        w["right"][h][j][row] if have_witness else None,
                        scalar_bound,
                    )
                    for row in range(params["d"])
                ])
            right_scalar_idx.append(h_idx)

    refresh_bits = int(params.get("refresh_bits", 0))
    refresh_r_idx: int | None = None
    refresh_z_idx: int | None = None
    refresh_ahat_idx: list[list[int]] = []
    refresh_ahat_lift_idx: list[int] = []
    refresh_ahat_values: list[list[list[list[int]]]] = []
    refresh_ahat_lifts: list[list[list[int]]] = []
    native_refresh = False
    if refresh_bits:
        r_bits = data.get("refresh_r_bits", [0] * refresh_bits)
        refresh_z = w.get("refresh_z", [])
        if (len(r_bits) != refresh_bits or
                any(bit not in (0, 1) for bit in r_bits)):
            raise ValueError("随机化刷新 r 不是规范的二进制向量")
        if have_witness and (len(refresh_z) != gadget_dim or
                             any(len(block) != params["d"]
                                 for block in refresh_z)):
            raise ValueError("随机化刷新 z 的维度错误")
        # One binary witness vector containing refresh_bits constant
        # polynomials.  A single native deg-0 constraint covers the complete
        # vector.  The old layout used one witness vector and one constraint
        # per bit, which inflated Labrador's vector metadata and Python/CFFI
        # call count without reducing the number of proved coefficients.
        refresh_r_idx = append_vector(
            [const_poly(n, bit) for bit in r_bits] if have_witness else None,
            refresh_bits,
            0,
        )
        z_flat = ([refresh_z[j][col] for j in range(gadget_dim)
                   for col in range(params["d"])] if have_witness else None)
        refresh_z_idx = append_vector(
            z_flat, gadget_dim * params["d"],
            gadget_dim * params["d"] * n)

        native_refresh = (have_witness and
                          os.environ.get(
                              "REVAC_NATIVE_REFRESH_WITNESS", "1") != "0" and
                          hasattr(data.get("refresh_A"), "coeffs") and
                          hasattr(lib, "labrador38_py_refresh_set_witness"))
        if have_witness and not native_refresh:
            refresh_a = data["refresh_A"]
            q_half = q // 2
            for row in range(params["d"]):
                row_chunks: list[list[list[int]]] = [
                    [] for _ in range(LINEAR_DECOMP_CHUNKS)
                ]
                row_lifts: list[list[int]] = []
                for gadget in range(gadget_dim):
                    for col in range(params["d"]):
                        poly = refresh_matrix_poly(
                            refresh_a, 0, gadget, row, col)
                        for bit, selected in enumerate(r_bits):
                            if selected:
                                poly = add(
                                    poly,
                                    refresh_matrix_poly(
                                        refresh_a, bit + 1, gadget, row, col))
                        centered = [
                            ((value + q_half) % q) - q_half for value in poly
                        ]
                        row_lifts.append([
                            (raw - reduced) // q
                            for raw, reduced in zip(poly, centered)
                        ])
                        for chunk, chunk_poly in enumerate(
                                signed_decompose_poly(centered)):
                            row_chunks[chunk].append(chunk_poly)
                refresh_ahat_values.append(row_chunks)
                refresh_ahat_lifts.append(row_lifts)
        for row in range(params["d"]):
            refresh_ahat_idx.append([
                append_vector(
                    (_NativeWitness() if native_refresh else
                     refresh_ahat_values[row][chunk]
                     if have_witness else None),
                    gadget_dim * params["d"],
                    gadget_dim * params["d"] * SMALL_CHUNK_NORM_SQ)
                for chunk in range(LINEAR_DECOMP_CHUNKS)
            ])
            refresh_ahat_lift_idx.append(append_vector(
                (_NativeWitness() if native_refresh else
                 refresh_ahat_lifts[row] if have_witness else None),
                gadget_dim * params["d"],
                gadget_dim * params["d"] * SMALL_CHUNK_NORM_SQ))
    native_exact_enabled = os.environ.get(
        "REVAC_NATIVE_EXACT_LIFT", "1") != "0"
    native_acc_lift = (have_witness and native_exact_enabled and
                       os.environ.get("REVAC_NATIVE_ACC_LIFT", "1") != "0" and
                       hasattr(lifts, "coeffs") and
                       hasattr(lib, "labrador38_py_decompose_set_witness"))
    lift_idx = [
        append_native_decomposed() if native_acc_lift else
        append_decomposed(lifts[i] if have_witness else zero_poly)
        for i in range(params["d"])
    ]

    link_r_idx: list[list[int]] = [[], []]
    link_r = w.get("link_r", [[], []])
    if have_witness and (len(link_r) != 2 or
                         any(len(part) != params["d"] for part in link_r)):
        raise ValueError("same-handle commitment opening has the wrong dimension")
    for part in range(2):
        link_r_idx[part] = [
            append_slot(link_r[part][row] if have_witness else None,
                        LINK_OPENING_POLY_NORM_SQ)
            for row in range(params["d"])
        ]
    link_lifts = data.get("link_lift") if have_witness else None
    if have_witness and link_lifts is None:
        link_lifts = compute_link_lift(data)
    if have_witness and (link_lifts is None or
                         len(link_lifts) != params["d"]):
        raise ValueError("same-handle exact lift is missing or malformed")
    native_link_lift = (have_witness and native_exact_enabled and
                        os.environ.get("REVAC_NATIVE_LINK_LIFT", "1") != "0" and
                        hasattr(link_lifts, "coeffs") and
                        hasattr(lib, "labrador38_py_decompose_set_witness"))
    link_lift_idx = [
        append_native_decomposed() if native_link_lift else
        append_decomposed(link_lifts[row] if have_witness else zero_poly)
        for row in range(params["d"])
    ]

    refresh_link_constraints = (params["d"] * gadget_dim * params["d"]
                                if refresh_bits else 0)
    num_constraints = (2 * params["d"] +
                       (handle_bits if grouped_right else 0) +
                       refresh_link_constraints)
    # The permutation/topology depends on the norm *type*, not the numerical
    # bound.  Reuse one native statement across deletion counts and refresh its
    # bounds before parameter selection; parameter objects remain separately
    # keyed by deletion count.
    statement_key = (
        tuple(witness_dimensions), tuple(int(bound == 0) for bound in norms), num_constraints,
        handle_bits + refresh_bits, bool(prover_zk), right_layout,
    )
    live_entry = (static_template.live_statements.get(statement_key)
                  if static_template is not None and grouped_right else None)
    new_statement = live_entry is None
    witness_cache_enabled = (
        os.environ.get("REVAC_LABRADOR_WITNESS_CACHE", "1") != "0" and
        static_template is not None and grouped_right and live_entry is not None)
    previous_keys = (static_template.live_witness_keys.get(statement_key)
                     if witness_cache_enabled else None)
    previous_norms = (static_template.live_witness_norms.get(statement_key)
                      if witness_cache_enabled else None)
    prepared_keys: list[Any] = [None] * len(witness_vectors)
    prepared_norms: list[int | None] = [None] * len(witness_vectors)
    unchanged_witness = [False] * len(witness_vectors)
    if have_witness:
        for idx, (polys, bound) in enumerate(zip(witness_vectors, norms)):
            assert polys is not None
            if isinstance(polys, _NativeWitness):
                continue
            exact_key = tuple(tuple(poly) for poly in polys)
            prepared_keys[idx] = exact_key
            unchanged = (
                previous_keys is not None and previous_norms is not None and
                idx < len(previous_keys) and idx < len(previous_norms) and
                previous_keys[idx] == exact_key)
            unchanged_witness[idx] = unchanged
            if unchanged:
                actual = previous_norms[idx]
                assert actual is not None
            elif bound == 0:
                if any(c not in (0, 1) for poly in polys for c in poly):
                    raise ValueError(
                        f"binary witness vector {idx} is not binary")
                actual = 0
            else:
                actual = sum(l2sq(poly) for poly in polys)
            prepared_norms[idx] = actual
            if bound != 0 and actual > bound:
                raise ValueError(
                    f"witness vector {idx} has norm² {actual}, "
                    f"exceeding the fixed proof bound {bound}"
                )

    print_timing("见证分解、向量分组与范数检查", step_started)
    print(f"Labrador right 布局：{right_layout}")
    print(f"Labrador 见证向量数：{len(witness_vectors)}")
    print(f"Labrador 见证多项式数：{sum(witness_dimensions)}")
    print(f"Labrador 约束数：{num_constraints}")
    step_started = time.perf_counter()
    if live_entry is None:
        ps = proof_statement([n] * len(witness_vectors), witness_dimensions, norms,
                             num_constraints=num_constraints,
                             primesize="38",
                             num_deg0_constraints=handle_bits + refresh_bits,
                             zk=prover_zk)
    else:
        ps = live_entry[0]
        ps.cur_witness_num = 0
        ps.witness_polys = 0
        native_norms = [0] * len(norms)
        for original_idx, (dimension, bound) in enumerate(
                zip(witness_dimensions, norms)):
            mapped_idx = ps.map[original_idx]
            native_norms[mapped_idx] = dimension * n if bound == 0 else bound
        output = lib.labrador38_py_statement_set_norms(
            ps.statement, len(native_norms), ffi.new("uint64_t[]", native_norms))
        if output != 0:
            raise ValueError(f"常驻 MemVer statement 范数更新失败：{output}")
    print(f"Labrador 零知识模式：{'启用' if prover_zk else '关闭（仅限降级测试）'}")
    print_timing("proof_statement 参数初始化/缓存命中", step_started)

    step_started = time.perf_counter()
    if have_witness:
        # Feed decoded integer coefficients directly into Labrador's native
        # witness storage.  The previous path constructed poly_t/polyvec_t
        # wrappers and immediately converted them back to int64 arrays.
        witness_copy_started = time.perf_counter()
        reused_vectors = 0
        for idx, polys in enumerate(witness_vectors):
            assert polys is not None
            if isinstance(polys, _NativeWitness):
                ps.cur_witness_num += 1
                continue
            if unchanged_witness[idx]:
                # append_witness_pointer only copies coefficients into the
                # already-allocated native slot.  Exact equality means that
                # slot is still current; preserve its accounting without the
                # redundant Python-to-C conversion and native copy.
                ps.witness_polys += len(polys) * n // 256
                ps.cur_witness_num += 1
                reused_vectors += 1
                continue
            coeffs = ffi.new(
                "int64_t[]",
                [coefficient for poly in polys for coefficient in poly],
            )
            ps.append_witness_pointer(coeffs, n, len(polys))
        if static_template is not None and grouped_right:
            static_template.live_witness_keys[statement_key] = prepared_keys
            static_template.live_witness_norms[statement_key] = prepared_norms
        if os.environ.get("REVAC_LABRADOR_PROFILE", "0") != "0":
            print(f"MemVer 普通 witness 缓存：复用 {reused_vectors}/"
                  f"{len(witness_vectors)} 个向量")
        print_timing("普通 witness vector 写入", witness_copy_started)
        native_refresh_started = time.perf_counter()
        if native_refresh:
            bank = data["refresh_A"]
            chunk_mapped = ffi.new(
                "size_t[]",
                [ps.map[idx] for row_indices in refresh_ahat_idx
                 for idx in row_indices])
            lift_mapped = ffi.new(
                "size_t[]", [ps.map[idx] for idx in refresh_ahat_lift_idx])
            r_native = ffi.new("uint8_t[]", data["refresh_r_bits"])
            output = lib.labrador38_py_refresh_set_witness(
                ps.witness, chunk_mapped, lift_mapped, bank.coeffs, r_native,
                refresh_bits, gadget_dim, params["d"], params["d"], n,
                q, LINEAR_DECOMP_BASE, LINEAR_DECOMP_CHUNKS)
            if output != 0:
                raise ValueError(
                    f"C 端 A-hat/exact-lift witness 构造失败：{output}")
        print_timing("随机刷新 witness 原生构造", native_refresh_started)
        exact_lift_started = time.perf_counter()
        for source, index_terms, label, enabled in (
                (lifts, lift_idx, "累加器 exact-lift", native_acc_lift),
                (link_lifts, link_lift_idx, "同句柄 exact-lift",
                 native_link_lift)):
            if not enabled:
                continue
            mapped = ffi.new(
                "size_t[]",
                [ps.map[idx] for poly_terms in index_terms
                 for _scale, idx in poly_terms])
            output = lib.labrador38_py_decompose_set_witness(
                ps.witness, mapped, source.coeffs, params["d"], n,
                LINEAR_DECOMP_BASE, LINEAR_DECOMP_CHUNKS)
            if output != 0:
                raise ValueError(f"C 端{label}分解失败：{output}")
        print_timing("exact-lift witness 原生分解", exact_lift_started)
    constraint_bind_started = time.perf_counter()
    if new_statement:
        for idx in x_idx:
            ps.append_deg0_statement(idx)
        if refresh_r_idx is not None:
            mapped_idx = ps.map[refresh_r_idx]
            output = lib.labrador38_py_append_deg0_vector_constraints(
                ps.statement, mapped_idx, refresh_bits)
            if output != 0:
                raise ValueError("无法添加紧凑 refresh-r deg0 约束")
            ps.cur_deg0_constraint_num += refresh_bits

    # Compile the immutable MemVer topology once.  The plan contains the
    # generic linear/quadratic/binary/degree wiring, but not request-specific
    # public RHS values or private witness values.
    relation_plan = (static_template.compiled_relation
                     if static_template is not None and grouped_right else None)
    layout_signature = tuple(witness_dimensions)
    if relation_plan is not None and relation_plan["layout"] != layout_signature:
        raise ValueError("缓存的 MemVer 约束骨架与当前见证布局不一致")
    if relation_plan is None:
        native_sparse_refresh = (
            refresh_bits > 0 and static_template is not None and
            os.environ.get("REVAC_NATIVE_SPARSE_REFRESH", "1") != "0" and
            hasattr(data.get("refresh_A"), "coeffs") and
            hasattr(lib, "labrador38_py_append_refresh_constraints"))
        bit_links: list[tuple[list[Any], list[int]]] = []
        if grouped_right:
            for group, group_range in enumerate(grouped_ranges):
                for local_h, h in enumerate(group_range):
                    selector = (static_template.selectors[group][local_h]
                                if static_template is not None else
                                polyvec_t(
                                    ring, len(group_range),
                                    [poly_t(ring, const_poly(
                                        n, int(index == local_h)))
                                     for index in range(len(group_range))]))
                    minus_one = (static_template.minus_one
                                 if static_template is not None else -one)
                    bit_links.append(([selector, minus_one],
                                      [x_vector_idx[group], x_idx[h]]))

        link_rows: list[tuple[list[Any], list[int]]] = []
        for row in range(params["d"]):
            stat_link = [static_template.one if static_template is not None else one]
            wit_link = [link_r_idx[0][row]]
            for col in range(params["d"]):
                stat_link.append(
                    static_template.link_a[row][col]
                    if static_template is not None else
                    poly_t(ring, data["link_A"][row][col]))
                wit_link.append(link_r_idx[1][col])
            for bit in range(handle_bits):
                stat_link.append(
                    static_template.link_d_bits[row][bit]
                    if static_template is not None else
                    poly_t(ring, negacyclic_monomial_mul(
                        data["link_D"][row][0], bit)))
                wit_link.append(x_idx[bit])
            for scale, idx in link_lift_idx[row]:
                stat_link.append(
                    static_template.q_scale_poly[scale]
                    if static_template is not None else
                    poly_t(ring, const_poly(n, q * scale)))
                wit_link.append(idx)
            link_rows.append((stat_link, wit_link))

        refresh_links: list[tuple[list[Any], list[int], Any]] = []
        if refresh_bits and not native_sparse_refresh:
            assert refresh_r_idx is not None
            for row in range(params["d"]):
                for entry in range(gadget_dim * params["d"]):
                    if static_template is not None:
                        coeff_vector = static_template.refresh_coeff[row][entry]
                        selector = static_template.refresh_selectors[row][entry]
                        rhs = -poly_t(
                            ring,
                            refresh_matrix_poly(
                                data["refresh_A"], 0,
                                entry // params["d"], row,
                                entry % params["d"]))
                    else:
                        coeff = [
                            poly_t(ring, refresh_matrix_poly(
                                data["refresh_A"], i + 1,
                                entry // params["d"], row,
                                entry % params["d"]))
                            for i in range(refresh_bits)
                        ]
                        coeff_vector = polyvec_t(ring, refresh_bits, coeff)
                        selector_values = [poly_t(ring, zero_poly)] * (
                            gadget_dim * params["d"])
                        selector_values[entry] = one
                        selector = polyvec_t(
                            ring, len(selector_values), selector_values)
                        rhs = -poly_t(
                            ring, refresh_matrix_poly(
                                data["refresh_A"], 0,
                                entry // params["d"], row,
                                entry % params["d"]))
                    ahat_stat = []
                    ahat_wit = []
                    scale = 1
                    for chunk in range(LINEAR_DECOMP_CHUNKS):
                        ahat_stat.append(-(scale * selector))
                        ahat_wit.append(refresh_ahat_idx[row][chunk])
                        scale *= LINEAR_DECOMP_BASE
                    ahat_stat.append(-(q * selector))
                    ahat_wit.append(refresh_ahat_lift_idx[row])
                    # The public refresh coefficients already form a vector
                    # paired with the compact r witness vector.
                    refresh_links.append(
                        ([coeff_vector] + ahat_stat,
                         [refresh_r_idx] + ahat_wit, rhs))

        accumulator_rows: list[tuple[list[Any], list[int], list[int],
                                    list[int], list[int]]] = []
        for row in range(params["d"]):
            stat_lin: list[Any] = []
            wit_lin: list[int] = []
            for scale, idx in left_idx[row]:
                stat_lin.append(static_template.scale_poly[scale]
                                if static_template is not None else
                                poly_t(ring, const_poly(n, scale)))
                wit_lin.append(idx)
            for col in range(params["d"]):
                for chunk, (scale, idx) in enumerate(v12_idx[col]):
                    stat_lin.append(static_template.atrap_scaled[chunk][row][col]
                                    if static_template is not None else
                                    poly_t(ring, scalar_mul(
                                        scale, data["Atrap"][row][col])))
                    wit_lin.append(idx)
            for i in range(trunc_k_l):
                for col in range(params["d"]):
                    stat_lin.append(static_template.trunc_neg[i][row][col]
                                    if static_template is not None else
                                    -poly_t(ring, data["trunc_B"][i][row][col]))
                    wit_lin.append(v2_idx[i][col])
            if grouped_right:
                for j in range(gadget_dim):
                    for group, group_range in enumerate(grouped_ranges):
                        for col in range(params["d"]):
                            stat_lin.append(
                                static_template.ba_grouped[j][group][row][col]
                                if static_template is not None else
                                polyvec_t(
                                    ring, len(group_range),
                                    [poly_t(ring, data["Ba"][h][j][row][col])
                                     for h in group_range]))
                            wit_lin.append(right_group_idx[j][group][col])
            else:
                for h in range(handle_bits):
                    for j in range(gadget_dim):
                        for col in range(params["d"]):
                            stat_lin.append(poly_t(
                                ring, data["Ba"][h][j][row][col]))
                            wit_lin.append(right_scalar_idx[h][j][col])
            for scale, idx in lift_idx[row]:
                stat_lin.append(static_template.q_scale_poly[scale]
                                if static_template is not None else
                                poly_t(ring, const_poly(n, q * scale)))
                wit_lin.append(idx)

            stat_quad: list[int] = []
            wit_quad1: list[int] = []
            wit_quad2: list[int] = []
            bexpi = 196
            for i in range(trunc_k_l):
                stat_quad.append(bexpi)
                wit_quad1.append(tag_idx)
                wit_quad2.append(v2_idx[i][row])
                bexpi *= 14
            if grouped_right:
                for j in range(gadget_dim):
                    for group in range(len(grouped_ranges)):
                        stat_quad.append(-(gadget_base**j))
                        wit_quad1.append(x_vector_idx[group])
                        wit_quad2.append(right_group_idx[j][group][row])
            else:
                for h in range(handle_bits):
                    for j in range(gadget_dim):
                        stat_quad.append(-(gadget_base**j))
                        wit_quad1.append(x_idx[h])
                        wit_quad2.append(right_scalar_idx[h][j][row])
            if refresh_bits:
                assert refresh_z_idx is not None
                scale = 1
                for chunk in range(LINEAR_DECOMP_CHUNKS):
                    stat_quad.append(scale)
                    wit_quad1.append(refresh_ahat_idx[row][chunk])
                    wit_quad2.append(refresh_z_idx)
                    scale *= LINEAR_DECOMP_BASE
            accumulator_rows.append(
                (stat_lin, wit_lin, stat_quad, wit_quad1, wit_quad2))

        relation_plan = {
            "layout": layout_signature,
            "native_sparse_refresh": native_sparse_refresh,
            "bit_links": bit_links,
            "link_rows": link_rows,
            "refresh_links": refresh_links,
            "accumulator_rows": accumulator_rows,
        }
        if static_template is not None and grouped_right:
            static_template.compiled_relation = relation_plan
        print("MemVer 通用约束骨架：首次编译")
    else:
        print("MemVer 通用约束骨架：命中缓存")

    zero_ring_poly = (static_template.zero if static_template is not None else
                      poly_t(ring, zero_poly))
    raw_relation = (static_template.compiled_raw_relation
                    if static_template is not None and grouped_right else None)
    map_signature = tuple(ps.map)
    if raw_relation is not None and raw_relation["map"] != map_signature:
        raise ValueError("缓存的原生 MemVer 约束映射与当前 statement 不一致")
    if raw_relation is None and static_template is not None and grouped_right:
        raw_relation = {
            "map": map_signature,
            "bit_links": [
                _compile_raw_linear(ps, stat_list, wit_list, zero_ring_poly)
                for stat_list, wit_list in relation_plan["bit_links"]
            ],
            "link_rows": [
                _compile_raw_linear(ps, stat_list, wit_list, zero_ring_poly)
                for stat_list, wit_list in relation_plan["link_rows"]
            ],
            "refresh_links": [
                _compile_raw_linear(ps, stat_list, wit_list, rhs)
                for stat_list, wit_list, rhs in
                relation_plan.get("refresh_links", [])
            ],
            "accumulator_rows": [
                _compile_raw_quadratic(ps, relation, zero_ring_poly)
                for relation in relation_plan["accumulator_rows"]
            ],
        }
        static_template.compiled_raw_relation = raw_relation
        print("MemVer CFFI 约束数组：首次扁平化并缓存")

    if raw_relation is not None and not new_statement:
        assert live_entry is not None
        static_hash, dynamic_first = live_entry[1], live_entry[2]
        dynamic_key = _dynamic_rhs_key(data)
        bind_hit = (_dynamic_bind_cache_enabled() and
                    static_template.live_dynamic_rhs.get(statement_key) ==
                    dynamic_key)
        if not bind_hit:
            lib.labrador38_py_statement_hash_reset(
                ps.statement, ffi.new("uint8_t[16]", static_hash))
            constraint_index = dynamic_first
            for row, raw in enumerate(raw_relation["link_rows"]):
                rhs = ffi.new("int64_t[]", data["link_commitment"][row])
                output = lib.labrador38_py_statement_rebind_constraint_raw(
                    ps.statement, constraint_index, raw["nvec"],
                    raw["lengths"], raw["degree"], raw["phi"], rhs,
                    0, ffi.NULL)
                if output != 0:
                    raise ValueError(
                        f"常驻 MemVer link statement 动态绑定失败：{output}")
                constraint_index += 1
            for row, raw in enumerate(raw_relation["accumulator_rows"]):
                rhs = ffi.new("int64_t[]", data["acc"][row])
                output = lib.labrador38_py_statement_rebind_constraint_raw(
                    ps.statement, constraint_index, raw["nlin"],
                    raw["lengths"], raw["degree"], raw["phi"], rhs,
                    raw["nprod"], raw["quad"])
                if output != 0:
                    raise ValueError(
                        f"常驻 MemVer accumulator statement 动态绑定失败：{output}")
                constraint_index += 1
            static_template.live_dynamic_rhs[statement_key] = dynamic_key
        if os.environ.get("REVAC_LABRADOR_PROFILE", "0") != "0":
            print("MemVer 动态 RHS 绑定缓存：" +
                  ("命中" if bind_hit else "未命中"))
        print("MemVer proof_statement 公共骨架：命中常驻缓存")
    elif raw_relation is not None:
        # Put every immutable constraint first.  The eight request-specific
        # commitment/accumulator RHS constraints form a mutable tail.
        for raw in raw_relation["bit_links"]:
            _append_raw_linear(ps, raw)
        if relation_plan.get("native_sparse_refresh", False):
            assert refresh_r_idx is not None
            bank = data["refresh_A"]
            chunk_mapped = ffi.new(
                "size_t[]",
                [ps.map[idx] for row_indices in refresh_ahat_idx
                 for idx in row_indices])
            lift_mapped = ffi.new(
                "size_t[]", [ps.map[idx] for idx in refresh_ahat_lift_idx])
            output = lib.labrador38_py_append_refresh_constraints(
                ps.statement, ps.map[refresh_r_idx], chunk_mapped,
                lift_mapped, bank.coeffs, refresh_bits, gadget_dim,
                params["d"], params["d"], n, q,
                LINEAR_DECOMP_BASE, LINEAR_DECOMP_CHUNKS)
            if output != 0:
                raise ValueError(
                    f"原生稀疏 A-hat 绑定约束构造失败：{output}")
            ps.cur_constraint_num += params["d"] * gadget_dim * params["d"]
        else:
            for raw in raw_relation["refresh_links"]:
                _append_raw_linear(ps, raw)
        dynamic_first = ps.cur_constraint_num
        static_hash_buf = ffi.new("uint8_t[16]")
        lib.labrador38_py_statement_hash(static_hash_buf, ps.statement)
        static_hash = bytes(ffi.buffer(static_hash_buf, 16))
        for row, raw in enumerate(raw_relation["link_rows"]):
            _append_raw_linear(ps, raw, data["link_commitment"][row])
        for row, raw in enumerate(raw_relation["accumulator_rows"]):
            _append_raw_quadratic(ps, raw, data["acc"][row])
        if static_template is not None:
            static_template.live_statements[statement_key] = (
                ps, static_hash, dynamic_first)
            static_template.live_dynamic_rhs[statement_key] = \
                _dynamic_rhs_key(data)
            static_template.seal_native_statement(layout_signature)
        print("MemVer proof_statement 公共骨架：首次构造并转为常驻缓存")
    else:
        for stat_list, wit_list in relation_plan["bit_links"]:
            ps.append_statement(stat_list, wit_list, zero_ring_poly)
        for row, (stat_list, wit_list) in enumerate(relation_plan["link_rows"]):
            ps.append_statement(stat_list, wit_list,
                                poly_t(ring, data["link_commitment"][row]))
        for stat_list, wit_list, rhs in relation_plan.get("refresh_links", []):
            ps.append_statement(stat_list, wit_list, rhs)
        for row, relation in enumerate(relation_plan["accumulator_rows"]):
            ps.append_quadratic_statement(
                *relation, poly_t(ring, data["acc"][row]))
    print_timing("动态约束构造/重绑定", constraint_bind_started)
    print_timing("填充见证与构造全部约束", step_started)

    if (have_witness and
            os.environ.get("REVAC_LABRADOR_SIMPLE_CHECK", "0") != "0"):
        simple_ok = int(ps.simple_verify(ps.statement, ps.witness)) == 1
        print(f"Labrador 原生关系预检查：{'通过' if simple_ok else '失败'}")
        if not simple_ok:
            raise ValueError("Labrador witness does not satisfy the compiled MemVer relation")

    if debug_hash:
        from lazer import ffi, lib
        st = ps.output_statement()
        h = ffi.new("uint8_t[16]")
        lib.labrador38_py_statement_hash(h, st)
        print("statement hash:", bytes(ffi.buffer(h, 16)).hex())
        if proof_in is None and proof_out is None:
            print_timing("Labrador statement 构造总计", total_started)
            return True

    if proof_in is not None or proof_bytes_in is not None:
        if proof_in is not None and proof_bytes_in is not None:
            raise ValueError("provide either proof_in or proof_bytes_in, not both")
        step_started = time.perf_counter()
        proof_bytes = proof_bytes_in if proof_bytes_in is not None else proof_in.read_bytes()
        expected_profile_id = memver_labrador_profile_id(data)
        generated_verify_params: LiveLabrador38Params | None = None
        effective_cached_params = cached_params_obj
        retain_generated = False
        try:
            if labrador38_proof_is_compact(proof_bytes):
                if labrador38_proof_profile_id(proof_bytes) != expected_profile_id:
                    raise ValueError("compact proof profile ID 与当前 MemVer statement 不一致")
                if effective_cached_params is None:
                    generated_verify_params = LiveLabrador38Params()
                    error = ps.gen_params(
                        generated_verify_params.params, ps.statement,
                        1 if ps.zk else 0, 1 if ps.debug else 0)
                    if error != 0:
                        raise ValueError(
                            f"GCS 本地生成 compact proof 参数失败：{error}")
                    effective_cached_params = generated_verify_params
                    print("GCS compact 参数：本地生成并完成 profile 绑定")
            params_from_proof, proof, refs = deserialize_labrador38_proof(
                proof_bytes, reset_comkey=not reuse_comkey,
                cached_params_obj=effective_cached_params,
                expected_profile_id=expected_profile_id)
            print_timing("读取并反序列化 Labrador proof", step_started)
            step_started = time.perf_counter()
            try:
                ok = pack_verify((params_from_proof, proof),
                                 ps.output_statement(), "38")
            finally:
                lib.labrador38_py_free_proof(proof)
                if effective_cached_params is None:
                    lib.labrador38_py_free_params(params_from_proof)
            print_timing("Pack Verify 核心计算", step_started)
            _ = refs
            if generated_verify_params is not None and params_obj_out is not None:
                params_obj_out.append(generated_verify_params)
                retain_generated = True
            print_timing("Labrador 验证脚本内部总计", total_started)
            return ok
        finally:
            if generated_verify_params is not None and not retain_generated:
                generated_verify_params.close()

    if not have_witness:
        raise ValueError("proving requires witness data")

    if cached_params_bytes is not None and cached_params_obj is not None:
        raise ValueError("不能同时提供序列化参数缓存和常驻原生参数缓存")
    cached_params = (cached_params_bytes is not None or
                     cached_params_obj is not None)
    generated_live_params: LiveLabrador38Params | None = None
    direct_params = None
    direct_proof = None
    direct_refs: list[Any] = []
    if cached_params_obj is not None:
        step_started = time.perf_counter()
        params_for_proof = cached_params_obj.checked_params(
            require_zk=prover_zk)
        direct_proof = ffi.new("labrador38_dch_pack_proof")
        proof_for_serialization = direct_proof
        print_timing("Pack 常驻参数缓存检查", step_started)
        print("Labrador Pack 参数：命中常驻原生 zk=True 缓存")
    elif cached_params_bytes is not None:
        step_started = time.perf_counter()
        direct_params, direct_refs = deserialize_labrador38_params(
            cached_params_bytes, require_zk=prover_zk,
            reset_comkey=not reuse_comkey)
        direct_proof = ffi.new("labrador38_dch_pack_proof")
        print_timing("Pack 参数缓存反序列化", step_started)
        params_for_proof = direct_params
        proof_for_serialization = direct_proof
        print("Labrador Pack 参数：命中经过 zk=True 检查的缓存")
    elif params_obj_out is not None:
        # Generate directly into a separately-owned native object.  The
        # proof_statement destructor then frees only this request's statement
        # and witness, while the worker retains the immutable public params.
        step_started = time.perf_counter()
        generated_live_params = LiveLabrador38Params()
        error = ps.gen_params(generated_live_params.params, ps.statement,
                              1 if ps.zk else 0,
                              1 if ps.debug else 0)
        print_timing("Pack 参数生成", step_started)
        if error != 0:
            generated_live_params.close()
            print(f"Labrador Pack 参数生成失败：{error}")
            return False
        params_for_proof = generated_live_params.params
        direct_proof = ffi.new("labrador38_dch_pack_proof")
        proof_for_serialization = direct_proof
        if params_bytes_out is not None:
            params_bytes_out.append(
                serialize_labrador38_params(params_for_proof,
                                             require_zk=prover_zk))
        print("Labrador Pack 参数：新生成并转为常驻原生缓存")
    else:
        step_started = time.perf_counter()
        error = ps.gen_params(ps.params, ps.statement, 1 if ps.zk else 0,
                              1 if ps.debug else 0)
        print_timing("Pack 参数生成", step_started)
        if error != 0:
            print(f"Labrador Pack 参数生成失败：{error}")
            return False
        params_for_proof = ps.params
        proof_for_serialization = ps.proof
        if params_bytes_out is not None:
            params_bytes_out.append(
                serialize_labrador38_params(params_for_proof,
                                             require_zk=prover_zk))
        print("Labrador Pack 参数：新生成并可缓存")

    if prover_zk and lib.labrador38_py_params_is_zk(params_for_proof) != 1:
        raise ValueError("Prover 请求了零知识模式，但生成参数未包含 ZK/LNP 阶段")
    step_started = time.perf_counter()
    ps.prove(proof_for_serialization, ps.statement, ps.witness,
             params_for_proof)
    print_timing("Dachshund + LaBRADOR + LNP Prove", step_started)
    if not cached_params and generated_live_params is None:
        ps.proven = True
    step_started = time.perf_counter()
    envelope_params_bytes = cached_params_bytes
    if not compact_params:
        if cached_params_obj is not None:
            envelope_params_bytes = cached_params_obj.serialized_bytes(
                require_zk=prover_zk)
        elif generated_live_params is not None:
            envelope_params_bytes = generated_live_params.serialized_bytes(
                require_zk=prover_zk)
    proof_bytes = serialize_labrador38_proof(
        params_for_proof, proof_for_serialization,
        require_zk=prover_zk,
        cached_params_bytes=envelope_params_bytes,
        compact_profile_id=(memver_labrador_profile_id(data)
                            if compact_params else None))
    print_timing("Labrador 参数与 proof 序列化", step_started)
    if direct_proof is not None:
        lib.labrador38_py_free_proof(direct_proof)
    if direct_params is not None:
        lib.labrador38_py_free_params(direct_params)
        _ = direct_refs
    if generated_live_params is not None:
        params_obj_out.append(generated_live_params)
    if proof_bytes_out is not None:
        proof_bytes_out.append(proof_bytes)
    if proof_out is not None:
        step_started = time.perf_counter()
        proof_out.write_bytes(proof_bytes)
        print_timing("proof 文件写入", step_started)
        print(f"序列化后的 Labrador proof 大小：{len(proof_bytes)} 字节")
    if not reuse_comkey:
        lib.labrador38_py_comkey_reset()
    if not self_verify:
        print("Labrador prover 内部自验证：已跳过（交给 verifier 验证）")
        print_timing("Labrador 证明脚本内部总计", total_started)
        return True
    step_started = time.perf_counter()
    roundtrip_params, roundtrip, refs = deserialize_labrador38_proof(
        proof_bytes, reset_comkey=not reuse_comkey)
    try:
        ok = pack_verify((roundtrip_params, roundtrip), ps.output_statement(), "38")
    finally:
        lib.labrador38_py_free_proof(roundtrip)
        lib.labrador38_py_free_params(roundtrip_params)
    print_timing("prover 内部自验证", step_started)
    _ = refs
    print_timing("Labrador 证明脚本内部总计", total_started)
    return ok


def main() -> int:
    script_started = time.perf_counter()
    parser = argparse.ArgumentParser()
    parser.add_argument("json_path")
    parser.add_argument(
        "--mode",
        choices=["check", "labrador-simple", "prove", "verify",
                 "statement-hash", "prove-non-zk-test"],
        default="check",
    )
    parser.add_argument("--proof-out", type=pathlib.Path)
    parser.add_argument("--proof-in", type=pathlib.Path)
    parser.add_argument("--debug-hash", action="store_true")
    parser.add_argument("--right-layout", choices=["scalar", "grouped"],
                        default="grouped")
    args = parser.parse_args()

    step_started = time.perf_counter()
    with open(args.json_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    print_timing("读取并解析 MemVer JSON", step_started)

    if args.mode == "verify":
        if args.proof_in is None:
            raise ValueError("--proof-in is required in verify mode")
        try:
            ok = build_labrador(data, None, proof_in=args.proof_in,
                                debug_hash=args.debug_hash,
                                right_layout=args.right_layout)
        except ValueError as exc:
            print(f"GCS 验证策略拒绝 proof：{exc}", file=sys.stderr,
                  flush=True)
            return 1
        print(f"Labrador proof 字节验证：{'通过' if ok else '失败'}")
        print_timing("Python 验证子进程总计", script_started)
        return 0 if ok else 1
    if args.mode == "statement-hash":
        if "witness" in data:
            lifts, _ = load_or_compute_lift(data)
        else:
            lifts = None
        ok = build_labrador(data, lifts, debug_hash=True,
                            right_layout=args.right_layout)
        print_timing("Python statement-hash 子进程总计", script_started)
        return 0 if ok else 1

    step_started = time.perf_counter()
    lifts, lift_source = load_or_compute_lift(data)
    print_timing("读取/计算整数关系 q-lift", step_started)
    print(f"lift 来源       ：{lift_source}")
    max_lift = max(abs(c) for row in lifts for c in row)
    print("真实 AccGenpp MemVer 导入")
    print("===========================")
    print(f"输入 JSON       ：{args.json_path}")
    if "x" in data:
        print(f"x               ：0x{data['x']:08x}")
    print(f"q               ：{data['params']['q']}")
    print(f"撤销次数        ：{data['params'].get('deletion_count', 0)}")
    print(f"right 分组界²   ：{right_group_norm_bound(data['params'])}")
    print("right 分组布局  ：" + ",".join(
        str(len(group))
        for group in right_group_ranges(data["params"])))
    print(f"lift 多项式数   ：{len(lifts)}")
    print(f"最大 |lift 系数|：{max_lift}")
    step_started = time.perf_counter()
    max_right_norm, max_right_group_norm, total_right_norm = right_witness_norms(data)
    print_timing("right 见证范数统计", step_started)
    print(f"最大 right 范数²：{max_right_norm}")
    print(f"最大 right 分组范数²：{max_right_group_norm}")
    print(f"right 总范数²   ：{total_right_norm}")
    print("整数 q-lifting 检查：通过")

    if args.mode in {"labrador-simple", "prove", "prove-non-zk-test"}:
        if args.mode in {"prove", "prove-non-zk-test"} and args.proof_out is None:
            raise ValueError("--proof-out is required in prove mode")
        prover_zk = args.mode != "prove-non-zk-test"
        ok = build_labrador(data, lifts, proof_out=args.proof_out,
                            debug_hash=args.debug_hash,
                            self_verify=args.mode == "labrador-simple",
                            right_layout=args.right_layout,
                            prover_zk=prover_zk)
        if args.mode == "prove-non-zk-test":
            print("已生成仅供降级回归测试使用的 non-ZK proof；生产 GCS 必须拒绝它")
        print(f"Labrador 导入/证明：{'通过' if ok else '失败'}")
        print_timing("Python 证明子进程总计", script_started)
        return 0 if ok else 1
    print("已跳过 Labrador 导入；如需执行，请在 SDE 下使用 --mode labrador-simple。")
    print_timing("Python 检查子进程总计", script_started)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
