#!/usr/bin/env python3
"""Persistent same-device worker for the q96 MemVer + same-handle proof.

The first request compiles the immutable statement: 96 exact rank-16/rank-8
quadratic-aggregation equalities, 48 accumulator relations, and four
same-handle relations.  A full prove request also installs the private witness;
later requests bound to the identical PP ID, sequence, handle and accumulator
reuse it while replacing the fresh same-handle opening.  Verification rebinds
only the 52 dynamic public RHS values.  The Unix socket is local and mode 0600;
private prover payloads are never an inter-device protocol.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import signal
import socket
import struct
import time
from collections import OrderedDict
from typing import Any

# Configure libgomp before importing the LaZer CFFI module.  Managed
# containers occasionally export the invalid value ``OMP_NUM_THREADS=0``;
# fixing it in main() would be too late because loading the shared library can
# initialize the OpenMP runtime during the imports below.
try:
    _omp_threads_valid = int(os.environ.get("OMP_NUM_THREADS", "")) > 0
except ValueError:
    _omp_threads_valid = False
if not _omp_threads_valid:
    os.environ["OMP_NUM_THREADS"] = "12"
# The q96 Pack prover reaches native reducers with shared scratch state.  On
# the pqac8 Xeon, every tested multi-threaded configuration eventually emitted
# a proof rejected by a later tail-reduction self-check.  Keep statement import
# parallel through OMP_NUM_THREADS, but serialize the transcript-sensitive
# reducers until their scratch ownership is made reentrant.  An explicit
# override remains available for diagnostics only.
os.environ.setdefault("REVAC_LABRADOR_THREADS", "1")
# q88 accumulator coefficients are lifted through three external CRT lanes.
# Use the exact six-prime LaZer lift by default; the approximate conversion can
# be selected explicitly for diagnostics, but must not be the production path.
os.environ.setdefault("REVAC_LABRADOR_EXACT_CRT", "1")
os.environ.setdefault("REVAC_LABRADOR_STATEMENT_THREADS", "8")
os.environ.setdefault("REVAC_LABRADOR_PYTHON_STATEMENT_THREADS", "8")
os.environ.setdefault("REVAC_LABRADOR_Q96_COMPACT_HASH", "1")

# The q88 production profile uses only kernels whose output partitions and
# scratch ownership have been audited.  Keep every setting overridable for
# profiling, but make the stable sub-10-second AutoDL profile the default.
os.environ.setdefault("REVAC_LABRADOR_QUADG_THREADS", "8")
os.environ.setdefault("REVAC_LABRADOR_GARBAGE_PARALLEL", "1")
os.environ.setdefault("REVAC_LABRADOR_POLXVEC_SCALAR_THREADS", "8")
os.environ.setdefault("REVAC_LABRADOR_FROMPOLY_THREADS", "8")
os.environ.setdefault("REVAC_LABRADOR_REFRESH_THREADS", "8")
os.environ.setdefault("REVAC_LABRADOR_EXACT_REFRESH_THREADS", "1")
os.environ.setdefault("REVAC_LABRADOR_LDD_PROJECTION_THREADS", "8")
os.environ.setdefault("REVAC_LABRADOR_LDD_PROJECTION_PARALLEL", "1")
os.environ.setdefault("REVAC_LABRADOR_LDD_LING_THREADS", "12")
os.environ.setdefault("REVAC_LABRADOR_LDD_LING_PARALLEL", "1")
os.environ.setdefault("REVAC_LABRADOR_LDD_COLLAPSE_JL_THREADS", "4")
os.environ.setdefault("REVAC_LABRADOR_LDD_COLLAPSE_JL_PARALLEL", "1")
os.environ.setdefault("REVAC_LABRADOR_LDD_ZQ_THREADS", "4")
os.environ.setdefault("REVAC_LABRADOR_LDD_ZQ_PARALLEL", "1")
os.environ.setdefault("REVAC_LABRADOR_LDD_ADDCHECK_THREADS", "1")
os.environ.setdefault("REVAC_LABRADOR_LDD_JL_SPROD_THREADS", "8")
os.environ.setdefault("REVAC_LABRADOR_LDD_JL_SPROD_PARALLEL", "1")
os.environ.setdefault("REVAC_LABRADOR_LDD_LIFT_THREADS", "4")
os.environ.setdefault("REVAC_LABRADOR_LDD_LIFT_PARALLEL", "1")

# These generic paths reach aliased views or transcript-sensitive reducers in
# the q88 statement, so production keeps them serialized.
os.environ.setdefault("REVAC_LABRADOR_LDD_PARALLEL", "0")
os.environ.setdefault("REVAC_LABRADOR_POLXVEC_PARALLEL", "0")
os.environ.setdefault("REVAC_LABRADOR_ZERO_PARALLEL", "0")
os.environ.setdefault("REVAC_LABRADOR_COMMIT_PARALLEL", "0")

from revac_acc96_memver_labrador import (
    Q96StatementCache, _aggregated_quadratics_enabled, build_statement,
)
from revac_acc96_wire import (
    INSTANCE_BYTES, LINK_BYTES, POLY_BYTES, PROVER_BYTES, InstancePayload,
    LinkPayload, ProverPayload, PublicLinkPayload, PublicParameters,
    ZeroProverPayload,
)
from revac_acc_memver_import_labrador import (
    deserialize_labrador38_proof, serialize_labrador38_proof,
)


HASH_BYTES = 32
MAX_FRAME = 64 * 1024 * 1024
REUSE_MAGIC = b"RV88RU01"
REUSE_PREFIX_BYTES = 8 + 4 + 4 + POLY_BYTES


class ReuseMiss(Exception):
    """The local worker has no witness matching a compact reuse request."""

def _bind_deployment_profile(entry: dict[str, Any],
                             context_hash: bytes) -> bytes:
    """Bind compact Pack proofs to public q96 PP and q19 link matrices.

    Pack parameter serialization contains process-local representation details
    and is therefore not a canonical cross-device identifier.  Both parties
    instead generate their own parameters from the same canonical zero-RHS
    statement and bind the compact envelope to the authenticated public PP,
    the seed that expands the q19 link matrices, and an explicit compiler ABI.
    """
    previous = entry.get("link_context_hash")
    if previous is not None and previous != context_hash:
        raise ValueError("q96 worker link-matrix context changed within PP cache")
    entry["link_context_hash"] = context_hash
    compiler_layout = (b"rank16x8-aggregate-v3-compact-hash" if
                       _aggregated_quadratics_enabled() else
                       b"direct-quadratics-v3-compact-hash")
    profile_id = hashlib.shake_256(
        b"REVAC-Q88-LABRADOR38-DEPLOYMENT-PROFILE-V1\x00" +
        compiler_layout + b"\x00" +
        entry["pp_id"] + context_hash
    ).digest(HASH_BYTES)
    for cache_name in ("prove_statement", "verify_statement"):
        cache: Q96StatementCache = entry[cache_name]
        if cache.profile_id is not None and cache.profile_id != profile_id:
            raise ValueError("q96 resident statement profile binding changed")
        cache.profile_id = profile_id
    return profile_id


def _recv_full(conn: socket.socket, size: int) -> bytes:
    chunks: list[bytes] = []
    remaining = size
    while remaining:
        chunk = conn.recv(remaining)
        if not chunk:
            raise ValueError("truncated q96 worker frame")
        chunks.append(chunk)
        remaining -= len(chunk)
    return b"".join(chunks)


def recv_frame(conn: socket.socket) -> bytes:
    length = struct.unpack("<I", _recv_full(conn, 4))[0]
    if length == 0 or length > MAX_FRAME:
        raise ValueError("invalid q96 worker frame length")
    return _recv_full(conn, length)


def send_frame(conn: socket.socket, payload: bytes) -> None:
    if not payload or len(payload) > MAX_FRAME:
        raise ValueError("invalid q96 worker response length")
    conn.sendall(struct.pack("<I", len(payload)) + payload)


def prove(dynamic: bytes, entry: dict[str, Any],
          reuse: bool = False) -> tuple[bytes, bool]:
    total_started = time.perf_counter()
    stage_started = time.perf_counter()
    if reuse:
        if (len(dynamic) != REUSE_PREFIX_BYTES + LINK_BYTES or
                dynamic[:8] != REUSE_MAGIC):
            raise ValueError("invalid q96 compact-reuse payload")
        prover = entry.get("last_prover")
        if prover is None:
            raise ReuseMiss
        sequence, handle = struct.unpack_from("<II", dynamic, 8)
        accumulator = dynamic[16:16 + POLY_BYTES]
        if (sequence != prover.sequence or handle != prover.handle or
                accumulator != prover.raw[16:16 + POLY_BYTES]):
            raise ReuseMiss
        link = LinkPayload.decode(dynamic[REUSE_PREFIX_BYTES:])
    else:
        if len(dynamic) != PROVER_BYTES + LINK_BYTES:
            raise ValueError("invalid q96 worker prove payload length")
        prover = ProverPayload(dynamic[:PROVER_BYTES])
        link = LinkPayload.decode(dynamic[PROVER_BYTES:])
    _bind_deployment_profile(entry, link.context_hash)
    print(f"[q96 worker timing] "
          f"{'reuse cached prover/link' if reuse else 'decode prover/link'}: "
          f"{time.perf_counter() - stage_started:.6f} s", flush=True)
    # revac_q96_link_wire_encode() accepts this local payload only after the C
    # side has checked the exact integer q19 opening and produced its lift.
    # The four Pack constraints enforce the opening again.  Recomputing the
    # same clear relation in Python therefore adds no verifier security and
    # used about 0.20 s per show.  Keep it as an opt-in compiler diagnostic.
    if os.environ.get("REVAC_LABRADOR_Q96_LINK_AUDIT", "0") != "0":
        stage_started = time.perf_counter()
        link.audit_relation(prover.handle_vector)
        print(f"[q96 worker timing] exact-link audit: "
              f"{time.perf_counter() - stage_started:.6f} s", flush=True)
    # Pack Prove mutates native working storage inside proof_statement.  Keep
    # the prover skeleton strictly separate from the verifier skeleton even
    # when this integration-test worker serves both roles in one process.
    resident_cache: Q96StatementCache = entry["prove_statement"]
    fresh_statement = os.environ.get(
        "REVAC_LABRADOR_Q96_FRESH_PROVER_STATEMENT", "0") != "0"
    if fresh_statement:
        # Diagnostic isolation mode.  Pack Prove is allowed to mutate its
        # request-local statement and witness, while the verifier keeps using
        # the canonical resident profile.  The compact envelope remains bound
        # to the same deployment profile and proof format.
        cache = Q96StatementCache()
        cache.profile_id = resident_cache.profile_id
        cache_hit = False
    else:
        cache = resident_cache
        cache_hit = cache.ps is not None
    # The C bridge has already run full q96 MemVer, constructed an exact
    # coefficient-wise lift, and checked the same-handle opening before this
    # local mode-0600 Unix-socket request.  Pack Prove enforces the compiled
    # relation again.  Labrador's clear evaluator is therefore an optional
    # compiler diagnostic instead of a per-show self-verification pass.
    clear_check = os.environ.get(
        "REVAC_LABRADOR_Q96_SIMPLE_CHECK", "0") != "0"
    stage_started = time.perf_counter()
    try:
        max_attempts = int(os.environ.get(
            "REVAC_LABRADOR_Q96_PROVE_ATTEMPTS", "8"))
    except ValueError as exc:
        raise ValueError(
            "REVAC_LABRADOR_Q96_PROVE_ATTEMPTS must be an integer") from exc
    if not 1 <= max_attempts <= 32:
        raise ValueError(
            "REVAC_LABRADOR_Q96_PROVE_ATTEMPTS must be in [1, 32]")
    from lazer import lib
    from labrados import pack_verify
    proof = None
    for attempt in range(1, max_attempts + 1):
        stage_started = time.perf_counter()
        ps, params, proof = build_statement(
            entry["public"], prover, link, True, check_witness=clear_check,
            cache=cache,
            reuse_accumulator_witness=(reuse and not fresh_statement))
        print(f"[q96 worker timing] statement/rebind + Pack Prove "
              f"(attempt {attempt}/{max_attempts}): "
              f"{time.perf_counter() - stage_started:.6f} s", flush=True)
        relation_valid = getattr(
            lib, "labrador38_py_last_prove_tail_valid", None)
        if relation_valid is not None and int(relation_valid()) != 1:
            lib.labrador38_py_free_proof(proof)
            proof = None
            print(f"[q96 worker] rejected invalid native relation attempt "
                  f"{attempt}/{max_attempts}",
                  flush=True)
            continue
        full_self_verify = os.environ.get(
            "REVAC_LABRADOR_Q96_SELF_VERIFY", "0") != "0"
        if relation_valid is not None and not full_self_verify:
            print("[q96 worker] native prove relation guard: passed",
                  flush=True)
            break
        stage_started = time.perf_counter()
        if pack_verify((params, proof), ps.output_statement(), "38"):
            print(f"[q96 worker timing] prover-side proof self-verify: "
                  f"{time.perf_counter() - stage_started:.6f} s", flush=True)
            break
        lib.labrador38_py_free_proof(proof)
        proof = None
        print(f"[q96 worker] rejected invalid native proof attempt "
              f"{attempt}/{max_attempts}", flush=True)
    if proof is None:
        raise ValueError(
            f"q96 proof generation guard failed after {max_attempts} attempts")
    stage_started = time.perf_counter()
    proof_bytes = serialize_labrador38_proof(
        params, proof, require_zk=True,
        compact_profile_id=cache.bind_profile_id())
    print(f"[q96 worker timing] compact proof serialization: "
          f"{time.perf_counter() - stage_started:.6f} s", flush=True)
    # The resident cache owns only ps and params.  Every returned proof is a
    # per-request native allocation, including the first (cold) request.
    # Retaining that first proof leaked several GiB in a long-lived worker.
    lib.labrador38_py_free_proof(proof)
    if not reuse:
        # Keep only the most recently proved local witness.  It is bound to
        # PP ID, sequence, handle and the exact packed accumulator before a
        # compact request may reuse it.  The proof itself is never reused.
        entry["last_prover"] = prover
    print(f"[q96 worker timing] prove request total: "
          f"{time.perf_counter() - total_started:.6f} s", flush=True)
    return proof_bytes, cache_hit


def prepare_prover(dynamic: bytes, entry: dict[str, Any]) -> bool:
    """Build the resident prover skeleton without producing a proof."""
    if len(dynamic) != PROVER_BYTES + LINK_BYTES:
        raise ValueError("invalid q96 prover-prewarm payload length")
    prover = ProverPayload(dynamic[:PROVER_BYTES])
    link = LinkPayload.decode(dynamic[PROVER_BYTES:])
    _bind_deployment_profile(entry, link.context_hash)
    cache: Q96StatementCache = entry["prove_statement"]
    cache_hit = cache.ps is not None
    started = time.perf_counter()
    build_statement(
        entry["public"], prover, link, False, check_witness=False,
        cache=cache, prepare_prover=True)
    entry["last_prover"] = prover
    print(f"[q96 worker timing] prover statement prewarm: "
          f"{time.perf_counter() - started:.6f} s", flush=True)
    return cache_hit


def prepare_verifier(dynamic: bytes, entry: dict[str, Any]) -> bool:
    """Build the public verifier skeleton without requiring proof bytes."""
    if len(dynamic) != INSTANCE_BYTES:
        raise ValueError("invalid q96 verifier-prewarm payload length")
    instance = InstancePayload.decode(dynamic)
    _bind_deployment_profile(entry, instance.context_hash)
    cache: Q96StatementCache = entry["verify_statement"]
    cache_hit = cache.ps is not None
    started = time.perf_counter()
    build_statement(
        entry["public"], ZeroProverPayload(instance),
        PublicLinkPayload(instance), False, check_witness=False, cache=cache)
    print(f"[q96 worker timing] verifier statement prewarm: "
          f"{time.perf_counter() - started:.6f} s", flush=True)
    return cache_hit


def verify(dynamic: bytes, entry: dict[str, Any]) -> tuple[bool, bool]:
    total_started = time.perf_counter()
    if len(dynamic) < INSTANCE_BYTES + 4:
        raise ValueError("truncated q96 worker verify payload")
    proof_len = struct.unpack_from("<I", dynamic, INSTANCE_BYTES)[0]
    if (proof_len == 0 or
            INSTANCE_BYTES + 4 + proof_len != len(dynamic)):
        raise ValueError("invalid q96 worker proof length")
    stage_started = time.perf_counter()
    instance = InstancePayload.decode(dynamic[:INSTANCE_BYTES])
    _bind_deployment_profile(entry, instance.context_hash)
    prover = ZeroProverPayload(instance)
    link = PublicLinkPayload(instance)
    cache: Q96StatementCache = entry["verify_statement"]
    cache_hit = cache.ps is not None
    ps, _, _ = build_statement(
        entry["public"], prover, link, False, check_witness=False,
        cache=cache)
    print(f"[q96 worker timing] decode/rebind verifier instance: "
          f"{time.perf_counter() - stage_started:.6f} s", flush=True)
    stage_started = time.perf_counter()
    params, proof, refs = deserialize_labrador38_proof(
        dynamic[INSTANCE_BYTES + 4:], require_zk=True,
        reset_comkey=False, cached_params_obj=cache,
        expected_profile_id=cache.bind_profile_id())
    print(f"[q96 worker timing] compact proof deserialization: "
          f"{time.perf_counter() - stage_started:.6f} s", flush=True)
    from lazer import lib
    from labrados import pack_verify
    stage_started = time.perf_counter()
    try:
        ok = bool(pack_verify((params, proof), ps.output_statement(), "38"))
    finally:
        lib.labrador38_py_free_proof(proof)
        refs.clear()
    print(f"[q96 worker timing] native Pack Verify: "
          f"{time.perf_counter() - stage_started:.6f} s", flush=True)
    print(f"[q96 worker timing] verify request total: "
          f"{time.perf_counter() - total_started:.6f} s", flush=True)
    return ok, cache_hit


def serve(path: str, cache_limit: int, max_requests: int = 0) -> int:
    entries: OrderedDict[bytes, dict[str, Any]] = OrderedDict()
    if os.path.exists(path):
        try:
            probe = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            probe.connect(path)
        except OSError:
            os.unlink(path)
        else:
            probe.close()
            raise RuntimeError(f"q96 worker already owns {path}")
    server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    server.bind(path)
    os.chmod(path, 0o600)
    server.listen(8)

    def cleanup(_signum: int, _frame: Any) -> None:
        server.close()
        try:
            os.unlink(path)
        except FileNotFoundError:
            pass
        raise SystemExit(0)

    signal.signal(signal.SIGTERM, cleanup)
    signal.signal(signal.SIGINT, cleanup)
    handled_requests = 0
    while True:
        conn, _ = server.accept()
        with conn:
            try:
                request = recv_frame(conn)
                if len(request) < 1 + HASH_BYTES + 4:
                    raise ValueError("truncated q96 worker request")
                command = request[:1]
                pp_id = request[1:1 + HASH_BYTES]
                dynamic_len = struct.unpack_from(
                    "<I", request, 1 + HASH_BYTES)[0]
                offset = 1 + HASH_BYTES + 4
                if dynamic_len == 0 or offset + dynamic_len > len(request):
                    raise ValueError("invalid q96 dynamic length")
                dynamic = request[offset:offset + dynamic_len]
                entry = entries.get(pp_id)
                if command in (b"P", b"R", b"V", b"U", b"G") and entry is None:
                    send_frame(conn, b"N")
                    continue
                if command in (b"I", b"J", b"K", b"L"):
                    public_offset = offset + dynamic_len
                    if public_offset + 4 > len(request):
                        raise ValueError("missing q96 public length")
                    public_len = struct.unpack_from(
                        "<I", request, public_offset)[0]
                    public_start = public_offset + 4
                    if (public_len == 0 or
                            public_start + public_len != len(request)):
                        raise ValueError("invalid q96 public length")
                    public_raw = request[public_start:]
                    if hashlib.shake_256(public_raw).digest(HASH_BYTES) != pp_id:
                        raise ValueError("q96 public parameter identifier mismatch")
                    entry = {
                        "pp_id": pp_id,
                        "public": PublicParameters.decode(public_raw),
                        "prove_statement": Q96StatementCache(),
                        "verify_statement": Q96StatementCache(),
                    }
                    entries[pp_id] = entry
                    entries.move_to_end(pp_id)
                    while len(entries) > cache_limit:
                        entries.popitem(last=False)
                    command = {
                        b"I": b"P", b"J": b"V",
                        b"K": b"U", b"L": b"G",
                    }[command]
                assert entry is not None
                entries.move_to_end(pp_id)
                if command == b"P":
                    proof, hit = prove(dynamic, entry)
                    send_frame(conn, b"S" + bytes((hit,)) + proof)
                elif command == b"R":
                    try:
                        proof, hit = prove(dynamic, entry, reuse=True)
                    except ReuseMiss:
                        send_frame(conn, b"N")
                    else:
                        send_frame(conn, b"S" + bytes((hit,)) + proof)
                elif command == b"V":
                    ok, hit = verify(dynamic, entry)
                    send_frame(conn, b"S" + bytes((hit, ok)))
                elif command == b"U":
                    hit = prepare_prover(dynamic, entry)
                    send_frame(conn, b"S" + bytes((hit,)))
                elif command == b"G":
                    hit = prepare_verifier(dynamic, entry)
                    send_frame(conn, b"S" + bytes((hit,)))
                else:
                    raise ValueError("unknown q96 worker command")
            except Exception as exc:
                send_frame(conn, b"F" +
                           str(exc).encode("utf-8", "replace")[:1024])
        handled_requests += 1
        if max_requests and handled_requests >= max_requests:
            # A normal process exit is important for profile-guided builds:
            # GCC writes its accumulated counters while unloading the shared
            # library.  Production leaves this at zero and remains resident.
            # Cache-miss negotiation frames are deliberately not counted.
            break
    server.close()
    try:
        os.unlink(path)
    except FileNotFoundError:
        pass
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--socket", required=True)
    parser.add_argument("--cache-limit", type=int, default=2)
    parser.add_argument(
        "--max-requests", type=int,
        default=int(os.environ.get("REVAC_LABRADOR_Q96_MAX_REQUESTS", "0")),
        help=("exit normally after this many handled prove/verify requests "
              "(zero keeps the worker resident)"))
    args = parser.parse_args()
    if args.max_requests < 0:
        parser.error("--max-requests must be non-negative")
    return serve(args.socket, args.cache_limit, args.max_requests)


if __name__ == "__main__":
    raise SystemExit(main())
