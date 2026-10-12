#!/usr/bin/env python3
"""Run the real q88-shaped isolated single-88-bit profile.

The input witness is a captured CRT show only to keep the bounded witness
distribution and 174-column layout identical.  A fresh target is derived in
the independent field; this is intentionally reported as a profile result,
not as a production accumulator proof.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import time

from q88_bridge_core import (SINGLE88_PRIME, build_grouped_pack_statement,
                             compile_single88_bridge_fast)
from q88_bridge_relation import build_single88_effective_relation
from q88_single_prime_profile import verify_modulus_certificate
from revac_acc96_wire import PROVER_BYTES, LinkPayload, ProverPayload, PublicParameters


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--public", default="/tmp/q88bridge_real/public.bin")
    parser.add_argument("--prover-link", default="/tmp/q88bridge_real/prover_link.bin")
    # 16 is the smallest tested grouping that reduced Pack reducer work while
    # keeping the exact norm gate and both serialized verification checks.
    parser.add_argument("--group-size", type=int, default=16)
    args = parser.parse_args()
    verify_modulus_certificate()
    public = PublicParameters.decode(Path(args.public).read_bytes())
    raw = Path(args.prover_link).read_bytes()
    prover = ProverPayload(raw[:PROVER_BYTES])
    link = LinkPayload.decode(raw[PROVER_BYTES:])
    link.audit_relation(prover.handle_vector)
    started = time.perf_counter()
    matrix, witness = build_single88_effective_relation(
        public, prover, modulus=SINGLE88_PRIME)
    bridge = compile_single88_bridge_fast(
        matrix, witness, base=128, witness_digits=9, quotient_digits=16)
    bridge_seconds = time.perf_counter() - started
    statement_started = time.perf_counter()
    statement, _ = build_grouped_pack_statement(bridge, group_size=args.group_size)
    statement_seconds = time.perf_counter() - statement_started
    if int(statement.simple_verify(statement.statement, statement.witness)) != 1:
        raise RuntimeError("single88 profile Clear Verify failed")
    prove_started = time.perf_counter()
    error, params, proof = statement.pack_prove()
    prove_seconds = time.perf_counter() - prove_started
    if error:
        raise RuntimeError(f"single88 profile Pack Prove failed: {error}")
    from labrados import ffi, lib, pack_verify
    verify_started = time.perf_counter()
    if not pack_verify((params, proof), statement.output_statement(), "38"):
        raise RuntimeError("single88 profile Pack Verify failed")
    verify_seconds = time.perf_counter() - verify_started
    size = lib.labrador38_py_proof_nbytes(params, proof)
    serialized = ffi.new("uint8_t[]", size)
    if lib.labrador38_py_proof_serialize(serialized, size, params, proof):
        raise RuntimeError("single88 profile serialization failed")
    decoded = ffi.new("labrador38_dch_pack_proof")
    if lib.labrador38_py_proof_deserialize(decoded, params, serialized, size):
        raise RuntimeError("single88 profile deserialization failed")
    if not pack_verify((params, decoded), statement.output_statement(), "38"):
        raise RuntimeError("single88 profile deserialized Pack Verify failed")
    print(f"single88 modulus={SINGLE88_PRIME}")
    print(f"columns={len(matrix)} witness_degree={len(witness[0])}")
    print(f"bridge={bridge_seconds:.3f}s statement={statement_seconds:.3f}s")
    print(f"clear=passed prove={prove_seconds:.3f}s verify={verify_seconds:.3f}s")
    print(f"serialized={int(size)} bytes ({int(size)/1024:.3f} KiB) reload=passed")


if __name__ == "__main__":
    main()
