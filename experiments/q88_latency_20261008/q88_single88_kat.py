#!/usr/bin/env python3
"""Independent single-88-bit-prime quotient/carry Pack KAT.

This is deliberately a bridge profile, not the production q88 CRT worker.
It exercises the exact relation, native degree-4096 Clear Verify, Pack Verify,
and proof serialization/deserialization before any full-profile integration.
"""

from __future__ import annotations

import time

from q88_bridge_core import (SINGLE88_PRIME, build_grouped_pack_statement,
                             compile_single88_bridge)


def main() -> None:
    degree = 4096
    matrix = [[0] * degree]
    witness = [[0] * degree]
    matrix[0][0], matrix[0][1] = 3, 1
    witness[0][0], witness[0][1] = 5, -2
    target = [0] * degree
    target[0], target[1], target[2] = 15, -1, -2
    target = [value % SINGLE88_PRIME for value in target]

    started = time.perf_counter()
    bridge = compile_single88_bridge(
        matrix, witness, target, base=128, witness_digits=8,
        quotient_digits=10)
    bridge_seconds = time.perf_counter() - started
    statement, _ = build_grouped_pack_statement(bridge, group_size=1)
    clear_started = time.perf_counter()
    if int(statement.simple_verify(statement.statement, statement.witness)) != 1:
        raise RuntimeError("single-88 native clear verification failed")
    clear_seconds = time.perf_counter() - clear_started

    error, params, proof = statement.pack_prove()
    if error:
        raise RuntimeError(f"single-88 Pack Prove failed: {error}")
    from labrados import ffi, lib, pack_verify
    verify_started = time.perf_counter()
    if not pack_verify((params, proof), statement.output_statement(), "38"):
        raise RuntimeError("single-88 Pack Verify failed")
    verify_seconds = time.perf_counter() - verify_started
    size = lib.labrador38_py_proof_nbytes(params, proof)
    serialized = ffi.new("uint8_t[]", size)
    if lib.labrador38_py_proof_serialize(serialized, size, params, proof):
        raise RuntimeError("single-88 proof serialization failed")
    decoded = ffi.new("labrador38_dch_pack_proof")
    if lib.labrador38_py_proof_deserialize(decoded, params, serialized, size):
        raise RuntimeError("single-88 proof deserialization failed")
    if not pack_verify((params, decoded), statement.output_statement(), "38"):
        raise RuntimeError("single-88 deserialized Pack Verify failed")
    print(f"single88 modulus={SINGLE88_PRIME}")
    print(f"bridge={bridge_seconds:.3f}s clear={clear_seconds:.3f}s")
    print(f"pack_verify={verify_seconds:.3f}s serialized={size} bytes")


if __name__ == "__main__":
    main()
