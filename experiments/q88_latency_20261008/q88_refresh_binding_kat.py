#!/usr/bin/env python3
"""Exact single-prime refresh-selector quotient/carry KAT.

This validates one Remark-4.3 refresh block in the independent field. It is
kept separate from the production CRT worker until all six blocks can share a
single statement without duplicating witness vectors.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import time

from q88_bridge_core import (SINGLE88_PRIME, build_grouped_pack_statement,
                             compile_single88_bridge_fast)
from q88_single_prime_profile import verify_modulus_certificate
from revac_acc96_wire import (PROVER_BYTES, ProverPayload, PublicParameters,
                              REFRESH_BITS, REFRESH_BLOCKS)


BASE = 128
BLOCK = 0


def center(value):
    value = int(value) % SINGLE88_PRIME
    return value - SINGLE88_PRIME if value > SINGLE88_PRIME // 2 else value


def digits(value, count):
    out = []
    remaining = int(value)
    for _ in range(count):
        digit = (remaining + BASE // 2) % BASE - BASE // 2
        out.append(digit)
        remaining = (remaining - digit) // BASE
    if remaining:
        raise ValueError("refresh radix capacity exhausted")
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--public", default="/tmp/q88bridge_real/public.bin")
    parser.add_argument("--prover", default="/tmp/q88bridge_real/prover_link.bin")
    args = parser.parse_args()
    verify_modulus_certificate()
    public = PublicParameters.decode(Path(args.public).read_bytes())
    prover = ProverPayload(Path(args.prover).read_bytes()[:PROVER_BYTES])
    degree = len(prover.witness_poly(0))
    one = tuple([1] + [0] * (degree - 1))
    rows = [public.refresh(index, BLOCK) for index in range(REFRESH_BITS + 1)]
    raw = [int(v) for v in rows[0]]
    for bit, selected in enumerate(prover.refresh_r_bits):
        if selected:
            raw = [left + right for left, right in zip(raw, rows[bit + 1])]
    centered = [center(v) for v in raw]
    lift = [(int(value) - int(c)) // SINGLE88_PRIME
            for value, c in zip(raw, centered)]
    ahat = [[0] * degree for _ in range(13)]
    for index, value in enumerate(centered):
        for digit, item in enumerate(digits(value, 13)):
            ahat[digit][index] = item
    selector_polys = [tuple([int(bit)] + [0] * (degree - 1))
                      for bit in prover.refresh_r_bits]
    witness = selector_polys + [tuple(layer) for layer in ahat] + [tuple(lift)]
    matrix = [tuple(int(v) % SINGLE88_PRIME for v in rows[index + 1])
              for index in range(REFRESH_BITS)]
    matrix += [tuple((-(BASE ** digit)) % SINGLE88_PRIME if i == 0 else 0
                     for i in range(degree)) for digit in range(13)]
    matrix.append(tuple((-SINGLE88_PRIME) % SINGLE88_PRIME
                        if i == 0 else 0 for i in range(degree)))
    target = tuple((-int(v)) % SINGLE88_PRIME for v in rows[0])
    bridge_started = time.perf_counter()
    bridge = compile_single88_bridge_fast(
        matrix, witness, target=target, base=BASE,
        witness_digits=9, quotient_digits=16)
    bridge_seconds = time.perf_counter() - bridge_started
    extra = []
    deg0 = []
    for index in range(REFRESH_BITS):
        ref = ("w", 0, index)
        one_poly = tuple(one)
        extra.append({"linear": [((tuple(-v for v in one_poly)), ref)],
                      "quadratic": [(1, ref, ref)]})
        deg0.append(ref)
    statement_started = time.perf_counter()
    statement, _ = build_grouped_pack_statement(
        bridge, group_size=1, extra_quadratics=extra, extra_deg0=deg0)
    statement_seconds = time.perf_counter() - statement_started
    from labrados import ffi, lib, pack_verify
    prove_started = time.perf_counter()
    error, params, proof = statement.pack_prove()
    prove_seconds = time.perf_counter() - prove_started
    if error:
        raise RuntimeError(f"refresh binding Pack Prove failed: {error}")
    verify_started = time.perf_counter()
    if not pack_verify((params, proof), statement.output_statement(), "38"):
        raise RuntimeError("refresh binding Pack Verify failed")
    verify_seconds = time.perf_counter() - verify_started
    size = lib.labrador38_py_proof_nbytes(params, proof)
    serialized = ffi.new("uint8_t[]", size)
    if lib.labrador38_py_proof_serialize(serialized, size, params, proof):
        raise RuntimeError("refresh binding serialization failed")
    decoded = ffi.new("labrador38_dch_pack_proof")
    if lib.labrador38_py_proof_deserialize(decoded, params, serialized, size):
        raise RuntimeError("refresh binding deserialization failed")
    if not pack_verify((params, decoded), statement.output_statement(), "38"):
        raise RuntimeError("refresh binding reload Verify failed")
    print(f"refresh_block={BLOCK} selectors={REFRESH_BITS} modulus={SINGLE88_PRIME}")
    print(f"bridge={bridge_seconds:.3f}s statement={statement_seconds:.3f}s")
    print(f"clear=passed prove={prove_seconds:.3f}s verify={verify_seconds:.3f}s")
    print(f"serialized={int(size)} bytes ({int(size)/1024:.3f} KiB) reload=passed")


if __name__ == "__main__":
    main()
