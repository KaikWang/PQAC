#!/usr/bin/env python3
"""Six-block refresh selector proof in one native Pack statement."""
from __future__ import annotations

from pathlib import Path
import time

from q88_bridge_core import (SINGLE88_PRIME, build_multi_bridge_pack_statement,
                             compile_single88_bridge_fast)
from q88_single_prime_profile import verify_modulus_certificate
from revac_acc96_wire import (PROVER_BYTES, ProverPayload, PublicParameters,
                              REFRESH_BITS, REFRESH_BLOCKS)

BASE = 128


def center(value):
    value = int(value) % SINGLE88_PRIME
    return value - SINGLE88_PRIME if value > SINGLE88_PRIME // 2 else value


def make_bridge(public, prover, block):
    n = len(prover.witness_poly(0))
    rows = [public.refresh(index, block) for index in range(REFRESH_BITS + 1)]
    raw = [int(v) for v in rows[0]]
    for bit, selected in enumerate(prover.refresh_r_bits):
        if selected:
            raw = [left + right for left, right in zip(raw, rows[bit + 1])]
    centered = [center(v) for v in raw]
    lift = [(int(v) - int(c)) // SINGLE88_PRIME
            for v, c in zip(raw, centered)]
    ahat = [[0] * n for _ in range(13)]
    for index, value in enumerate(centered):
        remaining = int(value)
        for digit in range(13):
            item = (remaining + BASE // 2) % BASE - BASE // 2
            ahat[digit][index] = item
            remaining = (remaining - item) // BASE
        if remaining:
            raise ValueError("refresh radix capacity exhausted")
    witness = [
        tuple([int(bit)] + [0] * (n - 1))
        for bit in prover.refresh_r_bits
    ] + [tuple(layer) for layer in ahat] + [tuple(lift)]
    matrix = [tuple(int(v) % SINGLE88_PRIME for v in rows[index + 1])
              for index in range(REFRESH_BITS)]
    matrix += [tuple((-(BASE ** digit)) % SINGLE88_PRIME if i == 0 else 0
                     for i in range(n)) for digit in range(13)]
    matrix.append(tuple((-SINGLE88_PRIME) % SINGLE88_PRIME
                        if i == 0 else 0 for i in range(n)))
    target = tuple((-int(v)) % SINGLE88_PRIME for v in rows[0])
    return compile_single88_bridge_fast(
        matrix, witness, target=target, base=BASE,
        witness_digits=9, quotient_digits=16)


def main():
    verify_modulus_certificate()
    public = PublicParameters.decode(Path("/tmp/q88bridge_real/public.bin").read_bytes())
    raw = Path("/tmp/q88bridge_real/prover_link.bin").read_bytes()
    prover = ProverPayload(raw[:PROVER_BYTES])
    started = time.perf_counter()
    bridges = [make_bridge(public, prover, block)
               for block in range(REFRESH_BLOCKS)]
    bridge_seconds = time.perf_counter() - started
    statement_started = time.perf_counter()
    statement, _ = build_multi_bridge_pack_statement(bridges, group_size=1)
    statement_seconds = time.perf_counter() - statement_started
    from labrados import ffi, lib, pack_verify
    prove_started = time.perf_counter()
    error, params, proof = statement.pack_prove()
    prove_seconds = time.perf_counter() - prove_started
    if error:
        raise RuntimeError(f"six-block refresh Pack Prove failed: {error}")
    verify_started = time.perf_counter()
    if not pack_verify((params, proof), statement.output_statement(), "38"):
        raise RuntimeError("six-block refresh Pack Verify failed")
    verify_seconds = time.perf_counter() - verify_started
    size = lib.labrador38_py_proof_nbytes(params, proof)
    serialized = ffi.new("uint8_t[]", size)
    if lib.labrador38_py_proof_serialize(serialized, size, params, proof):
        raise RuntimeError("six-block refresh serialization failed")
    decoded = ffi.new("labrador38_dch_pack_proof")
    if lib.labrador38_py_proof_deserialize(decoded, params, serialized, size):
        raise RuntimeError("six-block refresh deserialization failed")
    if not pack_verify((params, decoded), statement.output_statement(), "38"):
        raise RuntimeError("six-block refresh reload Verify failed")
    print(f"refresh_blocks={REFRESH_BLOCKS} selectors={REFRESH_BITS} modulus={SINGLE88_PRIME}")
    print(f"bridge={bridge_seconds:.3f}s statement={statement_seconds:.3f}s")
    print(f"clear=passed prove={prove_seconds:.3f}s verify={verify_seconds:.3f}s")
    print(f"serialized={int(size)} bytes ({int(size)/1024:.3f} KiB) reload=passed")


if __name__ == "__main__":
    main()
