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
from q88_bridge_relation import build_single88_handle_bound_relation
from q88_single_prime_profile import verify_modulus_certificate
from revac_acc96_wire import (PROVER_BYTES, LinkPayload, ProverPayload,
                              PublicParameters, Q19, Q19_D, Q19_N,
                              LINK_DECOMP_BASE, LINK_DECOMP_LIMBS)


def single88_binding_constraints(witness, handle_meta, link):
    extra_witness = []
    extra_quadratics = []
    extra_deg0 = []
    extra_linear = []
    extra_integers = []
    if handle_meta is not None:
        one = [1] + [0] * (len(witness[0]) - 1)
        for bit, value in enumerate(handle_meta["handle_bits"]):
            ref = ("bit", bit)
            extra_witness.append({"ref": ref,
                                  "values": [tuple([value] + [0] * (len(one) - 1))],
                                  "norm": 0})
            extra_deg0.append(ref)
            extra_quadratics.append({
                "linear": [((tuple(-v for v in one)), ref)],
                "quadratic": [(1, ref, ref)],
            })
        for digit, product_index in enumerate(handle_meta["product_refs"]):
            terms = [(1, ("bit", bit), ("w", 0,
                                           handle_meta["right_refs"][digit][bit]))
                     for bit in range(len(handle_meta["handle_bits"]))]
            extra_quadratics.append({
                "linear": [((tuple(-v for v in one)), ("w", 0, product_index))],
                "quadratic": terms,
            })
        # q19 link rows are degree-256 linear relations.  The q19 handle is a
        # separate scalar representation and is tied coefficient 0 to the
        # degree-4096 handle witness with integer constraints below.
        q19_one = tuple([1] + [0] * (Q19_N - 1))
        q19_bits = []
        for bit, value in enumerate(handle_meta["handle_bits"]):
            ref = ("q19bit", bit)
            q19_bits.append(ref)
            extra_witness.append({
                "ref": ref, "degree": Q19_N,
                "values": [tuple([value] + [0] * (Q19_N - 1))],
                "norm": 0,
            })
            extra_deg0.append(ref)
            coefficients = [0] * len(witness[0])
            coefficients[0] = 1
            full_ref = ("bit", bit)
            small_coefficients = [0] * Q19_N
            small_coefficients[0] = -1
            extra_integers.append({
                "coefficients": coefficients + small_coefficients,
                "refs": [full_ref, ref], "target": 0,
            })
        link_a, link_d = link.public_matrices()
        def q19_center(value):
            value = int(value) % Q19
            return value - Q19 if value > Q19 // 2 else value
        def shift(poly, amount):
            output = [0] * Q19_N
            for index, raw_value in enumerate(poly):
                value = q19_center(raw_value)
                destination = index + amount
                if destination < Q19_N:
                    output[destination] += value
                else:
                    output[destination - Q19_N] -= value
            return tuple(q19_center(value) for value in output)
        for part in range(2):
            if part == 0:
                for row in range(Q19_D):
                    ref = ("q19r0", row)
                    extra_witness.append({
                        "ref": ref, "degree": Q19_N,
                        "values": [link.randomness[0][row]], "norm": 0,
                    })
            else:
                for col in range(Q19_D):
                    ref = ("q19r1", col)
                    extra_witness.append({
                        "ref": ref, "degree": Q19_N,
                        "values": [link.randomness[1][col]], "norm": 0,
                    })
        for row in range(Q19_D):
            for limb in range(LINK_DECOMP_LIMBS):
                ref = ("q19lift", row, limb)
                extra_witness.append({
                    "ref": ref, "degree": Q19_N,
                    "values": [link.lift_limbs[row][limb]],
                    "norm": Q19_N * (LINK_DECOMP_BASE // 2) ** 2,
                })
            linear = []
            linear.append((q19_one, ("q19r0", row)))
            for col in range(Q19_D):
                linear.append((tuple(q19_center(v) for v in link_a[row][col]),
                               ("q19r1", col)))
            for bit, ref in enumerate(q19_bits):
                linear.append((shift(link_d[row][0], bit), ref))
            for limb in range(LINK_DECOMP_LIMBS):
                scalar = Q19 * (LINK_DECOMP_BASE ** limb)
                linear.append((tuple((scalar if i == 0 else 0)
                                     for i in range(Q19_N)),
                               ("q19lift", row, limb)))
            target = tuple(q19_center(v) for v in link.commitment[row])
            extra_linear.append({"degree": Q19_N, "linear": linear,
                                 "target": target})
    return extra_witness, extra_quadratics, extra_deg0, extra_linear, extra_integers


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--public", default="/tmp/q88bridge_real/public.bin")
    parser.add_argument("--prover-link", default="/tmp/q88bridge_real/prover_link.bin")
    # 16 is the smallest tested grouping that reduced Pack reducer work while
    # keeping the exact norm gate and both serialized verification checks.
    parser.add_argument("--group-size", type=int, default=16)
    parser.add_argument("--handle-bound", action="store_true",
                        help="bind handle bits and product witnesses in Pack")
    args = parser.parse_args()
    verify_modulus_certificate()
    public = PublicParameters.decode(Path(args.public).read_bytes())
    raw = Path(args.prover_link).read_bytes()
    prover = ProverPayload(raw[:PROVER_BYTES])
    link = LinkPayload.decode(raw[PROVER_BYTES:])
    link.audit_relation(prover.handle_vector)
    started = time.perf_counter()
    if args.handle_bound:
        matrix, witness, handle_meta = build_single88_handle_bound_relation(
            public, prover, modulus=SINGLE88_PRIME)
    else:
        matrix, witness = build_single88_effective_relation(
            public, prover, modulus=SINGLE88_PRIME)
        handle_meta = None
    bridge = compile_single88_bridge_fast(
        matrix, witness, base=128, witness_digits=9, quotient_digits=16)
    bridge_seconds = time.perf_counter() - started
    statement_started = time.perf_counter()
    (extra_witness, extra_quadratics, extra_deg0, extra_linear,
     extra_integers) = single88_binding_constraints(witness, handle_meta, link)
    statement, _ = build_grouped_pack_statement(
        bridge, group_size=(1 if handle_meta is not None else args.group_size),
        extra_witness=extra_witness, extra_quadratics=extra_quadratics,
        extra_deg0=extra_deg0, extra_linear=extra_linear,
        extra_integers=extra_integers)
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
    print(f"handle_binding={'enabled' if handle_meta is not None else 'effective-profile'}")
    print(f"columns={len(matrix)} witness_degree={len(witness[0])}")
    print(f"bridge={bridge_seconds:.3f}s statement={statement_seconds:.3f}s")
    print(f"clear=passed prove={prove_seconds:.3f}s verify={verify_seconds:.3f}s")
    print(f"serialized={int(size)} bytes ({int(size)/1024:.3f} KiB) reload=passed")


if __name__ == "__main__":
    main()
