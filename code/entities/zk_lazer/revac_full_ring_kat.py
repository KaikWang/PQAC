#!/usr/bin/env python3
"""Known-answer test for LaZer's degree-4096 extension-ring mapping."""

from __future__ import annotations

import argparse
import time

from lazer import polyring_t, poly_t
from labrados import pack_verify, proof_statement


PROOF_MODULUS = 2**38 - 107
FULL_DEGREE = 4096


def _sparse_poly(ring: polyring_t, terms: dict[int, int]) -> poly_t:
    coefficients = [0] * ring.deg
    for exponent, coefficient in terms.items():
        coefficients[exponent] = coefficient
    return poly_t(ring, coefficients)


def _check_bad_degree_zero(ring: polyring_t) -> None:
    bad_scalar = _sparse_poly(ring, {0: 1, 257: 1})
    statement = proof_statement(
        [ring.deg], [1], [0], num_constraints=0, primesize="38",
        num_deg0_constraints=1, zk=True)
    bad_index = statement.append_witness(bad_scalar)
    statement.append_deg0_statement(bad_index)
    if int(statement.simple_verify(statement.statement, statement.witness)) != 0:
        raise RuntimeError("degree-4096 degree-zero predicate accepted a high coefficient")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--skip-pack", action="store_true",
        help="run only clear-relation and degree-zero mapping checks")
    args = parser.parse_args()

    ring = polyring_t(FULL_DEGREE, PROOF_MODULUS)
    one = _sparse_poly(ring, {0: 1})
    left = _sparse_poly(ring, {
        0: 2, 1: -1, 255: 3, 256: 1, 2049: -2, 4095: 1,
    })
    right = _sparse_poly(ring, {
        0: -1, 17: 2, 511: 1, 2048: -1, 3073: 2,
    })
    public_multiplier = _sparse_poly(ring, {
        0: 3, 129: -1, 1024: 2, 4095: -1,
    })
    scalar = _sparse_poly(ring, {0: 1})

    linear_rhs = public_multiplier * left + right
    quadratic_rhs = left * right + scalar

    statement = proof_statement(
        [FULL_DEGREE] * 3,
        [1, 1, 1],
        [2**20, 2**20, 0],
        num_constraints=2,
        primesize="38",
        num_deg0_constraints=1,
        zk=True,
    )
    left_index = statement.append_witness(left)
    right_index = statement.append_witness(right)
    scalar_index = statement.append_witness(scalar)
    statement.append_deg0_statement(scalar_index)
    statement.append_statement(
        [public_multiplier, one], [left_index, right_index], linear_rhs)
    statement.append_quadratic_statement(
        [one], [scalar_index], [1], [left_index], [right_index],
        quadratic_rhs)

    started = time.perf_counter()
    if int(statement.simple_verify(statement.statement, statement.witness)) != 1:
        raise RuntimeError("degree-4096 linear/quadratic relation mapping mismatch")
    clear_seconds = time.perf_counter() - started
    _check_bad_degree_zero(ring)

    print(f"degree-4096 clear relation: passed ({clear_seconds:.3f} s)")
    print("degree-4096 degree-zero rejection: passed")
    print(f"backend witness polynomials: {statement.witness_polys}")

    if args.skip_pack:
        return 0

    started = time.perf_counter()
    error, params, proof = statement.pack_prove()
    prove_seconds = time.perf_counter() - started
    if error != 0:
        raise RuntimeError(f"degree-4096 Pack Prove failed: {error}")

    started = time.perf_counter()
    verified = pack_verify((params, proof), statement.output_statement(), "38")
    verify_seconds = time.perf_counter() - started
    if not verified:
        raise RuntimeError("degree-4096 Pack Verify rejected the KAT proof")

    print(f"degree-4096 Pack Prove: {prove_seconds:.3f} s")
    print(f"degree-4096 Pack Verify: {verify_seconds:.3f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
