"""End-to-end exact 100-bit congruence Pack proof on a small witness."""

import argparse
import time

from accumulator_pack_bridge import PACK_MODULUS, compile_bridge
from labrados import pack_verify, proof_statement
from lazer import polyring_t, poly_t

ACCUMULATOR_MODULUS = 1267650600228229401496702836737


def make_statement(bridge):
    ring = polyring_t(bridge.degree, PACK_MODULUS)
    vectors = []
    labels = []
    for col, digits in enumerate(bridge.w_digits):
        for layer, values in enumerate(digits):
            vectors.append((values, max(1, bridge.bounds["w"][col][layer])))
            labels.append(("w", col, layer))
    for layer, values in enumerate(bridge.quotient_digits):
        vectors.append((values, max(1, bridge.bounds["k"][layer])))
        labels.append(("k", layer))
    for layer, values in enumerate(bridge.carry_digits):
        vectors.append((values, max(1, bridge.bounds["c"][layer])))
        labels.append(("c", layer))
    statement = proof_statement(
        [bridge.degree] * len(vectors), [1] * len(vectors),
        [bound for _, bound in vectors], bridge.equation_count(),
        primesize="38", zk=True)
    indexes = {}
    for label, (values, _) in zip(labels, vectors):
        indexes[label] = statement.append_witness(poly_t(ring, values))
    one = poly_t(ring, {0: 1})
    minus_one = poly_t(ring, {0: -1})
    for layer in range(bridge.equation_count()):
        terms, witnesses = [], []
        for col, hdigits in enumerate(bridge.h_digits):
            for i, values in enumerate(hdigits):
                j = layer - i
                if 0 <= j < len(bridge.w_digits[col]):
                    terms.append(poly_t(ring, values))
                    witnesses.append(indexes[("w", col, j)])
        for i, digit in enumerate(bridge.q_digits):
            j = layer - i
            if digit and 0 <= j < len(bridge.quotient_digits):
                terms.append(poly_t(ring, {0: -digit}))
                witnesses.append(indexes[("k", j)])
        if layer:
            terms.append(one)
            witnesses.append(indexes[("c", layer - 1)])
        terms.append(poly_t(ring, {0: -bridge.base}))
        witnesses.append(indexes[("c", layer)])
        rhs = (bridge.target_digits[layer] if
               layer < len(bridge.target_digits) else [0] * bridge.degree)
        statement.append_statement(terms, witnesses, poly_t(ring, rhs))
    return statement


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree", type=int, default=256)
    args = parser.parse_args()
    n = args.degree
    if n < 256 or n & (n - 1):
        raise ValueError("Pack degree must be a power of two >= 256")
    public = [[0] * n]
    public[0][0] = ACCUMULATOR_MODULUS - 1
    private = [[0] * n]
    private[0][0] = -1
    target = [0] * n
    target[0] = 1
    bridge = compile_bridge(public, private, target, ACCUMULATOR_MODULUS,
                            witness_digits=1, quotient_digits=1)
    tampered = target.copy()
    tampered[0] = 2
    try:
        compile_bridge(public, private, tampered, ACCUMULATOR_MODULUS,
                       witness_digits=1, quotient_digits=1)
    except ValueError:
        pass
    else:
        raise AssertionError("tampered target was accepted")
    statement = make_statement(bridge)
    if int(statement.simple_verify(statement.statement, statement.witness)) != 1:
        raise RuntimeError("Pack clear relation failed")
    start = time.perf_counter()
    error, params, proof = statement.pack_prove()
    prove_seconds = time.perf_counter() - start
    if error:
        raise RuntimeError(f"Pack parameter generation failed: {error}")
    start = time.perf_counter()
    if not pack_verify((params, proof), statement.output_statement(), "38"):
        raise RuntimeError("Pack verification failed")
    verify_seconds = time.perf_counter() - start
    print({"result": "PASS", "accumulator_modulus_bits":
           ACCUMULATOR_MODULUS.bit_length(), "degree": n,
           "constraints": bridge.equation_count(),
           "pack_prove_seconds": prove_seconds,
           "pack_verify_seconds": verify_seconds})


if __name__ == "__main__":
    main()
