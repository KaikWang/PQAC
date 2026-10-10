#!/usr/bin/env python3
"""Small exact KAT for the q88 CRT quotient/carry bridge.

It intentionally runs on degree 8 so the reference convolution is easy to
inspect.  Production integration must feed the same compiler the q88 public
polynomials and auxiliary product witnesses from the full-ring relation.
"""
from q88_bridge_core import (CRT_MODULUS, LANE_MODULI,
                             compile_linear_crt_bridge, negacyclic_mul)


def main() -> None:
    witness = [[1, -2, 3, 0, 0, 0, 0, 0]]
    # Pick a canonical 88-bit coefficient close to M.  Multiplying by two
    # then wraps once, so the KAT exercises a nonzero quotient digit.
    matrix = [CRT_MODULUS - 5] + [0] * 7
    matrix = [matrix]
    product = negacyclic_mul(matrix[0], witness[0])
    lane_matrix = []
    lane_target = []
    for modulus in LANE_MODULI:
        lane_matrix.append([[value % modulus for value in matrix[0]]])
        lane_target.append([value % modulus for value in product])
    bridge = compile_linear_crt_bridge(
        lane_matrix, witness, lane_target,
        base=128, witness_digits=3, quotient_digits=3)
    bridge.verify()
    print("q88 CRT quotient/carry KAT: PASS")
    print(f"carry layers: {len(bridge.carry_digits)}")
    print(f"matrix radix layers: {len(bridge.matrix_digits[0])}")
    print(f"quotient coefficient: {bridge.quotient_digits[0][0]}")


if __name__ == "__main__":
    main()
