#!/usr/bin/env python3
"""Small exact KAT for the q88 CRT quotient/carry bridge.

It intentionally runs on degree 8 so the reference convolution is easy to
inspect.  Production integration must feed the same compiler the q88 public
polynomials and auxiliary product witnesses from the full-ring relation.
"""
from q88_bridge_core import (CRT_MODULUS, LANE_MODULI,
                             compile_linear_crt_bridge, negacyclic_mul)


def main() -> None:
    witness = [[2, -2, 3, 0, 0, 0, 0, 0]]
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

    # Auxiliary hidden-product checks are part of the bridge boundary.
    from q88_bridge_core import (aggregate_handle_witness,
                                 compose_refresh_witness,
                                 verify_auxiliary_equality)
    rights = [[1, 0, 0, 0, 0, 0, 0, 0] for _ in range(12)]
    aggregate = aggregate_handle_witness(rights, [1, 0],
                                         gadget_base=7, gadget_digits=6)
    verify_auxiliary_equality(
        aggregate, [sum(7 ** digit for digit in range(6)), 0, 0, 0,
                     0, 0, 0, 0], label="handle")
    refresh = compose_refresh_witness(
        [[3, 0, 0, 0, 0, 0, 0, 0],
         [5, 0, 0, 0, 0, 0, 0, 0],
         [11, 0, 0, 0, 0, 0, 0, 0]], [1, 0])
    verify_auxiliary_equality(refresh, [8, 0, 0, 0, 0, 0, 0, 0],
                              label="refresh")
    print("auxiliary hidden-product checks: PASS")


if __name__ == "__main__":
    main()
