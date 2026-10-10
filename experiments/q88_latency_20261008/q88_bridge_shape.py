#!/usr/bin/env python3
"""Shape audit for a q88 CRT-to-single-Pack quotient/carry bridge.

This is intentionally a sizing/compiler audit, not a production proof.  It
uses the q88 public modulus and the conservative 50,000-deletion witness gate
to calculate the radix-128 digit counts that a real CRT-combined full-ring
statement must support.  The production adapter must still construct the
combined public polynomials and prove every carry equation.
"""

from __future__ import annotations

import json
import math


LANE_MODULI = (1_073_692_673, 536_813_569, 536_690_689)
PACK_MODULUS = (1 << 38) - 107
BASE = 128
DEGREE = 4096
NON_TAG_WITNESS_POLYS = 174
GROUP = 12
GADGET_DIGITS = 6
HANDLE_BITS = 24
REFRESH_BLOCKS = 6


def digits_for_unsigned(bits: int) -> int:
    return (bits + 6) // 7


def digits_for_signed_bits(bits: int) -> int:
    # A balanced radix-B representation with d digits covers approximately
    # d*log2(B)-1 signed bits.  The +1 guard keeps the endpoint representable.
    return (bits + 7) // 7


def main() -> None:
    modulus = math.prod(LANE_MODULI)
    modulus_bits = modulus.bit_length()

    # This is the bound already enforced by OURS for the complete 50,000
    # deletion window: 50000 * 6 * N * 13003^2 + (2^30-1).
    witness_bound = (
        50_000 * 6 * DEGREE * 13_003 * 13_003 + ((1 << 30) - 1)
    )
    witness_bits = witness_bound.bit_length()
    public_digits = digits_for_unsigned(modulus_bits)
    witness_digits = digits_for_signed_bits(witness_bits)

    # Conservative quotient estimate for one matrix relation.  The exact
    # compiler must replace this with a norm-derived bound from every public
    # polynomial and every witness digit.  This estimate includes the
    # negacyclic sum and all q88 witness slots.
    quotient_bits = max(
        1,
        1 + witness_bits + modulus_bits +
        math.ceil(math.log2(DEGREE * NON_TAG_WITNESS_POLYS)) - modulus_bits,
    )
    quotient_digits = digits_for_signed_bits(quotient_bits)

    # Carry magnitude is profile dependent.  Reserve one extra radix layer
    # beyond the larger public/quotient convolution depth until exact bounds
    # are available.
    equation_count = max(
        public_digits + witness_digits - 1,
        public_digits + quotient_digits - 1,
    ) + 1
    carry_digits = equation_count

    witness_vectors = math.ceil(
        NON_TAG_WITNESS_POLYS * witness_digits / GROUP
    )
    quotient_vectors = quotient_digits
    carry_vectors = carry_digits
    estimated_vectors = witness_vectors + quotient_vectors + carry_vectors

    # A grouped full-ring statement would have one relation per digit layer;
    # refresh and same-handle constraints are added only after this audit.
    result = {
        "crt_modulus_bits": modulus_bits,
        "pack_modulus_bits": PACK_MODULUS.bit_length(),
        "radix": BASE,
        "witness_gate_bits": witness_bits,
        "public_digits": public_digits,
        "witness_digits": witness_digits,
        "quotient_bits_estimate": quotient_bits,
        "quotient_digits": quotient_digits,
        "carry_layers": carry_digits,
        "group_size": GROUP,
        "witness_vectors_estimate": witness_vectors,
        "quotient_vectors": quotient_vectors,
        "carry_vectors": carry_vectors,
        "single_relation_vectors_estimate": estimated_vectors,
        "equation_layers_estimate": equation_count,
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
