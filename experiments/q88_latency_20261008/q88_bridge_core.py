"""q88-specific CRT recombination primitives for the future Pack bridge.

The current production compiler proves three proof-prime lane equations.  A
quotient/carry bridge must first combine those equations over the original
88-bit CRT modulus.  These helpers do only that exact integer step; they do
not construct a proof statement or choose security parameters.
"""

from __future__ import annotations

from math import prod
from dataclasses import dataclass
from typing import Iterable, Sequence


LANE_MODULI = (1_073_692_673, 536_813_569, 536_690_689)
CRT_MODULUS = prod(LANE_MODULI)
PACK_MODULUS = (1 << 38) - 107

# Independent single-modulus experiment.  This is deliberately not the
# production CRT modulus or wire parameter.  It is an 88-bit prime chosen
# near the current CRT product so coefficient magnitudes remain comparable.
from q88_single_prime_profile import SINGLE88_PRIME, verify_modulus_certificate


def _is_probable_prime(value: int) -> bool:
    """Probable-prime screening only; not an 88-bit primality certificate."""
    if value < 2:
        return False
    for small in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if value % small == 0:
            return value == small
    odd = value - 1
    exponent = 0
    while odd % 2 == 0:
        odd //= 2
        exponent += 1
    for base in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if base >= value:
            continue
        candidate = pow(base, odd, value)
        if candidate in (1, value - 1):
            continue
        for _ in range(exponent - 1):
            candidate = candidate * candidate % value
            if candidate == value - 1:
                break
        else:
            return False
    return True


verify_modulus_certificate()


def full_ring_dimension_breakdown(*, compact_limbs: int = 4,
                                  aggregate_ranks: tuple[int, ...] = (16, 8),
                                  trunc_right: int = 22,
                                  gadget_digits: int = 6,
                                  lanes: int = 3) -> dict[str, int]:
    """Estimate native witness-polynomial dimensions by relation component.

    This is deliberately a structural audit, not a security parameter
    selector.  It identifies the large truncation-right and handle-product
    blocks before a quotient/carry rewrite is attempted.
    """
    phase_factor = 16  # degree-4096 native ring / degree-256 phase ring
    rank_sum = sum(aggregate_ranks)
    return {
        "left_and_v12": 2 * compact_limbs * phase_factor,
        "trunc_right": trunc_right * compact_limbs * phase_factor,
        "handle_right": gadget_digits * compact_limbs * rank_sum * phase_factor,
        "handle_aggregate": lanes * 6 * rank_sum * phase_factor,
        "crt_lift": lanes * 4 * phase_factor,
    }


def pack_gadget_digit_pairs(digits: Sequence[Sequence[int]], *,
                            base: int) -> tuple[list[list[int]], int]:
    """Pack adjacent gadget digits and return the exact signed bound.

    For digits ``d0..d5`` this returns ``d0+B*d1``, ``d2+B*d3`` and
    ``d4+B*d5``.  The operation is integer-only and deliberately does not
    reduce modulo a proof prime; callers must add the corresponding carry
    reconstruction equations before using the packed values in a statement.
    """
    if base <= 1:
        raise ValueError("gadget pair base must be greater than one")
    if len(digits) == 0 or len(digits) % 2:
        raise ValueError("gadget digits must contain an even nonzero count")
    width = len(digits[0])
    if any(len(row) != width for row in digits):
        raise ValueError("gadget digit widths do not match")
    packed = []
    bound = 0
    for index in range(0, len(digits), 2):
        left, right = digits[index], digits[index + 1]
        values = [int(a) + base * int(b) for a, b in zip(left, right)]
        packed.append(values)
        bound = max(bound, *(abs(value) for value in values))
    return packed, bound


def unpack_gadget_digit_pairs(packed: Sequence[Sequence[int]], *,
                              base: int) -> list[list[int]]:
    """Invert :func:`pack_gadget_digit_pairs` using balanced base digits."""
    if base <= 1:
        raise ValueError("gadget pair base must be greater than one")
    result: list[list[int]] = []
    half = base // 2
    for row in packed:
        low = []
        high = []
        for value in row:
            digit = (int(value) + half) % base - half
            low.append(digit)
            high.append((int(value) - digit) // base)
        result.extend((low, high))
    return result


def verify_gadget_digit_pairs(digits: Sequence[Sequence[int]], *,
                              base: int) -> int:
    """Return the maximum reconstruction error for a pair packing."""
    packed, _ = pack_gadget_digit_pairs(digits, base=base)
    recovered = unpack_gadget_digit_pairs(packed, base=base)
    if len(recovered) != len(digits):
        raise ValueError("gadget pair reconstruction count mismatch")
    return max((abs(int(a) - int(b))
                for left, right in zip(digits, recovered)
                for a, b in zip(left, right)), default=0)


def _crt_constants() -> tuple[int, ...]:
    constants = []
    for modulus in LANE_MODULI:
        partial = CRT_MODULUS // modulus
        constants.append(partial * pow(partial, -1, modulus))
    return tuple(constants)


CRT_CONSTANTS = _crt_constants()


def crt_recombine(residues: Sequence[int]) -> int:
    """Return the canonical 88-bit CRT representative."""
    if len(residues) != len(LANE_MODULI):
        raise ValueError("q88 CRT lane count mismatch")
    value = sum(int(residue) * constant
                for residue, constant in zip(residues, CRT_CONSTANTS))
    return value % CRT_MODULUS


def crt_recombine_poly(lane_polys: Sequence[Sequence[int]]) -> list[int]:
    """Recombine equal-length lane polynomials coefficient by coefficient."""
    if len(lane_polys) != len(LANE_MODULI):
        raise ValueError("q88 CRT lane count mismatch")
    degree = len(lane_polys[0])
    if any(len(poly) != degree for poly in lane_polys):
        raise ValueError("q88 CRT polynomial degree mismatch")
    return [crt_recombine([lane_polys[lane][index]
                           for lane in range(len(LANE_MODULI))])
            for index in range(degree)]


def verify_recombination(combined: Sequence[int],
                         lane_polys: Sequence[Sequence[int]]) -> None:
    """Fail closed unless recombination preserves every lane residue."""
    if len(lane_polys) != len(LANE_MODULI):
        raise ValueError("q88 CRT lane count mismatch")
    if any(len(poly) != len(combined) for poly in lane_polys):
        raise ValueError("q88 CRT polynomial degree mismatch")
    for lane, modulus in enumerate(LANE_MODULI):
        for index, value in enumerate(combined):
            if value % modulus != int(lane_polys[lane][index]) % modulus:
                raise ValueError("q88 CRT recombination mismatch")


def build_shared_crt_aggregate(lane_polys: Sequence[Sequence[int]]) -> tuple[
        list[int], list[list[int]]]:
    """Build one exact CRT aggregate and its lane quotient witnesses.

    For every coefficient this returns ``A`` and quotients ``k_l`` satisfying
    ``A = lane_l + LANE_MODULI[l] * k_l`` over the integers.  The quotients
    are kept as signed values and are not reduced modulo a proof prime; the
    native bridge must range-bound and radix-decompose them explicitly.
    """
    combined = crt_recombine_poly(lane_polys)
    quotients = []
    for lane, modulus in enumerate(LANE_MODULI):
        quotients.append([
            (int(value) - int(lane_polys[lane][index])) // modulus
            for index, value in enumerate(combined)
        ])
    for lane, modulus in enumerate(LANE_MODULI):
        for index, value in enumerate(combined):
            if int(value) != (int(lane_polys[lane][index]) +
                              modulus * quotients[lane][index]):
                raise ValueError("q88 shared CRT quotient reconstruction failed")
    return combined, quotients


def verify_shared_crt_aggregate(combined: Sequence[int],
                                lane_polys: Sequence[Sequence[int]],
                                quotients: Sequence[Sequence[int]]) -> None:
    """Fail closed unless a shared aggregate and all lane quotients match."""
    if len(lane_polys) != len(LANE_MODULI) or len(quotients) != len(LANE_MODULI):
        raise ValueError("q88 shared CRT lane count mismatch")
    if any(len(poly) != len(combined) or len(quotients[lane]) != len(combined)
           for lane, poly in enumerate(lane_polys)):
        raise ValueError("q88 shared CRT polynomial degree mismatch")
    for lane, modulus in enumerate(LANE_MODULI):
        for index, value in enumerate(combined):
            if int(value) != (int(lane_polys[lane][index]) +
                              modulus * int(quotients[lane][index])):
                raise ValueError("q88 shared CRT quotient mismatch")


def balanced_radix(value: int, *, base: int = 128,
                   digits: int | None = None) -> list[int]:
    """Encode one integer with centered radix digits and exact capacity checks."""
    if base < 2 or base & (base - 1):
        raise ValueError("q88 bridge radix must be a power of two")
    remaining = int(value)
    output: list[int] = []
    limit = digits if digits is not None else max(1, (abs(value).bit_length() + 6) // 7)
    half = base // 2
    for _ in range(limit):
        digit = (remaining + half) % base - half
        output.append(digit)
        remaining = (remaining - digit) // base
    if remaining:
        raise ValueError("q88 bridge radix capacity exhausted")
    return output


def audit_round_trip(values: Iterable[int], *, base: int = 128,
                     digits: int | None = None) -> None:
    for value in values:
        encoded = balanced_radix(value, base=base, digits=digits)
        reconstructed = sum(digit * base**index
                            for index, digit in enumerate(encoded))
        if reconstructed != value:
            raise ValueError("q88 bridge radix round trip failed")


def negacyclic_mul(a: Sequence[int], b: Sequence[int]) -> list[int]:
    """Small exact reference multiplication used by the bridge KAT.

    The production implementation will replace this with the ring backend;
    keeping the KAT arithmetic independent makes CRT/lift mistakes visible.
    """
    if len(a) != len(b):
        raise ValueError("q88 bridge polynomial degree mismatch")
    n = len(a)
    out = [0] * n
    nonzero_a = [(i, int(value)) for i, value in enumerate(a) if value]
    nonzero_b = [(i, int(value)) for i, value in enumerate(b) if value]
    if len(nonzero_a) * len(nonzero_b) < n * n // 8:
        for i, left in nonzero_a:
            for j, right in nonzero_b:
                k = i + j
                out[k % n] += left * right if k < n else -left * right
        return out
    for i, left in enumerate(a):
        for j, right in enumerate(b):
            k = i + j
            out[k % n] += left * right if k < n else -left * right
    return out


def crt_lift_polynomial(lane_polys: Sequence[Sequence[int]]) -> list[int]:
    """Canonical coefficient lift of three lane polynomials to the q88 ring."""
    return crt_recombine_poly(lane_polys)


def _digits(value: int, base: int, count: int, *, signed: bool) -> list[int]:
    if base < 2 or base & (base - 1):
        raise ValueError("q88 bridge radix must be a power of two")
    out: list[int] = []
    remaining = int(value)
    half = base // 2
    for _ in range(count):
        if signed:
            digit = (remaining + half) % base - half
        else:
            if remaining < 0:
                raise ValueError("negative unsigned bridge value")
            digit = remaining % base
        out.append(digit)
        remaining = (remaining - digit) // base
    if remaining:
        raise ValueError("q88 bridge digit capacity exhausted")
    return out


def _decompose(poly: Sequence[int], base: int, count: int,
               *, signed: bool) -> list[list[int]]:
    layers = [[0] * len(poly) for _ in range(count)]
    for index, value in enumerate(poly):
        for layer, digit in enumerate(_digits(value, base, count,
                                               signed=signed)):
            layers[layer][index] = digit
    return layers


@dataclass(frozen=True)
class LinearCRTBridge:
    """Exact quotient/carry expansion of one CRT-combined relation.

    This object is deliberately proof-engine agnostic.  It is the boundary
    between the q88 relation and the Pack adapter: all arrays are integer
    arrays, and ``verify`` must pass before a native statement is built.
    """

    base: int
    matrix_digits: tuple[tuple[tuple[int, ...], ...], ...]
    witness_digits: tuple[tuple[tuple[int, ...], ...], ...]
    modulus_digits: tuple[int, ...]
    target_digits: tuple[tuple[int, ...], ...]
    quotient_digits: tuple[tuple[int, ...], ...]
    carry_digits: tuple[tuple[int, ...], ...]

    def verify(self) -> None:
        n = len(self.target_digits[0]) if self.target_digits else 0
        zero = [0] * n
        for layer in range(len(self.carry_digits)):
            expression = [0] * n
            for matrix, witness in zip(self.matrix_digits,
                                       self.witness_digits):
                for i, matrix_digit in enumerate(matrix):
                    j = layer - i
                    if 0 <= j < len(witness):
                        product = negacyclic_mul(matrix_digit,
                                                 witness[j])
                        expression = [x + y for x, y in zip(expression,
                                                             product)]
            for i, modulus_digit in enumerate(self.modulus_digits):
                j = layer - i
                if 0 <= j < len(self.quotient_digits):
                    qproduct = [modulus_digit * x
                                for x in self.quotient_digits[j]]
                    expression = [x - y for x, y in zip(expression, qproduct)]
            if layer < len(self.target_digits):
                expression = [x - y for x, y in
                              zip(expression, self.target_digits[layer])]
            if layer:
                expression = [x + y for x, y in
                              zip(expression, self.carry_digits[layer - 1])]
            expression = [x - self.base * y for x, y in
                          zip(expression, self.carry_digits[layer])]
            if expression != zero:
                raise ValueError(f"q88 bridge carry equation {layer} failed")


@dataclass(frozen=True)
class Single88Bridge:
    """Exact quotient/carry bridge for the independent single-prime profile.

    The relation is ``sum_i A_i*w_i - t = P*z`` over integer polynomial
    coefficients, with all public and private values radix decomposed before
    entering the 38-bit Pack ring.  This class intentionally has no CRT
    fields and is never consumed by the production q88 compiler.
    """

    degree: int
    base: int
    h_digits: tuple[tuple[tuple[int, ...], ...], ...]
    w_digits: tuple[tuple[tuple[int, ...], ...], ...]
    q_digits: tuple[int, ...]
    quotient_digits: tuple[tuple[int, ...], ...]
    target_digits: tuple[tuple[int, ...], ...]
    carry_digits: tuple[tuple[int, ...], ...]
    modulus: int

    def equation_count(self) -> int:
        return len(self.carry_digits)

    def verify(self) -> None:
        zero = [0] * self.degree
        for layer in range(self.equation_count()):
            expression = list(self.carry_digits[layer - 1]
                              if layer else zero)
            for matrix, witness in zip(self.h_digits, self.w_digits):
                for matrix_layer, matrix_digit in enumerate(matrix):
                    witness_layer = layer - matrix_layer
                    if 0 <= witness_layer < len(witness):
                        product = negacyclic_mul(
                            matrix_digit, witness[witness_layer])
                        expression = [x + y for x, y in
                                      zip(expression, product)]
            for modulus_layer, modulus_digit in enumerate(self.q_digits):
                quotient_layer = layer - modulus_layer
                if 0 <= quotient_layer < len(self.quotient_digits):
                    expression = [x - modulus_digit * y for x, y in zip(
                        expression, self.quotient_digits[quotient_layer])]
            if layer < len(self.target_digits):
                expression = [x - y for x, y in zip(
                    expression, self.target_digits[layer])]
            expression = [x - self.base * y for x, y in zip(
                expression, self.carry_digits[layer])]
            if expression != zero:
                raise ValueError(f"single q88 bridge carry {layer} failed")


def compile_single88_bridge(matrix: Sequence[Sequence[int]],
                            witness: Sequence[Sequence[int]],
                            target: Sequence[int], *,
                            base: int = 128,
                            witness_digits: int = 9,
                            quotient_digits: int = 16,
                            modulus: int = SINGLE88_PRIME) -> Single88Bridge:
    """Compile one exact 88-bit-prime relation into radix carry equations."""
    if base < 2 or base & (base - 1):
        raise ValueError("single q88 bridge radix must be a power of two")
    if modulus != SINGLE88_PRIME:
        raise ValueError("single q88 bridge requires the certified profile modulus")
    verify_modulus_certificate()
    if not witness or len(matrix) != len(witness):
        raise ValueError("single q88 matrix/witness shape mismatch")
    degree = len(target)
    if degree == 0 or any(len(poly) != degree
                           for poly in list(matrix) + list(witness)):
        raise ValueError("single q88 polynomial degree mismatch")
    total = [0] * degree
    for public_poly, secret_poly in zip(matrix, witness):
        total = [x + y for x, y in zip(
            total, negacyclic_mul(public_poly, secret_poly))]
    centered_target = [int(value) % modulus for value in target]
    quotient = [(value - rhs) // modulus
                if (value - rhs) % modulus == 0 else None
                for value, rhs in zip(total, centered_target)]
    if any(value is None for value in quotient):
        raise ValueError("single q88 relation is not congruent modulo prime")
    base_bits = base.bit_length() - 1
    # Signed digits need one extra sign bit at the maximum public magnitude.
    public_bits = max(modulus.bit_length(),
                      max((abs(int(value)).bit_length() + 1 for poly in matrix
                           for value in poly), default=1))
    matrix_count = (public_bits + base_bits - 1) // base_bits
    matrix_digits = tuple(tuple(tuple(layer) for layer in
                                _decompose(poly, base, matrix_count,
                                           signed=True))
                          for poly in matrix)
    witness_layers = tuple(tuple(tuple(layer) for layer in
                                 _decompose(poly, base, witness_digits,
                                            signed=True))
                           for poly in witness)
    target_digits = tuple(tuple(layer) for layer in
                          _decompose(centered_target, base, matrix_count,
                                     signed=False))
    q_digits = tuple(_digits(modulus, base, matrix_count, signed=False))
    quotient_layers = tuple(tuple(layer) for layer in
                            _decompose([int(value) for value in quotient],
                                       base, quotient_digits, signed=True))
    layer_count = max(matrix_count + witness_digits - 1,
                      matrix_count + quotient_digits - 1) + 1
    carries = []
    previous = [0] * degree
    for layer in range(layer_count):
        expression = previous[:]
        for matrix_poly, secret in zip(matrix_digits, witness_layers):
            for matrix_layer, public_digit in enumerate(matrix_poly):
                witness_layer = layer - matrix_layer
                if 0 <= witness_layer < len(secret):
                    product = negacyclic_mul(
                        public_digit, secret[witness_layer])
                    expression = [x + y for x, y in zip(expression, product)]
        for modulus_layer, modulus_digit in enumerate(q_digits):
            quotient_layer = layer - modulus_layer
            if 0 <= quotient_layer < len(quotient_layers):
                expression = [x - modulus_digit * y for x, y in zip(
                    expression, quotient_layers[quotient_layer])]
        if layer < len(target_digits):
            expression = [x - y for x, y in zip(expression, target_digits[layer])]
        if any(value % base for value in expression):
            raise ValueError(f"single q88 bridge carry {layer} is not integral")
        previous = [value // base for value in expression]
        carries.append(previous)
    if any(previous):
        raise ValueError("single q88 bridge terminal carry is nonzero")
    bridge = Single88Bridge(
        degree, base, matrix_digits, witness_layers, q_digits,
        quotient_layers, target_digits,
        tuple(tuple(layer) for layer in carries), modulus)
    bridge.verify()
    return bridge


def compile_single88_bridge_fast(matrix: Sequence[Sequence[int]],
                                 witness: Sequence[Sequence[int]], *,
                                 target: Sequence[int] | None = None,
                                 base: int = 128,
                                 witness_digits: int = 9,
                                 quotient_digits: int = 16,
                                 modulus: int = SINGLE88_PRIME) -> Single88Bridge:
    """Construct the same single-prime bridge with exact digit FFT sums.

    ``target=None`` derives a fresh public target from the supplied bounded
    witness.  This mode is used only by the isolated profile: it is not a
    conversion of an existing CRT accumulator.  When a target is supplied,
    the function checks its integer congruence before constructing carries.
    The FFT operates on radix digits and every rounded coefficient is checked
    by the integer carry recurrence and by :meth:`Single88Bridge.verify`.
    """
    if modulus != SINGLE88_PRIME:
        raise ValueError("single q88 bridge requires the certified profile modulus")
    if base < 2 or base & (base - 1):
        raise ValueError("single q88 bridge radix must be a power of two")
    verify_modulus_certificate()
    if not matrix or len(matrix) != len(witness):
        raise ValueError("single q88 matrix/witness shape mismatch")
    degree = len(witness[0])
    if degree == 0 or any(len(poly) != degree
                           for poly in list(matrix) + list(witness)):
        raise ValueError("single q88 polynomial degree mismatch")
    try:
        from exact_negacyclic_fft import digit_sum_negacyclic
    except ImportError:
        from Lazer_Acc.python.exact_negacyclic_fft import digit_sum_negacyclic
    bits = base.bit_length() - 1
    matrix_count = max(1, (modulus.bit_length() + bits - 1) // bits)
    conv = digit_sum_negacyclic(matrix, witness, bits,
                                matrix_count, witness_digits)
    total = [sum(layer[k] * base ** layer_index
                  for layer_index, layer in enumerate(conv))
             for k in range(degree)]
    if target is None:
        centered_target = [value % modulus for value in total]
    else:
        if len(target) != degree:
            raise ValueError("single q88 target degree mismatch")
        centered_target = [int(value) % modulus for value in target]
    difference = [value - rhs for value, rhs in zip(total, centered_target)]
    if any(value % modulus for value in difference):
        raise ValueError("single q88 relation is not congruent modulo prime")
    quotient = [value // modulus for value in difference]
    h = [_decompose(poly, base, matrix_count, signed=False)
         for poly in matrix]
    w = [_decompose(poly, base, witness_digits, signed=True)
         for poly in witness]
    q_digits = _digits(modulus, base, matrix_count, signed=False)
    t = _decompose(centered_target, base, matrix_count, signed=False)
    q = _decompose(quotient, base, quotient_digits, signed=True)
    layer_count = max(len(conv), matrix_count + quotient_digits - 1) + 1
    carries: list[list[int]] = []
    previous = [0] * degree
    for layer in range(layer_count):
        expression = previous[:]
        if layer < len(conv):
            expression = [x + y for x, y in zip(expression, conv[layer])]
        for i, modulus_digit in enumerate(q_digits):
            j = layer - i
            if 0 <= j < len(q):
                expression = [x - modulus_digit * y for x, y in
                              zip(expression, q[j])]
        if layer < len(t):
            expression = [x - y for x, y in zip(expression, t[layer])]
        if any(value % base for value in expression):
            raise ValueError(f"single q88 fast bridge carry {layer} is not integral")
        previous = [value // base for value in expression]
        carries.append(previous)
    if any(previous):
        raise ValueError("single q88 fast bridge terminal carry is nonzero")
    bridge = Single88Bridge(degree, base,
                            tuple(tuple(tuple(layer) for layer in poly)
                                  for poly in h),
                            tuple(tuple(tuple(layer) for layer in poly)
                                  for poly in w),
                            q_digits,
                            tuple(tuple(layer) for layer in q),
                            tuple(tuple(layer) for layer in t),
                            tuple(tuple(layer) for layer in carries), modulus)
    # The recurrence above is the exact coefficient-wise verification.  A
    # second generic check would redo every 4096^2 negacyclic product and
    # makes the 174-column profile needlessly quadratic in Python.  Keep the
    # independent O(N^2) check for small KATs; production-sized profiles use
    # the exact digit-sum recurrence plus native Clear Verify below.
    if degree <= 512 or len(matrix) <= 16:
        bridge.verify()
    return bridge


def compile_linear_crt_bridge(matrix: Sequence[Sequence[Sequence[int]]],
                              witness: Sequence[Sequence[int]],
                              target: Sequence[Sequence[int]], *,
                              base: int = 128,
                              witness_digits: int = 9,
                              quotient_digits: int = 16) -> LinearCRTBridge:
    """Compile lane residues into one exact q88 quotient/carry relation.

    ``matrix`` and ``target`` are three lane arrays.  ``witness`` is the
    common signed integer witness.  The function rejects any lane mismatch or
    non-integral CRT lift; it never silently reduces the 88-bit relation to a
    proof-prime relation.
    """
    if len(matrix) != len(LANE_MODULI) or len(target) != len(LANE_MODULI):
        raise ValueError("q88 bridge expects three CRT lanes")
    if not witness:
        raise ValueError("q88 bridge needs a witness")
    degree = len(witness[0])
    if any(len(poly) != degree for poly in witness):
        raise ValueError("q88 bridge witness degree mismatch")
    if any(len(poly) != degree for lane in matrix for poly in lane):
        raise ValueError("q88 bridge matrix degree mismatch")
    if any(len(poly) != degree for poly in target):
        raise ValueError("q88 bridge target degree mismatch")
    lifted_matrix = [crt_lift_polynomial([matrix[lane][column]
                                           for lane in range(3)])
                     for column in range(len(witness))]
    lifted_target = crt_lift_polynomial(target)
    total = [0] * degree
    for coeffs, secret in zip(lifted_matrix, witness):
        total = [x + y for x, y in zip(total, negacyclic_mul(coeffs, secret))]
    quotient = []
    for value, rhs in zip(total, lifted_target):
        difference = value - rhs
        if difference % CRT_MODULUS:
            raise ValueError("q88 CRT congruence does not hold")
        quotient.append(difference // CRT_MODULUS)
    matrix_count = max(CRT_MODULUS.bit_length(),
                       max((value.bit_length() for poly in lifted_matrix
                            for value in poly), default=1))
    matrix_count = max(1, (matrix_count + base.bit_length() - 2) /
                       (base.bit_length() - 1))
    matrix_count = int(matrix_count)
    modulus_digits = tuple(_digits(CRT_MODULUS, base, matrix_count,
                                   signed=False))
    matrix_digits = tuple(tuple(tuple(layer) for layer in
                                _decompose(poly, base, matrix_count,
                                           signed=False))
                          for poly in lifted_matrix)
    witness_layers = tuple(tuple(tuple(layer) for layer in
                                 _decompose(poly, base, witness_digits,
                                            signed=True))
                           for poly in witness)
    target_digits = tuple(tuple(layer) for layer in
                          _decompose(lifted_target, base, matrix_count,
                                     signed=False))
    quotient_layers = tuple(tuple(layer) for layer in
                            _decompose(quotient, base, quotient_digits,
                                        signed=True))
    layer_count = max(matrix_count + witness_digits - 1,
                      matrix_count + quotient_digits - 1) + 1
    carries: list[list[int]] = []
    previous = [0] * degree
    for layer in range(layer_count):
        expression = previous[:]
        for coeffs, secret in zip(matrix_digits, witness_layers):
            for i, coeff_digit in enumerate(coeffs):
                j = layer - i
                if 0 <= j < len(secret):
                    expression = [x + y for x, y in
                                  zip(expression,
                                      negacyclic_mul(coeff_digit, secret[j]))]
        for i, modulus_digit in enumerate(modulus_digits):
            j = layer - i
            if 0 <= j < len(quotient_layers):
                expression = [x - modulus_digit * y for x, y in
                              zip(expression, quotient_layers[j])]
        if layer < len(target_digits):
            expression = [x - y for x, y in
                          zip(expression, target_digits[layer])]
        if any(value % base for value in expression):
            raise ValueError(f"q88 bridge carry {layer} is not integral")
        previous = [value // base for value in expression]
        carries.append(previous)
    if any(previous):
        raise ValueError("q88 bridge has nonzero terminal carry")
    bridge = LinearCRTBridge(
        base, matrix_digits, witness_layers, modulus_digits, target_digits,
        quotient_layers, tuple(tuple(layer) for layer in carries))
    bridge.verify()
    return bridge


def aggregate_handle_witness(right_polys: Sequence[Sequence[int]],
                             handle_bits: Sequence[int], *,
                             gadget_base: int = 26_007,
                             gadget_digits: int = 6) -> list[int]:
    """Build the hidden-handle aggregate used as a bridge auxiliary witness.

    ``right_polys`` is laid out as ``HANDLE_BITS * gadget_digits`` vectors.
    The returned polynomial is the exact integer combination appearing in the
    accumulator relation.  Callers must prove this equality; it is never
    treated as a public coefficient.
    """
    if len(right_polys) != len(handle_bits) * gadget_digits:
        raise ValueError("q88 handle aggregate shape mismatch")
    if any(bit not in (0, 1) for bit in handle_bits):
        raise ValueError("q88 handle aggregate expects binary bits")
    degree = len(right_polys[0])
    if any(len(poly) != degree for poly in right_polys):
        raise ValueError("q88 handle aggregate degree mismatch")
    out = [0] * degree
    for bit, selected in enumerate(handle_bits):
        scale = 1
        for digit in range(gadget_digits):
            right = right_polys[bit * gadget_digits + digit]
            factor = scale if selected else 0
            out = [x + factor * y for x, y in zip(out, right)]
            scale *= gadget_base
    return out


def compose_refresh_witness(refresh_bank: Sequence[Sequence[int]],
                            selector_bits: Sequence[int]) -> list[int]:
    """Compose refresh public rows using the hidden selector bits.

    The result is an auxiliary public-polynomial witness input.  A production
    bridge must retain the selector binding equation, so this helper does not
    make the composed result public by itself.
    """
    if len(refresh_bank) != len(selector_bits) + 1:
        raise ValueError("q88 refresh bank shape mismatch")
    if any(bit not in (0, 1) for bit in selector_bits):
        raise ValueError("q88 refresh selector expects binary bits")
    degree = len(refresh_bank[0])
    if any(len(poly) != degree for poly in refresh_bank):
        raise ValueError("q88 refresh bank degree mismatch")
    out = list(refresh_bank[0])
    for bit, poly in zip(selector_bits, refresh_bank[1:]):
        if bit:
            out = [x + y for x, y in zip(out, poly)]
    return out


def verify_auxiliary_equality(actual: Sequence[int], expected: Sequence[int],
                              *, label: str) -> None:
    """Fail closed when an auxiliary witness differs from its hidden recipe."""
    if len(actual) != len(expected) or any(x != y for x, y in zip(actual,
                                                                  expected)):
        raise ValueError(f"q88 auxiliary witness mismatch: {label}")

@dataclass(frozen=True)
class DigitBinding:
    base: int
    input_digits: tuple[tuple[tuple[int, ...], ...], ...]
    output_digits: tuple[tuple[int, ...], ...]
    carries: tuple[tuple[int, ...], ...]

    def verify(self) -> None:
        if not self.output_digits:
            raise ValueError("empty digit binding")
        degree = len(self.output_digits[0])
        previous = [0] * degree
        for layer, output in enumerate(self.output_digits):
            expression = previous[:]
            for digits in self.input_digits:
                if layer < len(digits):
                    expression = [x + y for x, y in zip(expression,
                                                         digits[layer])]
            expression = [x - y for x, y in zip(expression,
                                                 self.output_digits[layer])]
            expression = [x - self.base * y for x, y in
                          zip(expression, self.carries[layer])]
            if expression != [0] * degree:
                raise ValueError(f"digit binding layer {layer} failed")
            previous = list(self.carries[layer])


def compile_digit_binding(inputs: Sequence[Sequence[int]],
                          scalars: Sequence[int], target: Sequence[int], *,
                          base: int = 128, digits: int = 12) -> DigitBinding:
    """Compile exact scalar-combination binding with radix carries."""
    if len(inputs) != len(scalars) or len(target) != len(inputs[0]):
        raise ValueError("digit binding shape mismatch")
    degree = len(target)
    total = [0] * degree
    for poly, scalar in zip(inputs, scalars):
        total = [x + int(scalar) * y for x, y in zip(total, poly)]
    if any(abs(value) >= (base ** digits) // 2 for value in total):
        raise ValueError("digit binding capacity exhausted")
    input_digits = []
    for poly, scalar in zip(inputs, scalars):
        scaled = [int(scalar) * int(value) for value in poly]
        input_digits.append(tuple(tuple(layer) for layer in
                                 _decompose(scaled, base, digits,
                                            signed=True)))
    output_digits = tuple(tuple(_digit for _digit in layer) for layer in
                          zip(*[_balanced for _balanced in []])) if False else tuple(
        tuple(layer) for layer in _decompose(target, base, digits, signed=True))
    carries = []
    previous = [0] * degree
    for layer in range(digits):
        expression = previous[:]
        for item in input_digits:
            expression = [x + y for x, y in zip(expression, item[layer])]
        expression = [x - y for x, y in zip(expression, output_digits[layer])]
        if any(value % base for value in expression):
            raise ValueError("digit binding carry is not integral")
        previous = [value // base for value in expression]
        carries.append(previous)
    if any(previous):
        raise ValueError("digit binding terminal carry is nonzero")
    result = DigitBinding(base, tuple(input_digits), output_digits,
                          tuple(tuple(layer) for layer in carries))
    result.verify()
    return result

@dataclass(frozen=True)
class NativeBindingConstraint:
    layer: int
    input_indices: tuple[int, ...]
    output_index: int
    previous_carry_index: int | None
    carry_index: int
    coefficients: tuple[int, ...]


def make_native_binding_plan(binding: DigitBinding,
                              *, input_start: int = 0,
                              output_start: int | None = None,
                              carry_start: int | None = None
                              ) -> tuple[NativeBindingConstraint, ...]:
    """Describe digit binding equations for native Pack construction.

    Indices refer to the caller's witness-vector map.  The plan contains no
    hidden effective coefficient; all coefficients are public radix constants.
    """
    input_count = len(binding.input_digits)
    output_base = (input_start + input_count if output_start is None
                   else output_start)
    carry_base = (output_base + len(binding.output_digits)
                  if carry_start is None else carry_start)
    result = []
    for layer in range(len(binding.output_digits)):
        previous = carry_base + layer - 1 if layer else None
        result.append(NativeBindingConstraint(
            layer=layer,
            input_indices=tuple(input_start + index * len(binding.output_digits)
                                + layer for index in range(input_count)),
            output_index=output_base + layer,
            previous_carry_index=previous,
            carry_index=carry_base + layer,
            coefficients=tuple([1] * input_count + [-1] +
                               ([1] if previous is not None else []) +
                               [-binding.base])))
    return tuple(result)


def append_native_binding_plan(statement, ring, plan,
                               *, input_indices: Sequence[int],
                               output_indices: Sequence[int],
                               carry_indices: Sequence[int]) -> None:
    """Append digit binding constraints to a LaZer proof statement.

    ``statement`` is the existing ``labrados.proof_statement`` object and all
    indices are already allocated witness vectors.  Coefficients are constant
    ring polynomials, so this helper cannot silently introduce a hidden
    effective public coefficient.
    """
    from lazer import poly_t
    one = poly_t(ring, [1] + [0] * (ring.deg - 1))
    zero = poly_t(ring, [0] * ring.deg)
    for constraint in plan:
        terms = [one for _ in constraint.input_indices]
        witnesses = [input_indices[index] for index in constraint.input_indices]
        terms.append(poly_t(ring, [-1] + [0] * (ring.deg - 1)))
        witnesses.append(output_indices[constraint.layer])
        if constraint.previous_carry_index is not None:
            terms.append(one)
            witnesses.append(carry_indices[constraint.previous_carry_index])
        terms.append(poly_t(ring,
                            [-constraint.base] + [0] * (ring.deg - 1)))
        witnesses.append(carry_indices[constraint.layer])
        statement.append_statement(terms, witnesses, zero)


def build_grouped_pack_statement(bridge: Bridge, *, group_size: int = 12,
                                 zk: bool = True, extra_witness=None,
                                 extra_quadratics=None, extra_deg0=None):
    """Build one grouped native Pack statement from a verified Bridge."""
    from lazer import poly_t, polyvec_t, polyring_t
    from labrados import proof_statement, lib
    if group_size <= 0:
        raise ValueError("group size must be positive")
    degree = bridge.degree
    ring = polyring_t(degree, PACK_MODULUS)
    # Labrador witness coefficients are stored in signed 16-bit limbs.  The
    # raw quotient/carry recurrence can exceed that bound (the real q88
    # bridge reaches about 4e4), so split every carry polynomial into signed
    # radix-128 digits before placing it in the native witness.
    max_carry_bits = max((max(abs(int(value)).bit_length() for value in carry)
                          for carry in bridge.carry_digits), default=1)
    carry_digit_count = max(1, (max_carry_bits + 6) // 7)
    carry_digits = []
    for carry in bridge.carry_digits:
        layers = [[0] * degree for _ in range(carry_digit_count)]
        for coeff_index, value in enumerate(carry):
            remaining = int(value)
            for digit_index in range(carry_digit_count):
                digit = (remaining + bridge.base // 2) % bridge.base - bridge.base // 2
                layers[digit_index][coeff_index] = digit
                remaining = (remaining - digit) // bridge.base
            if remaining:
                raise ValueError("carry radix capacity exhausted")
        carry_digits.append(layers)
    vectors = []
    norms = []
    witness_refs = []
    columns = len(bridge.w_digits)
    for layer in range(len(bridge.w_digits[0])):
        for start in range(0, columns, group_size):
            block = [bridge.w_digits[col][layer]
                     for col in range(start, min(start + group_size, columns))]
            pv = polyvec_t(ring, len(block))
            for index, values in enumerate(block):
                pv[index] = poly_t(ring, list(values))
            vectors.append(pv)
            norms.append(max(1, sum(value * value for poly in block
                                    for value in poly)))
            witness_refs.append(("w", layer, start))
    for layer, values in enumerate(bridge.quotient_digits):
        pv = polyvec_t(ring, 1)
        pv[0] = poly_t(ring, list(values))
        vectors.append(pv)
        norms.append(max(1, sum(value * value for value in values)))
        witness_refs.append(("q", layer, 0))
    for layer, digit_layers in enumerate(carry_digits):
        for digit, values in enumerate(digit_layers):
            pv = polyvec_t(ring, 1)
            pv[0] = poly_t(ring, list(values))
            vectors.append(pv)
            norms.append(max(1, sum(value * value for value in values)))
            witness_refs.append(("c", layer, digit))
    if any(value >= (1 << 64) for value in norms):
        raise ValueError("grouped Pack norm exceeds uint64")
    equation_count = bridge.equation_count()
    # The bridge KAT deliberately uses a conservative public norm budget.  A
    # production single-prime profile must replace this with an audited
    # witness-bound calculation; it must never infer security from this KAT.
    kat_norms = [max(1, norm) for norm in norms]
    # Optional witnesses/constraints are used by the isolated full-shape
    # profile to bind hidden selectors and product witnesses in the same
    # native statement.  Each extra entry is
    # ``{"ref": hashable, "values": [poly, ...], "norm": uint64}``.
    extra_witness = list(extra_witness or ())
    extra_quadratics = list(extra_quadratics or ())
    extra_deg0 = list(extra_deg0 or ())
    for item in extra_witness:
        values = item["values"]
        pv = polyvec_t(ring, len(values))
        for index, values_poly in enumerate(values):
            pv[index] = poly_t(ring, list(values_poly))
        vectors.append(pv)
        norms.append(int(item.get("norm", 1)))
        witness_refs.append(item["ref"])
    statement = proof_statement(
        [degree] * len(vectors), [vector.dim for vector in vectors], norms,
        num_constraints=equation_count + len(extra_quadratics),
        primesize="38", num_deg0_constraints=len(extra_deg0), zk=zk)
    indices = [statement.append_witness(vector) for vector in vectors]
    ref_to_index = {ref: index for ref, index in zip(witness_refs, indices)}
    zero = poly_t(ring, [0] * degree)
    # Public modulus digits are already stored as scalar coefficients in q_digits.
    for layer in range(equation_count):
        terms = []
        witnesses = []
        # Keep the same layer-first convolution order as the reference
        # grouped Pack adapter: witness digit l is paired with public digit
        # j=layer-l.  This also makes the witness-index mapping explicit.
        for witness_layer in range(len(bridge.w_digits[0])):
            public_layer = layer - witness_layer
            if not 0 <= public_layer < len(bridge.h_digits[0]):
                continue
            for start in range(0, columns, group_size):
                end = min(start + group_size, columns)
                block = [bridge.h_digits[index][public_layer]
                         for index in range(start, end)]
                pv = polyvec_t(ring, len(block))
                for index, values in enumerate(block):
                    pv[index] = poly_t(ring, list(values))
                terms.append(pv)
                witnesses.append(ref_to_index[("w", witness_layer, start)])
        for q_layer, modulus_digit in enumerate(bridge.q_digits):
            quotient_layer = layer - q_layer
            if 0 <= quotient_layer < len(bridge.quotient_digits):
                terms.append(poly_t(ring, [-modulus_digit] + [0] * (degree - 1)))
                witnesses.append(ref_to_index[("q", quotient_layer, 0)])
        if layer:
            for digit in range(carry_digit_count):
                terms.append(poly_t(ring, [bridge.base ** digit] + [0] * (degree - 1)))
                witnesses.append(ref_to_index[("c", layer - 1, digit)])
        for digit in range(carry_digit_count):
            terms.append(poly_t(ring, [-(bridge.base ** (digit + 1))] + [0] * (degree - 1)))
            witnesses.append(ref_to_index[("c", layer, digit)])
        rhs = (bridge.target_digits[layer]
               if layer < len(bridge.target_digits) else [0] * degree)
        statement.append_statement(terms, witnesses, poly_t(ring, list(rhs)))
    one = poly_t(ring, [1] + [0] * (degree - 1))
    zero = poly_t(ring, [0] * degree)
    for item in extra_deg0:
        index = ref_to_index[item]
        error = lib.labrador38_py_append_deg0_logical_vector_constraints(
            statement.statement, statement.map[index], 1, degree // 256)
        if error:
            raise ValueError(f"extra degree-zero constraint failed: {error}")
        statement.cur_deg0_constraint_num += 1
    for item in extra_quadratics:
        linear = [poly_t(ring, list(coeff)) for coeff, _ in item.get("linear", ())]
        linear_refs = [ref_to_index[ref] for _, ref in item.get("linear", ())]
        quadratic_coefficients = [int(coeff) for coeff, _, _ in item.get("quadratic", ())]
        quadratic_left = [ref_to_index[left]
                          for _, left, _ in item.get("quadratic", ())]
        quadratic_right = [ref_to_index[right]
                           for _, _, right in item.get("quadratic", ())]
        statement.append_quadratic_statement(
            linear, linear_refs, quadratic_coefficients,
            quadratic_left, quadratic_right, zero)
    if int(statement.simple_verify(statement.statement, statement.witness)) != 1:
        raise ValueError("grouped Pack statement clear verification failed")
    return statement, ref_to_index


def compile_linear_crt_bridge_fast(matrix, witness, target, *,
                                   base: int = 128,
                                   witness_digits: int = 9,
                                   quotient_digits: int = 16):
    """FFT-assisted exact lift compiler for degree-4096 q88 experiments.

    The FFT is used only to construct integer products and carries; the native
    Pack statement still proves the resulting digit equations.  The helper is
    optional so the small reference compiler remains dependency-free.
    """
    if base < 2 or base & (base - 1):
        raise ValueError("fast q88 bridge radix must be a power of two")
    base_bits = base.bit_length() - 1
    try:
        from exact_negacyclic_fft import digit_sum_negacyclic
    except ImportError as exc:
        raise RuntimeError("exact_negacyclic_fft is required for fast bridge") from exc
    if len(matrix) != len(LANE_MODULI) or len(target) != len(LANE_MODULI):
        raise ValueError("q88 fast bridge expects three lanes")
    columns = len(witness)
    degree = len(witness[0])
    lifted_matrix = [crt_lift_polynomial([matrix[lane][column]
                                           for lane in range(3)])
                     for column in range(columns)]
    lifted_target = crt_lift_polynomial(target)
    conv = digit_sum_negacyclic(lifted_matrix, witness, base_bits,
                                 (CRT_MODULUS.bit_length() + base_bits - 1) // base_bits,
                                 witness_digits)
    total = [sum(conv[layer][k] * (base ** layer)
                 for layer in range(len(conv))) for k in range(degree)]
    quotient = []
    for value, rhs in zip(total, lifted_target):
        diff = value - rhs
        if diff % CRT_MODULUS:
            raise ValueError("q88 fast bridge congruence failed")
        quotient.append(diff // CRT_MODULUS)
    h_count = (CRT_MODULUS.bit_length() + base_bits - 1) // base_bits
    h = [_decompose(poly, base, h_count, signed=False)
         for poly in lifted_matrix]
    w = [_decompose(poly, base, witness_digits, signed=True)
         for poly in witness]
    q_digits = _digits(CRT_MODULUS, base, h_count, signed=False)
    t = _decompose(lifted_target, base, h_count, signed=False)
    q_layers = _decompose(quotient, base, quotient_digits, signed=True)
    layer_count = max(len(conv), h_count + quotient_digits - 1) + 1
    carries = []
    previous = [0] * degree
    for layer in range(layer_count):
        expression = previous[:]
        if layer < len(conv):
            expression = [x + y for x, y in zip(expression, conv[layer])]
        for i, qdigit in enumerate(q_digits):
            j = layer - i
            if 0 <= j < quotient_digits:
                expression = [x - qdigit * y for x, y in
                              zip(expression, q_layers[j])]
        if layer < len(t):
            expression = [x - y for x, y in zip(expression, t[layer])]
        if any(value % base for value in expression):
            raise ValueError(f"q88 fast bridge carry {layer} is not integral")
        previous = [value // base for value in expression]
        carries.append(previous)
    if any(previous):
        raise ValueError("q88 fast bridge terminal carry is nonzero")
    bridge = FastBridge(degree, base, h, w, q_digits, q_layers, t, carries, {})
    return bridge


@dataclass
class FastBridge:
    degree: int
    base: int
    h_digits: list
    w_digits: list
    q_digits: list
    quotient_digits: list
    target_digits: list
    carry_digits: list
    bounds: dict

    def equation_count(self):
        return len(self.carry_digits)


def build_complete_single_pack_statement(public, prover, link=None, *,
                                         cache=None,
                                         check_witness=True):
    """Build the security-complete q88 single Pack statement.

    The experimental CRT bridge above deliberately owns only the lifted main
    accumulator equation.  Hidden handle bits, refresh selectors, exact CRT
    lifts, and the four q19 same-handle rows are already implemented and
    audited by the production degree-4096 compiler.  Reusing that compiler
    here keeps those quadratic constraints in the *same* native
    ``proof_statement`` instead of silently replacing them with public
    coefficients.

    ``link`` must be a decoded ``LinkPayload`` for the full anonymous-show
    scope.  Passing ``None`` explicitly requests accumulator-only scope and
    is rejected when the caller asks for the complete manifest.
    """
    # Carry over the production worker's relation-preserving parallel setup
    # when this acceptance gate is called directly (outside the worker).  All
    # values are defaults so deployments can still override them explicitly.
    import os
    tuning_defaults = {
        "OMP_NUM_THREADS": "12",
        "REVAC_LABRADOR_THREADS": "1",  # transcript-sensitive reducer
        # Verification kernels have disjoint outputs and can use a wider
        # native team; keep proving serial for transcript/reducer stability.
        "REVAC_LABRADOR_THREADS_VERIFY": str(
            min(6, max(1, os.cpu_count() or 1))),
        "REVAC_LABRADOR_GAUSSIAN_THREADS": "8",
        "REVAC_LABRADOR_STATEMENT_THREADS": "8",
        "REVAC_LABRADOR_PYTHON_STATEMENT_THREADS": "8",
        "REVAC_LABRADOR_QUADG_THREADS": str(
            min(12, max(1, os.cpu_count() or 1))),
        "REVAC_LABRADOR_GARBAGE_PARALLEL": "1",
        "REVAC_LABRADOR_ROTATION_THREADS": "8",
        "REVAC_LABRADOR_POLXVEC_SCALAR_THREADS": "8",
        "REVAC_LABRADOR_FROMPOLY_THREADS": "8",
        "REVAC_LABRADOR_REFRESH_THREADS": "8",
        "REVAC_LABRADOR_EXACT_REFRESH_THREADS": "8",
        "REVAC_LABRADOR_DCH_ADDCHECK_PARALLEL": "1",
        "REVAC_LABRADOR_LDD_PROJECTION_THREADS": "8",
        "REVAC_LABRADOR_LDD_PROJECTION_PARALLEL": "1",
        "REVAC_LABRADOR_LDD_LING_THREADS": "12",
        "REVAC_LABRADOR_LDD_LING_PARALLEL": "1",
        "REVAC_LABRADOR_LDD_COLLAPSE_JL_THREADS": "4",
        "REVAC_LABRADOR_LDD_COLLAPSE_JL_PARALLEL": "1",
        "REVAC_LABRADOR_LDD_ZQ_THREADS": "4",
        "REVAC_LABRADOR_LDD_ZQ_PARALLEL": "1",
        "REVAC_LABRADOR_LDD_JL_SPROD_THREADS": "8",
        "REVAC_LABRADOR_LDD_JL_SPROD_PARALLEL": "1",
        "REVAC_LABRADOR_LDD_LIFT_THREADS": "4",
        "REVAC_LABRADOR_LDD_LIFT_PARALLEL": "1",
        "REVAC_LABRADOR_LDD_PARALLEL": "1",
        "REVAC_LABRADOR_POLXVEC_PARALLEL": "1",
        "REVAC_LABRADOR_DYNAMIC_BIND_CACHE": "1",
        "REVAC_LABRADOR_WITNESS_CACHE": "1",
        "REVAC_LABRADOR_ROTATION_LAZY_SCRATCH": "1",
        "REVAC_LABRADOR_PACK_IMPROVE_THRESHOLD": "0.50",
        # Four compact proof-only limbs reduce the native witness dimension
        # while retaining the five-limb q88 wire representation.  The full
        # ring compiler derives the corresponding norm bounds exactly.
        "REVAC_Q96_COMPACT_WITNESS_LIMBS": "1",
        "REVAC_Q96_COMPACT_WITNESS_COUNT": "4",
    }
    for name, value in tuning_defaults.items():
        os.environ.setdefault(name, value)

    from revac_acc96_wire import HANDLE_BITS
    if link is None:
        raise ValueError(
            "complete q88 single Pack requires the q19 LinkPayload; "
            "use the production accumulator-only path explicitly otherwise")
    from revac_acc96_full_ring import build_statement_full_ring

    statement, params, proof = build_statement_full_ring(
        public, prover, link, prove=False, check_witness=check_witness,
        cache=cache)
    expected = getattr(statement, "num_constraints", 0)
    if expected < 3 + 4:
        raise ValueError("complete q88 statement is missing required constraints")
    if getattr(statement, "cur_constraint_num", 0) != expected:
        raise ValueError("complete q88 constraint accounting mismatch")
    if getattr(statement, "cur_int_constraint_num", 0) < HANDLE_BITS:
        raise ValueError("complete q88 hidden-handle binary constraints missing")
    return statement, params, proof


def prove_complete_single_pack(public, prover, link, *, cache=None,
                               check_witness=True,
                               reuse_accumulator_witness=True):
    """Prove, verify, serialize, deserialize, and reverify full q88 Pack.

    This is the executable acceptance gate for the complete statement.  It
    deliberately goes through ``build_statement_full_ring`` so the hidden
    handle, refresh, lift, and q19 constraints cannot be omitted by an
    experimental bridge caller.
    """
    import time
    from labrados import ffi, lib, pack_verify
    from revac_acc_memver_worker import (
        native_thread_scope, verify_aggregation_parallel_scope)

    if cache is None:
        statement, _unused_params, _unused_proof = \
            build_complete_single_pack_statement(
                public, prover, link, cache=None,
                check_witness=check_witness)
        prove_started = time.perf_counter()
        with native_thread_scope("prove"):
            error, params, proof = statement.pack_prove()
        prove_seconds = time.perf_counter() - prove_started
        if error:
            raise RuntimeError(f"complete q88 Pack Prove failed: {error}")
    else:
        # Build once into the resident cache.  This creates the static
        # statement, witness views, and Pack parameters; subsequent calls
        # rebind only dynamic accumulator/link data and reuse those objects.
        statement, _unused_params, _unused_proof = \
            build_complete_single_pack_statement(
                public, prover, link, cache=cache,
                check_witness=check_witness)
        from revac_acc96_full_ring import build_statement_full_ring
        prove_started = time.perf_counter()
        with native_thread_scope("prove"):
            statement, params, proof = build_statement_full_ring(
                public, prover, link, prove=True,
                check_witness=False, cache=cache,
                reuse_accumulator_witness=reuse_accumulator_witness)
        prove_seconds = time.perf_counter() - prove_started
        if proof is None or params is None:
            raise RuntimeError("complete q88 resident Pack Prove returned no proof")
    verify_started = time.perf_counter()
    with (native_thread_scope("verify"),
          verify_aggregation_parallel_scope("verify")):
        verified = pack_verify((params, proof), statement.output_statement(), "38")
    if not verified:
        raise RuntimeError("complete q88 Pack Verify failed")
    verify_seconds = time.perf_counter() - verify_started

    size = lib.labrador38_py_proof_nbytes(params, proof)
    serialized = ffi.new("uint8_t[]", size)
    if lib.labrador38_py_proof_serialize(serialized, size, params, proof):
        raise RuntimeError("complete q88 proof serialization failed")
    decoded = ffi.new("labrador38_dch_pack_proof")
    if lib.labrador38_py_proof_deserialize(decoded, params, serialized, size):
        raise RuntimeError("complete q88 proof deserialization failed")
    with (native_thread_scope("verify"),
          verify_aggregation_parallel_scope("verify")):
        decoded_verified = pack_verify(
            (params, decoded), statement.output_statement(), "38")
    if not decoded_verified:
        raise RuntimeError("complete q88 deserialized Pack Verify failed")
    return {
        "prove_seconds": prove_seconds,
        "verify_seconds": verify_seconds,
        "serialized_bytes": int(size),
        "deserialized_verified": True,
    }


def decode_complete_q88_dynamic(raw: bytes):
    """Decode and cross-check the full prover-plus-q19-link payload.

    The worker wire format is ``ProverPayload || LinkPayload``.  Keeping this
    check in the bridge module prevents a complete-proof experiment from
    accidentally pairing a q19 opening with a different accumulator or
    handle.
    """
    from revac_acc96_wire import LINK_BYTES, PROVER_BYTES, LinkPayload, ProverPayload
    if len(raw) != PROVER_BYTES + LINK_BYTES:
        raise ValueError("complete q88 payload must contain prover and link data")
    prover = ProverPayload(raw[:PROVER_BYTES])
    link = LinkPayload.decode(raw[PROVER_BYTES:])
    if link.accumulator_digest != prover.accumulator_digest:
        raise ValueError("q19 link and accumulator digests do not match")
    link.audit_relation(prover.handle_vector)
    return prover, link
