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
    scale = 1
    for bit, selected in enumerate(handle_bits):
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
