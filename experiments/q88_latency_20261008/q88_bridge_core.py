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
                                 zk: bool = True):
    """Build one grouped native Pack statement from a verified Bridge."""
    from lazer import poly_t, polyvec_t, polyring_t
    from labrados import proof_statement
    if group_size <= 0:
        raise ValueError("group size must be positive")
    degree = bridge.degree
    ring = polyring_t(degree, PACK_MODULUS)
    # Labrador witness coefficients are stored in signed 16-bit limbs.  The
    # raw quotient/carry recurrence can exceed that bound (the real q88
    # bridge reaches about 4e4), so split every carry polynomial into signed
    # radix-128 digits before placing it in the native witness.
    carry_digit_count = 3
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
                pv[index] = poly_t(ring, values)
            vectors.append(pv)
            norms.append(max(1, sum(value * value for poly in block
                                    for value in poly)))
            witness_refs.append(("w", layer, start))
    for layer, values in enumerate(bridge.quotient_digits):
        pv = polyvec_t(ring, 1)
        pv[0] = poly_t(ring, values)
        vectors.append(pv)
        norms.append(max(1, sum(value * value for value in values)))
        witness_refs.append(("q", layer, 0))
    for layer, digit_layers in enumerate(carry_digits):
        for digit, values in enumerate(digit_layers):
            pv = polyvec_t(ring, 1)
            pv[0] = poly_t(ring, values)
            vectors.append(pv)
            norms.append(max(1, sum(value * value for value in values)))
            witness_refs.append(("c", layer, digit))
    if any(value >= (1 << 64) for value in norms):
        raise ValueError("grouped Pack norm exceeds uint64")
    equation_count = bridge.equation_count()
    statement = proof_statement(
        [degree] * len(vectors), [vector.dim for vector in vectors], norms,
        num_constraints=equation_count, primesize="38", zk=zk)
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
                    pv[index] = poly_t(ring, values)
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
        rhs = bridge.target_digits[layer] if layer < len(bridge.target_digits) else [0] * degree
        statement.append_statement(terms, witnesses, poly_t(ring, rhs))
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
