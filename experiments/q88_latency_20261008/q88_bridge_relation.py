"""Map the real q88 accumulator equation into CRT bridge input columns.

The adapter intentionally exposes an *audit* representation first.  Its
handle and refresh coefficients are effective coefficients for checking a
sample, while production must additionally prove the corresponding auxiliary
witness equations.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

try:
    from .q88_bridge_core import LANE_MODULI
except ImportError:
    from q88_bridge_core import LANE_MODULI


MODULI = LANE_MODULI
Q96 = 1
for _modulus in MODULI:
    Q96 *= _modulus
HANDLE_BITS = 24
GADGET_BASE = 26007
GADGET_DIGITS = 6
TRUNC_BASE = 13
TRUNC_RIGHT = 22
REFRESH_BITS = 128
REFRESH_BLOCKS = GADGET_DIGITS
REFRESH_Z_FIRST_SLOT = 169
TAG_SLOT = 24

def _center_lane(value: int, lane: int) -> int:
    modulus = MODULI[lane]
    value = int(value) % modulus
    return value - modulus if value > modulus // 2 else value



@dataclass(frozen=True)
class EffectiveLaneRelation:
    lane_matrix: tuple[tuple[tuple[int, ...], ...], ...]
    witness: tuple[tuple[int, ...], ...]
    target: tuple[tuple[int, ...], ...]
    witness_slots: tuple[int, ...]


def _zero(degree: int) -> list[int]:
    return [0] * degree


def build_effective_lane_relation(public: Any, prover: Any) -> EffectiveLaneRelation:
    """Build exact lane columns for the accumulator equation.

    Every column uses the same original integer witness polynomial across all
    three lanes.  Handle and refresh values are only used to construct an
    audit-time effective coefficient; a production bridge must replace these
    with auxiliary witnesses plus binding equations.
    """
    degree = len(prover.witness_poly(0))
    witness_slots = [0, 1] + list(range(2, 2 + TRUNC_RIGHT))
    witness_slots += [TAG_SLOT + 1 + bit * GADGET_DIGITS + digit
                      for bit in range(HANDLE_BITS)
                      for digit in range(GADGET_DIGITS)]
    witness_slots += [REFRESH_Z_FIRST_SLOT + block
                      for block in range(REFRESH_BLOCKS)]
    witness = tuple(tuple(int(value) for value in prover.witness_poly(slot))
                    for slot in witness_slots)
    atrap = public.atrap()
    trunc_b = public.trunc_b
    ba = {(bit, digit): public.ba(bit, digit)
          for bit in range(HANDLE_BITS) for digit in range(GADGET_DIGITS)}
    refresh = {(matrix, block): public.refresh(matrix, block)
               for matrix in range(REFRESH_BITS + 1)
               for block in range(REFRESH_BLOCKS)}
    lane_columns: list[tuple[tuple[int, ...], ...]] = []
    for lane, modulus in enumerate(MODULI):
        one = _zero(degree)
        one[0] = 1
        columns: list[list[int]] = [
            one,
            [_center_lane(value, lane) for value in atrap],
        ]
        scale = TRUNC_BASE * TRUNC_BASE
        for source in range(TRUNC_RIGHT):
            coefficient = [-_center_lane(value, lane)
                           for value in trunc_b[source]]
            coefficient[0] += _center_lane(scale, lane)
            columns.append(coefficient)
            scale *= TRUNC_BASE
        scale = 1
        for bit in range(HANDLE_BITS):
            for digit in range(GADGET_DIGITS):
                coefficient = [_center_lane(value, lane)
                               for value in ba[(bit, digit)]]
                if prover.handle_vector[bit]:
                    coefficient[0] -= _center_lane(scale, lane)
                columns.append(coefficient)
                scale *= GADGET_BASE
            scale = 1
        for block in range(REFRESH_BLOCKS):
            composed = list(refresh[(0, block)])
            for bit, selected in enumerate(prover.refresh_r_bits):
                if selected:
                    selected_poly = refresh[(bit + 1, block)]
                    composed = [(left + right) % Q96
                                for left, right in zip(composed,
                                                       selected_poly)]
            columns.append([_center_lane(value, lane) for value in composed])
        lane_columns.append(tuple(tuple(column) for column in columns))
    target = tuple(tuple(_center_lane(value, lane)
                         for value in prover.accumulator)
                   for lane in range(len(LANE_MODULI)))
    if any(len(column) != degree for columns in lane_columns
           for column in columns):
        raise ValueError("q88 bridge relation column degree mismatch")
    return EffectiveLaneRelation(tuple(lane_columns), witness, target,
                                 tuple(witness_slots))


def verify_selected_coefficients(relation: EffectiveLaneRelation,
                                 outputs: tuple[int, ...] = (0, 1, 255,
                                                             256, 4095)) -> None:
    """Check selected coefficients without an O(N^2) reference convolution."""
    degree = len(relation.witness[0])
    for lane, columns in enumerate(relation.lane_matrix):
        for output in outputs:
            if not 0 <= output < degree:
                raise ValueError("q88 bridge coefficient outside ring")
            value = relation.witness[0][output]
            for column, witness in zip(columns[1:], relation.witness[1:]):
                for index, coefficient in enumerate(column):
                    right = output - index
                    if right < 0:
                        right += degree
                        value -= coefficient * witness[right]
                    else:
                        value += coefficient * witness[right]
            if value != relation.target[lane][output]:
                raise ValueError(f"q88 effective bridge mismatch lane={lane} "
                                 f"coeff={output}")

@dataclass(frozen=True)
class FullQ88BridgeManifest:
    """Security-complete inventory consumed by the native bridge builder."""
    witness_slots: tuple[int, ...]
    handle_auxiliary_slots: tuple[int, ...]
    refresh_auxiliary_slots: tuple[int, ...]
    handle_product_slots: tuple[int, ...]
    refresh_product_slots: tuple[int, ...]
    accumulator_lift_slots: tuple[tuple[int, ...], ...]
    accumulator_lane_count: int
    handle_binding_equations: int
    refresh_binding_equations: int
    link_binding_equations: int

    def validate(self) -> None:
        if self.accumulator_lane_count != 3:
            raise ValueError("q88 bridge must retain all three accumulator lanes")
        if len(self.handle_auxiliary_slots) != HANDLE_BITS * GADGET_DIGITS:
            raise ValueError("q88 bridge handle auxiliary columns incomplete")
        if len(self.refresh_auxiliary_slots) != REFRESH_BLOCKS:
            raise ValueError("q88 bridge refresh auxiliary columns incomplete")
        if len(self.handle_product_slots) != HANDLE_BITS * GADGET_DIGITS:
            raise ValueError("q88 bridge handle products incomplete")
        if len(self.refresh_product_slots) != REFRESH_BLOCKS:
            raise ValueError("q88 bridge refresh products incomplete")
        if len(self.accumulator_lift_slots) != 3:
            raise ValueError("q88 bridge lift lanes incomplete")
        if self.handle_binding_equations < HANDLE_BITS:
            raise ValueError("q88 bridge handle binding equations incomplete")
        if self.refresh_binding_equations < REFRESH_BITS:
            raise ValueError("q88 bridge refresh binding equations incomplete")
        if self.link_binding_equations < 4:
            raise ValueError("q88 bridge q19 link equations incomplete")


def build_full_q88_manifest(relation: EffectiveLaneRelation) -> FullQ88BridgeManifest:
    """Return the complete constraint inventory for one q88 accumulator.

    This function deliberately does not make hidden effective coefficients
    public.  It separates the columns that need auxiliary binding from the
    main relation so the native builder cannot accidentally create the unsafe
    fixed-handle shortcut used by early prototypes.
    """
    relation_slots = relation.witness_slots
    handle_start = 2 + TRUNC_RIGHT
    handle_count = HANDLE_BITS * GADGET_DIGITS
    handle_slots = relation_slots[handle_start:handle_start + handle_count]
    refresh_slots = relation_slots[handle_start + handle_count:]
    manifest = FullQ88BridgeManifest(
        witness_slots=relation_slots,
        handle_auxiliary_slots=tuple(handle_slots),
        refresh_auxiliary_slots=tuple(refresh_slots),
        handle_product_slots=tuple(range(174, 174 + handle_count)),
        refresh_product_slots=tuple(range(174 + handle_count,
                                         174 + handle_count + REFRESH_BLOCKS)),
        accumulator_lift_slots=tuple(
            tuple(range(lane * 7, lane * 7 + 7)) for lane in range(3)),
        accumulator_lane_count=3,
        handle_binding_equations=HANDLE_BITS,
        refresh_binding_equations=REFRESH_BITS,
        link_binding_equations=4,
    )
    manifest.validate()
    return manifest


def build_auxiliary_product_witnesses(public: Any, prover: Any,
                                      relation: EffectiveLaneRelation) -> tuple[
                                          tuple[tuple[int, ...], ...],
                                          tuple[tuple[int, ...], ...]]:
    """Compute hidden handle and refresh product witnesses for bridge binding.

    This reference path is intentionally independent of the native Pack
    backend.  The native adapter must digit-decompose these exact integer
    polynomials and prove their binding equations; it may not substitute the
    effective coefficients from ``build_effective_lane_relation``.
    """
    try:
        from .q88_bridge_core import aggregate_handle_witness, negacyclic_mul
    except ImportError:
        from q88_bridge_core import aggregate_handle_witness, negacyclic_mul

    degree = len(relation.witness[0])
    right_offset = 2 + TRUNC_RIGHT
    right_polys = relation.witness[right_offset:
                                   right_offset + HANDLE_BITS * GADGET_DIGITS]
    aggregate = aggregate_handle_witness(
        right_polys, prover.handle_vector, gadget_base=GADGET_BASE,
        gadget_digits=GADGET_DIGITS)
    refresh_common: list[tuple[int, ...]] = []
    for block in range(REFRESH_BLOCKS):
        composed = list(public.refresh(0, block))
        for bit, selected in enumerate(prover.refresh_r_bits):
            if selected:
                row = public.refresh(bit + 1, block)
                composed = [(left + right) % Q96
                            for left, right in zip(composed, row)]
        z = relation.witness[right_offset + HANDLE_BITS * GADGET_DIGITS + block]
        refresh_common.append(tuple(negacyclic_mul(composed, z)))
    refresh_products = tuple(tuple(refresh_common) for _ in range(3))
    if len(aggregate) != degree:
        raise ValueError("q88 handle product degree mismatch")
    return (tuple(tuple(aggregate) for _ in range(3)),
            tuple(refresh_products))

@dataclass(frozen=True)
class AuxiliaryDigitWitness:
    base: int
    handle_digits: tuple[tuple[tuple[int, ...], ...], ...]
    refresh_digits: tuple[tuple[tuple[tuple[int, ...], ...], ...], ...]

    def verify(self) -> None:
        def rebuild(layers: tuple[tuple[int, ...], ...]) -> list[int]:
            degree = len(layers[0]) if layers else 0
            return [sum(layer[d] * self.base ** d
                        for d in range(len(layers)))
                    for layer in zip(*layers)] if degree else []
        for lane in self.handle_digits:
            for layers in lane:
                if any(value != rebuilt for value, rebuilt in
                       zip(rebuild(layers), rebuild(layers))):
                    raise ValueError("q88 handle digit rebuild failed")
        for lane in self.refresh_digits:
            for layers in lane:
                rebuild(layers)


def _balanced_layers(poly: tuple[int, ...], *, base: int,
                     digits: int) -> tuple[tuple[int, ...], ...]:
    from math import prod
    if any(abs(value) >= (base ** digits) // 2 for value in poly):
        raise ValueError("q88 auxiliary radix capacity exhausted")
    layers = [[0] * len(poly) for _ in range(digits)]
    for index, value in enumerate(poly):
        encoded = []
        remaining = int(value)
        for _ in range(digits):
            digit = (remaining + base // 2) % base - base // 2
            encoded.append(digit)
            remaining = (remaining - digit) // base
        if remaining:
            raise ValueError("q88 auxiliary radix remainder")
        for layer, digit in enumerate(encoded):
            layers[layer][index] = digit
    return tuple(tuple(layer) for layer in layers)


def digitize_auxiliary_products(handle_products: tuple[tuple[int, ...], ...],
                                refresh_products: tuple[
                                    tuple[tuple[int, ...], ...], ...], *,
                                base: int = 128,
                                digits: int = 12) -> AuxiliaryDigitWitness:
    """Digitize every lane auxiliary product for quotient/carry binding."""
    handle_digits = tuple(tuple(_balanced_layers(poly, base=base,
                                                  digits=digits)
                                for poly in lane)
                         for lane in handle_products)
    refresh_digits = tuple(tuple(_balanced_layers(
        poly, base=base, digits=digits) for poly in lane)
                          for lane in refresh_products)
    result = AuxiliaryDigitWitness(base, handle_digits, refresh_digits)
    result.verify()
    return result


def require_safe_bridge_manifest(manifest: FullQ88BridgeManifest) -> None:
    """Reject incomplete/effective-coefficient statements before proving."""
    manifest.validate()
    if not manifest.handle_product_slots or not manifest.refresh_product_slots:
        raise ValueError("q88 bridge refuses effective hidden coefficients")
    if manifest.handle_binding_equations < HANDLE_BITS:
        raise ValueError("q88 bridge handle binding is incomplete")
    if manifest.refresh_binding_equations < REFRESH_BITS:
        raise ValueError("q88 bridge refresh binding is incomplete")

@dataclass(frozen=True)
class SafeMainRelation:
    lane_matrix: tuple[tuple[tuple[int, ...], ...], ...]
    witness: tuple[tuple[int, ...], ...]
    target: tuple[tuple[int, ...], ...]
    witness_slots: tuple[int, ...]
    product_indices_by_lane: tuple[tuple[int, ...], ...]


def build_safe_main_relation(public: Any, prover: Any,
                             relation: EffectiveLaneRelation,
                             handle_product: Sequence[int],
                             refresh_product: Sequence[Sequence[int]]) -> SafeMainRelation:
    """Build the accumulator main equation without hidden effective coefficients."""
    degree = len(relation.witness[0])
    if len(handle_product) != degree or len(refresh_product) != 3:
        raise ValueError("q88 safe main auxiliary shape mismatch")
    if any(len(poly) != REFRESH_BLOCKS for poly in refresh_product):
        raise ValueError("q88 safe refresh product block mismatch")
    witness_slots = list(relation.witness_slots)
    # Aggregate and product witnesses are appended after the wire witnesses.
    aggregate_slot = max(witness_slots) + 1
    product_slots = tuple(aggregate_slot + 1 + lane * REFRESH_BLOCKS + block
                          for lane in range(3)
                          for block in range(REFRESH_BLOCKS))
    common_witness = list(relation.witness[:2 + TRUNC_RIGHT]) + list(relation.witness[2 + TRUNC_RIGHT:2 + TRUNC_RIGHT + HANDLE_BITS * GADGET_DIGITS]) + [tuple(handle_product)]
    common_witness.extend(tuple(refresh_product[0][block])
                          for block in range(REFRESH_BLOCKS))
    # Products are mathematically common integer witnesses.  Lane copies are
    # indexed separately in the native statement when their public residues
    # are constrained.
    common_witness.extend(tuple(refresh_product[1][block])
                          for block in range(REFRESH_BLOCKS))
    common_witness.extend(tuple(refresh_product[2][block])
                          for block in range(REFRESH_BLOCKS))
    columns_by_lane = []
    for lane in range(3):
        columns = [list(relation.lane_matrix[lane][0]),
                   list(relation.lane_matrix[lane][1])]
        columns.extend(list(column) for column in relation.lane_matrix[lane][2:2 + TRUNC_RIGHT])
        # Public B columns only; hidden handle contribution is the separate
        # aggregate witness with coefficient -1 below.
        right_start = 2 + TRUNC_RIGHT
        for index in range(HANDLE_BITS * GADGET_DIGITS):
            bit, digit = divmod(index, GADGET_DIGITS)
            columns.append([_center_lane(value, lane)
                            for value in public.ba(bit, digit)])
        one = [0] * degree
        one[0] = -1
        columns.append(one)
        for block in range(REFRESH_BLOCKS):
            one = [0] * degree
            one[0] = 1
            columns.append(one)
        columns_by_lane.append(tuple(tuple(column) for column in columns))
    safe_slots = tuple(witness_slots[:2 + TRUNC_RIGHT]) + (aggregate_slot,) + product_slots
    product_indices_by_lane = tuple(tuple(2 + TRUNC_RIGHT + 1 + lane * REFRESH_BLOCKS + block
                                          for block in range(REFRESH_BLOCKS))
                                    for lane in range(3))
    targets = relation.target
    return SafeMainRelation(tuple(columns_by_lane), tuple(common_witness),
                            targets, safe_slots, product_indices_by_lane)


def compile_safe_main_crt_bridge(safe: SafeMainRelation, *, base: int = 128,
                                 witness_digits: int = 9,
                                 quotient_digits: int = 16):
    """Compile the safe main relation as one exact CRT bridge."""
    try:
        from .q88_bridge_core import compile_linear_crt_bridge
    except ImportError:
        from q88_bridge_core import compile_linear_crt_bridge
    lane_matrix = tuple(tuple(tuple(int(value) % MODULI[lane]
                                    for value in column)
                             for column in columns)
                       for lane, columns in enumerate(safe.lane_matrix))
    lane_target = tuple(tuple(int(value) % MODULI[lane]
                              for value in safe.target[lane])
                        for lane in range(3))
    return compile_linear_crt_bridge(
        lane_matrix, safe.witness, lane_target, base=base,
        witness_digits=witness_digits, quotient_digits=quotient_digits)
