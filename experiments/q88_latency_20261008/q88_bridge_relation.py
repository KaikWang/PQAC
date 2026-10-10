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

from code.entities.zk_lazer.revac_acc96_wire import (
    GADGET_BASE, GADGET_DIGITS, HANDLE_BITS, MODULI, REFRESH_BITS,
    REFRESH_BLOCKS, REFRESH_Z_FIRST_SLOT, TAG_SLOT, TRUNC_BASE,
    TRUNC_RIGHT, _center_lane,
)


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
                    composed = [(left + right) % (1 << 88)
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
