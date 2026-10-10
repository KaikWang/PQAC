"""q88-specific CRT recombination primitives for the future Pack bridge.

The current production compiler proves three proof-prime lane equations.  A
quotient/carry bridge must first combine those equations over the original
88-bit CRT modulus.  These helpers do only that exact integer step; they do
not construct a proof statement or choose security parameters.
"""

from __future__ import annotations

from math import prod
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
