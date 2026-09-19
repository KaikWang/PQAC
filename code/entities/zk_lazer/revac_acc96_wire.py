#!/usr/bin/env python3
"""Decode and independently audit the q88 MemVer Labrador wire format.

The prover payload is private and is only suitable for a same-device Unix
socket.  This module deliberately has no network code.  The public instance
contains the accumulator state and q19 commitment needed by the production
joint MemVer/same-handle relation, but never the q88 witness or link opening.
"""

from __future__ import annotations

import argparse
import hashlib
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

try:
    import numpy as _np
except ImportError:  # Preserve the dependency-free audit-tool path.
    _np = None


N = 4096
PHASES = 16
PHASE_N = 256
LANES = 3
MODULI = (1073692673, 536813569, 536690689)
Q96 = MODULI[0] * MODULI[1] * MODULI[2]
HANDLE_BITS = 24
GADGET_BASE = 26007
GADGET_DIGITS = 6
TRUNC_BASE = 13
TRUNC_RIGHT = 22
WITNESS_POLYS = 169
TAG_SLOT = 24
NON_TAG_POLYS = 168
WITNESS_LIMBS = 5
LIFT_LIMBS = 7
# This proof radix is deliberately below the int16_t storage limit.  The
# smaller digits keep LaBRADOR's conservative JL norm parameterizable.
LIMB_RADIX = 4095
LIMB_BOUND = 2047
PUBLIC_MAGIC = b"RV88PP02"
PROVER_MAGIC = b"RV88PV02"
COEFF_BYTES = 11
POLY_BYTES = N * COEFF_BYTES
PUBLIC_BYTES = 8 + 32 + TRUNC_RIGHT * POLY_BYTES
PROVER_BYTES = (8 + 8 + POLY_BYTES + HANDLE_BITS +
                NON_TAG_POLYS * WITNESS_LIMBS * N * 4 +
                LANES * LIFT_LIMBS * N * 4)
Q19 = 405769
Q19_D = 4
Q19_N = 256
LINK_DECOMP_BASE = 256
LINK_DECOMP_LIMBS = 4
LINK_MAGIC = b"RV88LK01"
LINK_BYTES = (8 + 32 +
              (Q19_D + 2 * Q19_D + Q19_D * LINK_DECOMP_LIMBS) *
              Q19_N * 4)
INSTANCE_MAGIC = b"RV88IN01"
INSTANCE_BYTES = 8 + 4 + POLY_BYTES + 32 + Q19_D * Q19_N * 4


def _i32_view(raw: bytes, offset: int, count: int) -> memoryview:
    end = offset + count * 4
    if end > len(raw):
        raise ValueError("truncated q96 int32 limb payload")
    return memoryview(raw)[offset:end].cast("i")


def _i32_within_bound(values: memoryview, bound: int) -> bool:
    """Constant-result range gate with a native vectorized fast path."""
    if not values:
        return True
    if _np is not None:
        array = _np.frombuffer(values, dtype=_np.int32)
        # Comparing min/max directly also handles INT32_MIN without the
        # signed-overflow corner case of abs(INT32_MIN).
        return int(array.min()) >= -bound and int(array.max()) <= bound
    return all(-bound <= int(value) <= bound for value in values)


def _canonical_coefficients(raw: bytes | memoryview) -> Iterator[int]:
    if len(raw) % COEFF_BYTES:
        raise ValueError("q88 packed coefficient data is not coefficient aligned")
    for offset in range(0, len(raw), COEFF_BYTES):
        value = int.from_bytes(raw[offset:offset + COEFF_BYTES], "little")
        if value >= Q96:
            raise ValueError("non-canonical q88 coefficient")
        yield value


def _center_lane(canonical: int, lane: int) -> int:
    modulus = MODULI[lane]
    value = canonical % modulus
    return value - modulus if value > modulus // 2 else value


def _uniform_poly(seed: bytes, domain: int, first: int,
                  second: int) -> list[int]:
    absorbed = seed + struct.pack("<III", domain, first, second)
    output_length = COEFF_BYTES * (N + 16)
    while True:
        stream = hashlib.shake_256(absorbed).digest(output_length)
        values: list[int] = []
        for offset in range(0, len(stream) - COEFF_BYTES + 1, COEFF_BYTES):
            candidate = int.from_bytes(
                stream[offset:offset + COEFF_BYTES], "little")
            if candidate < Q96:
                values.append(candidate)
                if len(values) == N:
                    return values
        output_length *= 2


@dataclass(frozen=True)
class PublicParameters:
    seed: bytes
    trunc_b: tuple[tuple[int, ...], ...]

    @classmethod
    def decode(cls, raw: bytes) -> "PublicParameters":
        if len(raw) != PUBLIC_BYTES or raw[:8] != PUBLIC_MAGIC:
            raise ValueError("invalid q96 public-parameter payload")
        seed = raw[8:40]
        offset = 40
        trunc_b = []
        for _ in range(TRUNC_RIGHT):
            trunc_b.append(tuple(_canonical_coefficients(
                raw[offset:offset + POLY_BYTES])))
            offset += POLY_BYTES
        return cls(seed, tuple(trunc_b))

    def atrap(self) -> list[int]:
        return _uniform_poly(self.seed, 0x41393631, 0, 0)

    def ba(self, bit: int, digit: int) -> list[int]:
        if not 0 <= bit < HANDLE_BITS or not 0 <= digit < GADGET_DIGITS:
            raise IndexError("q96 Ba index out of range")
        return _uniform_poly(self.seed, 0x44393631, bit, digit)


class ProverPayload:
    def __init__(self, raw: bytes):
        if len(raw) != PROVER_BYTES or raw[:8] != PROVER_MAGIC:
            raise ValueError("invalid q96 prover payload")
        self.raw = raw
        self.sequence, self.handle = struct.unpack_from("<II", raw, 8)
        offset = 16
        self.accumulator = tuple(_canonical_coefficients(
            raw[offset:offset + POLY_BYTES]))
        offset += POLY_BYTES
        self.handle_vector = tuple(raw[offset:offset + HANDLE_BITS])
        offset += HANDLE_BITS
        expected_bits = tuple((self.handle >> bit) & 1
                              for bit in range(HANDLE_BITS))
        if self.handle_vector != expected_bits:
            raise ValueError("q96 prover handle bits are not canonical")
        limb_count = NON_TAG_POLYS * WITNESS_LIMBS * N
        self.witness_limbs = _i32_view(raw, offset, limb_count)
        offset += limb_count * 4
        self.lift_limbs = _i32_view(raw, offset, LANES * LIFT_LIMBS * N)
        offset += LANES * LIFT_LIMBS * N * 4
        if offset != len(raw):
            raise ValueError("trailing q96 prover payload data")
        if not _i32_within_bound(self.witness_limbs, LIMB_BOUND):
            raise ValueError("q96 witness limb exceeds LaBRADOR int16 gate")
        if not _i32_within_bound(self.lift_limbs, LIMB_BOUND):
            raise ValueError("q96 exact-lift limb exceeds LaBRADOR int16 gate")

    @staticmethod
    def _wire_slot(slot: int) -> int:
        if not 0 <= slot < WITNESS_POLYS or slot == TAG_SLOT:
            raise IndexError("q96 witness slot is not transmitted")
        return slot if slot < TAG_SLOT else slot - 1

    def witness_phase(self, slot: int, limb: int,
                      phase: int) -> memoryview:
        wire_slot = self._wire_slot(slot)
        if not 0 <= limb < WITNESS_LIMBS or not 0 <= phase < PHASES:
            raise IndexError("q96 witness limb/phase index out of range")
        start = ((wire_slot * WITNESS_LIMBS + limb) * N +
                 phase * PHASE_N)
        return self.witness_limbs[start:start + PHASE_N]

    def witness_poly(self, slot: int) -> list[int]:
        phases = [[int(value) for value in self.witness_phase(slot, limb, phase)]
                  for limb in range(WITNESS_LIMBS)
                  for phase in range(PHASES)]
        out = [0] * N
        for phase in range(PHASES):
            for coefficient in range(PHASE_N):
                out[PHASES * coefficient + phase] = sum(
                    phases[limb * PHASES + phase][coefficient] *
                    LIMB_RADIX ** limb for limb in range(WITNESS_LIMBS))
        return out

    def lift_poly(self, lane: int) -> list[int]:
        if not 0 <= lane < LANES:
            raise IndexError("q96 lift lane out of range")
        out = [0] * N
        lane_start = lane * LIFT_LIMBS * N
        for phase in range(PHASES):
            for coefficient in range(PHASE_N):
                phased = phase * PHASE_N + coefficient
                out[PHASES * coefficient + phase] = sum(
                    int(self.lift_limbs[lane_start + limb * N + phased]) *
                    LIMB_RADIX ** limb for limb in range(LIFT_LIMBS))
        return out

    def lift_phase(self, lane: int, limb: int, phase: int) -> memoryview:
        if (not 0 <= lane < LANES or not 0 <= limb < LIFT_LIMBS or
                not 0 <= phase < PHASES):
            raise IndexError("q96 lift lane/limb/phase index out of range")
        start = ((lane * LIFT_LIMBS + limb) * N + phase * PHASE_N)
        return self.lift_limbs[start:start + PHASE_N]


@dataclass(frozen=True)
class InstancePayload:
    """Public per-show input consumed by a GCS verifier."""

    sequence: int
    accumulator: tuple[int, ...]
    context_hash: bytes
    commitment: tuple[tuple[int, ...], ...]

    @classmethod
    def decode(cls, raw: bytes) -> "InstancePayload":
        if len(raw) != INSTANCE_BYTES or raw[:8] != INSTANCE_MAGIC:
            raise ValueError("invalid q96 public verifier instance")
        sequence = struct.unpack_from("<I", raw, 8)[0]
        offset = 12
        accumulator = tuple(_canonical_coefficients(
            raw[offset:offset + POLY_BYTES]))
        offset += POLY_BYTES
        context_hash = raw[offset:offset + 32]
        offset += 32
        commitment = []
        for _ in range(Q19_D):
            values = struct.unpack_from(f"<{Q19_N}I", raw, offset)
            offset += Q19_N * 4
            if any(value >= Q19 for value in values):
                raise ValueError("non-canonical q19 commitment in q96 instance")
            commitment.append(tuple(values))
        if offset != len(raw):
            raise ValueError("trailing q96 public verifier instance data")
        return cls(sequence, accumulator, context_hash, tuple(commitment))


class ZeroProverPayload:
    """Shape-only witness used to construct a public verifier statement."""

    def __init__(self, instance: InstancePayload):
        self.sequence = instance.sequence
        self.handle = 0
        self.accumulator = instance.accumulator
        self.handle_vector = (0,) * HANDLE_BITS
        self._zero = (0,) * PHASE_N

    def witness_phase(self, slot: int, limb: int,
                      phase: int) -> tuple[int, ...]:
        ProverPayload._wire_slot(slot)
        if not 0 <= limb < WITNESS_LIMBS or not 0 <= phase < PHASES:
            raise IndexError("q96 witness limb/phase index out of range")
        return self._zero

    def lift_phase(self, lane: int, limb: int,
                   phase: int) -> tuple[int, ...]:
        if (not 0 <= lane < LANES or not 0 <= limb < LIFT_LIMBS or
                not 0 <= phase < PHASES):
            raise IndexError("q96 lift lane/limb/phase index out of range")
        return self._zero


def _q19_uniform_poly(seed: bytes, domain: int, row: int,
                      col: int) -> tuple[int, ...]:
    """Match poly_q_uniform for the q19 same-handle public matrices."""
    absorbed = seed + struct.pack("<IQQ", domain, row, col)
    output_length = 2048
    while True:
        stream = hashlib.shake_128(absorbed).digest(output_length)
        coefficients: list[int] = []
        for offset in range(0, len(stream) - 4, 5):
            word = int.from_bytes(stream[offset:offset + 5], "little")
            low = word & ((1 << 19) - 1)
            high = (word >> 20) & ((1 << 19) - 1)
            if low < Q19:
                coefficients.append(low)
            if len(coefficients) < Q19_N and high < Q19:
                coefficients.append(high)
            if len(coefficients) == Q19_N:
                return tuple(coefficients)
        output_length *= 2


def _q19_center(value: int) -> int:
    value %= Q19
    return value - Q19 if value > Q19 // 2 else value


@dataclass(frozen=True)
class LinkPayload:
    context_hash: bytes
    commitment: tuple[tuple[int, ...], ...]
    randomness: tuple[tuple[tuple[int, ...], ...], ...]
    lift_limbs: tuple[tuple[tuple[int, ...], ...], ...]

    @classmethod
    def decode(cls, raw: bytes) -> "LinkPayload":
        if len(raw) != LINK_BYTES or raw[:8] != LINK_MAGIC:
            raise ValueError("invalid q96 same-handle link payload")
        offset = 8
        context_hash = raw[offset:offset + 32]
        offset += 32

        def read_i32_poly() -> tuple[int, ...]:
            nonlocal offset
            values = struct.unpack_from(f"<{Q19_N}i", raw, offset)
            offset += Q19_N * 4
            return tuple(values)

        commitment = tuple(read_i32_poly() for _ in range(Q19_D))
        if any(not 0 <= value < Q19
               for poly in commitment for value in poly):
            raise ValueError("non-canonical q19 link commitment")
        randomness = tuple(
            tuple(read_i32_poly() for _ in range(Q19_D))
            for _ in range(2))
        if any(value not in (0, 1)
               for part in randomness for poly in part for value in poly):
            raise ValueError("q19 link opening is not binary")
        lift_limbs = tuple(
            tuple(read_i32_poly() for _ in range(LINK_DECOMP_LIMBS))
            for _ in range(Q19_D))
        if any(abs(value) > LINK_DECOMP_BASE // 2
               for row in lift_limbs for poly in row for value in poly):
            raise ValueError("q19 link lift digit exceeds radix-256 bound")
        if offset != len(raw):
            raise ValueError("trailing q96 same-handle link payload data")
        return cls(context_hash, commitment, randomness, lift_limbs)

    def public_matrices(self) -> tuple[
            tuple[tuple[tuple[int, ...], ...], ...],
            tuple[tuple[tuple[int, ...], ...], ...]]:
        link_a = tuple(tuple(
            _q19_uniform_poly(self.context_hash, 29, row, col)
            for col in range(Q19_D)) for row in range(Q19_D))
        link_d = tuple(tuple(
            _q19_uniform_poly(self.context_hash, 30, row, col)
            for col in range(Q19_D)) for row in range(Q19_D))
        return link_a, link_d

    def audit_relation(self, handle_bits: tuple[int, ...]) -> None:
        """Independently check the exact q19 commitment opening equation."""
        if len(handle_bits) != HANDLE_BITS or any(
                bit not in (0, 1) for bit in handle_bits):
            raise ValueError("invalid common handle bits for q19 link audit")
        link_a, link_d = self.public_matrices()
        handle_poly = tuple(handle_bits) + (0,) * (Q19_N - HANDLE_BITS)

        def mul_coefficient(left: tuple[int, ...], right: tuple[int, ...],
                            output: int) -> int:
            result = 0
            for index, value in enumerate(left):
                right_index = output - index
                if right_index < 0:
                    result -= value * right[right_index + Q19_N]
                else:
                    result += value * right[right_index]
            return result

        for row in range(Q19_D):
            centered_a = [tuple(_q19_center(value) for value in poly)
                          for poly in link_a[row]]
            centered_d = tuple(_q19_center(value)
                               for value in link_d[row][0])
            for coefficient in range(Q19_N):
                expression = self.randomness[0][row][coefficient]
                for col in range(Q19_D):
                    expression += mul_coefficient(
                        centered_a[col], self.randomness[1][col],
                        coefficient)
                expression += mul_coefficient(
                    centered_d, handle_poly, coefficient)
                lift = sum(
                    self.lift_limbs[row][limb][coefficient] *
                    LINK_DECOMP_BASE ** limb
                    for limb in range(LINK_DECOMP_LIMBS))
                target = _q19_center(self.commitment[row][coefficient])
                if expression + Q19 * lift != target:
                    raise ValueError(
                        "q19 same-handle exact relation mismatch at "
                        f"row={row}, coefficient={coefficient}")


class PublicLinkPayload:
    """Same-handle statement data without its private commitment opening."""

    def __init__(self, instance: InstancePayload):
        self.context_hash = instance.context_hash
        self.commitment = instance.commitment
        zero = (0,) * Q19_N
        self.randomness = tuple(tuple(zero for _ in range(Q19_D))
                                for _ in range(2))
        self.lift_limbs = tuple(
            tuple(zero for _ in range(LINK_DECOMP_LIMBS))
            for _ in range(Q19_D))

    def public_matrices(self):
        link_a = tuple(tuple(
            _q19_uniform_poly(self.context_hash, 29, row, col)
            for col in range(Q19_D)) for row in range(Q19_D))
        link_d = tuple(tuple(
            _q19_uniform_poly(self.context_hash, 30, row, col)
            for col in range(Q19_D)) for row in range(Q19_D))
        return link_a, link_d


def _negacyclic_coefficient(left: list[int] | tuple[int, ...],
                            right: list[int], output: int) -> int:
    result = 0
    for index, value in enumerate(left):
        right_index = output - index
        if right_index < 0:
            right_index += N
            result -= value * right[right_index]
        else:
            result += value * right[right_index]
    return result


def audit_integer_relation(public: PublicParameters, prover: ProverPayload,
                           coefficients: tuple[int, ...] =
                           (0, 1, 255, 256, 4095)) -> None:
    """Cross-check selected coefficients of all three exact CRT equations."""
    witness = {slot: prover.witness_poly(slot)
               for slot in range(WITNESS_POLYS) if slot != TAG_SLOT}
    lifts = [prover.lift_poly(lane) for lane in range(LANES)]
    atrap = public.atrap()
    ba = {(bit, digit): public.ba(bit, digit)
          for bit in range(HANDLE_BITS)
          for digit in range(GADGET_DIGITS)}
    for lane, modulus in enumerate(MODULI):
        atrap_lane = [_center_lane(value, lane) for value in atrap]
        trunc_lane = [[_center_lane(value, lane) for value in poly]
                      for poly in public.trunc_b]
        ba_lane = {key: [_center_lane(value, lane) for value in poly]
                   for key, poly in ba.items()}
        for output in coefficients:
            expression = witness[0][output]
            expression += _negacyclic_coefficient(
                atrap_lane, witness[1], output)
            scale = TRUNC_BASE * TRUNC_BASE
            for index in range(TRUNC_RIGHT):
                v2 = witness[2 + index]
                expression -= _negacyclic_coefficient(
                    trunc_lane[index], v2, output)
                expression += _center_lane(scale, lane) * v2[output]
                scale *= TRUNC_BASE
            scale = 1
            for bit in range(HANDLE_BITS):
                for digit in range(GADGET_DIGITS):
                    right = witness[TAG_SLOT + 1 + bit * GADGET_DIGITS + digit]
                    expression += _negacyclic_coefficient(
                        ba_lane[(bit, digit)], right, output)
                    if prover.handle_vector[bit]:
                        expression -= _center_lane(scale, lane) * right[output]
                    scale *= GADGET_BASE
                scale = 1
            accumulator = _center_lane(prover.accumulator[output], lane)
            if expression + modulus * lifts[lane][output] != accumulator:
                raise ValueError(
                    f"q96 exact relation mismatch at lane={lane}, coeff={output}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--public", required=True, type=Path)
    parser.add_argument("--prover", required=True, type=Path)
    parser.add_argument("--full-audit", action="store_true",
                        help="audit selected exact integer coefficients")
    args = parser.parse_args()
    public = PublicParameters.decode(args.public.read_bytes())
    prover = ProverPayload(args.prover.read_bytes())
    if args.full_audit:
        audit_integer_relation(public, prover)
    print("q96 wire decode passed: sequence=%d, handle=0x%08x, "
          "public=%.2f MiB, private-local=%.2f MiB%s" % (
              prover.sequence, prover.handle,
              PUBLIC_BYTES / (1024 * 1024), PROVER_BYTES / (1024 * 1024),
              ", exact integer relation audited" if args.full_audit else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
