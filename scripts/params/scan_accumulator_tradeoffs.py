#!/usr/bin/env python3
"""Scan l=20 accumulator profiles under the PKC'26 Theorem 4.3 bound.

This is a parameter-selection aid, not a security certificate.  It deliberately
ignores the Core-SVP acceptance threshold at the user's request, but still
requires distinct NTT-prime CRT factors and computes the truncated-sampler
width used by audit_accumulator_security.py.  The reported gamma_max is the
largest integer for which that audit has beta_1 < q.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

from audit_accumulator_security import audit, is_prime_64, ntt_prime_below


RING_DEGREE = 4096
TRUNC_BASE = 13
TRUNC_L = 2
SPECTRAL_SLACK = 32.0


@dataclass(frozen=True)
class Profile:
    label: str
    factor_bits: tuple[int, ...]
    digits: int
    base: int | None = None


PROFILES_20 = (
    Profile("q72/k8, compact accumulator", (24, 24, 24), 8),
    Profile("q72/k10, 10k edge", (24, 24, 24), 10),
    Profile("q75/k7, compact 1k", (25, 25, 25), 7),
    Profile("q80/k6, compact 1k", (27, 27, 26), 6),
    Profile("q80/k7, compact 10k", (27, 27, 26), 7),
    Profile("q85/k6, balanced 10k", (29, 28, 28), 6),
    Profile("q88/k5, compact witness at 1k", (30, 29, 29), 5),
    Profile("q88/k6, balanced 10k", (30, 29, 29), 6),
    Profile("q94/k5, compact witness at 10k", (32, 31, 31), 5),
    Profile(
        "q96/k5, current arithmetic",
        (32, 32, 32),
        5,
        602249,
    ),
    Profile("q100/k5, high capacity", (25, 25, 25, 25), 5),
)

PROFILES_24 = (
    Profile("q72/k11, minimum accumulator", (24, 24, 24), 11),
    Profile("q80/k7, compact wire", (27, 27, 26), 7),
    Profile("q85/k6, balanced size", (29, 28, 28), 6),
    Profile("q88/k6, larger theorem margin", (30, 29, 29), 6),
    Profile("q94/k5, compact witness", (32, 31, 31), 5),
    Profile(
        "q96/k5, current arithmetic",
        (32, 32, 32),
        5,
        602249,
    ),
)


def ntt_factors(bits: tuple[int, ...]) -> tuple[int, ...]:
    factors: list[int] = []
    used: set[int] = set()
    for width in bits:
        factor = ntt_prime_below(width, RING_DEGREE)
        while factor in used:
            factor -= 2 * RING_DEGREE
            while not is_prime_64(factor):
                factor -= 2 * RING_DEGREE
        factors.append(factor)
        used.add(factor)
    return tuple(factors)


def minimum_odd_base(q: int, digits: int) -> int:
    low = 3
    high = 3
    while high**digits < q:
        high *= 2
    while low < high:
        middle = (low + high) // 2
        if middle**digits >= q:
            high = middle
        else:
            low = middle + 1
    if low % 2 == 0:
        low += 1
    while low > 3 and (low - 2) ** digits >= q:
        low -= 2
    return low


def run_audit(
    q: int,
    factors: tuple[int, ...],
    base: int,
    gamma: int,
    handle_bits: int,
):
    return audit(
        q,
        ring_degree=RING_DEGREE,
        module_rows=1,
        deletion_bound=gamma,
        handle_bits=handle_bits,
        eval_base=base,
        eval_digit_bound=(base - 1) // 2,
        trunc_base=TRUNC_BASE,
        trunc_l=TRUNC_L,
        spectral_slack=SPECTRAL_SLACK,
        q_factors=factors,
        target_quantum_bits=0,
    )


def maximum_deletions(
    q: int, factors: tuple[int, ...], base: int, handle_bits: int
) -> int:
    low = 0
    high = 1
    while run_audit(q, factors, base, high, handle_bits).beta1_below_q:
        low = high
        high *= 2
    while low + 1 < high:
        middle = (low + high) // 2
        if run_audit(q, factors, base, middle, handle_bits).beta1_below_q:
            low = middle
        else:
            high = middle
    return low


def delta4_worst_kib(
    digits: int, base: int, accumulated_updates: int = 1
) -> tuple[int, float]:
    if accumulated_updates < 1:
        raise ValueError("accumulated_updates must be positive")
    digit_bound = (base - 1) // 2
    coefficient_bound = (
        accumulated_updates
        * digits
        * RING_DEGREE
        * digit_bound
        * digit_bound
    )
    zigzag_bits = (2 * coefficient_bound).bit_length()
    payload = 4 * digits * RING_DEGREE * zigzag_bits / 8 / 1024
    return zigzag_bits, payload


def print_table(targets: tuple[int, ...], handle_bits: int) -> None:
    print(
        "| profile | CRT factors | q bits | base | digits | gamma_max | "
        + " | ".join(f"margin@{target}" for target in targets)
        + " | accumulator tight/fixed KiB | witness tight/fixed MiB | "
        "Ba tight/fixed MiB | "
        f"DELTA4 single/cumulative@{max(targets)} KiB |"
    )
    print(
        "|---|---|---:|---:|---:|---:|"
        + "---:|" * len(targets)
        + "---:|---:|---:|---:|"
    )
    profiles = PROFILES_24 if handle_bits == 24 else PROFILES_20
    for profile in profiles:
        factors = ntt_factors(profile.factor_bits)
        q = math.prod(factors)
        base = profile.base or minimum_odd_base(q, profile.digits)
        if base**profile.digits < q:
            raise ValueError(f"{profile.label}: gadget does not cover q")
        maximum = maximum_deletions(q, factors, base, handle_bits)
        audits = [
            run_audit(q, factors, base, target, handle_bits)
            for target in targets
        ]
        reference = audits[0]
        _, delta_kib = delta4_worst_kib(profile.digits, base)
        _, cumulative_delta_kib = delta4_worst_kib(
            profile.digits, base, max(targets)
        )
        byte_scale = math.ceil(q.bit_length() / 8) * 8 / q.bit_length()
        margins = " | ".join(
            f"{item.beta1_margin_log2:.3f}" for item in audits
        )
        print(
            f"| {profile.label} | {' x '.join(map(str, factors))} | "
            f"{q.bit_length()} | {base} | {profile.digits} | {maximum:,} | "
            f"{margins} | {reference.accumulator_kib:.2f}/"
            f"{reference.accumulator_kib * byte_scale:.2f} | "
            f"{reference.full_witness_mib:.3f}/"
            f"{reference.full_witness_mib * byte_scale:.3f} | "
            f"{reference.ba_mib:.3f}/"
            f"{reference.ba_mib * byte_scale:.3f} | "
            f"{delta_kib:.1f}/{cumulative_delta_kib:.1f} |"
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--handle-bits",
        type=int,
        choices=(20, 24),
        default=20,
        help="revocation-handle length",
    )
    parser.add_argument(
        "--targets",
        type=int,
        nargs="+",
        default=(1000, 10000),
        help="deletion counts for which beta_1 margins are printed",
    )
    args = parser.parse_args()
    if any(target <= 0 for target in args.targets):
        parser.error("targets must be positive")
    print_table(tuple(args.targets), args.handle_bits)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
