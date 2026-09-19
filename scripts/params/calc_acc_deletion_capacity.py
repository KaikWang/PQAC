#!/usr/bin/env python3
"""Reproducible balanced-base-3 accumulator deletion-capacity calculator.

This script does not modify the production parameter profile.  It reports:

1. the correctness-only ceiling obtained from beta_0 = B_0 + gamma*C_delta;
2. a stricter PKC'26-Theorem-4.3 *necessary* ceiling, obtained by setting the
   security theorem's Gaussian width to its smallest deletion-dependent value
   s=T and evaluating beta_1 < q.

The second number is not a final security claim.  A final profile must also
derive the truncated sampler covariance/tail bound, include every lower bound
on s, account for the proof system's relaxed extraction, and pass an
M-SIS/M-LWE estimator.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Capacity:
    q: int
    q_bits: int
    eval_digits: int
    trunc_digits: int
    delta_slope: float
    correctness_ceiling: int
    theorem43_necessary_ceiling: int
    accumulator_bytes: int
    full_witness_bytes: int
    ba_bytes: int
    relation_polynomials_projection: int


def ceil_log(value: int, base: int) -> int:
    if value <= 1 or base <= 1:
        raise ValueError("value and base must be greater than one")
    power = 1
    digits = 0
    while power < value:
        power *= base
        digits += 1
    return digits


def capacity(
    q: int,
    *,
    ring_n: int = 256,
    module_d: int = 4,
    handle_bits: int = 32,
    eval_base: int = 3,
    digit_bound: int = 1,
    trunc_base: int = 14,
    trunc_l: int = 2,
    fresh_bound: float = 0.0,
) -> Capacity:
    q_bits = q.bit_length()
    eval_digits = ceil_log(q, eval_base)
    trunc_digits = ceil_log(q, trunc_base)
    if trunc_digits <= trunc_l:
        raise ValueError("truncation parameter l must be smaller than k")

    m_prime = module_d * eval_digits
    left_ring_coordinates = module_d * (2 + trunc_digits - trunc_l)
    gadget_coefficients = ring_n * m_prime

    # PKC'26 Lemma 4.1 / Theorem 4.2 with the actual signed-digit bound delta.
    delta_slope = (
        digit_bound
        * math.sqrt(gadget_coefficients)
        * math.sqrt(
            1.0
            + (
                digit_bound
                * ring_n
                * m_prime
                * math.sqrt(handle_bits - 1)
            )
            ** 2
        )
    )
    correctness_ceiling = max(
        0, math.floor((q - 1.0 - fresh_bound) / delta_slope)
    )

    # Translate the dimensions in PKC'26 Theorem 4.3:
    # paper d -> ring_n, paper m' -> m_prime, paper m -> left coordinates.
    theta2 = math.sqrt(2.0 * ring_n / 3.0) * (
        math.sqrt(left_ring_coordinates + handle_bits * m_prime)
        + math.log(math.sqrt(ring_n))
    )
    full_coefficients = ring_n * (
        left_ring_coordinates + handle_bits * m_prime
    )

    # With the smallest deletion-dependent choice s=T, beta_1 is affine in
    # T=gamma*C_delta.  Any larger sampler/security lower bound on s can only
    # reduce this ceiling.
    beta1_multiplier = (
        math.sqrt(left_ring_coordinates * ring_n)
        + theta2
        * (math.sqrt(ring_n * handle_bits * m_prime) + 1.0)
        + (math.sqrt(full_coefficients) + 1.0)
        * math.sqrt(1.0 + theta2 * theta2)
    )
    beta1_constant = math.sqrt(left_ring_coordinates * ring_n)
    theorem43_necessary_ceiling = max(
        0,
        math.floor(
            (q - 1.0 - beta1_constant)
            / (delta_slope * beta1_multiplier)
        ),
    )

    vector_count = (
        2 + (trunc_digits - trunc_l) + handle_bits * eval_digits
    )
    accumulator_bytes = (
        module_d * ring_n * q_bits + 7
    ) // 8
    full_witness_bytes = (
        vector_count * module_d * ring_n * q_bits + ring_n + 7
    ) // 8
    ba_bytes = (
        handle_bits
        * eval_digits
        * module_d
        * module_d
        * ring_n
        * q_bits
        + 7
    ) // 8

    # Current grouped relation has 1685 polynomials at 12 eval digits, of
    # which 32*12*4 belong to the right witness.  This is a dimension
    # projection, not a serialized-proof-size prediction.
    fixed_relation_polynomials = 1685 - 32 * 12 * 4
    relation_polynomials_projection = (
        fixed_relation_polynomials
        + handle_bits * eval_digits * module_d
    )

    return Capacity(
        q=q,
        q_bits=q_bits,
        eval_digits=eval_digits,
        trunc_digits=trunc_digits,
        delta_slope=delta_slope,
        correctness_ceiling=correctness_ceiling,
        theorem43_necessary_ceiling=theorem43_necessary_ceiling,
        accumulator_bytes=accumulator_bytes,
        full_witness_bytes=full_witness_bytes,
        ba_bytes=ba_bytes,
        relation_polynomials_projection=relation_polynomials_projection,
    )


def minimum_q_for_target(target: int, model: str) -> Capacity:
    if target < 1:
        raise ValueError("target must be positive")
    # The selected truncated sampler requires k>l.
    high = 14**2 + 1
    field = (
        "correctness_ceiling"
        if model == "correctness"
        else "theorem43_necessary_ceiling"
    )
    while getattr(capacity(high), field) < target:
        high *= 2
    low = 14**2
    while low + 1 < high:
        mid = (low + high) // 2
        if getattr(capacity(mid), field) >= target:
            high = mid
        else:
            low = mid
    return capacity(high)


def mib(value: int) -> float:
    return value / (1024.0 * 1024.0)


def print_capacity(label: str, item: Capacity) -> None:
    print(
        f"{label}: q≈2^{math.log2(item.q):.3f} ({item.q_bits} bit), "
        f"eval digits={item.eval_digits}, trunc digits={item.trunc_digits}"
    )
    print(f"  C_delta                         {item.delta_slope:,.3f}")
    print(
        "  correctness-only gamma ceiling "
        f"{item.correctness_ceiling:,}"
    )
    print(
        "  Theorem-4.3 necessary ceiling  "
        f"{item.theorem43_necessary_ceiling:,}"
    )
    print(
        f"  accumulator                    {item.accumulator_bytes / 1024:.3f} KiB"
    )
    print(f"  full witness                   {mib(item.full_witness_bytes):.3f} MiB")
    print(f"  Ba                             {mib(item.ba_bytes):.3f} MiB")
    print(
        "  grouped relation projection    "
        f"{item.relation_polynomials_projection:,} polynomials"
    )


def print_rss_diagnostic(
    q: int,
    fresh_bound: float,
    per_update_rms: float,
    target: int,
) -> None:
    """Print a diagnostic square-sum projection, never a security claim.

    This model is useful for deciding whether a centered/subgaussian proof is
    worth pursuing.  The PKC'26 theorem currently used by the project applies
    the triangle inequality instead, so this projection must not be reported
    as a certified deletion capacity.
    """
    projected = math.sqrt(
        fresh_bound * fresh_bound
        + target * per_update_rms * per_update_rms
    )
    available_rms = math.sqrt(
        max(0.0, (q * q - fresh_bound * fresh_bound) / target)
    )
    print("centered square-sum research diagnostic (not certified)")
    print(f"  target deletions                  {target:,}")
    print(f"  supplied per-update RMS norm      {per_update_rms:,.3f}")
    print(f"  projected total norm              {projected:,.3f}")
    print(f"  projected total / q               {projected / q:.6f}")
    print(f"  maximum RMS fitting below q       {available_rms:,.3f}")
    print(
        "  proof obligation                  centered/subgaussian increments, "
        "adaptive transcript control, and a concrete tail probability"
    )


def print_group_schedule_diagnostic(
    measured_max_group_norm_sq: int,
    target: int,
    per_deletion: int,
    maximum: int,
) -> None:
    """Compare a measured trace with a research-only Labrador norm schedule."""
    selected = min(maximum, target * per_deletion)
    print("Labrador grouped-right research schedule (not certified)")
    print(f"  target deletions                  {target:,}")
    print(f"  per-deletion norm² budget         {per_deletion:,}")
    print(f"  schedule cap norm²                {maximum:,}")
    print(f"  selected public norm² bound       {selected:,}")
    print(f"  measured maximum group norm²      {measured_max_group_norm_sq:,}")
    print(
        "  trace result                       " +
        ("pass" if measured_max_group_norm_sq <= selected else "FAIL")
    )
    print(
        "  measured slack                     "
        f"{selected / measured_max_group_norm_sq:.3f}x"
    )
    print(
        "  warning                            trace acceptance is a ZK-layout "
        "measurement, not an M-SIS security proof"
    )


def print_rademacher_window_diagnostic(
    fresh_bound: float,
    per_update_rms: float,
    max_update_norm: float,
    target: int,
    tail_bits: int,
    protected_members: int,
) -> None:
    """Evaluate a signed-update martingale bound without calling it security.

    If the deletion target is fixed before an independently sampled, signed
    Rademacher bit is committed to the authenticated update, the signed
    witness increments form a Hilbert-space martingale difference sequence.
    Pinelis--Bernstein gives

      Pr[max_k ||sum_{t<=k} X_t|| >= z]
        <= 2 exp(-z^2 / (2(V + L z / 3))).

    V and L still need analytical, distribution-independent bounds.
    Supplying benchmark estimates here is useful only for feasibility
    screening; it is not a replacement for those proof obligations.
    """
    if max_update_norm <= 0.0:
        return
    variance = target * per_update_rms * per_update_rms
    log_factor = (
        (tail_bits * math.log(2.0))
        + math.log(2.0 * protected_members)
    )
    linear = log_factor * max_update_norm / 3.0
    update_bound = linear + math.sqrt(
        linear * linear + 2.0 * log_factor * variance
    )
    total_bound = fresh_bound + update_bound
    print("signed-update martingale diagnostic (conditional, not certified)")
    print(f"  target deletions                  {target:,}")
    print(f"  protected members                 {protected_members:,}")
    print(f"  joint tail budget                 <= 2^-{tail_bits}")
    print(f"  supplied per-update RMS norm      {per_update_rms:,.3f}")
    print(f"  supplied maximum update norm      {max_update_norm:,.3f}")
    print(f"  quadratic variation V             {variance:,.3f}")
    print(f"  Pinelis--Bernstein update bound    {update_bound:,.3f}")
    print(f"  fresh + update bound              {total_bound:,.3f}")
    print(
        "  minimum bit length above bound    "
        f"{max(1, math.ceil(math.log2(total_bound + 1.0)))}"
    )
    print(
        "  proof obligations                  unpredictable non-grindable "
        "signs; analytical V,L bounds; adaptive deletion handling; "
        "Theorem-4.3 reduction and lattice estimator"
    )


def signed_window_bound(
    fresh_bound: float,
    per_update_rms: float,
    max_update_norm: float,
    target: int,
    tail_bits: int,
    protected_members: int,
) -> float:
    """Return the conditional Pinelis--Bernstein fresh+update envelope."""
    variance = target * per_update_rms * per_update_rms
    log_factor = tail_bits * math.log(2.0) + math.log(
        2.0 * protected_members
    )
    linear = log_factor * max_update_norm / 3.0
    return fresh_bound + linear + math.sqrt(
        linear * linear + 2.0 * log_factor * variance
    )


def print_balanced_power_scan(args: argparse.Namespace) -> None:
    """Screen q=3^k candidates without changing the production profile."""
    envelope = signed_window_bound(
        args.fresh_bound,
        args.per_update_rms,
        args.max_update_norm,
        args.target,
        args.tail_bits,
        args.protected_members,
    )
    print("balanced-power candidate scan (conditional, not certified)")
    print(
        "  k   q=3^k       bits  acc(KiB)  witness(MiB)  "
        "relation-polys  envelope/q"
    )
    for exponent in range(12, 17):
        item = capacity(3**exponent, fresh_bound=args.fresh_bound)
        print(
            f" {exponent:2d}  {item.q:10d}  {item.q_bits:4d}  "
            f"{item.accumulator_bytes / 1024:8.3f}  "
            f"{mib(item.full_witness_bytes):12.3f}  "
            f"{item.relation_polynomials_projection:14,d}  "
            f"{envelope / item.q:10.4f}"
        )
    print(
        "  interpretation: envelope/q < 1 is only the signed-update "
        "correctness screen; M-SIS/M-LWE security and the integrated "
        "registered-revocation reduction remain separate obligations."
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--q-bits",
        type=int,
        default=61,
        help="evaluate q=2^bits-1 (default: 61)",
    )
    parser.add_argument(
        "--target",
        type=int,
        default=128,
        help="target deletion count for fixed-point searches",
    )
    parser.add_argument(
        "--per-update-rms",
        type=float,
        default=0.0,
        help=("optional measured per-update RMS norm for a non-certified "
              "square-sum research diagnostic"),
    )
    parser.add_argument(
        "--max-update-norm",
        type=float,
        default=0.0,
        help=("optional maximum update norm for the conditional signed-update "
              "martingale diagnostic"),
    )
    parser.add_argument(
        "--tail-bits",
        type=int,
        default=128,
        help="joint overflow-tail target for the martingale diagnostic",
    )
    parser.add_argument(
        "--protected-members",
        type=int,
        default=128,
        help="members covered by the martingale diagnostic union bound",
    )
    parser.add_argument(
        "--fresh-bound",
        type=float,
        default=305738.704,
        help=("fresh-witness norm used by the correctness/RSS diagnostics "
              "(default: legacy analytical envelope 305738.704)"),
    )
    parser.add_argument(
        "--measured-max-group-norm-sq",
        type=int,
        default=0,
        help="optional trace maximum for the grouped-right ZK diagnostic",
    )
    parser.add_argument(
        "--group-per-deletion-norm-sq",
        type=int,
        default=1 << 23,
        help="research-only grouped-right norm² slope (default: 2^23)",
    )
    parser.add_argument(
        "--group-max-norm-sq",
        type=int,
        default=1 << 31,
        help="research-only grouped-right norm² cap (default: 2^31)",
    )
    parser.add_argument(
        "--scan-balanced-powers",
        action="store_true",
        help=("also screen q=3^12,...,3^16 under the supplied signed-window "
              "diagnostic inputs; does not modify production parameters"),
    )
    args = parser.parse_args()
    if args.q_bits < 2:
        parser.error("--q-bits must be at least 2")

    if args.fresh_bound < 0.0:
        parser.error("--fresh-bound cannot be negative")
    if (args.max_update_norm < 0.0 or args.tail_bits < 1 or
            args.protected_members < 1):
        parser.error("signed-update diagnostic values are invalid")
    if (args.measured_max_group_norm_sq < 0 or
            args.group_per_deletion_norm_sq <= 0 or
            args.group_max_norm_sq <= 0):
        parser.error("grouped-right diagnostic values are invalid")
    current = capacity(405769, fresh_bound=args.fresh_bound)
    candidate = capacity((1 << args.q_bits) - 1)
    correctness = minimum_q_for_target(args.target, "correctness")
    theorem43 = minimum_q_for_target(args.target, "theorem43")

    print_capacity("current candidate (fresh-bound envelope included)", current)
    print()
    print_capacity(f"{args.q_bits}-bit accumulator-modulus candidate", candidate)
    print()
    print_capacity(
        f"minimum update-correctness q for gamma={args.target}", correctness
    )
    print()
    print_capacity(
        f"minimum Theorem-4.3 necessary q for gamma={args.target}", theorem43
    )
    if args.per_update_rms > 0.0:
        print()
        print_rss_diagnostic(
            current.q, args.fresh_bound, args.per_update_rms, args.target
        )
        if args.max_update_norm > 0.0:
            print()
            print_rademacher_window_diagnostic(
                args.fresh_bound,
                args.per_update_rms,
                args.max_update_norm,
                args.target,
                args.tail_bits,
                args.protected_members,
            )
    if args.measured_max_group_norm_sq > 0:
        print()
        print_group_schedule_diagnostic(
            args.measured_max_group_norm_sq,
            args.target,
            args.group_per_deletion_norm_sq,
            args.group_max_norm_sq,
        )
    if (args.scan_balanced_powers and args.per_update_rms > 0.0 and
            args.max_update_norm > 0.0):
        print()
        print_balanced_power_scan(args)
    print()
    print(
        "WARNING: the Theorem-4.3 result is a necessary parameter floor, "
        "not a final security certification."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
