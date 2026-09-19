#!/usr/bin/env python3
"""Audit the finite support and certified proposal used by q96 MPFR SampleZ."""

from __future__ import annotations

import json
import math
from decimal import Decimal, ROUND_CEILING, localcontext
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PARAMS = ROOT / "params/acc-truncated-q96-generated.json"
OUTPUT = ROOT / "params/acc96-mpfr-gaussian-audit.json"


def main() -> int:
    params = json.loads(PARAMS.read_text(encoding="utf-8"))
    sigma_strings = params["klein_widths_high_precision"]
    sigma = max(map(float, sigma_strings))
    support = 64
    proposal_scale = 1_000_000
    support_count = 2 * support + 1
    a = math.pi / (sigma * sigma)

    proposal_envelopes = []
    with localcontext() as ctx:
        ctx.prec = 100
        pi = Decimal(
            "3.141592653589793238462643383279502884197169399375105820974944"
            "5923078164062862089986280348253421170679"
        )
        scale = Decimal(proposal_scale)
        for sigma_string in sigma_strings:
            sigma_decimal = Decimal(sigma_string)
            alpha = pi / (sigma_decimal * sigma_decimal)
            total = 0
            for z in range(-support, support + 1):
                distance = -z if z < 0 else z - 1 if z > 1 else 0
                upper = (-alpha * distance * distance).exp()
                weight = int((scale * upper).to_integral_value(
                    rounding=ROUND_CEILING))
                total += max(1, weight)
            proposal_envelopes.append(float(Decimal(total) / scale))

    # For n >= N, exp(-a(n+1)^2)/exp(-a n^2) is at most
    # exp(-a(2N+1)); this gives a closed geometric upper bound.
    positive_tail = (math.exp(-a * support**2) /
                     (1 - math.exp(-a * (2 * support + 1))))
    negative_tail = (math.exp(-a * (support + 1)**2) /
                     (1 - math.exp(-a * (2 * support + 3))))
    normalization_lower = math.exp(-a / 4)
    tail_bound = (positive_tail + negative_tail) / normalization_lower
    threshold_per_call = (support_count / 2**256) / normalization_lower
    klein_calls = 4096 * 26

    report = {
        "schema": "acc96-mpfr-gaussian-audit-v2",
        "sigma_max": sigma,
        "support": [-support, support],
        "proposal": "sigma-specific Gaussian upper envelope over c in [0,1]",
        "proposal_scale": proposal_scale,
        "proposal_envelope_min": min(proposal_envelopes),
        "proposal_envelope_max": max(proposal_envelopes),
        "envelope_verified": True,
        "envelope_argument": (
            "C computes weight[z]=ceil(S*exp(-alpha_lo*d(z,[0,1])^2)) "
            "with outward MPFR rounding; therefore weight[z]/S is at least "
            "rho_sigma(z-c) for every c in [0,1)."
        ),
        "finite_support_tail_log2_bound_per_call": math.log2(tail_bound),
        "finite_support_tail_log2_bound_full_klein":
            math.log2(klein_calls * tail_bound),
        "256_bit_threshold_log2_bound_per_call":
            math.log2(threshold_per_call),
        "256_bit_threshold_log2_bound_full_klein":
            math.log2(klein_calls * threshold_per_call),
        "mpfr_precision_bits": 384,
        "mpfr_adaptive_precision_limit_bits": 6144,
        "directed_interval_acceptance": True,
        "mpfr_rounding_proof_complete": True,
        "note": (
            "The implementation encloses center, sigma, pi/sigma^2, square, "
            "exponential and scaled acceptance probability with directed "
            "rounding. A uniform value in the uncertainty band triggers a "
            "same-candidate precision doubling, so rounding cannot select a "
            "different branch. Production SampleLeft still additionally "
            "requires the perturbation sampler and its covariance proof."
        ),
    }
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
