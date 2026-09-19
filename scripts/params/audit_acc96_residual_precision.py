#!/usr/bin/env python3
"""Conservative precision budget for q96 Klein residual centers.

This bounds only the distribution shift caused by rounding c in D_{Z,s,c}; it
does not certify the internal arithmetic or randomness of a SampleZ backend.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", type=Path, required=True)
    ap.add_argument("--center-bits", type=int, default=53)
    ap.add_argument("--target-bits", type=int, default=128)
    ap.add_argument("--composition-margin-bits", type=int, default=32)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    data = json.loads(args.params.read_text(encoding="utf-8"))
    calls = int(data["degree"]) * int(data["rank"]) * int(data["k"])
    sigma_min = min(float(x) for x in data["klein_widths"])

    # For p_c(z) proportional to exp(-pi(z-c)^2/s^2), integration of the
    # score function gives TV(D_c,D_c') <= L |c-c'|.  The deliberately loose
    # E|Z-c| <= s+1 bound gives L <= 4*pi*(s+1)/s^2.
    lipschitz = 4 * math.pi * (sigma_min + 1) / (sigma_min * sigma_min)
    log2_composed_error = (
        math.log2(calls) + math.log2(lipschitz) - args.center_bits
    )
    required_for_target = math.ceil(
        args.target_bits + math.log2(calls * lipschitz)
    )
    required_with_margin = math.ceil(
        args.target_bits + args.composition_margin_bits
        + math.log2(calls * lipschitz)
    )
    recommended = 64 * math.ceil(required_with_margin / 64)
    report = {
        "scope": "Klein residual-center rounding only",
        "klein_calls_per_sampleleft": calls,
        "sigma_min": sigma_min,
        "center_lipschitz_upper_bound": lipschitz,
        "tested_center_precision_bits": args.center_bits,
        "composed_rounding_error_log2_upper_bound": log2_composed_error,
        "target_statistical_bits": args.target_bits,
        "required_center_bits_for_target": required_for_target,
        "composition_margin_bits": args.composition_margin_bits,
        "required_center_bits_with_margin": required_with_margin,
        "recommended_backend_precision_bits": recommended,
        "double_center_is_sufficient":
            log2_composed_error <= -args.target_bits,
        "limitations": [
            "does not cover probability-evaluation rounding inside SampleZ",
            "does not cover PRNG bias or perturbation sampler arithmetic",
            "uses a conservative score-function Lipschitz bound",
        ],
    }
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
