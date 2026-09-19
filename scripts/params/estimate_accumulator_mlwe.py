#!/usr/bin/env python3
"""Run the official lattice-estimator on the two PKC'26 M-LWE instances.

Theorem 4.3 uses M-LWE_(m-n,m,q,U(S_inf_1)) and
M-LWE_(m-n,m,q,U(B_1)).  Ring/module instances are flattened coefficient-wise
for this concrete estimator screen.  This is intentionally reported separately
from the theorem/reduction audit.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


def serialise_cost(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return str(value)


def log2_cost(value) -> float:
    # Sage RealField values retain exponents far beyond IEEE double range.
    # Convert only after taking the logarithm.
    try:
        return float(value.log(2))
    except (AttributeError, TypeError, ValueError):
        return math.log2(float(value))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--estimator", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--q", type=int, default=1208925396278498255577089)
    parser.add_argument("--ring-degree", type=int, default=4096)
    parser.add_argument("--module-rows", type=int, default=1)
    parser.add_argument("--trapdoor-cols", type=int, default=22)
    args = parser.parse_args()

    sys.path.insert(0, str(args.estimator.resolve()))
    from estimator import LWE, ND  # pylint: disable=import-error,import-outside-toplevel

    module_secret = args.trapdoor_cols - args.module_rows
    flat_secret = module_secret * args.ring_degree
    flat_samples = args.trapdoor_cols * args.ring_degree
    distributions = {
        "S_inf_1": ND.Uniform(-1, 1),
        "B_1": ND.Uniform(0, 1),
    }
    instances = {}
    for name, distribution in distributions.items():
        result = LWE.estimate.rough(
            LWE.Parameters(
                n=flat_secret,
                q=args.q,
                Xs=distribution,
                Xe=distribution,
                m=flat_samples,
            )
        )
        attacks = {
            key: {
                "rop_log2": log2_cost(cost["rop"]),
                "beta": serialise_cost(cost.get("beta", "n/a")),
            }
            for key, cost in result.items()
        }
        best_name, best = min(attacks.items(), key=lambda item: item[1]["rop_log2"])
        instances[name] = {
            "distribution": str(distribution),
            "attacks": attacks,
            "best_attack": best_name,
            "best_rop_log2": best["rop_log2"],
            "best_beta": best["beta"],
        }

    report = {
        "schema": "pqac-accumulator-mlwe-screen-v1",
        "warning": "rough lattice-estimator screen, not an end-to-end security claim",
        "q": args.q,
        "ring_degree": args.ring_degree,
        "module_rows": args.module_rows,
        "trapdoor_cols": args.trapdoor_cols,
        "flattened_secret_dimension": flat_secret,
        "flattened_samples": flat_samples,
        "instances": instances,
        "release_ready": False,
    }
    rendered = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
