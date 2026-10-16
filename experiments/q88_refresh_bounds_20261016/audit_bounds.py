#!/usr/bin/env python3
"""Exact arithmetic audit of witness-age bounds; never a security certificate.

Run from any directory. Production source and parameters are read only.
The fresh coefficient cap is checked at production admission; this script
reports enforcement locations but does not certify the sampler's tail law.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from decimal import Decimal, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PI = Decimal("3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679")


def norm_budget_holds(age: int, fresh_sq: int, delta_sq: int, beta_sq: int) -> bool:
    # sqrt(fresh_sq) + age*sqrt(delta_sq) <= sqrt(beta_sq), exactly.
    slack = beta_sq - fresh_sq - age * age * delta_sq
    return slack >= 0 and slack * slack >= 4 * fresh_sq * age * age * delta_sq


def exact_capacity(fresh_sq: int, delta_sq: int, beta_sq: int) -> int:
    if fresh_sq > beta_sq:
        return -1
    lo, hi = 0, 1
    while norm_budget_holds(hi, fresh_sq, delta_sq, beta_sq):
        hi *= 2
    while lo + 1 < hi:
        middle = (lo + hi) // 2
        if norm_budget_holds(middle, fresh_sq, delta_sq, beta_sq):
            lo = middle
        else:
            hi = middle
    return lo


def audit() -> dict:
    files = {
        "sampler": ROOT / "params/acc-truncated-q88-generated.json",
        "profile": ROOT / "params/revac-q88-window50000-candidate-v2.json",
        "header": ROOT / "code/entities/common/core/include/acc96_truncated_params_generated.h",
        "eval": ROOT / "code/entities/common/core/src/acc80_eval.c",
        "layout": ROOT / "code/entities/common/core/include/acc96_proof_layout.h",
        "sampleleft": ROOT / "code/entities/common/core/src/acc96_sampleleft.c",
        "relation": ROOT / "code/entities/common/core/src/acc96_relation.c",
        "relation_header": ROOT / "code/entities/common/core/include/acc96_relation.h",
        "admission": ROOT / "code/entities/common/core/src/acc96_accumulator.c",
        "ta": ROOT / "code/entities/ta/revac_ta.c",
        "sync": ROOT / "code/entities/common/revac_sync.c",
    }
    p = json.loads(files["sampler"].read_text())
    profile = json.loads(files["profile"].read_text())
    header = files["header"].read_text()
    norm_exp = int(re.search(r"#define ACC96_PUBLIC_WITNESS_NORM_SQ_BOUND 0x1p(\d+)L", header)[1])
    # Freeze the dimensions used in the accompanying source-level proof.
    assert (p["degree"], p["rank"], p["base"], p["k"], p["trunc_l"]) == (4096, 1, 13, 24, 2)
    assert profile["accumulator"]["gadget"]["digits"] == 6
    assert profile["accumulator"]["gadget"]["digit_bound"] == 13003
    assert profile["accumulator"]["handle_bits"] == 24
    n, k, ell, digit, bits = 4096, 6, 24, 13003, 128
    fresh_bits = int(re.search(r"#define ACC96_FRESH_COEFF_BOUND .*UINT64_C\(1\) << (\d+)", files["relation_header"].read_text())[1])
    fresh_cap = (1 << fresh_bits) - 1
    assert fresh_bits == 30
    per_coeff = k * n * digit**2
    # At most ell-1 convolution layers; final layer is G^-1(u), not a product.
    delta_sq = k * n * ((ell - 1) * per_coeff**2 + digit**2)
    # 24 Gaussian head polynomials, 144+6 ternary polynomials, tag=1, r bits.
    fresh_sq = 24 * n * fresh_cap**2 + (ell * k + k) * n + 1 + bits
    beta_sq = 1 << norm_exp
    proof_age = profile["window"]["runtime_deletion_limit"]
    wire_cap = fresh_cap + proof_age * per_coeff
    capacity = exact_capacity(fresh_sq, delta_sq, beta_sq)
    assert norm_budget_holds(capacity, fresh_sq, delta_sq, beta_sq)
    assert not norm_budget_holds(capacity + 1, fresh_sq, delta_sq, beta_sq)
    # Reproduce the old screen and expose its theorem identity explicitly.
    sys.path.insert(0, str(ROOT / "scripts/params"))
    from audit_accumulator_security import audit as old_audit
    historical = old_audit(p["q"], ring_degree=n, module_rows=1,
        deletion_bound=proof_age, handle_bits=ell, eval_base=26007,
        eval_digit_bound=digit, trunc_base=13, trunc_l=2,
        spectral_slack=32, q_factors=(1073692673, 536813569, 536690689))
    # Appendix E formulas with a formal dimension mapping only. This does NOT
    # establish that structured truncated/ternary/CRT sampling meets its games.
    m = 24
    theta1 = math.sqrt(2*n/3) * (math.sqrt(m)+math.sqrt(k)+math.log(math.sqrt(n)))
    theta2 = math.sqrt(2*n/3) * (math.sqrt(m)+math.sqrt(ell*k)+math.log(math.sqrt(n)))
    factor_e = math.sqrt(1+theta2**2+(bits+1)**2*theta1**2)
    e_rows = []
    for lifetime in (50000, 100000):
        t = lifetime * math.sqrt(delta_sq)
        # Set s=T as a diagnostic, not as a production sampler parameter.
        s = t
        beta0 = s*math.sqrt(n*(m+k*(ell+1)))+t
        e1 = 2*math.sqrt(m*n)+theta2*t+beta0*factor_e
        e2 = ((s+1)*math.sqrt(m*n)+theta2*t
              +s*math.sqrt(n*k*(ell+1))*math.sqrt(factor_e**2-1)
              +beta0*factor_e)
        e_rows.append({"lifetime_deletions": lifetime,
            "log2_T": math.log2(t), "log2_beta1_case1": math.log2(e1),
            "log2_beta1_case2": math.log2(e2),
            "log2_q": math.log2(p["q"]),
            "case2_s_ge_T_for_production_max_width": p["s_11"] >= t})
    with localcontext() as ctx:
        ctx.prec = 90
        delta = Decimal(delta_sq).sqrt()
        fresh = Decimal(fresh_sq).sqrt()
        beta = Decimal(beta_sq).sqrt()
        # Conditional diagnostic ONLY: assumes a uniform coordinate tail for
        # the *output* head, which covariance construction alone does not prove.
        gaussian_coordinates = 24 * n
        tail_u = ((Decimal(128) * Decimal(2).ln() + Decimal(2 * gaussian_coordinates).ln()) / PI).sqrt()
        rho = Decimal(p["sampler_scalars_high_precision"]["spectral_norm_bound_sq"]).sqrt()
        factor = ((Decimal(2 * n * p["k"]).ln() + Decimal(128) * Decimal(2).ln()) / PI).sqrt() * (Decimal(13) + Decimal(1) / 13)
        sigma_max = factor * (Decimal(4) * Decimal(13**4 - 1) / (13**2 - 1) + 3 * rho**2).sqrt()
        tail_cap = int((tail_u * sigma_max).to_integral_value(rounding="ROUND_CEILING"))
        tail_fresh_sq = 24 * n * tail_cap**2 + (ell * k + k) * n + 1 + bits
        ages = []
        for age in (0, 50000, 56245, 100000):
            bound = fresh + age * delta
            ages.append({
                "age_since_last_refresh": age,
                "norm_upper_bound": str(bound),
                "log2_norm_upper_bound": float(bound.ln() / Decimal(2).ln()),
                "runtime_l2_budget_holds": norm_budget_holds(age, fresh_sq, delta_sq, beta_sq),
                "right_coefficient_upper_bound": 1 + age * per_coeff,
                "proof_coefficient_budget_holds": max(fresh_cap, 1 + age * per_coeff) <= wire_cap,
            })
        return {
            "status": "conditional_correctness_bound_and_security_reduction_audit_only",
            "security_certified": False,
            "production_changed": False,
            "q": p["q"], "ring_degree": n, "handle_bits": ell,
            "fresh_condition": "all 24 head polynomial coefficients <= 2^30-1 in absolute value; right/z ternary; tag=1",
            "fresh_condition_enforced_by_sampler": "ok = acc96_relation_fresh_witness_check(out, ctx);" in files["sampleleft"].read_text(),
            "fresh_gate_source_locations": [str(files[k].relative_to(ROOT)) for k in ("sampleleft", "admission", "ta", "sync")],
            "fresh_gate_failure_policy": "one attempt, explicit failure, no implicit resampling or clipping",
            "fresh_coefficient_cap": fresh_cap,
            "delta_coefficient_bound": per_coeff,
            "delta_norm_squared_exact": delta_sq,
            "delta_norm": str(delta),
            "fresh_norm_squared_exact_under_condition": fresh_sq,
            "fresh_norm_under_condition": str(fresh),
            "runtime_beta_squared_exact": beta_sq,
            "runtime_beta": str(beta),
            "norm_exponents": {
                "runtime_header": norm_exp,
                "sampler_json": p["public_witness_norm_sq_exp"],
                "profile_manifest": profile["window"]["public_witness_norm_sq_bound_log2"],
                "consistent": norm_exp == p["public_witness_norm_sq_exp"],
            },
            "l2_only_capacity_under_fresh_condition": capacity,
            "l2_only_capacity_using_sampler_json_gate": exact_capacity(fresh_sq, delta_sq, 1 << p["public_witness_norm_sq_exp"]),
            "proof_coefficient_cap": wire_cap,
            "proof_coefficient_only_capacity_under_fresh_condition": (wire_cap - 1) // per_coeff,
            "current_runtime_deletion_counter_cap": proof_age,
            "ages": ages,
            "conditional_sampler_tail_diagnostic": {
                "certified": False,
                "caution": "JS25 Theorem 4.2 gives a K-pushforward distribution, not automatically a Gaussian on the collapsed output; this tail needs a separate argument",
                "assumption": "uniformly for every residual/history, head coordinate tail <= 2 exp(-pi*u^2) with output width <= sigma_max; implementation distribution error added separately",
                "sigma_max": str(sigma_max),
                "head_coordinates": gaussian_coordinates,
                "tail_u_for_union_bound_2_neg128": str(tail_u),
                "head_coefficient_cap": tail_cap,
                "fresh_norm": str(Decimal(tail_fresh_sq).sqrt()),
                "fresh_gate_tail_exponent_without_implementation_error": str(PI * (Decimal(fresh_cap) / sigma_max)**2 / Decimal(2).ln() - Decimal(2 * gaussian_coordinates).ln() / Decimal(2).ln()),
            },
            "historical_screen": {
                "theorem": "Theorem 4.3, not randomized Appendix E",
                "log2_beta0": historical.beta0_log2,
                "log2_beta1": historical.beta1_log2,
                "log2_q_over_beta1": historical.beta1_margin_log2,
                "release_ready": historical.release_ready,
            },
            "appendix_e_formal_dimension_diagnostic": {
                "valid_instantiation_proved": False,
                "s_choice": "T for arithmetic diagnostic only; not actual sampler widths",
                "theta1": theta1, "theta2": theta2, "extraction_factor": factor_e,
                "rows": e_rows,
            },
            "source_sha256": {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in files.items()},
        }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(audit(), indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.write_text(result)
    else:
        print(result, end="")
