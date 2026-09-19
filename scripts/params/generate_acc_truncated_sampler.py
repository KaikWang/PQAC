#!/usr/bin/env python3
"""Generate accumulator-only truncated-gadget sampler constants.

This deliberately does not import the anonymous-signature parameter module.
The formulas are the ring specialization of Jeudy--Sanders, Theorem 4.2
(IACR ePrint 2024/1952).  The generated Klein tables are derived from an
integer basis of Lambda_q^perp(g), rather than copied from another modulus.
"""

from __future__ import annotations

import argparse
import json
import math
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path


def ceil_log(q: int, base: int) -> int:
    value, digits = 1, 0
    while value < q:
        value *= base
        digits += 1
    return digits


def base_digits(value: int, base: int, count: int) -> list[int]:
    out = []
    for _ in range(count):
        out.append(value % base)
        value //= base
    if value:
        raise ValueError("gadget dimension does not cover q")
    return out


def dot(a: list[Fraction], b: list[Fraction]) -> Fraction:
    return sum((x * y for x, y in zip(a, b)), Fraction(0))


def derive(q: int, degree: int, rank: int, base: int, trunc_l: int,
           tag_norm: int, statistical_bits: int, spectral_slack: float) -> dict:
    k = ceil_log(q, base)
    if not 0 < trunc_l < k:
        raise ValueError("need 0 < trunc_l < ceil(log_base(q))")

    # The old generator evaluated security-critical Gaussian widths in binary64.
    # Keep binary64 compatibility fields, but derive a 256-bit source of truth
    # for the production q96 sampler and its finite-precision error proof.
    with localcontext() as decimal_ctx:
        decimal_ctx.prec = 100
        d2 = Decimal(2)
        pi = Decimal(
            "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679"
        )
        eta_hp = (((Decimal(2 * degree * rank * k).ln() +
                    Decimal(statistical_bits) * d2.ln()) / pi).sqrt())
        norm_r_hp = Decimal(3) / Decimal(4) * (
            Decimal(rank * degree).sqrt()
            + Decimal((k - trunc_l) * rank * degree).sqrt()
            + Decimal(str(spectral_slack))
        )
        factor_hp = eta_hp * (Decimal(base) + Decimal(1) / Decimal(base))
        s_g_hp = eta_hp * Decimal(base * base + 1).sqrt()
        s_1_hp = factor_hp * (
            Decimal(4 * tag_norm * tag_norm)
            + Decimal(3) * norm_r_hp * norm_r_hp).sqrt()
        s_2_hp = factor_hp * Decimal(2 * tag_norm)
        s_3_hp = factor_hp * Decimal(3).sqrt() * norm_r_hp
        s_4_hp = factor_hp * Decimal(3).sqrt()
        geometric_hp = (Decimal(base ** (2 * trunc_l) - 1) /
                        Decimal(base * base - 1))
        s_11_hp = factor_hp * (
            geometric_hp * Decimal(4 * tag_norm * tag_norm)
            + Decimal(3) * norm_r_hp * norm_r_hp).sqrt()

    eta = float(eta_hp)
    norm_r = float(norm_r_hp)
    factor = float(factor_hp)
    s_g = float(s_g_hp)
    s_1 = float(s_1_hp)
    s_2 = float(s_2_hp)
    s_3 = float(s_3_hp)
    s_4 = float(s_4_hp)
    s_11 = float(s_11_hp)

    # Column basis [b*I with -1 below diagonal | digits_base(q)].
    q_digits = base_digits(q, base, k)
    basis = [[0 for _ in range(k)] for _ in range(k)]
    for col in range(k - 1):
        basis[col][col] = base
        basis[col + 1][col] = -1
    for row, digit in enumerate(q_digits):
        basis[row][k - 1] = digit

    columns = [[Fraction(basis[row][col]) for row in range(k)]
               for col in range(k)]
    gso: list[list[Fraction]] = []
    for column in columns:
        orthogonal = list(column)
        for previous in gso:
            scale = dot(orthogonal, previous) / dot(previous, previous)
            orthogonal = [x - scale * y
                          for x, y in zip(orthogonal, previous)]
        gso.append(orthogonal)
    neg_scaled = [[0.0 for _ in range(k)] for _ in range(k)]
    neg_scaled_exact = [[{"numerator": 0, "denominator": 1}
                         for _ in range(k)] for _ in range(k)]
    gso_norms = []
    gso_norms_hp = []
    for col, vector in enumerate(gso):
        norm_sq = dot(vector, vector)
        gso_norms.append(math.sqrt(float(norm_sq)))
        with localcontext() as decimal_ctx:
            decimal_ctx.prec = 100
            gso_norms_hp.append(
                (Decimal(norm_sq.numerator) /
                 Decimal(norm_sq.denominator)).sqrt())
        for row, value in enumerate(vector):
            exact = -value / norm_sq
            neg_scaled[row][col] = float(exact)
            neg_scaled_exact[row][col] = {
                "numerator": exact.numerator,
                "denominator": exact.denominator,
            }
    klein_widths = [s_g / value for value in gso_norms]
    with localcontext() as decimal_ctx:
        decimal_ctx.prec = 100
        klein_widths_hp = [s_g_hp / value for value in gso_norms_hp]
    exact_column_denominators = []
    exact_column_numerators = [[0 for _ in range(k)] for _ in range(k)]
    for col in range(k):
        fractions = [Fraction(
            neg_scaled_exact[row][col]["numerator"],
            neg_scaled_exact[row][col]["denominator"])
            for row in range(k)]
        common_denominator = math.lcm(*(x.denominator for x in fractions))
        exact_column_denominators.append(common_denominator)
        for row, value in enumerate(fractions):
            exact_column_numerators[row][col] = (
                value.numerator * (common_denominator // value.denominator)
            )

    with localcontext() as decimal_ctx:
        decimal_ctx.prec = 100
        s4sq_minus_sgsq_hp = s_4_hp * s_4_hp - s_g_hp * s_g_hp
        sampler_scalars_hp = {
            "spectral_norm_bound_sq": str(norm_r_hp * norm_r_hp),
            "s1_sq": str(s_1_hp * s_1_hp),
            "s2_sq": str(s_2_hp * s_2_hp),
            "s3_sq": str(s_3_hp * s_3_hp),
            "sqrt_s4sq_minus_sgsq": str(s4sq_minus_sgsq_hp.sqrt()),
            "neg_inv_sgsq_minus_inv_s4sq": str(
                -Decimal(1) /
                (Decimal(1) / (s_g_hp * s_g_hp)
                 - Decimal(1) / (s_4_hp * s_4_hp))),
            "neg_sgsq_div_s4sq_minus_sgsq": str(
                -(s_g_hp * s_g_hp) / s4sq_minus_sgsq_hp),
            "neg_sgsq": str(-(s_g_hp * s_g_hp)),
        }
    s4sq_minus_sgsq = float(s4sq_minus_sgsq_hp)
    if s4sq_minus_sgsq <= 0:
        raise ValueError("invalid covariance: s4^2 <= sG^2")

    return {
        "q": q,
        "degree": degree,
        "rank": rank,
        "base": base,
        "k": k,
        "trunc_l": trunc_l,
        "tag_norm_bound": tag_norm,
        "statistical_bits": statistical_bits,
        "spectral_slack": spectral_slack,
        "spectral_norm_bound_sq": norm_r * norm_r,
        "eta": eta,
        "s_g": s_g,
        "s_1": s_1,
        "s_2": s_2,
        "s_3": s_3,
        "s_4": s_4,
        "s_11": s_11,
        "s1_sq": s_1 * s_1,
        "s2_sq": s_2 * s_2,
        "s3_sq": s_3 * s_3,
        "sqrt_s4sq_minus_sgsq": math.sqrt(s4sq_minus_sgsq),
        "neg_inv_sgsq_minus_inv_s4sq": -1 / (1 / (s_g * s_g) - 1 / (s_4 * s_4)),
        "neg_sgsq_div_s4sq_minus_sgsq": -(s_g * s_g) / s4sq_minus_sgsq,
        "neg_sgsq": -(s_g * s_g),
        "sampler_scalars_high_precision": sampler_scalars_hp,
        "q_base_digits": q_digits,
        "gadget_basis": basis,
        "klein_widths": klein_widths,
        "high_precision_decimal_digits": 100,
        "klein_widths_high_precision": [str(x) for x in klein_widths_hp],
        "neg_scaled_gadget_gso": neg_scaled,
        "neg_scaled_gadget_gso_exact": neg_scaled_exact,
        "neg_scaled_gadget_gso_column_denominators":
            exact_column_denominators,
        "neg_scaled_gadget_gso_column_numerators": exact_column_numerators,
        "basis_kernel_check": all(
            sum((base ** row) * basis[row][col] for row in range(k)) % q == 0
            for col in range(k)
        ),
    }


def c_float(value: float) -> str:
    return f"{value:.20f}"


def emit_header(p: dict) -> str:
    k = p["k"]
    q = p["q"]
    q_bits = q.bit_length()
    q_lo = q & ((1 << 64) - 1)
    q_hi = q >> 64
    q_legacy = (f"#define ACC_TRUNC_PARAM_Q UINT64_C({q})\n"
                if q_bits <= 64 else
                "/* ACC_TRUNC_PARAM_Q is intentionally unavailable: the "
                "modulus does not fit uint64_t. */\n")
    rows_f = lambda rows: ",\n".join(
        "  {" + ", ".join(c_float(x) for x in row) + "}" for row in rows
    )
    rows_i = lambda rows: ",\n".join(
        "  {" + ", ".join(str(x) for x in row) + "}" for row in rows
    )
    column_den = [str(x) for x in
                  p["neg_scaled_gadget_gso_column_denominators"]]
    widths_hp_c = ", ".join(
        '"' + x + '"' for x in p["klein_widths_high_precision"])
    column_den_c = ", ".join('"' + x + '"' for x in column_den)
    hp = p["sampler_scalars_high_precision"]
    witness_norm = p.get("public_witness_norm_sq_exp")
    witness_norm_c = (
        f"#define ACC96_PUBLIC_WITNESS_NORM_SQ_BOUND 0x1p{witness_norm}L\n"
        if witness_norm is not None else "")
    column_num = [[str(x) for x in row] for row in
                  p["neg_scaled_gadget_gso_column_numerators"]]
    rows_s = lambda rows: ",\n".join(
        "  {" + ", ".join(f'\"{x}\"' for x in row) + "}" for row in rows
    )
    return f"""/* Generated by scripts/params/generate_acc_truncated_sampler.py.
 * Source formulas: Jeudy--Sanders, ePrint 2024/1952, Theorem 4.2.
 * Do not replace these with anonymous-signature sampler constants. */
#ifndef REVAC_ACC_TRUNCATED_PARAMS_Q{q_bits}_H
#define REVAC_ACC_TRUNCATED_PARAMS_Q{q_bits}_H

#include <stdint.h>

/* Exact, width-safe encoding.  Q = (Q_HI << 64) | Q_LO.  Keeping the
 * split form prevents a >64-bit literal from being silently truncated by C. */
#define ACC_TRUNC_PARAM_Q_BITS {q_bits}
#define ACC_TRUNC_PARAM_Q_LO UINT64_C({q_lo})
#define ACC_TRUNC_PARAM_Q_HI UINT64_C({q_hi})
{q_legacy.rstrip()}
#define ACC_TRUNC_PARAM_N {p['degree']}
#define ACC_TRUNC_PARAM_D {p['rank']}
#define ACC_TRUNC_PARAM_B {p['base']}
#define ACC_TRUNC_PARAM_K {k}
#define ACC_TRUNC_PARAM_L {p['trunc_l']}
#define ACC_TRUNC_PARAM_TAG_NORM {p['tag_norm_bound']}
#define ACC_TRUNC_PARAM_STATISTICAL_BITS {p['statistical_bits']}
#define ACC_TRUNC_R_MAX_SQ_SPECTRAL_NORM {c_float(p['spectral_norm_bound_sq'])}
#define ACC_TRUNC_S1SQ {c_float(p['s1_sq'])}
#define ACC_TRUNC_S2SQ {c_float(p['s2_sq'])}
#define ACC_TRUNC_S3SQ {c_float(p['s3_sq'])}
#define ACC_TRUNC_SQRT_S4SQ_SGSQ {c_float(p['sqrt_s4sq_minus_sgsq'])}
#define ACC_TRUNC_SGINVSQ_S4INVSQ {c_float(p['neg_inv_sgsq_minus_inv_s4sq'])}
#define ACC_TRUNC_NEGSGSQ_DIV_S4SQ_SGSQ {c_float(p['neg_sgsq_div_s4sq_minus_sgsq'])}
#define ACC_TRUNC_NEG_SGSQ {c_float(p['neg_sgsq'])}
{witness_norm_c.rstrip()}

/* High-precision source strings for perturbation covariance construction. */
#define ACC_TRUNC_R_MAX_SQ_SPECTRAL_NORM_HP "{hp['spectral_norm_bound_sq']}"
#define ACC_TRUNC_S1SQ_HP "{hp['s1_sq']}"
#define ACC_TRUNC_S2SQ_HP "{hp['s2_sq']}"
#define ACC_TRUNC_S3SQ_HP "{hp['s3_sq']}"
#define ACC_TRUNC_SQRT_S4SQ_SGSQ_HP "{hp['sqrt_s4sq_minus_sgsq']}"
#define ACC_TRUNC_SGINVSQ_S4INVSQ_HP "{hp['neg_inv_sgsq_minus_inv_s4sq']}"
#define ACC_TRUNC_NEGSGSQ_DIV_S4SQ_SGSQ_HP "{hp['neg_sgsq_div_s4sq_minus_sgsq']}"
#define ACC_TRUNC_NEG_SGSQ_HP "{hp['neg_sgsq']}"

static const double acc_trunc_klein_widths[ACC_TRUNC_PARAM_K] = {{
  {', '.join(c_float(x) for x in p['klein_widths'])}
}};

/* 100-decimal-digit values consumed by the >=192-bit residual sampler.
 * The double table above remains only for compatibility/KAT timing. */
#define ACC_TRUNC_KLEIN_WIDTH_DECIMAL_DIGITS {p['high_precision_decimal_digits']}
static const char *const
acc_trunc_klein_widths_high_precision[ACC_TRUNC_PARAM_K] = {{
  {widths_hp_c}
}};

static const double
acc_trunc_neg_scaled_gadget_gso[ACC_TRUNC_PARAM_K][ACC_TRUNC_PARAM_K] = {{
{rows_f(p['neg_scaled_gadget_gso'])}
}};

/* Faster equivalent representation: every column uses one positive common
 * denominator, turning center formation into a wide-integer dot product. */
static const char *const
acc_trunc_neg_scaled_gadget_gso_column_den[ACC_TRUNC_PARAM_K] = {{
  {column_den_c}
}};
static const char *const
acc_trunc_neg_scaled_gadget_gso_column_num[ACC_TRUNC_PARAM_K][ACC_TRUNC_PARAM_K] = {{
{rows_s(column_num)}
}};

static const int64_t
acc_trunc_gadget_basis[ACC_TRUNC_PARAM_K][ACC_TRUNC_PARAM_K] = {{
{rows_i(p['gadget_basis'])}
}};

#endif
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--degree", type=int, required=True)
    ap.add_argument("--rank", type=int, required=True)
    ap.add_argument("--base", type=int, required=True)
    ap.add_argument("--trunc-l", type=int, required=True)
    ap.add_argument("--tag-norm", type=int, default=1)
    ap.add_argument("--statistical-bits", type=int, default=128)
    ap.add_argument("--spectral-slack", type=float, default=6.0)
    ap.add_argument("--public-witness-norm-sq-exp", type=int)
    ap.add_argument("--header", type=Path)
    ap.add_argument("--json", type=Path)
    args = ap.parse_args()
    params = derive(args.q, args.degree, args.rank, args.base, args.trunc_l,
                    args.tag_norm, args.statistical_bits, args.spectral_slack)
    params["public_witness_norm_sq_exp"] = args.public_witness_norm_sq_exp
    if not params["basis_kernel_check"]:
        raise SystemExit("generated gadget basis is not in the q-ary kernel")
    if args.header:
        args.header.write_text(emit_header(params), encoding="utf-8")
    if args.json:
        args.json.write_text(json.dumps(params, indent=2) + "\n", encoding="utf-8")
    if not args.header and not args.json:
        print(json.dumps(params, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
