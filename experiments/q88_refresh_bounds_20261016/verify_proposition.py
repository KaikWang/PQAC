#!/usr/bin/env python3
"""Integer certificates for the numerical inequalities in the proposition.

This checks arithmetic only, not cryptographic security or a formal proof.
"""
import json
from math import isqrt
from pathlib import Path

from audit_bounds import audit


def main():
    a = audit()
    n, ell, k, d, fresh, gamma = 4096, 24, 6, 13003, (1 << 30)-1, 50000
    p = k*n*d*d
    b2 = 24*n*fresh*fresh + 150*n + 129
    c2 = k*n*((ell-1)*p*p+d*d)
    assert b2 == a["fresh_norm_squared_exact_under_condition"]
    assert c2 == a["delta_norm_squared_exact"]
    bceil = isqrt(b2)+1
    cceil = isqrt(c2)+1
    assert (bceil-1)**2 < b2 <= bceil**2
    assert (cceil-1)**2 < c2 <= cceil**2
    upper = bceil+gamma*cceil
    w = fresh+gamma*p
    capacity = (4095**5-1)//2
    assert upper < 2**68
    assert upper**2 < a["runtime_beta_squared_exact"] == 2**153
    assert 1+gamma*p <= w < capacity < a["q"]//2
    certificate = {
        "B_fresh_squared": b2,
        "C_delta_squared": c2,
        "ceil_B_fresh": bceil,
        "ceil_C_delta": cceil,
        "integer_norm_upper_bound_at_50000": upper,
        "norm_squared_bound_strictly_below_2_pow_153": True,
        "coefficient_update_bound": p,
        "coefficient_encoding_gate_W": w,
        "right_coefficient_bound_at_50000": 1+gamma*p,
        "balanced_radix_4095_five_digit_capacity": capacity,
        "centered_modulus_boundary": a["q"]//2,
        "accepted_fresh_witnesses_only": True,
        "lifetime_security_claim": False,
    }
    out = Path(__file__).with_name("proposition_certificate.json")
    out.write_text(json.dumps(certificate, indent=2)+"\n")
    print("PROPOSITION_INTEGER_INEQUALITIES_PASS")


if __name__ == "__main__":
    main()
