#!/usr/bin/env python3
"""Boundary validation for the actual mixed radix encoder; no private input."""
import numpy as np
import revac_acc96_full_ring as compiler

q = compiler.PROOF_MODULUS
b = 3779
public_top = (q // 2 + (b*b-1)//2)//(b*b)
# Include both field endpoints and all balanced-digit transition neighborhoods.
points = {0, -q//2 + 1, q//2}
for scale in (b, b*b):
    for k in range(-min(q//(2*scale), 100), min(q//(2*scale),100)+1):
        for delta in (-scale//2-1, -scale//2, -scale//2+1, 0, scale//2):
            value=k*scale+delta
            if abs(value)<=q//2: points.add(value)
rng = np.random.default_rng(4096)
values=np.concatenate((np.array(sorted(points),dtype=np.int64),
                      rng.integers(-q//2+1,q//2+1,size=100000,dtype=np.int64)))
digits=compiler._aggregate_mixed_limbs(values,b)
assert np.array_equal(values,digits[0]+b*digits[1]+b*b*digits[2])
assert all(np.max(np.abs(d))<=b//2 for d in digits[:2])
assert np.max(np.abs(digits[2]))<=public_top<=16383
print(f"PASS {len(values)} exact reconstructions; public top bound {public_top}")
