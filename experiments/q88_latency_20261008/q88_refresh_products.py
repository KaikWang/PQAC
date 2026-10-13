"""Bind each refresh radix digit to its refresh-z product in one Pack.

The auxiliary product is reconstructed from three radix-128 layers. All
bounds and layer counts are public profile constants. The quadratic equality
has no wrap in the 38-bit field: Cauchy-Schwarz bounds a convolution
coefficient by n*64 for ||z||^2 <= n and ||a_digit||^2 <= n*64^2;
the reconstructed product coefficient is at most
sqrt(n)*64*(1+128+128^2). Their sum is < PACK_MODULUS for n=4096.
"""
from math import isqrt
from q88_bridge_core import PACK_MODULUS, _decompose
from revac_acc96_wire import REFRESH_BITS, REFRESH_Z_FIRST_SLOT


def refresh_product_constraints(bridges, prover):
    from exact_negacyclic_fft import exact_sum_negacyclic
    witnesses, constraints = [], []
    n, base = bridges[0].degree, bridges[0].base
    if n != 4096 or base != 128:
        raise ValueError("refresh product layout requires degree-4096 radix-128")
    upper = n*64 + 64*isqrt(n)*(1+128+128**2)
    if upper >= PACK_MODULUS:
        raise ValueError("refresh product no-wrap bound failed")
    for block, bridge in enumerate(bridges):
        if len(bridge.w_digits) != REFRESH_BITS + 13:
            raise ValueError("refresh products require the compact bridge layout")
        z = tuple(int(v) for v in prover.witness_poly(REFRESH_Z_FIRST_SLOT + block))
        if any(abs(v) > 1 for v in z):
            raise ValueError("captured refresh-z is outside this ternary profile")
        zref = ("refresh-z", block)
        witnesses.append({"ref": zref, "values": [z], "norm": n})
        for digit in range(13):
            layers = bridge.w_digits[REFRESH_BITS + digit]
            if len(layers) != 1:
                raise ValueError("refresh digit must have exactly one radix layer")
            ahat = layers[0]
            # Small unsigned +/- components allow the existing integer FFT
            # constructor to build signed products. Native Clear Verify
            # independently checks every polynomial quadratic equality.
            product = exact_sum_negacyclic(
                [[max(v, 0) for v in ahat], [max(-v, 0) for v in ahat]],
                [list(z), [-v for v in z]], bits=7,
                left_count=1, right_count=1)
            linear = []
            for layer, values in enumerate(_decompose(product, base, 3, signed=True)):
                ref = ("refresh-product", block, digit, layer)
                witnesses.append({"ref": ref, "values": [tuple(values)],
                                  "norm": n*64**2})
                coeff = tuple([-(base**layer)] + [0]*(n-1))
                linear.append((coeff, ref))
            constraints.append({
                "linear": linear,
                "quadratic": [(1, (block, "w", 0, REFRESH_BITS+digit), zref)],
            })
    return witnesses, constraints
