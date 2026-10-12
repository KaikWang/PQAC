"""Parameter certificate for the isolated single-88-bit experiment.

This file is not imported by the production CRT compiler.  The modulus is
selected by the same NTT condition needed by the degree-4096 negacyclic ring:
``p = 1 (mod 8192)``.  The Pocklington certificate below checks primality and
the NTT congruence at runtime; it is deliberately separate from the q88 wire
manifest and does not change production parameters.
"""

from __future__ import annotations

from math import gcd

SINGLE88_PRIME = 309333912899981927159480321
SINGLE88_RING_DEGREE = 4096
SINGLE88_NTT_ORDER = 2 * SINGLE88_RING_DEGREE
_FACTORS = ((2, 14, 3), (5, 1, 3), (7, 1, 2),
            (57399227, 1, 2), (9397958172739, 1, 2))


def _prime_64(n: int) -> bool:
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while not d & 1:
        d //= 2
        s += 1
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def verify_modulus_certificate() -> None:
    p = SINGLE88_PRIME
    if p.bit_length() != 88 or p % SINGLE88_NTT_ORDER != 1:
        raise ValueError("single q88 modulus fails 88-bit/NTT requirements")
    factored = 1
    for factor, exponent, witness in _FACTORS:
        if not _prime_64(factor):
            raise ValueError("single q88 factor primality check failed")
        factored *= factor ** exponent
        if pow(witness, p - 1, p) != 1:
            raise ValueError("single q88 Pocklington exponent check failed")
        if gcd(pow(witness, (p - 1) // factor, p) - 1, p) != 1:
            raise ValueError("single q88 primitive-factor check failed")
    if factored != p - 1:
        raise ValueError("single q88 factorization certificate mismatch")


verify_modulus_certificate()
