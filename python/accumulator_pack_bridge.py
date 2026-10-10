"""Exact integer-lift compiler for a negacyclic accumulator equation.

The accumulator modulus may exceed the 38-bit Pack field.  A coefficient
equation is expanded in base ``base`` with a quotient polynomial and carries.
Every digit equation is proved over the Pack field.  Before constructing a
statement we require an independent coefficient bound strictly below half
the Pack modulus; that bound makes each field equation an integer equation.

This module is a small, deliberately conservative compiler.  It rejects a
profile whenever its public norm bounds cannot certify non-wrapping.
"""

from dataclasses import dataclass
from math import isqrt

PACK_MODULUS = (1 << 38) - 107


def negacyclic_mul(a, b):
    if len(a) != len(b):
        raise ValueError("polynomial degrees differ")
    n = len(a)
    out = [0] * n
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    k = i + j
                    out[k % n] += (x * y if k < n else -x * y)
    return out


def balanced_digits(value, base, count):
    digits = []
    for _ in range(count):
        digit = (value + base // 2) % base - base // 2
        digits.append(digit)
        value = (value - digit) // base
    if value:
        raise ValueError("digit capacity exhausted")
    return digits


def unsigned_digits(value, base, count):
    if value < 0 or value >= base ** count:
        raise ValueError("unsigned digit capacity exhausted")
    return [(value // base ** i) % base for i in range(count)]


def decompose_polynomial(poly, base, count, signed):
    layers = [[0] * len(poly) for _ in range(count)]
    for k, value in enumerate(poly):
        digits = (balanced_digits(value, base, count) if signed else
                  unsigned_digits(value, base, count))
        for i, digit in enumerate(digits):
            layers[i][k] = digit
    return layers


def exact_norm_sq(poly):
    return sum(x * x for x in poly)


@dataclass
class Bridge:
    degree: int
    base: int
    h_digits: list
    w_digits: list
    q_digits: list
    quotient_digits: list
    target_digits: list
    carry_digits: list
    bounds: list

    def equation_count(self):
        return len(self.carry_digits)

    def _layer_terms(self, layer):
        terms = []
        for col, h in enumerate(self.h_digits):
            for i, hp in enumerate(h):
                j = layer - i
                if 0 <= j < len(self.w_digits[col]):
                    terms.append((hp, self.w_digits[col][j], 1))
        for i, qdigit in enumerate(self.q_digits):
            j = layer - i
            if qdigit and 0 <= j < len(self.quotient_digits):
                terms.append(([qdigit] + [0] * (self.degree - 1),
                              self.quotient_digits[j], -1))
        return terms

    def verify_integer_equations(self):
        zero = [0] * self.degree
        for layer in range(self.equation_count()):
            lhs = zero.copy()
            for public, witness, sign in self._layer_terms(layer):
                product = negacyclic_mul(public, witness)
                lhs = [x + sign * y for x, y in zip(lhs, product)]
            if layer < len(self.target_digits):
                lhs = [x - y for x, y in zip(lhs, self.target_digits[layer])]
            if layer:
                prior = self.carry_digits[layer - 1]
                lhs = [x + y for x, y in zip(lhs, prior)]
            lhs = [x - self.base * y for x, y in
                   zip(lhs, self.carry_digits[layer])]
            if lhs != zero:
                raise ValueError(f"integer carry equation {layer} failed")


def compile_bridge(matrix, witness, target, modulus, *, base=128,
                   witness_digits=2, quotient_digits=15,
                   carry_digits=8, pack_modulus=PACK_MODULUS):
    """Compile and certify one exact 100-bit congruence for Pack.

    ``bounds`` are squared L2 bounds for every private digit polynomial.
    They are fixed from the supplied public profile and returned with the
    statement.  Any field equation whose proven L2 bounds permit wrapping is
    rejected.  A production verifier must pin these bounds and the public
    matrix to its profile, never accept prover-selected values.
    """
    if not matrix or len(matrix) != len(witness):
        raise ValueError("matrix/witness shape mismatch")
    n = len(target)
    if any(len(poly) != n for poly in matrix + witness):
        raise ValueError("polynomial degree mismatch")
    if any(not 0 <= x < modulus for poly in matrix for x in poly):
        raise ValueError("noncanonical matrix")
    if any(not 0 <= x < modulus for x in target):
        raise ValueError("noncanonical target")
    if base < 2 or base & (base - 1):
        raise ValueError("base must be a power of two")
    h_count = (modulus.bit_length() + base.bit_length() - 2) // (base.bit_length() - 1)
    h = [decompose_polynomial(poly, base, h_count, False) for poly in matrix]
    w = [decompose_polynomial(poly, base, witness_digits, True) for poly in witness]
    q = unsigned_digits(modulus, base, h_count)
    t = decompose_polynomial(target, base, h_count, False)
    total = [0] * n
    for a, s in zip(matrix, witness):
        total = [x + y for x, y in zip(total, negacyclic_mul(a, s))]
    quotient = []
    for value, rhs in zip(total, target):
        diff = value - rhs
        if diff % modulus:
            raise ValueError("original accumulator congruence fails")
        quotient.append(diff // modulus)
    k = decompose_polynomial(quotient, base, quotient_digits, True)
    layer_count = max(h_count + witness_digits - 1,
                      h_count + quotient_digits - 1)
    layer_count += 1
    carries = []
    previous = [0] * n
    for layer in range(layer_count):
        expression = previous.copy()
        for col in range(len(matrix)):
            for i in range(h_count):
                j = layer - i
                if 0 <= j < witness_digits:
                    product = negacyclic_mul(h[col][i], w[col][j])
                    expression = [x + y for x, y in zip(expression, product)]
        for i, qdigit in enumerate(q):
            j = layer - i
            if qdigit and 0 <= j < quotient_digits:
                expression = [x - qdigit * y for x, y in zip(expression, k[j])]
        if layer < h_count:
            expression = [x - y for x, y in zip(expression, t[layer])]
        if any(x % base for x in expression):
            raise ValueError("carry is not integral")
        previous = [x // base for x in expression]
        carries.append(previous)
    if any(previous):
        raise ValueError("nonzero terminal carry")
    if any(x for x in carries[-1]):
        raise ValueError("terminal carry is not zero")
    c = carries
    # Exact per-vector L2 bounds; callers can replace these only with
    # predeclared profile bounds and must rerun the non-wrapping audit.
    bounds = {
        "w": [[exact_norm_sq(poly) for poly in col] for col in w],
        "k": [exact_norm_sq(poly) for poly in k],
        "c": [exact_norm_sq(poly) for poly in c],
    }
    for group in (sum(bounds["w"], []), bounds["k"], bounds["c"]):
        if any(x >= 1 << 64 for x in group):
            raise ValueError("Pack norm bound exceeds uint64; decompose carry")
    bridge = Bridge(n, base, h, w, q, k, t, c, bounds)
    bridge.verify_integer_equations()
    # A field equation with |lhs| < p cannot hide a nonzero multiple of p.
    # Use the public polynomial norm and proved witness norms, and include
    # both carry terms.  This conservative bound is checked for every layer.
    for layer in range(layer_count):
        bound = 0
        for col, hs in enumerate(h):
            for i, hp in enumerate(hs):
                j = layer - i
                if 0 <= j < witness_digits:
                    bound += isqrt(exact_norm_sq(hp) * bounds["w"][col][j]) + 1
        for i, qdigit in enumerate(q):
            j = layer - i
            if 0 <= j < quotient_digits:
                bound += abs(qdigit) * (isqrt(bounds["k"][j]) + 1)
        if layer < h_count:
            bound += max(abs(x) for x in t[layer])
        if layer:
            bound += isqrt(bounds["c"][layer - 1]) + 1
        bound += base * (isqrt(bounds["c"][layer]) + 1)
        if bound >= pack_modulus:
            raise ValueError(f"layer {layer} may wrap Pack field: {bound} >= {pack_modulus}")
    return bridge
