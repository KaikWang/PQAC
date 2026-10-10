"""Exact integer negacyclic convolution by small-base FFT accumulation.

The inputs are decomposed into base-256 digits.  Each digit convolution is
small enough for float64 FFT rounding; the rounded integer coefficients are
then carried in Python integers.  This is used only to construct the 100-bit
lift quotient, never as a proof-side arithmetic primitive.
"""

from __future__ import annotations

import numpy as np


def _digits(values, signed, bits, count):
    base = 1 << bits
    half = base // 2
    out = []
    for value in values:
        value = int(value)
        row = []
        if signed:
            for _ in range(count):
                digit = (value + half) % base - half
                row.append(digit)
                value = (value - digit) // base
        else:
            if value < 0:
                raise ValueError("unsigned input is negative")
            for _ in range(count):
                row.append(value & (base - 1))
                value >>= bits
        if value:
            raise ValueError("base-256 digit capacity exhausted")
        out.append(row)
    return np.asarray(out, dtype=np.float64).T


def exact_sum_negacyclic(matrix, witness, bits=8, left_count=13,
                         right_count=13):
    """Return integer coefficients of sum_i matrix[i]*witness[i]."""
    if len(matrix) != len(witness) or not matrix:
        raise ValueError("matrix/witness shape mismatch")
    n = len(matrix[0])
    if n & (n - 1):
        raise ValueError("degree must be a power of two")
    if any(len(x) != n for x in matrix + witness):
        raise ValueError("degree mismatch")
    size = 2 * n
    result = digit_sum_negacyclic(matrix, witness, bits, left_count,
                                  right_count)
    return [sum(result[layer][k] * ((1 << bits) ** layer)
                for layer in range(len(result))) for k in range(n)]


def digit_sum_negacyclic(matrix, witness, bits=8, left_count=13,
                         right_count=13):
    """Return convolution coefficients grouped by base digit layer."""
    n = len(matrix[0])
    size = 2 * n
    result = [[0] * n for _ in range(left_count + right_count)]
    # Accumulate each base-256 digit pair.  Coefficients are below 2^53,
    # hence inverse FFT rounding is exact for the supported profile.
    left_digits = [_digits(left, False, bits, left_count) for left in matrix]
    right_digits = [_digits(right, True, bits, right_count) for right in witness]
    left_fft = [[np.fft.rfft(np.pad(d[i], (0, n)))
                 for i in range(left_count)] for d in left_digits]
    right_fft = [[np.fft.rfft(np.pad(d[j], (0, n)))
                  for j in range(right_count)] for d in right_digits]
    for i in range(left_count):
        for j in range(right_count):
            spectrum = sum((left_fft[col][i] * right_fft[col][j]
                            for col in range(len(matrix))))
            conv = np.rint(np.fft.irfft(spectrum, size)).astype(np.int64)
            folded = conv[:n].astype(object) - conv[n:].astype(object)
            for k, value in enumerate(folded):
                result[i + j][k] += int(value)
    return result
