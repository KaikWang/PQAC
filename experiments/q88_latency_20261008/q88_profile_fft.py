"""Batched preparation FFT for independent q88 profiles only.

Scalar factors are multiplied as integers. Dense small-digit products use
bounded FFT batches; native integer bridge checks remain the acceptance gate.
Nothing in the transcript, witness layout or bound selection is changed.
"""
import numpy as np


class DigitConvolution(list):
    """Prepared digits live only for this invocation, not across witnesses."""


def layers(polys, bits, count, signed):
    base = 1 << bits
    remaining = np.asarray(polys, dtype=object)
    result = []
    for _ in range(count):
        digit = (remaining + base//2) % base - base//2 if signed else remaining % base
        result.append(np.asarray(digit, dtype=np.int64))
        remaining = (remaining-digit)//base
    if np.any(remaining):
        raise ValueError("batched FFT digit capacity exhausted")
    return np.stack(result, axis=1)


def digit_sum_negacyclic(matrix, witness, bits=7, left_count=13,
                         right_count=3, *, public_cache=None):
    if not matrix or len(matrix) != len(witness):
        raise ValueError("batched FFT matrix/witness mismatch")
    n = len(matrix[0])
    if n & (n-1) or any(len(p) != n for p in list(matrix)+list(witness)):
        raise ValueError("batched FFT degree mismatch")
    # Conservative integer overflow / rounding budget for a complete layer.
    if len(matrix)*n*(1 << bits)**2*min(left_count, right_count) >= 1 << 42:
        raise ValueError("batched FFT exceeds rounding budget")
    result = np.zeros((left_count+right_count, n), dtype=np.int64)
    h_digits, w_digits = [], []
    for start in range(0, len(matrix), 16):
        columns = matrix[start:start+16]
        h = (layers(columns, bits, left_count, False) if public_cache is None
             else public_cache.digits(columns, bits, left_count))
        w = layers(witness[start:start+16], bits, right_count, True)
        h_digits.extend(h.tolist())
        w_digits.extend(w.tolist())
        dense = []
        for col in range(len(h)):
            if not np.any(h[col, :, 1:]):
                for i in range(left_count):
                    for j in range(right_count):
                        result[i+j] += h[col, i, 0]*w[col, j]
            elif not np.any(w[col, :, 1:]):
                for i in range(left_count):
                    for j in range(right_count):
                        result[i+j] += w[col, j, 0]*h[col, i]
            else:
                dense.append(col)
        if not dense:
            continue
        hf = np.fft.rfft(h[dense], n=2*n, axis=-1)
        wf = np.fft.rfft(w[dense], n=2*n, axis=-1)
        for layer in range(left_count+right_count-1):
            spectrum = np.zeros(n+1, dtype=np.complex128)
            for i in range(max(0, layer-right_count+1), min(left_count, layer+1)):
                j = layer-i
                spectrum += np.sum(hf[:, i]*wf[:, j], axis=0)
            conv = np.fft.irfft(spectrum, n=2*n)
            rounded = np.rint(conv)
            if np.max(np.abs(conv-rounded)) > 0.125:
                raise ValueError("batched FFT rounding margin failed")
            ints = rounded.astype(np.int64)
            result[layer] += ints[:n]-ints[n:]
    output = DigitConvolution(result.tolist())
    output.h_digits = h_digits
    output.w_digits = w_digits
    return output


class PublicDigitCache:
    """Bounded public-only cache; no witness digits, masks or responses."""
    def __init__(self, byte_limit=128 << 20):
        from collections import OrderedDict
        self.entries = OrderedDict()
        self.byte_limit, self.bytes, self.hits = byte_limit, 0, 0

    def digits(self, columns, bits, count):
        key = (bits, count, tuple(tuple(p) for p in columns))
        if key in self.entries:
            self.hits += 1
            self.entries.move_to_end(key)
            return self.entries[key][0]
        values = layers(columns, bits, count, False)
        values.flags.writeable = False
        # Include the key's Python integer/tuple footprint in the budget.
        footprint = values.nbytes + sum(48*len(p)+64 for p in columns)
        if footprint <= self.byte_limit:
            while self.entries and self.bytes+footprint > self.byte_limit:
                _, (_, previous_size) = self.entries.popitem(last=False)
                self.bytes -= previous_size
            self.entries[key] = values, footprint
            self.bytes += footprint
        return values

    def __call__(self, matrix, witness, bits=7, left_count=13, right_count=3):
        return digit_sum_negacyclic(matrix, witness, bits, left_count, right_count,
                                    public_cache=self)
