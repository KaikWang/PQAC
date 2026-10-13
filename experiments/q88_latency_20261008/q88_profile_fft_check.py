"""Independent reference checks for preparation arithmetic and digit reuse."""
import random

from q88_profile_fft import digit_sum_negacyclic as fast, PublicDigitCache
from q88_bridge_core import compile_single88_bridge_fast, _decompose
from exact_negacyclic_fft import digit_sum_negacyclic as reference


def main():
    random_source = random.Random(20261013)
    for n in (16, 64, 256, 4096):
        matrix = [[random_source.randrange(128**3) for _ in range(n)] for _ in range(19)]
        witness = [[random_source.randrange(-6000, 6000) for _ in range(n)] for _ in range(19)]
        matrix[0] = [73]+[0]*(n-1)
        witness[1] = [1]+[0]*(n-1)
        matrix[2] = [0]*n
        result = fast(matrix, witness, 7, 3, 2)
        assert result == reference(matrix, witness, 7, 3, 2), f"convolution n={n}"
        assert result.h_digits == [_decompose(p, 128, 3, signed=False) for p in matrix]
        assert result.w_digits == [_decompose(p, 128, 2, signed=True) for p in witness]
        cache = PublicDigitCache(32 << 20)
        assert cache(matrix, witness, 7, 3, 2) == result
        changed = [p[:] for p in witness]
        changed[0][0] += 1
        assert cache(matrix, changed, 7, 3, 2) == reference(matrix, changed, 7, 3, 2)
        assert cache.hits == 2, "public cache did not reuse both chunks"
        if n == 16:
            # Both bridges also perform the independent integer O(n^2)
            # polynomial recurrence check, beyond FFT-to-FFT comparison.
            assert compile_single88_bridge_fast(
                matrix, witness, witness_digits=2, convolution=fast) == compile_single88_bridge_fast(
                    matrix, witness, witness_digits=2)
    print("FFT, prepared digits and integer bridge equality: passed")


if __name__ == "__main__":
    main()
