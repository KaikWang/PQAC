#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>

#define TEST_BASE UINT64_C(1025)
#define TEST_DIGITS 8
#define TEST_BOUND INT64_C(512)
#define TEST_RANDOM_CASES 100000

typedef __int128 i128;
typedef unsigned __int128 u128;

static uint64_t splitmix64(uint64_t *state) {
  uint64_t value = (*state += UINT64_C(0x9e3779b97f4a7c15));
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

static int decompose(int64_t digits[TEST_DIGITS], i128 value) {
  for (size_t index = 0; index < TEST_DIGITS; index++) {
    i128 remainder = value % (i128)TEST_BASE;
    if (remainder < 0) remainder += (i128)TEST_BASE;
    if (remainder > TEST_BOUND) remainder -= (i128)TEST_BASE;
    digits[index] = (int64_t)remainder;
    value = (value - remainder) / (i128)TEST_BASE;
  }
  return value == 0;
}

static int check_value(i128 input) {
  int64_t digits[TEST_DIGITS];
  i128 reconstructed = 0;
  i128 power = 1;

  if (!decompose(digits, input)) return 0;
  for (size_t index = 0; index < TEST_DIGITS; index++) {
    if (digits[index] < -TEST_BOUND || digits[index] > TEST_BOUND) return 0;
    reconstructed += (i128)digits[index] * power;
    power *= (i128)TEST_BASE;
  }
  return reconstructed == input;
}

int main(void) {
  const u128 p0 = UINT64_C(1099511480321);
  const u128 p1 = UINT64_C(1099511390209);
  const u128 q = p0 * p1;
  const i128 half = (i128)((q - 1) / 2);
  const i128 boundaries[] = {-half, -half + 1, -1, 0, 1, half - 1, half};
  uint64_t state = UINT64_C(0x243f6a8885a308d3);

  if (p0 % UINT64_C(8192) != 1 || p1 % UINT64_C(8192) != 1) {
    fprintf(stderr, "CRT factors are not degree-4096 NTT friendly.\n");
    return 1;
  }
  for (size_t index = 0; index < sizeof(boundaries) / sizeof(boundaries[0]);
       index++) {
    if (!check_value(boundaries[index])) {
      fprintf(stderr, "balanced odd-base boundary KAT failed.\n");
      return 1;
    }
  }
  for (size_t index = 0; index < TEST_RANDOM_CASES; index++) {
    const u128 raw = (((u128)splitmix64(&state) << 64) |
                      (u128)splitmix64(&state)) % q;
    const i128 centered = raw > (u128)half ? (i128)(raw - q) : (i128)raw;
    if (!check_value(centered)) {
      fprintf(stderr, "balanced odd-base random KAT failed at %zu.\n", index);
      return 1;
    }
  }
  printf("balanced odd-base KAT passed: base=%" PRIu64
         ", digits=%d, bound=%" PRId64 ", CRT bits=80, cases=%d.\n",
         TEST_BASE, TEST_DIGITS, TEST_BOUND,
         TEST_RANDOM_CASES + (int)(sizeof(boundaries) / sizeof(boundaries[0])));
  return 0;
}
