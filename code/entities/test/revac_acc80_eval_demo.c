#define _POSIX_C_SOURCE 200809L

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "acc80_eval.h"

static uint64_t splitmix64(uint64_t *state) {
  uint64_t value = (*state += UINT64_C(0x9e3779b97f4a7c15));
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

static double seconds(void) {
  struct timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);
  return now.tv_sec + now.tv_nsec / 1e9;
}

static acc80_u128 digit_norm2(
    const int32_t digits[ACC80_GADGET_DIGITS]) {
  acc80_u128 norm = 0;
  for (size_t i = 0; i < ACC80_GADGET_DIGITS; i++) {
    const acc80_i128 digit = digits[i];
    const acc80_u128 magnitude =
        digit < 0 ? (acc80_u128)(-digit) : (acc80_u128)digit;
    norm += magnitude * magnitude;
  }
  return norm;
}

static int decomposition_reconstructs(
    const int32_t digits[ACC80_GADGET_DIGITS], acc80_u128 expected) {
  acc80_i128 reconstructed = 0;
  acc80_i128 power = 1;
  const acc80_i128 modulus = (acc80_i128)acc80_modulus();
  for (size_t i = 0; i < ACC80_GADGET_DIGITS; i++) {
    if (digits[i] < -ACC80_GADGET_BOUND ||
        digits[i] > ACC80_GADGET_BOUND)
      return 0;
    reconstructed += (acc80_i128)digits[i] * power;
    power *= (acc80_i128)ACC80_GADGET_BASE;
  }
  reconstructed %= modulus;
  if (reconstructed < 0) reconstructed += modulus;
  return (acc80_u128)reconstructed == expected;
}

static int check_min_norm_decomposition(const acc80_ctx_t *ctx) {
  static const acc80_u128 small_cases[] = {0, 1, 2};
  const acc80_u128 modulus = acc80_modulus();
  const acc80_u128 half = modulus / 2;
  uint64_t state = UINT64_C(0xd1b54a32d192ed03);
  size_t strict_improvements = 0;
  size_t checked = 0;

  for (size_t i = 0; i < sizeof(small_cases) / sizeof(small_cases[0]); i++) {
    const acc80_u128 values[] = {
        small_cases[i], modulus - 1 - small_cases[i],
        half - small_cases[i], half + small_cases[i]};
    for (size_t j = 0; j < sizeof(values) / sizeof(values[0]); j++) {
      acc80_coeff_t value;
      int32_t canonical[ACC80_GADGET_DIGITS];
      int32_t minimum[ACC80_GADGET_DIGITS];
      int32_t repeated[ACC80_GADGET_DIGITS];
      acc80_coeff_from_u128(&value, values[j]);
      if (!acc80_balanced_decompose_canonical(canonical, &value, ctx) ||
          !acc80_balanced_decompose_min_norm(minimum, &value, ctx) ||
          !acc80_balanced_decompose_min_norm(repeated, &value, ctx) ||
          !decomposition_reconstructs(canonical, values[j]) ||
          !decomposition_reconstructs(minimum, values[j]) ||
          memcmp(minimum, repeated, sizeof(minimum)) != 0 ||
          digit_norm2(minimum) > digit_norm2(canonical))
        return 0;
      strict_improvements += digit_norm2(minimum) < digit_norm2(canonical);
      checked++;
    }
  }

  /* The representable interval is only slightly wider than q, so strict
   * improvements occur near the two centered-lift boundaries. */
  for (acc80_u128 offset = 0; offset <= UINT64_C(1000000000000);
       offset = offset == 0 ? 1 : offset * 1000) {
    const acc80_u128 values[] = {half - offset, half + 1 + offset};
    for (size_t j = 0; j < sizeof(values) / sizeof(values[0]); j++) {
      acc80_coeff_t value;
      int32_t canonical[ACC80_GADGET_DIGITS];
      int32_t minimum[ACC80_GADGET_DIGITS];
      acc80_coeff_from_u128(&value, values[j]);
      if (!acc80_balanced_decompose_canonical(canonical, &value, ctx) ||
          !acc80_balanced_decompose_min_norm(minimum, &value, ctx) ||
          !decomposition_reconstructs(minimum, values[j]) ||
          digit_norm2(minimum) > digit_norm2(canonical))
        return 0;
      strict_improvements += digit_norm2(minimum) < digit_norm2(canonical);
      checked++;
    }
  }

  for (size_t trial = 0; trial < 1024; trial++) {
    const acc80_u128 raw =
        (((acc80_u128)splitmix64(&state) << 64) | splitmix64(&state)) %
        modulus;
    acc80_coeff_t value;
    int32_t canonical[ACC80_GADGET_DIGITS];
    int32_t minimum[ACC80_GADGET_DIGITS];
    acc80_coeff_from_u128(&value, raw);
    if (!acc80_balanced_decompose_canonical(canonical, &value, ctx) ||
        !acc80_balanced_decompose_min_norm(minimum, &value, ctx) ||
        !decomposition_reconstructs(canonical, raw) ||
        !decomposition_reconstructs(minimum, raw) ||
        digit_norm2(minimum) > digit_norm2(canonical))
      return 0;
    strict_improvements += digit_norm2(minimum) < digit_norm2(canonical);
    checked++;
  }
  printf("q88 congruence-class G^-1 KAT passed: %zu coefficients, "
         "%zu strict l2 improvements.\n", checked, strict_improvements);
  return 1;
}

int main(void) {
  acc80_ctx_t ctx;
  acc80_ba_t *ba = calloc(1, sizeof(*ba));
  acc80_gmat_t *result = calloc(1, sizeof(*result));
  acc80_poly_t *value = calloc(1, sizeof(*value));
  acc80_poly_t *acc_delta = calloc(1, sizeof(*acc_delta));
  acc80_poly_t *witness_delta = calloc(1, sizeof(*witness_delta));
  acc80_evalfx_delta_t *delta = calloc(1, sizeof(*delta));
  uint64_t state = UINT64_C(0x13198a2e03707344);
  if (ba == NULL || result == NULL || value == NULL || acc_delta == NULL ||
      witness_delta == NULL || delta == NULL || !acc80_ctx_init(&ctx))
    return 1;

  if (!check_min_norm_decomposition(&ctx)) {
    fprintf(stderr, "q88 minimum-norm congruence-class G^-1 failed.\n");
    return 1;
  }

  for (size_t bit = 0; bit < ACC80_HANDLE_BITS; bit++) {
    for (size_t digit = 0; digit < ACC80_GADGET_DIGITS; digit++) {
      for (size_t coeff = 0; coeff < ACC80_N; coeff++) {
        const acc80_u128 raw =
            (((acc80_u128)splitmix64(&state) << 64) | splitmix64(&state)) %
            acc80_modulus();
        acc80_coeff_from_u128(&ba->bit[bit].block[digit].coeff[coeff], raw);
      }
    }
  }
  if (!acc80_ginv_reconstructs(&ba->bit[0], &ctx)) {
    fprintf(stderr, "q88 G^-1 reconstruction failed.\n");
    return 1;
  }
  const double started = seconds();
  acc80_evalf(result, ba, UINT32_C(0x00c39e71), &ctx);
  const double elapsed = seconds() - started;
  for (size_t coeff = 0; coeff < ACC80_N; coeff++) {
    const acc80_u128 raw =
        (((acc80_u128)splitmix64(&state) << 64) | splitmix64(&state)) %
        acc80_modulus();
    acc80_coeff_from_u128(&value->coeff[coeff], raw);
  }
  acc80_evalf_apply(acc_delta, result, value, &ctx);
  double fx_elapsed = 0.0;
  for (unsigned differing_bit = 0; differing_bit < ACC80_HANDLE_BITS;
       differing_bit++) {
    const uint32_t member =
        UINT32_C(0x00c39e71) ^ (UINT32_C(1) << differing_bit);
    const double fx_started = seconds();
    acc80_evalfx(delta, ba, UINT32_C(0x00c39e71), member, value, &ctx);
    fx_elapsed += seconds() - fx_started;
    acc80_member_matrix_apply(witness_delta, ba, member, delta, &ctx);
    if (!acc80_poly_equal(acc_delta, witness_delta)) {
      fprintf(stderr, "q88 EvalFx indicator relation failed at bit %u.\n",
              differing_bit);
      return 1;
    }
  }
  fx_elapsed /= ACC80_HANDLE_BITS;
  const acc80_u128 checksum =
      acc80_coeff_to_u128(&result->block[0].coeff[0], &ctx);
  printf("q88 EvalF/EvalFx passed: %.3f/%.3f s, "
         "checksum=%016lx%016lx, Ba memory=%.2f MiB, all %u first-difference "
         "positions checked.\n",
         elapsed, fx_elapsed, (unsigned long)(checksum >> 64), (unsigned long)checksum,
         sizeof(*ba) / 1048576.0, ACC80_HANDLE_BITS);
  free(delta);
  free(witness_delta);
  free(acc_delta);
  free(value);
  free(result);
  free(ba);
  return 0;
}
