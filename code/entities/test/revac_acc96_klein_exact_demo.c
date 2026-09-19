#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "acc96_klein_exact.h"

static double seconds(void) {
  struct timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);
  return now.tv_sec + now.tv_nsec / 1e9;
}

static int64_t nearest_residual(double center, double sigma, void *opaque) {
  (void)sigma;
  (void)opaque;
  return center < 0.5 ? 0 : 1;
}

static int gadget_relation_holds(
    const acc96_poly_t y[ACC_TRUNC_PARAM_K], const acc96_poly_t *w) {
  static const uint64_t moduli[ACC96_LIMBS] = {ACC96_P0, ACC96_P1, ACC96_P2};
  for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
    for (size_t lane = 0; lane < ACC96_LIMBS; lane++) {
      uint64_t sum = 0;
      uint64_t power = 1;
      for (size_t row = 0; row < ACC_TRUNC_PARAM_K; row++) {
        const uint64_t term =
            (uint64_t)(((acc96_u128)y[row].coeff[coefficient].limb[lane] *
                        power) % moduli[lane]);
        sum += term;
        if (sum >= moduli[lane]) sum -= moduli[lane];
        power = (uint64_t)(((acc96_u128)power * ACC_TRUNC_PARAM_B) %
                           moduli[lane]);
      }
      if (sum != w->coeff[coefficient].limb[lane]) return 0;
    }
  }
  return 1;
}

int main(void) {
  acc96_ctx_t ctx;
  acc96_i128 v[ACC_TRUNC_PARAM_K] = {0};
  acc96_coeff_t integer_part, recombined;
  acc96_klein_exact_table_t *table;
  double fraction;
  mpz_t residual_numerator, residual_denominator;
  uint64_t checksum = 0;
  acc96_poly_t *syndrome = calloc(1, sizeof(*syndrome));
  acc96_poly_t *y = calloc(ACC_TRUNC_PARAM_K, sizeof(*y));
  if (syndrome == NULL || y == NULL || !acc96_ctx_init(&ctx)) return 1;
  mpz_inits(residual_numerator, residual_denominator, NULL);
  table = acc96_klein_exact_table_new();
  if (table == NULL) return 1;

  /* Deliberately exceed the 53-bit exact range of double. */
  const acc96_i128 w = ((acc96_i128)1 << 86) + 12345;
  v[0] = -54321;
  v[7] = ((acc96_i128)1 << 82) + 17;
  v[23] = -((acc96_i128)1 << 84) + 91;
  for (size_t column = 0; column < ACC_TRUNC_PARAM_K; column++) {
    if (!acc96_klein_center_split(table, column, w, v, &integer_part,
                                  &fraction) ||
        !isfinite(fraction) || fraction < 0.0 || fraction >= 1.0)
      return 1;
    acc96_klein_center_recombine(&recombined, &integer_part, -7, &ctx);
    const acc96_u128 got = acc96_coeff_to_u128(&recombined, &ctx);
    const acc96_u128 base = acc96_coeff_to_u128(&integer_part, &ctx);
    const acc96_u128 expected = base >= 7 ? base - 7 : acc96_modulus() - (7 - base);
    if (got != expected) return 1;
    if (!acc96_klein_center_split_rational(
            table, column, w, v, &integer_part, residual_numerator,
            residual_denominator) ||
        mpz_sgn(residual_numerator) < 0 ||
        mpz_cmp(residual_numerator, residual_denominator) >= 0)
      return 1;
  }
  if (acc96_klein_center_split(table, ACC_TRUNC_PARAM_K, w, v,
                               &integer_part, &fraction))
    return 1;
  const double started = seconds();
  for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
    v[coefficient % ACC_TRUNC_PARAM_K] +=
        (acc96_i128)(coefficient * UINT64_C(0x9e3779b1) + 1);
    for (size_t column = 0; column < ACC_TRUNC_PARAM_K; column++) {
      if (!acc96_klein_center_split(table, column, w + coefficient, v,
                                    &integer_part, &fraction))
        return 1;
      checksum ^= integer_part.limb[column % ACC96_LIMBS] +
                  (uint64_t)(fraction * 4294967296.0);
      checksum = (checksum << 7) | (checksum >> (64 - 7));
    }
  }
  const double all_centers_time = seconds() - started;
  for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++)
    acc96_coeff_from_u128(&syndrome->coeff[coefficient],
                          (((acc96_u128)1 << 86) + coefficient * 7919) %
                              acc96_modulus());
  const double klein_started = seconds();
  if (!acc96_klein_sample_exact(y, syndrome, table, nearest_residual, NULL,
                                &ctx) ||
      !gadget_relation_holds(y, syndrome))
    return 1;
  const double klein_time = seconds() - klein_started;
  acc96_klein_exact_table_free(table);
  mpz_clears(residual_denominator, residual_numerator, NULL);
  printf("q88 exact Klein center passed: 88-bit integer centers remain exact; "
         "only a [0,1) residual is converted to double; signed recombination "
         "was checked for all 24 columns. Full 4096x24 center pass %.3f s; "
         "checksum 0x%016llx. Wide-center-safe Klein traversal %.3f s and gadget "
         "preimage relation passed.\n",
         all_centers_time, (unsigned long long)checksum, klein_time);
  free(y);
  free(syndrome);
  return 0;
}
