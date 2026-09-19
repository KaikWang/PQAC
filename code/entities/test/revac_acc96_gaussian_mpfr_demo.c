#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "acc96_gaussian_mpfr.h"
#include "acc96_klein_exact.h"
#include "acc96_relation.h"

static double seconds(void) {
  struct timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);
  return (double)now.tv_sec + (double)now.tv_nsec / 1e9;
}

static int gadget_relation_holds(const acc96_poly_t y[ACC_TRUNC_PARAM_K],
                                 const acc96_poly_t *w) {
  const uint64_t moduli[ACC96_LIMBS] = {ACC96_P0, ACC96_P1, ACC96_P2};

  for (size_t lane = 0; lane < ACC96_LIMBS; lane++) {
    uint64_t powers[ACC_TRUNC_PARAM_K];
    powers[0] = 1;
    for (size_t i = 1; i < ACC_TRUNC_PARAM_K; i++)
      powers[i] = (powers[i - 1] * ACC_TRUNC_PARAM_B) % moduli[lane];
    for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
      uint64_t value = 0;
      for (size_t i = 0; i < ACC_TRUNC_PARAM_K; i++)
        value = (uint64_t)(((acc96_u128)y[i].coeff[coefficient].limb[lane] *
                            powers[i] + value) % moduli[lane]);
      if (value != w->coeff[coefficient].limb[lane]) return 0;
    }
  }
  return 1;
}

int main(void) {
  const uint8_t seed[ACC96_SEED_BYTES] = {0x96, 0x47, 0x4d};
  acc96_gaussian_mpfr_stats_t stats = {0};
  acc96_gaussian_mpfr_stats_t general_stats = {0};
  acc96_gaussian_mpfr_stats_t parallel_stats[8] = {{0}};
  void *parallel_opaque[8];
  acc96_ctx_t ctx;
  acc96_klein_exact_table_t *table = NULL;
  acc96_poly_t *w = NULL;
  acc96_poly_t *y = NULL;
  mpz_t numerator, denominator;
  mpfr_t general_center, general_variance;
  long double mean = 0.0L;
  double started;
  int ok = 0;

  mpz_inits(numerator, denominator, NULL);
  mpfr_init2(general_center, ACC96_GAUSSIAN_MPFR_PRECISION_BITS);
  mpfr_init2(general_variance, ACC96_GAUSSIAN_MPFR_PRECISION_BITS);
  mpz_set_ui(numerator, 1);
  mpz_set_ui(denominator, 3);
  for (size_t i = 0; i < 1024; i++) {
    int64_t sample;
    if (!acc96_gaussian_mpfr_sample_rational(
            &sample, numerator, denominator,
            acc_trunc_klein_widths_high_precision[0], &stats))
      goto cleanup;
    if (sample < -ACC96_GAUSSIAN_SUPPORT ||
        sample > ACC96_GAUSSIAN_SUPPORT)
      goto cleanup;
    mean += sample;
  }
  mean /= 1024.0L;
  if (mean < 0.0L || mean > 0.7L || stats.proposals <= stats.calls)
    goto cleanup;

  mpfr_set_d(general_center, 123.25, MPFR_RNDN);
  mpfr_set_str(general_variance,
               "11026.890136725950367197279247715404090967458444314", 10,
               MPFR_RNDN);
  for (size_t i = 0; i < 128; i++) {
    int64_t general_sample;
    if (!acc96_gaussian_mpfr_sample_center_variance(
            &general_sample, general_center, general_variance,
            &general_stats) ||
        general_sample < -2000 || general_sample > 2250)
      goto cleanup;
  }
  if (general_stats.calls != 128 ||
      general_stats.certified_accepts != general_stats.calls ||
      general_stats.certified_accepts + general_stats.certified_rejects !=
          general_stats.proposals)
    goto cleanup;
  printf("q96 general center/variance Gaussian passed: 128 samples, "
         "proposals/call=%.2f, refinements=%llu.\n",
         (double)general_stats.proposals / (double)general_stats.calls,
         (unsigned long long)general_stats.interval_refinements);

  if (!acc96_ctx_init(&ctx)) goto cleanup;
  table = acc96_klein_exact_table_new();
  w = calloc(1, sizeof(*w));
  y = calloc(ACC_TRUNC_PARAM_K, sizeof(*y));
  if (table == NULL || w == NULL || y == NULL) goto cleanup;
  acc96_relation_uniform_poly(w, seed, UINT32_C(0x4d504652), 0, 0);
  started = seconds();
  if (!acc96_klein_sample_exact_rational(
          y, w, table, acc96_gaussian_mpfr_sample_rational, &stats, &ctx) ||
      !gadget_relation_holds(y, w) ||
      stats.certified_accepts != stats.calls ||
      stats.certified_accepts + stats.certified_rejects != stats.proposals ||
      stats.maximum_precision_bits < ACC96_GAUSSIAN_MPFR_PRECISION_BITS)
    goto cleanup;
  printf("q96 interval-certified MPFR residual Gaussian + exact-rational "
         "Klein passed; "
         "sequential traversal %.3f s, total calls=%llu, "
         "proposals/call=%.2f, refinements=%llu, max precision=%u bits, "
         "D_Z mean smoke test=%.3Lf.\n",
         seconds() - started, (unsigned long long)stats.calls,
         (double)stats.proposals / (double)stats.calls,
         (unsigned long long)stats.interval_refinements,
         stats.maximum_precision_bits, mean);
  for (size_t i = 0; i < 8; i++) parallel_opaque[i] = &parallel_stats[i];
  started = seconds();
  if (!acc96_klein_sample_exact_rational_parallel(
          y, w, table, acc96_gaussian_mpfr_sample_rational,
          parallel_opaque, 8, &ctx) || !gadget_relation_holds(y, w))
    goto cleanup;
  {
    uint64_t calls = 0;
    uint64_t proposals = 0;
    for (size_t i = 0; i < 8; i++) {
      calls += parallel_stats[i].calls;
      proposals += parallel_stats[i].proposals;
      if (parallel_stats[i].certified_accepts != parallel_stats[i].calls ||
          parallel_stats[i].certified_accepts +
                  parallel_stats[i].certified_rejects !=
              parallel_stats[i].proposals ||
          parallel_stats[i].maximum_precision_bits <
              ACC96_GAUSSIAN_MPFR_PRECISION_BITS)
        goto cleanup;
    }
    printf("q96 parallel exact-rational Klein passed: 8 workers, %.3f s, "
           "proposals/call=%.2f.\n", seconds() - started,
           (double)proposals / (double)calls);
  }
  ok = 1;

cleanup:
  for (size_t i = 0; i < 8; i++)
    acc96_gaussian_mpfr_stats_clear(&parallel_stats[i]);
  acc96_gaussian_mpfr_stats_clear(&general_stats);
  acc96_gaussian_mpfr_stats_clear(&stats);
  free(y);
  free(w);
  acc96_klein_exact_table_free(table);
  mpfr_clear(general_variance);
  mpfr_clear(general_center);
  mpz_clears(denominator, numerator, NULL);
  return ok ? 0 : 1;
}
