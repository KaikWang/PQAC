#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdio.h>
#include <time.h>

#include "flint/arb_poly.h"

#include "acc96_samplefz.h"

static double seconds(void) {
  struct timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);
  return (double)now.tv_sec + (double)now.tv_nsec / 1e9;
}

int main(void) {
  acc96_samplefz_ctx_t *sampler = acc96_samplefz_ctx_new(384);
  arb_poly_t covariance, center;
  arb_t value;
  acc96_poly_t sample;
  acc96_ctx_t ctx;
  long double norm_square = 0.0L;
  double started;
  int ok = 0;
  if (sampler == NULL || !acc96_ctx_init(&ctx)) return 1;
  arb_poly_init(covariance);
  arb_poly_init(center);
  arb_init(value);
  arb_set_str(value,
              "11026.890136725950367197279247715404090967458444314", 384);
  arb_poly_set_coeff_arb(covariance, 0, value);
  arb_set_str(value, "0.25", 384);
  arb_poly_set_coeff_arb(center, 0, value);

  started = seconds();
  if (!acc96_samplefz(&sample, covariance, center, sampler)) goto cleanup;
  for (size_t i = 0; i < ACC96_N; i++) {
    const acc96_i128 coefficient =
        acc96_coeff_to_centered(&sample.coeff[i], &ctx);
    const long double converted = (long double)coefficient;
    norm_square += converted * converted;
  }
  if (!(norm_square > 0.0L) ||
      !isfinite(acc96_samplefz_maximum_input_radius(sampler)))
    goto cleanup;
  {
    const acc96_gaussian_mpfr_stats_t *stats =
        acc96_samplefz_gaussian_stats(sampler);
    if (stats == NULL || stats->calls != ACC96_N ||
        stats->certified_accepts != stats->calls)
      goto cleanup;
    printf("q96 SampleFz N=4096 passed: %.3f s, norm^2 %.3Le, "
           "proposals/call %.2f, max Arb midpoint radius %.3e.\n",
           seconds() - started, norm_square,
           (double)stats->proposals / (double)stats->calls,
           acc96_samplefz_maximum_input_radius(sampler));
  }

  /* Exercise the nonconstant/self-adjoint covariance path.  In
   * Z[x]/(x^N+1), conjugate(x)=-x^(N-1), hence 8x-8x^(N-1) is
   * self-adjoint.  Its canonical values differ from the positive constant
   * by at most 16, leaving a very large positivity margin. */
  arb_poly_zero(covariance);
  arb_set_si(value, 11027);
  arb_poly_set_coeff_arb(covariance, 0, value);
  arb_set_si(value, 8);
  arb_poly_set_coeff_arb(covariance, 1, value);
  arb_set_si(value, -8);
  arb_poly_set_coeff_arb(covariance, ACC96_N - 1, value);
  arb_poly_zero(center);
  arb_set_str(value, "-17.375", 384);
  arb_poly_set_coeff_arb(center, 1, value);
  started = seconds();
  if (!acc96_samplefz(&sample, covariance, center, sampler)) goto cleanup;
  {
    const acc96_gaussian_mpfr_stats_t *stats =
        acc96_samplefz_gaussian_stats(sampler);
    if (stats == NULL || stats->calls != 2 * ACC96_N ||
        stats->certified_accepts != stats->calls)
      goto cleanup;
    printf("q96 SampleFz nonconstant covariance passed: %.3f s, "
           "max Arb midpoint radius %.3e.\n",
           seconds() - started,
           acc96_samplefz_maximum_input_radius(sampler));
  }

  /* A covariance that is not positive definite must fail closed. */
  arb_poly_zero(covariance);
  arb_set_si(value, -1);
  arb_poly_set_coeff_arb(covariance, 0, value);
  arb_poly_zero(center);
  if (acc96_samplefz(&sample, covariance, center, sampler)) goto cleanup;
  puts("q96 SampleFz rejected a non-positive covariance.");
  ok = 1;

cleanup:
  arb_clear(value);
  arb_poly_clear(center);
  arb_poly_clear(covariance);
  acc96_samplefz_ctx_free(sampler);
  return ok ? 0 : 1;
}
