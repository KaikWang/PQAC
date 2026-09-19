#ifndef REVAC_ACC96_GAUSSIAN_MPFR_H
#define REVAC_ACC96_GAUSSIAN_MPFR_H

#include <stdint.h>
#include <gmp.h>
#include <mpfr.h>

#define ACC96_GAUSSIAN_MPFR_PRECISION_BITS 384U
#define ACC96_GAUSSIAN_UNIFORM_BITS 256U
#define ACC96_GAUSSIAN_SUPPORT 64
#define ACC96_GAUSSIAN_PROPOSAL_SCALE UINT32_C(1000000)
#define ACC96_GAUSSIAN_PROPOSAL_CACHE 32U

typedef struct {
  uint64_t calls;
  uint64_t proposals;
  uint64_t certified_accepts;
  uint64_t certified_rejects;
  uint64_t interval_refinements;
  uint32_t maximum_precision_bits;
  void *workspace;
} acc96_gaussian_mpfr_stats_t;

void acc96_gaussian_mpfr_stats_clear(acc96_gaussian_mpfr_stats_t *stats);

/* Sample the finite-support approximation to D_{Z,sigma,c}, where c is the
 * exact rational numerator/denominator.  Every transcendental acceptance
 * probability is enclosed by directed-rounding MPFR lower/upper bounds.  If
 * the 256-bit uniform variate falls inside the rounding uncertainty interval,
 * the same decision is recomputed at twice the precision.  Consequently an
 * MPFR rounding error cannot change an acceptance decision.  The fixed +/-64
 * support is selected for the generated q96 widths (sigma <= 6.605); its
 * omitted tail and the 256-bit Bernoulli discretization are accounted for by
 * the profile's statistical-distance budget. */
int acc96_gaussian_mpfr_sample_rational(
    int64_t *sample, const mpz_t center_numerator,
    const mpz_t center_denominator, const char *sigma_decimal,
    void *opaque);

/* General high-precision entry point used by the q96 SampleFz recursion.
 * center and variance are treated as exact MPFR dyadic values; square roots
 * and acceptance probabilities are enclosed with directed rounding.  The
 * proposal covers at least +/-12 standard deviations, giving a per-call
 * omitted-tail bound below 2^-640.  Callers separately account for the
 * distance between their exact theorem parameters and the supplied dyadics. */
int acc96_gaussian_mpfr_sample_center_variance(
    int64_t *sample, mpfr_srcptr center, mpfr_srcptr variance,
    void *opaque);

#endif
