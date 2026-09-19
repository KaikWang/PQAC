#ifndef REVAC_ACC96_SAMPLEFZ_H
#define REVAC_ACC96_SAMPLEFZ_H

/* FLINT exposes arf_get_mpfr only when mpfr.h has already been included. */
#include <mpfr.h>
#include <stddef.h>

#include "flint/arb_poly.h"

#include "acc96_crt.h"
#include "acc96_gaussian_mpfr.h"

typedef struct acc96_samplefz_ctx acc96_samplefz_ctx_t;

acc96_samplefz_ctx_t *acc96_samplefz_ctx_new(long precision_bits);
void acc96_samplefz_ctx_free(acc96_samplefz_ctx_t *sampler);

/* GM18 SampleFz specialized to Z[x]/(x^4096+1).  Covariance and center are
 * Arb polynomials.  Recursive ring arithmetic remains interval-valued; each
 * scalar leaf is sampled from the exact dyadic midpoint and the largest Arb
 * input radius is exposed for the statistical-distance audit. */
int acc96_samplefz(acc96_poly_t *out, const arb_poly_t covariance,
                   const arb_poly_t center, acc96_samplefz_ctx_t *sampler);

/* Prepare and use an immutable covariance factorization.  The fixed plans
 * are intended for SampleLeft, where the three covariance polynomials are
 * reused for every issuance.  All stored factors remain Arb intervals; this
 * is a setup-time optimization and does not replace the certified arithmetic
 * with midpoint floating point values. */
int acc96_samplefz_prepare_fixed(acc96_samplefz_ctx_t *sampler,
                                 const arb_poly_t covariance, size_t slot);
int acc96_samplefz_fixed(acc96_poly_t *out, const arb_poly_t center,
                         acc96_samplefz_ctx_t *sampler, size_t slot);

/* Outward-rounded arithmetic in the same cyclotomic ring and precision as
 * SampleFz.  These helpers let SampleLeft build and cache the Schur-complement
 * covariance once, instead of duplicating the delicate DFT/inversion code. */
void acc96_samplefz_ring_mul(arb_poly_t out, const arb_poly_t left,
                             const arb_poly_t right,
                             acc96_samplefz_ctx_t *sampler);
int acc96_samplefz_ring_invert(arb_poly_t out, const arb_poly_t input,
                              acc96_samplefz_ctx_t *sampler);
void acc96_samplefz_ring_conjugate(arb_poly_t out,
                                   const arb_poly_t input,
                                   acc96_samplefz_ctx_t *sampler);

const acc96_gaussian_mpfr_stats_t *acc96_samplefz_gaussian_stats(
    const acc96_samplefz_ctx_t *sampler);
double acc96_samplefz_maximum_input_radius(const acc96_samplefz_ctx_t *sampler);

#endif
