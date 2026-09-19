#ifndef REVAC_ACC96_KLEIN_EXACT_H
#define REVAC_ACC96_KLEIN_EXACT_H

#include <stddef.h>
#include <stdint.h>
#include <gmp.h>

#include "acc96_crt.h"
#include "acc96_truncated_params_generated.h"

typedef struct acc96_klein_exact_table acc96_klein_exact_table_t;

/* Cache the exact rational -b*_j/<b*_j,b*_j> table used by Klein. */
acc96_klein_exact_table_t *acc96_klein_exact_table_new(void);
void acc96_klein_exact_table_free(acc96_klein_exact_table_t *table);

/* Form one wide-modulus Klein center exactly.  The returned representation is
 * floor(center) mod q plus a fractional center in [0,1).  A caller samples
 * SampleZ(fraction, width[column]) and adds that small signed result with
 * acc96_klein_center_recombine().  No 88-bit integer is converted to double;
 * the residual still has the numerical precision of the selected sampler. */
int acc96_klein_center_split(
    const acc96_klein_exact_table_t *table, size_t column, acc96_i128 w,
    const acc96_i128 v[ACC_TRUNC_PARAM_K], acc96_coeff_t *integer_mod_q,
    double *fraction);

/* Strict backend entry point. residual_numerator and residual_denominator
 * must be initialized mpz objects.  It returns the exact residual with
 * 0 <= numerator < denominator, without any floating-point conversion. */
int acc96_klein_center_split_rational(
    const acc96_klein_exact_table_t *table, size_t column, acc96_i128 w,
    const acc96_i128 v[ACC_TRUNC_PARAM_K], acc96_coeff_t *integer_mod_q,
    mpz_t residual_numerator, mpz_t residual_denominator);

void acc96_klein_center_recombine(acc96_coeff_t *out,
                                  const acc96_coeff_t *integer_mod_q,
                                  int64_t centered_sample,
                                  const acc96_ctx_t *ctx);

typedef int64_t (*acc96_residual_gaussian_fn)(double fractional_center,
                                               double sigma, void *opaque);

typedef int (*acc96_residual_gaussian_rational_fn)(
    int64_t *sample, const mpz_t center_numerator,
    const mpz_t center_denominator, const char *sigma_decimal,
    void *opaque);

/* Wide-center-safe Klein traversal.  The callback only ever sees a center in
 * [0,1), while the returned sample is recombined with the exact wide integer
 * part.  On success, sum_i 13^i y_i = w in R_q.  A security claim must also
 * account for the callback's finite-precision statistical error. */
int acc96_klein_sample_exact(
    acc96_poly_t y[ACC_TRUNC_PARAM_K], const acc96_poly_t *w,
    const acc96_klein_exact_table_t *table,
    acc96_residual_gaussian_fn residual_sampler, void *opaque,
    const acc96_ctx_t *ctx);

/* Production-oriented traversal: neither the rational residual center nor
 * the generated Gaussian width is rounded through binary64. */
int acc96_klein_sample_exact_rational(
    acc96_poly_t y[ACC_TRUNC_PARAM_K], const acc96_poly_t *w,
    const acc96_klein_exact_table_t *table,
    acc96_residual_gaussian_rational_fn residual_sampler, void *opaque,
    const acc96_ctx_t *ctx);

/* Coefficients are independent Klein chains. Each worker receives its own
 * sampler state through thread_opaque[worker], so no PRNG/MPFR object is
 * shared. */
int acc96_klein_sample_exact_rational_parallel(
    acc96_poly_t y[ACC_TRUNC_PARAM_K], const acc96_poly_t *w,
    const acc96_klein_exact_table_t *table,
    acc96_residual_gaussian_rational_fn residual_sampler,
    void *const thread_opaque[], size_t thread_count,
    const acc96_ctx_t *ctx);

#endif
