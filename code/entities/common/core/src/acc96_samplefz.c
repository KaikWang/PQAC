#include "acc96_samplefz.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

#include "flint/acb_dft.h"
#include "flint/fmpz.h"

#define ACC96_SAMPLEFZ_LOG_N 12U
#define ACC96_SAMPLEFZ_MIN_PRECISION 256
#define ACC96_SAMPLEFZ_WORK_ITEMS 7U
#define ACC96_SAMPLEFZ_FIXED_PLANS 3U

typedef struct acc96_samplefz_plan_node {
  arb_poly_t multiplier;
  arb_poly_t variance;
  struct acc96_samplefz_plan_node *odd;
  struct acc96_samplefz_plan_node *even;
} acc96_samplefz_plan_node_t;

struct acc96_samplefz_ctx {
  slong precision;
  acb_dft_pre_struct pre[ACC96_SAMPLEFZ_LOG_N];
  /* One set of temporary polynomials per recursion level.  The covariance
   * tree is traversed depth-first, so a level's seven temporaries remain
   * live while its child level is active and can be reused by the next
   * SampleFz invocation. */
  arb_poly_t work[ACC96_SAMPLEFZ_LOG_N][ACC96_SAMPLEFZ_WORK_ITEMS];
  acc96_samplefz_plan_node_t *fixed[ACC96_SAMPLEFZ_FIXED_PLANS];
  arb_t temporary;
  arb_t auxiliary;
  mpfr_t leaf_center;
  mpfr_t leaf_variance;
  acc96_gaussian_mpfr_stats_t gaussian_stats;
  double maximum_input_radius;
};

static void record_radius(acc96_samplefz_ctx_t *sampler, const arb_t value) {
  const double radius = mag_get_d(arb_radref(value));
  if (radius > sampler->maximum_input_radius)
    sampler->maximum_input_radius = radius;
}

acc96_samplefz_ctx_t *acc96_samplefz_ctx_new(long precision_bits) {
  acc96_samplefz_ctx_t *sampler = calloc(1, sizeof(*sampler));
  if (sampler == NULL) return NULL;
  sampler->precision = precision_bits >= ACC96_SAMPLEFZ_MIN_PRECISION
                           ? (slong)precision_bits
                           : ACC96_SAMPLEFZ_MIN_PRECISION;
  for (size_t level = 1; level <= ACC96_SAMPLEFZ_LOG_N; level++)
    acb_dft_precomp_init(&sampler->pre[level - 1], (slong)2U << level,
                         sampler->precision);
  for (size_t level = 0; level < ACC96_SAMPLEFZ_LOG_N; level++)
    for (size_t item = 0; item < ACC96_SAMPLEFZ_WORK_ITEMS; item++)
      arb_poly_init(sampler->work[level][item]);
  arb_init(sampler->temporary);
  arb_init(sampler->auxiliary);
  mpfr_init2(sampler->leaf_center, (mpfr_prec_t)sampler->precision);
  mpfr_init2(sampler->leaf_variance, (mpfr_prec_t)sampler->precision);
  return sampler;
}

static void plan_free(acc96_samplefz_plan_node_t *node) {
  if (node == NULL) return;
  plan_free(node->odd);
  plan_free(node->even);
  arb_poly_clear(node->variance);
  arb_poly_clear(node->multiplier);
  free(node);
}

void acc96_samplefz_ctx_free(acc96_samplefz_ctx_t *sampler) {
  if (sampler == NULL) return;
  acc96_gaussian_mpfr_stats_clear(&sampler->gaussian_stats);
  mpfr_clear(sampler->leaf_variance);
  mpfr_clear(sampler->leaf_center);
  arb_clear(sampler->auxiliary);
  arb_clear(sampler->temporary);
  for (size_t level = 0; level < ACC96_SAMPLEFZ_LOG_N; level++)
    acb_dft_precomp_clear(&sampler->pre[level]);
  for (size_t level = 0; level < ACC96_SAMPLEFZ_LOG_N; level++)
    for (size_t item = 0; item < ACC96_SAMPLEFZ_WORK_ITEMS; item++)
      arb_poly_clear(sampler->work[level][item]);
  for (size_t slot = 0; slot < ACC96_SAMPLEFZ_FIXED_PLANS; slot++)
    plan_free(sampler->fixed[slot]);
  free(sampler);
}

static void multiply(arb_poly_t out, const arb_poly_t left,
                     const arb_poly_t right, size_t log_dimension,
                     acc96_samplefz_ctx_t *sampler) {
  const slong dimension = (slong)1U << log_dimension;
  if (log_dimension == 0) {
    arb_poly_get_coeff_arb(sampler->temporary, left, 0);
    arb_poly_get_coeff_arb(sampler->auxiliary, right, 0);
    arb_mul(sampler->temporary, sampler->temporary, sampler->auxiliary,
            sampler->precision);
    arb_poly_zero(out);
    arb_poly_set_coeff_arb(out, 0, sampler->temporary);
    return;
  }
  arb_poly_mul(out, left, right, sampler->precision);
  for (slong i = 0; i < dimension; i++) {
    arb_poly_get_coeff_arb(sampler->temporary, out, i);
    arb_poly_get_coeff_arb(sampler->auxiliary, out, i + dimension);
    arb_sub(sampler->temporary, sampler->temporary, sampler->auxiliary,
            sampler->precision);
    arb_poly_set_coeff_arb(out, i, sampler->temporary);
    arb_poly_set_coeff_si(out, i + dimension, 0);
  }
}

static int invert(arb_poly_t out, const arb_poly_t input,
                  size_t log_dimension, acc96_samplefz_ctx_t *sampler) {
  if (log_dimension == 0) {
    arb_poly_get_coeff_arb(sampler->temporary, input, 0);
    if (!arb_is_positive(sampler->temporary)) return 0;
    arb_inv(sampler->temporary, sampler->temporary, sampler->precision);
    arb_poly_zero(out);
    arb_poly_set_coeff_arb(out, 0, sampler->temporary);
    return arb_is_finite(sampler->temporary);
  }
  const slong dimension = (slong)1U << log_dimension;
  acb_ptr values = _acb_vec_init(2 * dimension);
  if (values == NULL) return 0;
  for (slong i = 0; i < dimension; i++) {
    arb_poly_get_coeff_arb(sampler->temporary, input, i);
    acb_set_arb(&values[i], sampler->temporary);
  }
  acb_dft_precomp(values, values, &sampler->pre[log_dimension - 1],
                  sampler->precision);
  for (slong i = 0; i < 2 * dimension; i++) {
    if (acb_contains_zero(&values[i])) {
      _acb_vec_clear(values, 2 * dimension);
      return 0;
    }
    acb_inv(&values[i], &values[i], sampler->precision);
  }
  acb_dft_inverse_precomp(values, values,
                          &sampler->pre[log_dimension - 1],
                          sampler->precision);
  arb_poly_zero(out);
  for (slong i = 0; i < dimension; i++) {
    acb_get_real(sampler->temporary, &values[i]);
    acb_get_real(sampler->auxiliary, &values[i + dimension]);
    arb_sub(sampler->temporary, sampler->temporary, sampler->auxiliary,
            sampler->precision);
    arb_poly_set_coeff_arb(out, i, sampler->temporary);
  }
  _acb_vec_clear(values, 2 * dimension);
  return 1;
}

static void conjugate(arb_poly_t polynomial, size_t log_dimension,
                      acc96_samplefz_ctx_t *sampler) {
  if (log_dimension == 0) return;
  const slong dimension = (slong)1U << log_dimension;
  for (slong i = 1; i < dimension / 2; i++) {
    arb_poly_get_coeff_arb(sampler->temporary, polynomial, i);
    arb_poly_get_coeff_arb(sampler->auxiliary, polynomial, dimension - i);
    arb_neg(sampler->temporary, sampler->temporary);
    arb_neg(sampler->auxiliary, sampler->auxiliary);
    arb_poly_set_coeff_arb(polynomial, i, sampler->auxiliary);
    arb_poly_set_coeff_arb(polynomial, dimension - i,
                           sampler->temporary);
  }
  arb_poly_get_coeff_arb(sampler->temporary, polynomial, dimension / 2);
  arb_neg(sampler->temporary, sampler->temporary);
  arb_poly_set_coeff_arb(polynomial, dimension / 2, sampler->temporary);
}

static void split(arb_poly_t even, arb_poly_t odd,
                  const arb_poly_t polynomial,
                  acc96_samplefz_ctx_t *sampler) {
  arb_poly_zero(even);
  arb_poly_zero(odd);
  for (slong i = 0; i < arb_poly_length(polynomial); i++) {
    arb_poly_get_coeff_arb(sampler->temporary, polynomial, i);
    arb_poly_set_coeff_arb((i & 1) == 0 ? even : odd, i / 2,
                           sampler->temporary);
  }
}

/* Build the fixed conditional-covariance tree used by SampleFz.  At a
 * non-leaf node the odd child uses f_even and the even child uses
 * f_even - (f_even^{-1} f_odd) f_odd^*.  Only the multiplier and leaf
 * variance are retained; all stored values remain outward-rounded Arb balls. */
static acc96_samplefz_plan_node_t *plan_build(
    const arb_poly_t covariance, size_t log_dimension,
    acc96_samplefz_ctx_t *sampler) {
  acc96_samplefz_plan_node_t *node = calloc(1, sizeof(*node));
  arb_poly_t f_even, f_odd, conditional, temporary;
  if (node == NULL) return NULL;
  arb_poly_init(node->multiplier);
  arb_poly_init(node->variance);
  if (log_dimension == 0) {
    arb_poly_set(node->variance, covariance);
    return node;
  }

  arb_poly_init(f_even);
  arb_poly_init(f_odd);
  arb_poly_init(conditional);
  arb_poly_init(temporary);
  split(f_even, f_odd, covariance, sampler);
  if (!invert(temporary, f_even, log_dimension - 1, sampler))
    goto fail;
  multiply(node->multiplier, temporary, f_odd, log_dimension - 1, sampler);
  conjugate(f_odd, log_dimension - 1, sampler);
  multiply(temporary, node->multiplier, f_odd, log_dimension - 1, sampler);
  arb_poly_sub(conditional, f_even, temporary, sampler->precision);

  node->odd = plan_build(f_even, log_dimension - 1, sampler);
  node->even = plan_build(conditional, log_dimension - 1, sampler);
  if (node->odd != NULL && node->even != NULL) {
    arb_poly_clear(temporary);
    arb_poly_clear(conditional);
    arb_poly_clear(f_odd);
    arb_poly_clear(f_even);
    return node;
  }

fail:
  plan_free(node->odd);
  plan_free(node->even);
  node->odd = NULL;
  node->even = NULL;
  arb_poly_clear(temporary);
  arb_poly_clear(conditional);
  arb_poly_clear(f_odd);
  arb_poly_clear(f_even);
  plan_free(node);
  return NULL;
}

void acc96_samplefz_ring_mul(arb_poly_t out, const arb_poly_t left,
                             const arb_poly_t right,
                             acc96_samplefz_ctx_t *sampler) {
  if (sampler != NULL)
    multiply(out, left, right, ACC96_SAMPLEFZ_LOG_N, sampler);
}

int acc96_samplefz_ring_invert(arb_poly_t out, const arb_poly_t input,
                              acc96_samplefz_ctx_t *sampler) {
  return sampler != NULL &&
         invert(out, input, ACC96_SAMPLEFZ_LOG_N, sampler);
}

void acc96_samplefz_ring_conjugate(arb_poly_t out,
                                   const arb_poly_t input,
                                   acc96_samplefz_ctx_t *sampler) {
  if (sampler == NULL) return;
  arb_poly_set(out, input);
  conjugate(out, ACC96_SAMPLEFZ_LOG_N, sampler);
}

static int sample_recursive(arb_poly_t result, const arb_poly_t covariance,
                            const arb_poly_t center, size_t log_dimension,
                            size_t depth, acc96_samplefz_ctx_t *sampler) {
  if (log_dimension == 0) {
    int64_t sample;
    arb_poly_get_coeff_arb(sampler->temporary, center, 0);
    arb_poly_get_coeff_arb(sampler->auxiliary, covariance, 0);
    if (!arb_is_finite(sampler->temporary) ||
        !arb_is_positive(sampler->auxiliary))
      return 0;
    record_radius(sampler, sampler->temporary);
    record_radius(sampler, sampler->auxiliary);
    arf_get_mpfr(sampler->leaf_center, arb_midref(sampler->temporary),
                 MPFR_RNDN);
    arf_get_mpfr(sampler->leaf_variance, arb_midref(sampler->auxiliary),
                 MPFR_RNDN);
    if (!acc96_gaussian_mpfr_sample_center_variance(
            &sample, sampler->leaf_center, sampler->leaf_variance,
            &sampler->gaussian_stats))
      return 0;
    arb_poly_zero(result);
    arb_poly_set_coeff_si(result, 0, (slong)sample);
    return 1;
  }

  arb_poly_struct *f_even = sampler->work[depth][0];
  arb_poly_struct *f_odd = sampler->work[depth][1];
  arb_poly_struct *c_even = sampler->work[depth][2];
  arb_poly_struct *c_odd = sampler->work[depth][3];
  arb_poly_struct *q_even = sampler->work[depth][4];
  arb_poly_struct *q_odd = sampler->work[depth][5];
  arb_poly_struct *temporary = sampler->work[depth][6];
  int ok = 0;

  split(f_even, f_odd, covariance, sampler);
  split(c_even, c_odd, center, sampler);
  if (!sample_recursive(q_odd, f_even, c_odd, log_dimension - 1,
                        depth + 1, sampler) ||
      !invert(temporary, f_even, log_dimension - 1, sampler))
    goto cleanup;
  multiply(temporary, temporary, f_odd, log_dimension - 1, sampler);
  arb_poly_sub(c_odd, q_odd, c_odd, sampler->precision);
  multiply(c_odd, temporary, c_odd, log_dimension - 1, sampler);
  arb_poly_add(c_even, c_even, c_odd, sampler->precision);

  conjugate(f_odd, log_dimension - 1, sampler);
  multiply(temporary, temporary, f_odd, log_dimension - 1, sampler);
  arb_poly_sub(temporary, f_even, temporary, sampler->precision);
  if (!sample_recursive(q_even, temporary, c_even, log_dimension - 1,
                        depth + 1, sampler))
    goto cleanup;

  arb_poly_zero(result);
  for (slong i = 0; i < ((slong)1U << log_dimension); i++) {
    arb_poly_get_coeff_arb(sampler->temporary,
                           (i & 1) == 0 ? q_even : q_odd, i / 2);
    arb_poly_set_coeff_arb(result, i, sampler->temporary);
  }
  ok = 1;
cleanup:
  return ok;
}

static int sample_fixed_recursive(arb_poly_t result, const arb_poly_t center,
                                  const acc96_samplefz_plan_node_t *node,
                                  size_t log_dimension, size_t depth,
                                  acc96_samplefz_ctx_t *sampler) {
  if (log_dimension == 0) {
    int64_t sample;
    arb_poly_get_coeff_arb(sampler->temporary, center, 0);
    arb_poly_get_coeff_arb(sampler->auxiliary, node->variance, 0);
    if (!arb_is_finite(sampler->temporary) ||
        !arb_is_positive(sampler->auxiliary))
      return 0;
    record_radius(sampler, sampler->temporary);
    record_radius(sampler, sampler->auxiliary);
    arf_get_mpfr(sampler->leaf_center, arb_midref(sampler->temporary),
                 MPFR_RNDN);
    arf_get_mpfr(sampler->leaf_variance, arb_midref(sampler->auxiliary),
                 MPFR_RNDN);
    if (!acc96_gaussian_mpfr_sample_center_variance(
            &sample, sampler->leaf_center, sampler->leaf_variance,
            &sampler->gaussian_stats))
      return 0;
    arb_poly_zero(result);
    arb_poly_set_coeff_si(result, 0, (slong)sample);
    return 1;
  }

  arb_poly_struct *c_even = sampler->work[depth][0];
  arb_poly_struct *c_odd = sampler->work[depth][1];
  arb_poly_struct *q_even = sampler->work[depth][2];
  arb_poly_struct *q_odd = sampler->work[depth][3];
  arb_poly_struct *temporary = sampler->work[depth][4];
  split(c_even, c_odd, center, sampler);
  if (!sample_fixed_recursive(q_odd, c_odd, node->odd, log_dimension - 1,
                              depth + 1, sampler))
    return 0;
  arb_poly_sub(c_odd, q_odd, c_odd, sampler->precision);
  multiply(temporary, node->multiplier, c_odd, log_dimension - 1, sampler);
  arb_poly_add(c_even, c_even, temporary, sampler->precision);
  if (!sample_fixed_recursive(q_even, c_even, node->even, log_dimension - 1,
                              depth + 1, sampler))
    return 0;

  arb_poly_zero(result);
  for (slong i = 0; i < ((slong)1U << log_dimension); i++) {
    arb_poly_get_coeff_arb(sampler->temporary,
                           (i & 1) == 0 ? q_even : q_odd, i / 2);
    arb_poly_set_coeff_arb(result, i, sampler->temporary);
  }
  return 1;
}

static void set_signed_coefficient(acc96_coeff_t *out, int64_t value) {
  if (value < 0) {
    const uint64_t magnitude = (uint64_t)(-(value + 1)) + 1;
    acc96_coeff_from_u128(out, acc96_modulus() - magnitude);
  } else {
    acc96_coeff_from_u128(out, (uint64_t)value);
  }
}

int acc96_samplefz(acc96_poly_t *out, const arb_poly_t covariance,
                   const arb_poly_t center, acc96_samplefz_ctx_t *sampler) {
  arb_poly_t sampled;
  fmpz_t integer;
  int ok = 0;
  if (out == NULL || center == NULL || sampler == NULL)
    return 0;
  sampler->maximum_input_radius = 0.0;
  arb_poly_init(sampled);
  fmpz_init(integer);
  if (!sample_recursive(sampled, covariance, center, ACC96_SAMPLEFZ_LOG_N,
                        0, sampler))
    goto cleanup;
  acc96_poly_zero(out);
  for (size_t i = 0; i < ACC96_N; i++) {
    arb_poly_get_coeff_arb(sampler->temporary, sampled, (slong)i);
    if (!arb_get_unique_fmpz(integer, sampler->temporary) ||
        !fmpz_fits_si(integer))
      goto cleanup;
    set_signed_coefficient(&out->coeff[i], (int64_t)fmpz_get_si(integer));
  }
  ok = 1;
cleanup:
  fmpz_clear(integer);
  arb_poly_clear(sampled);
  return ok;
}

int acc96_samplefz_prepare_fixed(acc96_samplefz_ctx_t *sampler,
                                 const arb_poly_t covariance, size_t slot) {
  acc96_samplefz_plan_node_t *plan;
  if (sampler == NULL || covariance == NULL ||
      slot >= ACC96_SAMPLEFZ_FIXED_PLANS)
    return 0;
  plan = plan_build(covariance, ACC96_SAMPLEFZ_LOG_N, sampler);
  if (plan == NULL) return 0;
  plan_free(sampler->fixed[slot]);
  sampler->fixed[slot] = plan;
  return 1;
}

int acc96_samplefz_fixed(acc96_poly_t *out, const arb_poly_t center,
                         acc96_samplefz_ctx_t *sampler, size_t slot) {
  arb_poly_t sampled;
  fmpz_t integer;
  int ok = 0;
  if (out == NULL || center == NULL || sampler == NULL ||
      slot >= ACC96_SAMPLEFZ_FIXED_PLANS || sampler->fixed[slot] == NULL)
    return 0;
  sampler->maximum_input_radius = 0.0;
  arb_poly_init(sampled);
  fmpz_init(integer);
  if (!sample_fixed_recursive(sampled, center, sampler->fixed[slot],
                              ACC96_SAMPLEFZ_LOG_N, 0, sampler))
    goto cleanup;
  acc96_poly_zero(out);
  for (size_t i = 0; i < ACC96_N; i++) {
    arb_poly_get_coeff_arb(sampler->temporary, sampled, (slong)i);
    if (!arb_get_unique_fmpz(integer, sampler->temporary) ||
        !fmpz_fits_si(integer))
      goto cleanup;
    set_signed_coefficient(&out->coeff[i], (int64_t)fmpz_get_si(integer));
  }
  ok = 1;
cleanup:
  fmpz_clear(integer);
  arb_poly_clear(sampled);
  return ok;
}

const acc96_gaussian_mpfr_stats_t *acc96_samplefz_gaussian_stats(
    const acc96_samplefz_ctx_t *sampler) {
  return sampler == NULL ? NULL : &sampler->gaussian_stats;
}

double acc96_samplefz_maximum_input_radius(
    const acc96_samplefz_ctx_t *sampler) {
  return sampler == NULL ? INFINITY : sampler->maximum_input_radius;
}
