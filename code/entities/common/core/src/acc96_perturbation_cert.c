#include "acc96_perturbation_cert.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

#include "flint/acb.h"
#include "flint/acb_dft.h"
#include "flint/arb.h"
#include "acc96_truncated_params_generated.h"

static void update_minimum(arf_t minimum, int *initialized,
                           const arb_t value, slong precision) {
  arf_t lower;
  arf_init(lower);
  arb_get_lbound_arf(lower, value, precision);
  if (!*initialized || arf_cmp(lower, minimum) < 0) {
    arf_set(minimum, lower);
    *initialized = 1;
  }
  arf_clear(lower);
}

static void update_maximum(arf_t maximum, const arb_t value,
                           slong precision) {
  arf_t upper;
  arf_init(upper);
  arb_get_ubound_arf(upper, value, precision);
  if (arf_cmp(upper, maximum) > 0) arf_set(maximum, upper);
  arf_clear(upper);
}

static int set_decimal(arb_t output, const char *decimal, slong precision) {
  return arb_set_str(output, decimal, precision) == 0;
}

int acc96_perturbation_covariance_certify_acb(
    const acc96_trapdoor_t *trapdoor, const acc96_ctx_t *ctx,
    long precision_bits, acc96_perturbation_covariance_cert_t *certificate) {
  acb_ptr row0 = NULL;
  acb_ptr row1 = NULL;
  arb_ptr rr00 = NULL;
  arb_ptr rr11 = NULL;
  acb_ptr rr01 = NULL;
  acb_t twist, twist_step, conjugate, product, c01;
  arb_t norm_square, scalar, s1_square, s2_square, s3_square;
  arb_t negative_sg_square, c00, c11, cross_square, quotient;
  arb_t adjusted_schur, independent_variance, row_energy;
  arf_t minimum_pivot, minimum_schur, minimum_independent;
  arf_t maximum_row_energy;
  int minimum_pivot_initialized = 0;
  int minimum_schur_initialized = 0;
  int minimum_independent_initialized = 0;
  int certified = 1;
  const slong precision = precision_bits >= 128 ? (slong)precision_bits : 128;

  if (trapdoor == NULL || ctx == NULL || !ctx->initialized ||
      certificate == NULL)
    return 0;
  memset(certificate, 0, sizeof(*certificate));

  row0 = _acb_vec_init(ACC96_N);
  row1 = _acb_vec_init(ACC96_N);
  rr00 = _arb_vec_init(ACC96_N);
  rr11 = _arb_vec_init(ACC96_N);
  rr01 = _acb_vec_init(ACC96_N);
  if (row0 == NULL || row1 == NULL || rr00 == NULL || rr11 == NULL ||
      rr01 == NULL) {
    certified = 0;
    goto vector_cleanup;
  }

  acb_init(twist);
  acb_init(twist_step);
  acb_init(conjugate);
  acb_init(product);
  acb_init(c01);
  arb_init(norm_square);
  arb_init(scalar);
  arb_init(s1_square);
  arb_init(s2_square);
  arb_init(s3_square);
  arb_init(negative_sg_square);
  arb_init(c00);
  arb_init(c11);
  arb_init(cross_square);
  arb_init(quotient);
  arb_init(adjusted_schur);
  arb_init(independent_variance);
  arb_init(row_energy);
  arf_init(minimum_pivot);
  arf_init(minimum_schur);
  arf_init(minimum_independent);
  arf_init(maximum_row_energy);
  arf_zero(maximum_row_energy);

  if (!set_decimal(scalar, ACC_TRUNC_SGINVSQ_S4INVSQ_HP, precision) ||
      !set_decimal(s1_square, ACC_TRUNC_S1SQ_HP, precision) ||
      !set_decimal(s2_square, ACC_TRUNC_S2SQ_HP, precision) ||
      !set_decimal(s3_square, ACC_TRUNC_S3SQ_HP, precision) ||
      !set_decimal(negative_sg_square, ACC_TRUNC_NEG_SGSQ_HP, precision)) {
    certified = 0;
    goto scalar_cleanup;
  }

  for (size_t frequency = 0; frequency < ACC96_N; frequency++) {
    arb_zero(&rr00[frequency]);
    arb_zero(&rr11[frequency]);
    acb_zero(&rr01[frequency]);
  }
  acb_unit_root(twist_step, 2U * ACC96_N, precision);

  for (size_t column = 0; column < ACC96_TRUNC_RIGHT; column++) {
    acb_one(twist);
    for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
      const acc96_i128 wide0 = acc96_coeff_to_centered(
          &trapdoor->r[0][column].coeff[coefficient], ctx);
      const acc96_i128 wide1 = acc96_coeff_to_centered(
          &trapdoor->r[1][column].coeff[coefficient], ctx);
      if (wide0 < -1 || wide0 > 1 || wide1 < -1 || wide1 > 1) {
        certified = 0;
        goto scalar_cleanup;
      }
      acb_mul_si(&row0[coefficient], twist, (slong)wide0, precision);
      acb_mul_si(&row1[coefficient], twist, (slong)wide1, precision);
      acb_mul(twist, twist, twist_step, precision);
    }
    acb_dft_rad2_inplace(row0, 12, precision);
    acb_dft_rad2_inplace(row1, 12, precision);
    for (size_t frequency = 0; frequency < ACC96_N; frequency++) {
      arb_mul(norm_square, acb_realref(&row0[frequency]),
              acb_realref(&row0[frequency]), precision);
      arb_addmul(norm_square, acb_imagref(&row0[frequency]),
                 acb_imagref(&row0[frequency]), precision);
      arb_add(&rr00[frequency], &rr00[frequency], norm_square, precision);

      arb_mul(norm_square, acb_realref(&row1[frequency]),
              acb_realref(&row1[frequency]), precision);
      arb_addmul(norm_square, acb_imagref(&row1[frequency]),
                 acb_imagref(&row1[frequency]), precision);
      arb_add(&rr11[frequency], &rr11[frequency], norm_square, precision);

      acb_conj(conjugate, &row1[frequency]);
      acb_mul(product, &row0[frequency], conjugate, precision);
      acb_add(&rr01[frequency], &rr01[frequency], product, precision);
    }
  }

  arb_add(independent_variance, s2_square, negative_sg_square, precision);
  if (!arb_is_positive(independent_variance)) certified = 0;
  update_minimum(minimum_independent, &minimum_independent_initialized,
                 independent_variance, precision);

  for (size_t frequency = 0; frequency < ACC96_N; frequency++) {
    arb_mul(c00, scalar, &rr00[frequency], precision);
    arb_add(c00, c00, s1_square, precision);
    arb_mul(c11, scalar, &rr11[frequency], precision);
    arb_add(c11, c11, s3_square, precision);
    acb_mul_arb(c01, &rr01[frequency], scalar, precision);
    arb_mul(cross_square, acb_realref(c01), acb_realref(c01), precision);
    arb_addmul(cross_square, acb_imagref(c01), acb_imagref(c01),
               precision);

    if (!arb_is_positive(c11)) {
      certified = 0;
      continue;
    }
    arb_div(quotient, cross_square, c11, precision);
    arb_add(adjusted_schur, c00, negative_sg_square, precision);
    arb_sub(adjusted_schur, adjusted_schur, quotient, precision);
    if (!arb_is_positive(adjusted_schur)) certified = 0;

    update_minimum(minimum_pivot, &minimum_pivot_initialized, c11,
                   precision);
    update_minimum(minimum_schur, &minimum_schur_initialized,
                   adjusted_schur, precision);
    arb_max(row_energy, &rr00[frequency], &rr11[frequency], precision);
    update_maximum(maximum_row_energy, row_energy, precision);
  }

  certificate->precision_bits = (long)precision;
  certificate->minimum_second_pivot_lower_bound =
      arf_get_d(minimum_pivot, ARF_RND_FLOOR);
  certificate->minimum_adjusted_schur_lower_bound =
      arf_get_d(minimum_schur, ARF_RND_FLOOR);
  certificate->independent_variance_lower_bound =
      arf_get_d(minimum_independent, ARF_RND_FLOOR);
  certificate->maximum_row_energy_upper_bound =
      arf_get_d(maximum_row_energy, ARF_RND_CEIL);
  if (!minimum_pivot_initialized || !minimum_schur_initialized ||
      !minimum_independent_initialized ||
      !isfinite(certificate->minimum_second_pivot_lower_bound) ||
      !isfinite(certificate->minimum_adjusted_schur_lower_bound) ||
      !isfinite(certificate->independent_variance_lower_bound) ||
      !isfinite(certificate->maximum_row_energy_upper_bound))
    certified = 0;

scalar_cleanup:
  arf_clear(maximum_row_energy);
  arf_clear(minimum_independent);
  arf_clear(minimum_schur);
  arf_clear(minimum_pivot);
  arb_clear(row_energy);
  arb_clear(independent_variance);
  arb_clear(adjusted_schur);
  arb_clear(quotient);
  arb_clear(cross_square);
  arb_clear(c11);
  arb_clear(c00);
  arb_clear(negative_sg_square);
  arb_clear(s3_square);
  arb_clear(s2_square);
  arb_clear(s1_square);
  arb_clear(scalar);
  arb_clear(norm_square);
  acb_clear(c01);
  acb_clear(product);
  acb_clear(conjugate);
  acb_clear(twist_step);
  acb_clear(twist);
vector_cleanup:
  if (rr01 != NULL) _acb_vec_clear(rr01, ACC96_N);
  if (rr11 != NULL) _arb_vec_clear(rr11, ACC96_N);
  if (rr00 != NULL) _arb_vec_clear(rr00, ACC96_N);
  if (row1 != NULL) _acb_vec_clear(row1, ACC96_N);
  if (row0 != NULL) _acb_vec_clear(row0, ACC96_N);
  return certified;
}
