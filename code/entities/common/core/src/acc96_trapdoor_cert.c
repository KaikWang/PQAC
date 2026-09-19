#include "acc96_trapdoor_cert.h"

#include <math.h>
#include <stddef.h>

#include "flint/acb.h"
#include "flint/acb_dft.h"
#include "flint/arb.h"
#include "acc96_truncated_params_generated.h"

int acc96_trapdoor_spectral_certify_acb(
    const acc96_trapdoor_t *trapdoor, const acc96_ctx_t *ctx,
    long precision_bits, double *certified_upper_bound) {
  acb_ptr values = NULL;
  arb_ptr energy = NULL;
  acb_t twist, twist_step;
  arb_t squared, theorem_bound;
  arf_t upper, maximum_upper;
  int certified = 1;
  const slong precision = precision_bits >= 96 ? (slong)precision_bits : 96;
  if (trapdoor == NULL || ctx == NULL || !ctx->initialized ||
      certified_upper_bound == NULL)
    return 0;
  values = _acb_vec_init(ACC96_N);
  energy = _arb_vec_init(ACC96_N);
  if (values == NULL || energy == NULL) {
    if (energy != NULL) _arb_vec_clear(energy, ACC96_N);
    if (values != NULL) _acb_vec_clear(values, ACC96_N);
    return 0;
  }
  acb_init(twist);
  acb_init(twist_step);
  arb_init(squared);
  arb_init(theorem_bound);
  arf_init(upper);
  arf_init(maximum_upper);
  arf_zero(maximum_upper);
  arb_set_d(theorem_bound, ACC_TRUNC_R_MAX_SQ_SPECTRAL_NORM);
  acb_unit_root(twist_step, 2U * ACC96_N, precision);

  for (size_t row = 0; row < 2; row++) {
    for (size_t frequency = 0; frequency < ACC96_N; frequency++)
      arb_zero(&energy[frequency]);
    for (size_t column = 0; column < ACC96_TRUNC_RIGHT; column++) {
      acb_one(twist);
      for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
        const acc96_i128 wide = acc96_coeff_to_centered(
            &trapdoor->r[row][column].coeff[coefficient], ctx);
        if (wide < -1 || wide > 1) {
          certified = 0;
          goto cleanup;
        }
        acb_mul_si(&values[coefficient], twist, (slong)wide, precision);
        acb_mul(twist, twist, twist_step, precision);
      }
      acb_dft_rad2_inplace(values, 12, precision);
      for (size_t frequency = 0; frequency < ACC96_N; frequency++) {
        arb_mul(squared, acb_realref(&values[frequency]),
                acb_realref(&values[frequency]), precision);
        arb_addmul(squared, acb_imagref(&values[frequency]),
                   acb_imagref(&values[frequency]), precision);
        arb_add(&energy[frequency], &energy[frequency], squared, precision);
      }
    }
    for (size_t frequency = 0; frequency < ACC96_N; frequency++) {
      if (!arb_lt(&energy[frequency], theorem_bound)) certified = 0;
      arb_get_ubound_arf(upper, &energy[frequency], precision);
      if (arf_cmp(upper, maximum_upper) > 0) arf_set(maximum_upper, upper);
    }
  }
  *certified_upper_bound = arf_get_d(maximum_upper, ARF_RND_CEIL);
  if (!isfinite(*certified_upper_bound)) certified = 0;

cleanup:
  arf_clear(maximum_upper);
  arf_clear(upper);
  arb_clear(theorem_bound);
  arb_clear(squared);
  acb_clear(twist_step);
  acb_clear(twist);
  _arb_vec_clear(energy, ACC96_N);
  _acb_vec_clear(values, ACC96_N);
  return certified;
}
