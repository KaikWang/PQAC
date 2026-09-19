#include "acc96_exact_lift.h"

#include <limits.h>
#include <stdlib.h>

#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

static const uint64_t lane_modulus[ACC96_LIMBS] = {
    ACC96_P0, ACC96_P1, ACC96_P2};

static void fmpz_set_i128(fmpz_t out, acc96_i128 value) {
  const int negative = value < 0;
  acc96_u128 magnitude = negative
      ? (acc96_u128)(-(value + 1)) + 1U : (acc96_u128)value;
  const ulong low = (ulong)(uint64_t)magnitude;
  const ulong high = (ulong)(uint64_t)(magnitude >> 64);
  if (negative)
    fmpz_neg_uiui(out, high, low);
  else
    fmpz_set_uiui(out, high, low);
}

static int fmpz_get_i128_checked(acc96_i128 *out, const fmpz_t value) {
  fmpz_t magnitude;
  ulong high;
  ulong low;
  int negative;
  int ok = 0;
  fmpz_init(magnitude);
  negative = fmpz_sgn(value) < 0;
  fmpz_abs(magnitude, value);
  if (fmpz_bits(magnitude) <= 127) {
    fmpz_get_uiui(&high, &low, magnitude);
    acc96_u128 assembled = ((acc96_u128)(uint64_t)high << 64) |
                           (acc96_u128)(uint64_t)low;
    *out = negative ? -(acc96_i128)assembled : (acc96_i128)assembled;
    ok = 1;
  }
  fmpz_clear(magnitude);
  return ok;
}

static int64_t lane_centered(const acc96_coeff_t *coefficient, size_t lane) {
  const uint64_t p = lane_modulus[lane];
  const uint64_t value = coefficient->limb[lane];
  return value > p / 2 ? (int64_t)(value - p) : (int64_t)value;
}

static void poly_from_lane(fmpz_poly_t out, const acc96_poly_t *in,
                           size_t lane) {
  fmpz_poly_zero(out);
  for (size_t i = 0; i < ACC96_N; i++)
    fmpz_poly_set_coeff_si(out, (slong)i,
                           (slong)lane_centered(&in->coeff[i], lane));
}

static void poly_from_common(fmpz_poly_t out, const acc96_poly_t *in,
                             const acc96_ctx_t *ctx) {
  fmpz_t value;
  fmpz_init(value);
  fmpz_poly_zero(out);
  for (size_t i = 0; i < ACC96_N; i++) {
    fmpz_set_i128(value, acc96_coeff_to_centered(&in->coeff[i], ctx));
    fmpz_poly_set_coeff_fmpz(out, (slong)i, value);
  }
  fmpz_clear(value);
}

static void negacyclic_muladd(fmpz_poly_t expression,
                              fmpz_poly_t lhs, fmpz_poly_t rhs,
                              fmpz_poly_t product, fmpz_t low,
                              fmpz_t high, fmpz_t term, int subtract) {
  fmpz_poly_mul(product, lhs, rhs);
  for (size_t k = 0; k < ACC96_N; k++) {
    fmpz_poly_get_coeff_fmpz(low, product, (slong)k);
    fmpz_poly_get_coeff_fmpz(high, product, (slong)(k + ACC96_N));
    fmpz_sub(term, low, high);
    if (subtract) fmpz_neg(term, term);
    fmpz_poly_get_coeff_fmpz(low, expression, (slong)k);
    fmpz_add(term, low, term);
    fmpz_poly_set_coeff_fmpz(expression, (slong)k, term);
  }
}

static int tag_is_one(const acc96_poly_t *tag) {
  for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
    const uint64_t expected = coefficient == 0 ? 1U : 0U;
    for (size_t lane = 0; lane < ACC96_LIMBS; lane++)
      if (tag->coeff[coefficient].limb[lane] != expected) return 0;
  }
  return 1;
}

static void scalar_common_poly_add(fmpz_poly_t expression,
                                   const acc96_poly_t *polynomial,
                                   int64_t scalar,
                                   const acc96_ctx_t *ctx,
                                   fmpz_t coefficient,
                                   fmpz_t term) {
  for (size_t k = 0; k < ACC96_N; k++) {
    fmpz_set_i128(coefficient,
                  acc96_coeff_to_centered(&polynomial->coeff[k], ctx));
    fmpz_mul_si(term, coefficient, (slong)scalar);
    fmpz_poly_get_coeff_fmpz(coefficient, expression, (slong)k);
    fmpz_add(term, coefficient, term);
    fmpz_poly_set_coeff_fmpz(expression, (slong)k, term);
  }
}

static int memver_exact_lift_lane(
    acc96_exact_lift_t *lift, unsigned *maximum_abs_bits, size_t lane,
    const acc96_relation_pp_t *pp, const acc96_poly_t *accumulator,
    uint32_t handle, const acc96_relation_witness_t *witness,
    const acc96_ctx_t *ctx) {
  fmpz_poly_t expression, lhs, rhs, product;
  fmpz_t low, high, term, difference, quotient;
  acc96_u128 scale;
  unsigned max_bits = 0;
  int ok = 1;
  fmpz_poly_init(expression);
  fmpz_poly_init(lhs);
  fmpz_poly_init(rhs);
  fmpz_poly_init(product);
  fmpz_init(low);
  fmpz_init(high);
  fmpz_init(term);
  fmpz_init(difference);
  fmpz_init(quotient);

  poly_from_common(expression, &witness->left, ctx);
  poly_from_lane(lhs, &pp->atrap, lane);
  poly_from_common(rhs, &witness->v12, ctx);
  negacyclic_muladd(expression, lhs, rhs, product, low, high, term, 0);

  scale = (acc96_u128)ACC96_TRUNC_BASE * ACC96_TRUNC_BASE;
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++) {
    poly_from_lane(lhs, &pp->trunc_b[i], lane);
    poly_from_common(rhs, &witness->v2[i], ctx);
    negacyclic_muladd(expression, lhs, rhs, product, low, high, term, 1);

    {
      const uint64_t lane_scale = (uint64_t)(scale % lane_modulus[lane]);
      const int64_t centered_scale = lane_scale > lane_modulus[lane] / 2
          ? (int64_t)(lane_scale - lane_modulus[lane])
          : (int64_t)lane_scale;
      scalar_common_poly_add(expression, &witness->v2[i], centered_scale,
                             ctx, low, term);
    }
    scale *= ACC96_TRUNC_BASE;
  }

  for (size_t bit = 0; bit < ACC96_HANDLE_BITS; bit++) {
    scale = 1;
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++) {
      poly_from_lane(lhs, &pp->ba.bit[bit].block[digit], lane);
      poly_from_common(rhs, &witness->right.right[bit][digit], ctx);
      negacyclic_muladd(expression, lhs, rhs, product, low, high, term, 0);
      if ((handle >> bit) & 1U) {
        const uint64_t lane_scale =
            (uint64_t)(scale % lane_modulus[lane]);
        const int64_t centered_scale =
            lane_scale > lane_modulus[lane] / 2
                ? (int64_t)(lane_scale - lane_modulus[lane])
                : (int64_t)lane_scale;
        /* Keep the private handle factor separate from the public Ba
         * coefficient.  This is congruent to centering (Ba-xG) first, but
         * gives LaBRADOR the clean quadratic term
         *   -center(G_j) * x_bit * right[bit][j]. */
        scalar_common_poly_add(expression,
                               &witness->right.right[bit][digit],
                               -centered_scale, ctx, low, term);
      }
      scale *= ACC96_GADGET_BASE;
    }
  }

  for (size_t k = 0; k < ACC96_N && ok; k++) {
    fmpz_poly_get_coeff_fmpz(term, expression, (slong)k);
    fmpz_set_si(difference,
                (slong)lane_centered(&accumulator->coeff[k], lane));
    fmpz_sub(difference, difference, term);
    if (!fmpz_divisible_si(difference, (slong)lane_modulus[lane])) {
      ok = 0;
      break;
    }
    fmpz_divexact_ui(quotient, difference, (ulong)lane_modulus[lane]);
    {
      const unsigned bits = (unsigned)fmpz_bits(quotient);
      if (bits > max_bits) max_bits = bits;
    }
    if (!fmpz_get_i128_checked(&lift->coeff[lane][k], quotient)) {
      ok = 0;
      break;
    }
  }

  if (maximum_abs_bits != NULL) *maximum_abs_bits = max_bits;
  fmpz_clear(quotient);
  fmpz_clear(difference);
  fmpz_clear(term);
  fmpz_clear(high);
  fmpz_clear(low);
  fmpz_poly_clear(product);
  fmpz_poly_clear(rhs);
  fmpz_poly_clear(lhs);
  fmpz_poly_clear(expression);
  return ok;
}

int acc96_memver_exact_lift(
    acc96_exact_lift_t *lift, unsigned *maximum_abs_bits,
    const acc96_relation_pp_t *pp, const acc96_poly_t *accumulator,
    uint32_t handle, const acc96_relation_witness_t *witness,
    const acc96_ctx_t *ctx) {
  unsigned max_bits = 0;
  int ok = lift != NULL && pp != NULL && accumulator != NULL &&
           witness != NULL && ctx != NULL && ctx->initialized &&
           acc96_window_proof_witness_check(witness, ctx, NULL) &&
           tag_is_one(&witness->tag);
  if (!ok) return 0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(ACC96_LIMBS) schedule(static) \
    reduction(&:ok) reduction(max:max_bits)
#endif
  for (size_t lane = 0; lane < ACC96_LIMBS; lane++) {
    unsigned lane_bits = 0;
    const int lane_ok = memver_exact_lift_lane(
        lift, &lane_bits, lane, pp, accumulator, handle, witness, ctx);
    ok &= lane_ok;
    if (lane_bits > max_bits) max_bits = lane_bits;
  }
  if (maximum_abs_bits != NULL) *maximum_abs_bits = max_bits;
  return ok;
}
