#include "acc96_exact_lift.h"

#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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
                              fmpz_t high, fmpz_t term, int subtract,
                              fmpz_poly_t upper, int bulk) {
  fmpz_poly_t reference;
  if (bulk == 2) { fmpz_poly_init(reference); fmpz_poly_set(reference, expression); }
  fmpz_poly_mul(product, lhs, rhs);
  if (bulk) {
    /* deg(lhs), deg(rhs) < N: reduction modulo X^N+1 is low-high.
       All operations below are exact integer operations, with reusable storage. */
    fmpz_poly_shift_right(upper, product, ACC96_N);
    fmpz_poly_truncate(product, ACC96_N);
    fmpz_poly_sub(product, product, upper);
    if (subtract) fmpz_poly_sub(expression, expression, product);
    else fmpz_poly_add(expression, expression, product);
    if (bulk != 2) return;
    /* Recompute the original product for the coefficient-loop audit. */
    fmpz_poly_mul(product, lhs, rhs);
  }
  for (size_t k = 0; k < ACC96_N; k++) {
    fmpz_poly_get_coeff_fmpz(low, product, (slong)k);
    fmpz_poly_get_coeff_fmpz(high, product, (slong)(k + ACC96_N));
    fmpz_sub(term, low, high);
    if (subtract) fmpz_neg(term, term);
    fmpz_poly_get_coeff_fmpz(low, bulk == 2 ? reference : expression, (slong)k);
    fmpz_add(term, low, term);
    fmpz_poly_set_coeff_fmpz(bulk == 2 ? reference : expression, (slong)k, term);
  }
  if (bulk == 2) {
    if (!fmpz_poly_equal(reference, expression)) abort();
    fmpz_poly_clear(reference);
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

static void scalar_cached_poly_add(fmpz_poly_t expression,
    const acc96_poly_t *polynomial, const fmpz_poly_t converted,
    int64_t scalar, const acc96_ctx_t *ctx, fmpz_t coefficient, fmpz_t term,
    fmpz_poly_t scratch, int bulk) {
  if (!bulk) {
    scalar_common_poly_add(expression, polynomial, scalar, ctx, coefficient, term);
    return;
  }
  fmpz_poly_t reference;
  if (bulk == 2) {
    fmpz_poly_init(reference); fmpz_poly_set(reference, expression);
    scalar_common_poly_add(reference, polynomial, scalar, ctx, coefficient, term);
  }
  fmpz_poly_scalar_mul_si(scratch, converted, (slong)scalar);
  fmpz_poly_add(expression, expression, scratch);
  if (bulk == 2) {
    if (!fmpz_poly_equal(reference, expression)) abort();
    fmpz_poly_clear(reference);
  }
}

#define LIFT_COMMON_RIGHT_START (2U + ACC96_TRUNC_RIGHT)
#define LIFT_COMMON_REFRESH_START (LIFT_COMMON_RIGHT_START + ACC96_HANDLE_BITS * ACC96_GADGET_DIGITS)
#define LIFT_COMMON_COUNT (LIFT_COMMON_REFRESH_START + ACC96_REFRESH_BLOCKS)

static const acc96_poly_t *lift_common_source(const acc96_relation_witness_t *w, size_t index) {
  if (index == 0) return &w->left;
  if (index == 1) return &w->v12;
  if (index < LIFT_COMMON_RIGHT_START) return &w->v2[index-2];
  if (index < LIFT_COMMON_REFRESH_START) {
    index -= LIFT_COMMON_RIGHT_START;
    return &w->right.right[index/ACC96_GADGET_DIGITS][index%ACC96_GADGET_DIGITS];
  }
  return &w->refresh_z[index-LIFT_COMMON_REFRESH_START];
}

static void lift_common_load(fmpz_poly_t out, const acc96_relation_witness_t *w,
    size_t index, const acc96_ctx_t *ctx, const fmpz_poly_struct *common) {
  if (common) fmpz_poly_set(out, &common[index]);
  else poly_from_common(out, lift_common_source(w,index), ctx);
}

static int memver_exact_lift_lane(
    acc96_exact_lift_t *lift, unsigned *maximum_abs_bits, size_t lane,
    const acc96_relation_pp_t *pp, const acc96_poly_t *accumulator,
    uint32_t handle, const acc96_relation_witness_t *witness,
    const acc96_ctx_t *ctx, const acc96_poly_t *shared_refresh,
    const fmpz_poly_struct *common) {
  fmpz_poly_t expression, lhs, rhs, product, upper;
  fmpz_t low, high, term, difference, quotient;
  acc96_poly_t *refresh_matrix = shared_refresh ? NULL : malloc(sizeof(*refresh_matrix));
  acc96_poly_t *refresh_expanded = shared_refresh ? NULL : malloc(sizeof(*refresh_expanded));
  acc96_u128 scale;
  unsigned max_bits = 0;
  const char *bulk_value = getenv("REVAC_Q88_BULK_LIFT");
  const int bulk = bulk_value && bulk_value[0] == '2' ? 2 :
      (!bulk_value || bulk_value[0] == '1' ? 1 : 0);
  int ok = shared_refresh || (refresh_matrix != NULL && refresh_expanded != NULL);
  if (!ok) {
    free(refresh_expanded);
    free(refresh_matrix);
    return 0;
  }
  fmpz_poly_init(expression);
  fmpz_poly_init(lhs);
  fmpz_poly_init(rhs);
  fmpz_poly_init(product);
  fmpz_poly_init(upper);
  fmpz_init(low);
  fmpz_init(high);
  fmpz_init(term);
  fmpz_init(difference);
  fmpz_init(quotient);

  lift_common_load(expression,witness,0,ctx,common);
  poly_from_lane(lhs, &pp->atrap, lane);
  lift_common_load(rhs,witness,1,ctx,common);
  negacyclic_muladd(expression, lhs, rhs, product, low, high, term, 0, upper, bulk);

  scale = (acc96_u128)ACC96_TRUNC_BASE * ACC96_TRUNC_BASE;
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++) {
    poly_from_lane(lhs, &pp->trunc_b[i], lane);
    lift_common_load(rhs,witness,2+i,ctx,common);
    negacyclic_muladd(expression, lhs, rhs, product, low, high, term, 1, upper, bulk);

    {
      const uint64_t lane_scale = (uint64_t)(scale % lane_modulus[lane]);
      const int64_t centered_scale = lane_scale > lane_modulus[lane] / 2
          ? (int64_t)(lane_scale - lane_modulus[lane])
          : (int64_t)lane_scale;
      scalar_cached_poly_add(expression, &witness->v2[i], rhs, centered_scale,
                             ctx, low, term, upper, bulk);
    }
    scale *= ACC96_TRUNC_BASE;
  }

  for (size_t bit = 0; bit < ACC96_HANDLE_BITS; bit++) {
    scale = 1;
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++) {
      poly_from_lane(lhs, &pp->ba.bit[bit].block[digit], lane);
      lift_common_load(rhs,witness,LIFT_COMMON_RIGHT_START+bit*ACC96_GADGET_DIGITS+digit,ctx,common);
      negacyclic_muladd(expression, lhs, rhs, product, low, high, term, 0, upper, bulk);
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
        scalar_cached_poly_add(expression,
                               &witness->right.right[bit][digit], rhs,
                               -centered_scale, ctx, low, term, upper, bulk);
      }
      scale *= ACC96_GADGET_BASE;
    }
  }

  for (size_t block = 0; block < ACC96_REFRESH_BLOCKS; block++) {
    if (shared_refresh) {
      poly_from_lane(lhs, &shared_refresh[block], lane);
    } else {
      acc96_relation_refresh_matrix_poly(refresh_matrix, pp, 0, block);
      for (size_t bit = 0; bit < ACC96_REFRESH_BITS; bit++) {
        if (((witness->refresh_r[bit >> 3] >> (bit & 7U)) & 1U) == 0U) continue;
        acc96_relation_refresh_matrix_poly(refresh_expanded, pp, bit + 1U, block);
        acc96_poly_add(refresh_matrix, refresh_matrix, refresh_expanded);
      }
      poly_from_lane(lhs, refresh_matrix, lane);
    }
    lift_common_load(rhs,witness,LIFT_COMMON_REFRESH_START+block,ctx,common);
    negacyclic_muladd(expression, lhs, rhs, product, low, high, term, 0, upper, bulk);
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
  if (bulk == 2 && ok) fprintf(stderr,"[q88 bulk lift audit] lane=%zu PASS\n",lane);
  fmpz_clear(quotient);
  fmpz_clear(difference);
  fmpz_clear(term);
  fmpz_clear(high);
  fmpz_clear(low);
  fmpz_poly_clear(upper);
  fmpz_poly_clear(product);
  fmpz_poly_clear(rhs);
  fmpz_poly_clear(lhs);
  fmpz_poly_clear(expression);
  free(refresh_expanded);
  free(refresh_matrix);
  return ok;
}

static int acc96_memver_exact_lift_impl(
    acc96_exact_lift_t *lift, unsigned *maximum_abs_bits,
    const acc96_relation_pp_t *pp, const acc96_poly_t *accumulator,
    uint32_t handle, const acc96_relation_witness_t *witness,
    const acc96_ctx_t *ctx, int shared, int shared_common) {
  unsigned max_bits = 0;
  int ok = lift != NULL && pp != NULL && accumulator != NULL &&
           witness != NULL && ctx != NULL && ctx->initialized &&
           acc96_window_proof_witness_check(witness, ctx, NULL) &&
           tag_is_one(&witness->tag);
  if (!ok) return 0;
  acc96_poly_t *refresh = NULL;
  if (shared) {
    refresh = malloc(ACC96_REFRESH_BLOCKS * sizeof(*refresh));
    if (!refresh) return 0;
    int expanded_ok = 1;
#ifdef _OPENMP
#pragma omp parallel for num_threads(ACC96_LIMBS) schedule(static) reduction(&:expanded_ok)
#endif
    for (size_t block=0; block<ACC96_REFRESH_BLOCKS; block++) {
      acc96_poly_t *term = malloc(sizeof(*term));
      if (!term) { expanded_ok=0; continue; }
      acc96_relation_refresh_matrix_poly(&refresh[block], pp, 0, block);
      for (size_t bit=0; bit<ACC96_REFRESH_BITS; bit++) {
        if (((witness->refresh_r[bit>>3] >> (bit & 7U)) & 1U) == 0U) continue;
        acc96_relation_refresh_matrix_poly(term, pp, bit+1U, block);
        acc96_poly_add(&refresh[block], &refresh[block], term);
      }
      free(term);
    }
    if (!expanded_ok) { free(refresh); return 0; }
  }

  fmpz_poly_struct *common=NULL;
  if (shared_common) {
    common=calloc(LIFT_COMMON_COUNT,sizeof(*common));
    if (!common) { free(refresh); return 0; }
    /* Per-call only: private converted values are cleared before return.
       All three CRT lanes read the same bounded common witness. */
#ifdef _OPENMP
#pragma omp parallel for num_threads(ACC96_LIMBS) schedule(static)
#endif
    for (size_t i=0;i<LIFT_COMMON_COUNT;i++) {
      fmpz_poly_init(&common[i]);
      poly_from_common(&common[i],lift_common_source(witness,i),ctx);
    }
  }

#ifdef _OPENMP
#pragma omp parallel for num_threads(ACC96_LIMBS) schedule(static) \
    reduction(&:ok) reduction(max:max_bits)
#endif
  for (size_t lane = 0; lane < ACC96_LIMBS; lane++) {
    unsigned lane_bits = 0;
    const int lane_ok = memver_exact_lift_lane(
        lift, &lane_bits, lane, pp, accumulator, handle, witness, ctx, refresh, common);
    ok &= lane_ok;
    if (lane_bits > max_bits) max_bits = lane_bits;
  }
  if (maximum_abs_bits != NULL) *maximum_abs_bits = max_bits;
  if (common) {
    for (size_t i=0;i<LIFT_COMMON_COUNT;i++) fmpz_poly_clear(&common[i]);
    free(common);
  }
  free(refresh);
  return ok;
}

int acc96_memver_exact_lift(
    acc96_exact_lift_t *lift, unsigned *maximum_abs_bits,
    const acc96_relation_pp_t *pp, const acc96_poly_t *accumulator,
    uint32_t handle, const acc96_relation_witness_t *witness,
    const acc96_ctx_t *ctx) {
  const char *mode=getenv("REVAC_Q88_LIFT_SHARED_REFRESH");
  /* Share only deterministic per-call preparation across CRT lanes.  Mode 0
     retains the reference path; mode 4 audits the default against it. */
  if (!mode || !mode[0]) mode="3";
  const int shared=mode && mode[0]>='1' && mode[0]<='4';
  const int shared_common=mode && (mode[0]=='3' || mode[0]=='4');
  unsigned bits=0;
  const int ok=acc96_memver_exact_lift_impl(lift,&bits,pp,accumulator,handle,witness,ctx,shared,shared_common);
  if (maximum_abs_bits) *maximum_abs_bits=bits;
  if (mode && (mode[0]=='2' || mode[0]=='4')) {
    acc96_exact_lift_t *reference=malloc(sizeof(*reference));
    if (!reference) return 0;
    unsigned reference_bits=0;
    const int reference_ok=acc96_memver_exact_lift_impl(reference,&reference_bits,pp,
        accumulator,handle,witness,ctx,0,0);
    if (reference_ok!=ok || (ok && (reference_bits!=bits ||
        memcmp(reference,lift,sizeof(*reference))))) {
      fprintf(stderr,"q88 shared-refresh exact lift mismatch\n");abort();
    }
    free(reference);
    fprintf(stderr,"[q88 shared-refresh audit] result=%d PASS\n",ok);
  }
  return ok;
}
