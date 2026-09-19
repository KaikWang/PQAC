#include "acc96_sampleleft.h"

#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "flint/arb_poly.h"
#include "flint/fmpz.h"

#include "acc96_gaussian_mpfr.h"
#include "acc96_klein_exact.h"
#include "acc96_perturbation_cert.h"
#include "acc96_samplefz.h"
#include "acc96_truncated_params_generated.h"
#include "randombytes.h"

struct acc96_sampleleft_ctx {
  const acc96_trapdoor_t *trapdoor;
  acc96_samplefz_ctx_t *samplefz;
  acc96_klein_exact_table_t *klein;
  arb_poly_t covariance_left;
  arb_poly_t covariance_second;
  arb_poly_t covariance_independent;
  arb_poly_t schur_multiplier;
  long precision;
  size_t threads;
  acc96_gaussian_mpfr_stats_t *p2_stats;
  acc96_gaussian_mpfr_stats_t *klein_stats;
};

static void set_signed(acc96_coeff_t *out, int64_t value) {
  static const uint64_t modulus[ACC96_LIMBS] = {
      ACC96_P0, ACC96_P1, ACC96_P2};
  if (value < 0) {
    const uint64_t magnitude = (uint64_t)(-(value + 1)) + 1;
    for (size_t lane = 0; lane < ACC96_LIMBS; lane++) {
      const uint64_t reduced = magnitude < modulus[lane]
          ? magnitude : magnitude % modulus[lane];
      out->limb[lane] = reduced == 0 ? 0 : modulus[lane] - reduced;
    }
  } else {
    const uint64_t magnitude = (uint64_t)value;
    for (size_t lane = 0; lane < ACC96_LIMBS; lane++)
      out->limb[lane] = magnitude < modulus[lane]
          ? magnitude : magnitude % modulus[lane];
  }
}

static void fmpz_set_i128_local(fmpz_t out, acc96_i128 value) {
  const int negative = value < 0;
  const acc96_u128 magnitude = negative
      ? (acc96_u128)(-(value + 1)) + 1 : (acc96_u128)value;
  fmpz_set_ui(out, (uint64_t)(magnitude >> 64));
  fmpz_mul_2exp(out, out, 64);
  fmpz_add_ui(out, out, (uint64_t)magnitude);
  if (negative) fmpz_neg(out, out);
}

static void arb_poly_from_q(arb_poly_t out, const acc96_poly_t *input,
                            const acc96_ctx_t *ctx) {
  fmpz_t integer;
  arb_t value;
  fmpz_init(integer);
  arb_init(value);
  arb_poly_zero(out);
  for (size_t i = 0; i < ACC96_N; i++) {
    fmpz_set_i128_local(integer,
        acc96_coeff_to_centered(&input->coeff[i], ctx));
    arb_set_fmpz(value, integer);
    arb_poly_set_coeff_arb(out, (slong)i, value);
  }
  arb_clear(value);
  fmpz_clear(integer);
}

static void arb_add_constant(arb_poly_t polynomial, const arb_t value,
                             long precision) {
  arb_t coefficient;
  arb_init(coefficient);
  arb_poly_get_coeff_arb(coefficient, polynomial, 0);
  arb_add(coefficient, coefficient, value, precision);
  arb_poly_set_coeff_arb(polynomial, 0, coefficient);
  arb_clear(coefficient);
}

static void q_conjugate(acc96_poly_t *out, const acc96_poly_t *input) {
  acc96_poly_zero(out);
  out->coeff[0] = input->coeff[0];
  for (size_t i = 1; i < ACC96_N; i++) {
    const acc96_coeff_t *source = &input->coeff[i];
    acc96_coeff_t *target = &out->coeff[ACC96_N - i];
    for (size_t lane = 0; lane < ACC96_LIMBS; lane++) {
      static const uint64_t modulus[ACC96_LIMBS] = {
          ACC96_P0, ACC96_P1, ACC96_P2};
      target->limb[lane] = source->limb[lane] == 0
          ? 0 : modulus[lane] - source->limb[lane];
    }
  }
}

static void q_muladd(acc96_poly_t *accumulator, const acc96_poly_t *left,
                     const acc96_poly_t *right, const acc96_ctx_t *ctx) {
  /* acc96_poly_add is coefficient-wise and therefore supports an in-place
   * destination.  The old implementation allocated two ~96 KiB objects for
   * every multiply-add (one product and one sum).  Keeping only the product
   * on the stack removes allocator traffic without changing any arithmetic. */
  uint64_t ntt_left[ACC96_N];
  uint64_t ntt_right[ACC96_N];
  acc96_poly_t product;
  acc96_poly_mul_workspace(&product, left, right, ctx, ntt_left, ntt_right);
  acc96_poly_add(accumulator, accumulator, &product);
}

static void q_mulsub(acc96_poly_t *accumulator, const acc96_poly_t *left,
                     const acc96_poly_t *right, const acc96_ctx_t *ctx) {
  /* As above, subtraction is also safe in place because each coefficient is
   * read before it is overwritten. */
  uint64_t ntt_left[ACC96_N];
  uint64_t ntt_right[ACC96_N];
  acc96_poly_t product;
  acc96_poly_mul_workspace(&product, left, right, ctx, ntt_left, ntt_right);
  acc96_poly_sub(accumulator, accumulator, &product);
}

static void q_scale(acc96_poly_t *out, const acc96_poly_t *input,
                    acc96_u128 scalar) {
  static const uint64_t modulus[ACC96_LIMBS] = {
      ACC96_P0, ACC96_P1, ACC96_P2};
  acc96_coeff_t reduced;
  acc96_coeff_from_u128(&reduced, scalar);
  for (size_t lane = 0; lane < ACC96_LIMBS; lane++)
    for (size_t i = 0; i < ACC96_N; i++)
      out->coeff[i].limb[lane] =
          (input->coeff[i].limb[lane] * reduced.limb[lane]) % modulus[lane];
}

static int build_covariances(acc96_sampleleft_ctx_t *sampler,
                             const acc96_trapdoor_t *trapdoor,
                             const acc96_ctx_t *ctx) {
  acc96_poly_t *rr00 = calloc(1, sizeof(*rr00));
  acc96_poly_t *rr11 = calloc(1, sizeof(*rr11));
  acc96_poly_t *rr01 = calloc(1, sizeof(*rr01));
  acc96_poly_t *conjugate = calloc(1, sizeof(*conjugate));
  arb_poly_t c00, c01, c10, inverse, temporary;
  arb_t scalar, s1_square, s2_square, s3_square, negative_sg_square;
  int ok = 0;
  if (rr00 == NULL || rr11 == NULL || rr01 == NULL || conjugate == NULL)
    goto allocation_cleanup;
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++) {
    q_conjugate(conjugate, &trapdoor->r[0][i]);
    q_muladd(rr00, &trapdoor->r[0][i], conjugate, ctx);
    q_conjugate(conjugate, &trapdoor->r[1][i]);
    q_muladd(rr11, &trapdoor->r[1][i], conjugate, ctx);
    q_muladd(rr01, &trapdoor->r[0][i], conjugate, ctx);
  }

  arb_poly_init(c00);
  arb_poly_init(c01);
  arb_poly_init(c10);
  arb_poly_init(inverse);
  arb_poly_init(temporary);
  arb_init(scalar);
  arb_init(s1_square);
  arb_init(s2_square);
  arb_init(s3_square);
  arb_init(negative_sg_square);
  if (arb_set_str(scalar, ACC_TRUNC_SGINVSQ_S4INVSQ_HP,
                  sampler->precision) != 0 ||
      arb_set_str(s1_square, ACC_TRUNC_S1SQ_HP, sampler->precision) != 0 ||
      arb_set_str(s2_square, ACC_TRUNC_S2SQ_HP, sampler->precision) != 0 ||
      arb_set_str(s3_square, ACC_TRUNC_S3SQ_HP, sampler->precision) != 0 ||
      arb_set_str(negative_sg_square, ACC_TRUNC_NEG_SGSQ_HP,
                  sampler->precision) != 0)
    goto arb_cleanup;

  arb_poly_from_q(c00, rr00, ctx);
  arb_poly_from_q(c01, rr01, ctx);
  arb_poly_from_q(sampler->covariance_second, rr11, ctx);
  arb_poly_scalar_mul(c00, c00, scalar, sampler->precision);
  arb_poly_scalar_mul(c01, c01, scalar, sampler->precision);
  arb_poly_scalar_mul(sampler->covariance_second,
                      sampler->covariance_second, scalar,
                      sampler->precision);
  arb_add_constant(c00, s1_square, sampler->precision);
  arb_add_constant(sampler->covariance_second, s3_square,
                   sampler->precision);

  if (!acc96_samplefz_ring_invert(inverse, sampler->covariance_second,
                                  sampler->samplefz))
    goto arb_cleanup;
  acc96_samplefz_ring_mul(sampler->schur_multiplier, c01, inverse,
                          sampler->samplefz);
  acc96_samplefz_ring_conjugate(c10, c01, sampler->samplefz);
  acc96_samplefz_ring_mul(temporary, sampler->schur_multiplier, c10,
                          sampler->samplefz);
  arb_poly_sub(sampler->covariance_left, c00, temporary,
               sampler->precision);
  arb_add_constant(sampler->covariance_left, negative_sg_square,
                   sampler->precision);
  arb_poly_zero(sampler->covariance_independent);
  arb_add(s2_square, s2_square, negative_sg_square, sampler->precision);
  arb_poly_set_coeff_arb(sampler->covariance_independent, 0, s2_square);
  ok = arb_is_positive(s2_square);

arb_cleanup:
  arb_clear(negative_sg_square);
  arb_clear(s3_square);
  arb_clear(s2_square);
  arb_clear(s1_square);
  arb_clear(scalar);
  arb_poly_clear(temporary);
  arb_poly_clear(inverse);
  arb_poly_clear(c10);
  arb_poly_clear(c01);
  arb_poly_clear(c00);
allocation_cleanup:
  free(conjugate);
  free(rr01);
  free(rr11);
  free(rr00);
  return ok;
}

acc96_sampleleft_ctx_t *acc96_sampleleft_ctx_new(
    const acc96_trapdoor_t *trapdoor, const acc96_ctx_t *ctx,
    long precision_bits, size_t thread_count) {
  acc96_sampleleft_ctx_t *sampler;
  acc96_perturbation_covariance_cert_t certificate;
  if (trapdoor == NULL || ctx == NULL || !ctx->initialized) return NULL;
  sampler = calloc(1, sizeof(*sampler));
  if (sampler == NULL) return NULL;
  sampler->precision = precision_bits >= 256 ? precision_bits : 256;
  sampler->trapdoor = trapdoor;
#ifdef _OPENMP
  if (thread_count == 0) thread_count = (size_t)omp_get_max_threads();
#else
  (void)thread_count;
  thread_count = 1;
#endif
  if (thread_count == 0) thread_count = 1;
  if (thread_count > 16) thread_count = 16;
  sampler->threads = thread_count;
  sampler->samplefz = acc96_samplefz_ctx_new(sampler->precision);
  sampler->klein = acc96_klein_exact_table_new();
  sampler->p2_stats = calloc(thread_count, sizeof(*sampler->p2_stats));
  sampler->klein_stats = calloc(thread_count, sizeof(*sampler->klein_stats));
  arb_poly_init(sampler->covariance_left);
  arb_poly_init(sampler->covariance_second);
  arb_poly_init(sampler->covariance_independent);
  arb_poly_init(sampler->schur_multiplier);
  if (sampler->samplefz == NULL || sampler->klein == NULL ||
      sampler->p2_stats == NULL || sampler->klein_stats == NULL ||
      !acc96_perturbation_covariance_certify_acb(
          trapdoor, ctx, sampler->precision, &certificate) ||
      !build_covariances(sampler, trapdoor, ctx) ||
      !acc96_samplefz_prepare_fixed(sampler->samplefz,
                                    sampler->covariance_independent, 0) ||
      !acc96_samplefz_prepare_fixed(sampler->samplefz,
                                    sampler->covariance_second, 1) ||
      !acc96_samplefz_prepare_fixed(sampler->samplefz,
                                    sampler->covariance_left, 2)) {
    acc96_sampleleft_ctx_free(sampler);
    return NULL;
  }
  return sampler;
}

void acc96_sampleleft_ctx_free(acc96_sampleleft_ctx_t *sampler) {
  if (sampler == NULL) return;
  if (sampler->p2_stats != NULL)
    for (size_t i = 0; i < sampler->threads; i++)
      acc96_gaussian_mpfr_stats_clear(&sampler->p2_stats[i]);
  if (sampler->klein_stats != NULL)
    for (size_t i = 0; i < sampler->threads; i++)
      acc96_gaussian_mpfr_stats_clear(&sampler->klein_stats[i]);
  free(sampler->klein_stats);
  free(sampler->p2_stats);
  acc96_klein_exact_table_free(sampler->klein);
  acc96_samplefz_ctx_free(sampler->samplefz);
  arb_poly_clear(sampler->schur_multiplier);
  arb_poly_clear(sampler->covariance_independent);
  arb_poly_clear(sampler->covariance_second);
  arb_poly_clear(sampler->covariance_left);
  free(sampler);
}

static int sample_right(acc96_evalfx_delta_t *right) {
  uint8_t buffer[4096];
  size_t offset = sizeof(buffer);
  for (size_t bit = 0; bit < ACC96_HANDLE_BITS; bit++) {
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++) {
      for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
        uint8_t sample;
        do {
          if (offset == sizeof(buffer)) {
            randombytes(buffer, sizeof(buffer));
            offset = 0;
          }
          sample = buffer[offset++] & 3U;
        } while (sample == 3U);
        set_signed(&right->right[bit][digit].coeff[coefficient],
                   (int64_t)sample - 1);
      }
    }
  }
  return 1;
}

static int sample_p2(acc96_poly_t p2[ACC96_TRUNC_RIGHT],
                     acc96_sampleleft_ctx_t *sampler) {
  int all_ok = 1;
#ifdef _OPENMP
#pragma omp parallel num_threads(sampler->threads) reduction(&:all_ok)
#endif
  {
    size_t worker = 0;
    mpfr_t center, variance;
#ifdef _OPENMP
    worker = (size_t)omp_get_thread_num();
#endif
    mpfr_init2(center, (mpfr_prec_t)sampler->precision);
    mpfr_init2(variance, (mpfr_prec_t)sampler->precision);
    mpfr_set_zero(center, 0);
    if (mpfr_set_str(variance, ACC_TRUNC_SQRT_S4SQ_SGSQ_HP, 10,
                     MPFR_RNDN) != 0)
      all_ok = 0;
    mpfr_sqr(variance, variance, MPFR_RNDN);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (size_t index = 0; index < ACC96_TRUNC_RIGHT * ACC96_N; index++) {
      int64_t value;
      if (!acc96_gaussian_mpfr_sample_center_variance(
              &value, center, variance, &sampler->p2_stats[worker])) {
        all_ok = 0;
      } else {
        set_signed(&p2[index / ACC96_N].coeff[index % ACC96_N], value);
      }
    }
    mpfr_clear(variance);
    mpfr_clear(center);
  }
  return all_ok;
}

static void center_from_q_scaled(arb_poly_t out, const acc96_poly_t *input,
                                 const arb_t scale,
                                 acc96_sampleleft_ctx_t *sampler,
                                 const acc96_ctx_t *ctx) {
  arb_poly_from_q(out, input, ctx);
  arb_poly_scalar_mul(out, out, scale, sampler->precision);
}

int acc96_sampleleft_secure(acc96_relation_witness_t *out,
                            const acc96_relation_pp_t *pp,
                            const acc96_trapdoor_t *trapdoor,
                            const acc96_poly_t *accumulator,
                            uint32_t handle, void *opaque,
                            const acc96_ctx_t *ctx) {
  acc96_sampleleft_ctx_t *sampler = opaque;
  acc96_poly_t *p_left = calloc(ACC96_TRUNC_L, sizeof(*p_left));
  acc96_poly_t *p12 = calloc(1, sizeof(*p12));
  acc96_poly_t *p2 = calloc(ACC96_TRUNC_RIGHT, sizeof(*p2));
  acc96_poly_t *y = calloc(ACC96_TRUNC_K, sizeof(*y));
  acc96_poly_t *w = calloc(1, sizeof(*w));
  acc96_poly_t *temporary = calloc(1, sizeof(*temporary));
  acc96_poly_t *right_contribution = calloc(1, sizeof(*right_contribution));
  arb_poly_t center0, center12, sampled12, difference, correction, zero;
  arb_t center_scale;
  acc96_u128 power = ACC96_TRUNC_BASE;
  int ok = 0;
  if (out == NULL || pp == NULL || trapdoor == NULL || accumulator == NULL ||
      sampler == NULL || ctx == NULL || !ctx->initialized ||
      sampler->trapdoor != trapdoor || p_left == NULL || p12 == NULL ||
      p2 == NULL || y == NULL || w == NULL || temporary == NULL ||
      right_contribution == NULL)
    goto raw_cleanup;
  memset(out, 0, sizeof(*out));
  if (!sample_right(&out->right)) goto raw_cleanup;
  acc96_member_matrix_apply(right_contribution, &pp->ba, handle,
                            &out->right, ctx);
  acc96_poly_sub(w, accumulator, right_contribution);

  arb_poly_init(center0);
  arb_poly_init(center12);
  arb_poly_init(sampled12);
  arb_poly_init(difference);
  arb_poly_init(correction);
  arb_poly_init(zero);
  arb_init(center_scale);
  arb_poly_zero(zero);
  if (arb_set_str(center_scale, ACC_TRUNC_NEGSGSQ_DIV_S4SQ_SGSQ_HP,
                  sampler->precision) != 0 ||
      !acc96_samplefz_fixed(&p_left[1], zero, sampler->samplefz, 0) ||
      !sample_p2(p2, sampler))
    goto arb_cleanup;

  acc96_poly_zero(&p_left[0]);
  acc96_poly_zero(p12);
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++) {
    q_muladd(&p_left[0], &trapdoor->r[0][i], &p2[i], ctx);
    q_muladd(p12, &trapdoor->r[1][i], &p2[i], ctx);
  }
  center_from_q_scaled(center0, &p_left[0], center_scale, sampler, ctx);
  center_from_q_scaled(center12, p12, center_scale, sampler, ctx);
  if (!acc96_samplefz_fixed(p12, center12, sampler->samplefz, 1))
    goto arb_cleanup;
  arb_poly_from_q(sampled12, p12, ctx);
  arb_poly_sub(difference, sampled12, center12, sampler->precision);
  acc96_samplefz_ring_mul(correction, difference,
                          sampler->schur_multiplier, sampler->samplefz);
  arb_poly_add(center0, center0, correction, sampler->precision);
  if (!acc96_samplefz_fixed(&p_left[0], center0, sampler->samplefz, 2))
    goto arb_cleanup;

  acc96_poly_sub(w, w, &p_left[0]);
  q_scale(temporary, &p_left[1], power);
  acc96_poly_sub(w, w, temporary);
  q_mulsub(w, &pp->atrap, p12, ctx);
  power *= ACC96_TRUNC_BASE;
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++) {
    q_muladd(w, &pp->trunc_b[i], &p2[i], ctx);
    q_scale(temporary, &p2[i], power);
    acc96_poly_sub(w, w, temporary);
    power *= ACC96_TRUNC_BASE;
  }

  {
    void *thread_opaque[16];
    for (size_t i = 0; i < sampler->threads; i++)
      thread_opaque[i] = &sampler->klein_stats[i];
    if (!acc96_klein_sample_exact_rational_parallel(
            y, w, sampler->klein,
            acc96_gaussian_mpfr_sample_rational, thread_opaque,
            sampler->threads, ctx))
      goto arb_cleanup;
  }

  acc96_poly_set(&out->left, &p_left[0]);
  q_scale(temporary, &p_left[1], ACC96_TRUNC_BASE);
  acc96_poly_add(&out->left, &out->left, temporary);
  acc96_poly_add(&out->left, &out->left, &y[0]);
  q_scale(temporary, &y[1], ACC96_TRUNC_BASE);
  acc96_poly_add(&out->left, &out->left, temporary);
  acc96_poly_set(&out->v12, p12);
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++) {
    q_muladd(&out->left, &trapdoor->r[0][i], &y[ACC96_TRUNC_L + i], ctx);
    q_muladd(&out->v12, &trapdoor->r[1][i],
             &y[ACC96_TRUNC_L + i], ctx);
    acc96_poly_add(&out->v2[i], &p2[i], &y[ACC96_TRUNC_L + i]);
  }
  acc96_poly_zero(&out->tag);
  acc96_coeff_from_u128(&out->tag.coeff[0], 1);
  ok = 1;

arb_cleanup:
  arb_clear(center_scale);
  arb_poly_clear(zero);
  arb_poly_clear(correction);
  arb_poly_clear(difference);
  arb_poly_clear(sampled12);
  arb_poly_clear(center12);
  arb_poly_clear(center0);
raw_cleanup:
  free(right_contribution);
  free(temporary);
  free(w);
  free(y);
  free(p2);
  free(p12);
  free(p_left);
  if (!ok && out != NULL) memset(out, 0, sizeof(*out));
  return ok;
}

size_t acc96_sampleleft_thread_count(const acc96_sampleleft_ctx_t *sampler) {
  return sampler == NULL ? 0 : sampler->threads;
}

double acc96_sampleleft_maximum_arb_radius(
    const acc96_sampleleft_ctx_t *sampler) {
  return sampler == NULL ? 0.0 :
      acc96_samplefz_maximum_input_radius(sampler->samplefz);
}
