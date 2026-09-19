#include "acc96_klein_exact.h"

#include <gmp.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

struct acc96_klein_exact_table {
  mpz_t numerator[ACC_TRUNC_PARAM_K][ACC_TRUNC_PARAM_K];
  mpz_t denominator[ACC_TRUNC_PARAM_K];
  size_t initialized_numerators;
  size_t initialized_denominators;
};

/* The rational center split is on the hot path of the exact Klein sampler.
 * Keeping these GMP temporaries per worker avoids initializing and clearing
 * five big integers for every coordinate/column step.  The objects are never
 * shared between workers, so this is only an allocation optimization; the
 * arithmetic and the resulting distribution are unchanged. */
typedef struct {
  mpz_t center_numerator;
  mpz_t value;
  mpz_t floor_value;
  mpz_t floor_mod;
  mpz_t modulus;
  mpz_t residual_numerator;
  mpz_t residual_denominator;
} klein_rational_scratch_t;

static void klein_rational_scratch_init(klein_rational_scratch_t *scratch) {
  mpz_inits(scratch->center_numerator, scratch->value,
            scratch->floor_value, scratch->floor_mod, scratch->modulus,
            scratch->residual_numerator, scratch->residual_denominator,
            NULL);
}

static void klein_rational_scratch_clear(klein_rational_scratch_t *scratch) {
  mpz_clears(scratch->residual_denominator,
             scratch->residual_numerator, scratch->modulus,
             scratch->floor_mod, scratch->floor_value, scratch->value,
             scratch->center_numerator, NULL);
}

static void mpz_set_u128(mpz_t out, acc96_u128 value) {
  const uint64_t low = (uint64_t)value;
  const uint64_t high = (uint64_t)(value >> 64);
  mpz_set_ui(out, high);
  mpz_mul_2exp(out, out, 64);
  mpz_add_ui(out, out, low);
}

static void mpz_set_i128(mpz_t out, acc96_i128 value) {
  if (value < 0) {
    /* Avoid negating the minimum signed value in its signed type. */
    const acc96_u128 magnitude = (acc96_u128)(-(value + 1)) + 1;
    mpz_set_u128(out, magnitude);
    mpz_neg(out, out);
  } else {
    mpz_set_u128(out, (acc96_u128)value);
  }
}

static acc96_u128 mpz_get_u128(const mpz_t value) {
  uint64_t words[2] = {0, 0};
  size_t count = 0;
  mpz_export(words, &count, -1, sizeof(words[0]), 0, 0, value);
  return (acc96_u128)words[0] |
         (count > 1 ? (acc96_u128)words[1] << 64 : 0);
}

acc96_klein_exact_table_t *acc96_klein_exact_table_new(void) {
  acc96_klein_exact_table_t *table = calloc(1, sizeof(*table));
  if (table == NULL) return NULL;
  for (size_t column = 0; column < ACC_TRUNC_PARAM_K; column++) {
    mpz_init(table->denominator[column]);
    table->initialized_denominators++;
    if (mpz_set_str(table->denominator[column],
                    acc_trunc_neg_scaled_gadget_gso_column_den[column], 10) !=
            0 || mpz_sgn(table->denominator[column]) <= 0) {
      acc96_klein_exact_table_free(table);
      return NULL;
    }
  }
  for (size_t row = 0; row < ACC_TRUNC_PARAM_K; row++) {
    for (size_t column = 0; column < ACC_TRUNC_PARAM_K; column++) {
      mpz_init(table->numerator[row][column]);
      table->initialized_numerators++;
      if (mpz_set_str(table->numerator[row][column],
                      acc_trunc_neg_scaled_gadget_gso_column_num[row][column],
                      10) != 0) {
        acc96_klein_exact_table_free(table);
        return NULL;
      }
    }
  }
  return table;
}

void acc96_klein_exact_table_free(acc96_klein_exact_table_t *table) {
  if (table == NULL) return;
  for (size_t i = 0; i < table->initialized_numerators; i++)
    mpz_clear(table->numerator[i / ACC_TRUNC_PARAM_K]
                              [i % ACC_TRUNC_PARAM_K]);
  for (size_t i = 0; i < table->initialized_denominators; i++)
    mpz_clear(table->denominator[i]);
  free(table);
}

int acc96_klein_center_split(
    const acc96_klein_exact_table_t *table, size_t column, acc96_i128 w,
    const acc96_i128 v[ACC_TRUNC_PARAM_K], acc96_coeff_t *integer_mod_q,
    double *fraction) {
  mpz_t residual_numerator, residual_denominator;
  int ok;
  if (fraction == NULL) return 0;
  mpz_inits(residual_numerator, residual_denominator, NULL);
  ok = acc96_klein_center_split_rational(
      table, column, w, v, integer_mod_q, residual_numerator,
      residual_denominator);
  if (ok) {
    *fraction = mpz_get_d(residual_numerator) /
                mpz_get_d(residual_denominator);
    ok = *fraction >= 0.0 && *fraction < 1.0;
  }
  mpz_clears(residual_denominator, residual_numerator, NULL);
  return ok;
}

static int klein_center_split_rational_scratch(
    const acc96_klein_exact_table_t *table, size_t column, acc96_i128 w,
    const acc96_i128 v[ACC_TRUNC_PARAM_K], acc96_coeff_t *integer_mod_q,
    mpz_t residual_numerator, mpz_t residual_denominator,
    klein_rational_scratch_t *scratch) {
  if (table == NULL || column >= ACC_TRUNC_PARAM_K || v == NULL ||
      integer_mod_q == NULL || scratch == NULL)
    return 0;
  mpz_t *center_numerator = &scratch->center_numerator;
  mpz_t *value = &scratch->value;
  mpz_t *floor_value = &scratch->floor_value;
  mpz_t *floor_mod = &scratch->floor_mod;
  mpz_t *modulus = &scratch->modulus;
  mpz_set_ui(*center_numerator, 0);
  for (size_t row = 0; row < ACC_TRUNC_PARAM_K; row++) {
    const acc96_i128 input = row == 0 ? w + v[0] : v[row];
    mpz_set_i128(*value, input);
    mpz_addmul(*center_numerator, *value, table->numerator[row][column]);
  }
  mpz_fdiv_qr(*floor_value, residual_numerator, *center_numerator,
              table->denominator[column]);
  mpz_set(residual_denominator, table->denominator[column]);

  mpz_set_u128(*modulus, acc96_modulus());
  mpz_fdiv_r(*floor_mod, *floor_value, *modulus);
  acc96_coeff_from_u128(integer_mod_q, mpz_get_u128(*floor_mod));
  return mpz_sgn(residual_numerator) >= 0 &&
         mpz_cmp(residual_numerator, residual_denominator) < 0;
}

int acc96_klein_center_split_rational(
    const acc96_klein_exact_table_t *table, size_t column, acc96_i128 w,
    const acc96_i128 v[ACC_TRUNC_PARAM_K], acc96_coeff_t *integer_mod_q,
    mpz_t residual_numerator, mpz_t residual_denominator) {
  klein_rational_scratch_t scratch;
  int ok;
  if (integer_mod_q == NULL || residual_numerator == NULL ||
      residual_denominator == NULL)
    return 0;
  klein_rational_scratch_init(&scratch);
  ok = klein_center_split_rational_scratch(
      table, column, w, v, integer_mod_q, residual_numerator,
      residual_denominator, &scratch);
  klein_rational_scratch_clear(&scratch);
  return ok;
}

void acc96_klein_center_recombine(acc96_coeff_t *out,
                                  const acc96_coeff_t *integer_mod_q,
                                  int64_t centered_sample,
                                  const acc96_ctx_t *ctx) {
  const acc96_u128 modulus = acc96_modulus();
  acc96_u128 value = acc96_coeff_to_u128(integer_mod_q, ctx);
  if (centered_sample < 0) {
    const uint64_t magnitude = (uint64_t)(-(centered_sample + 1)) + 1;
    /* magnitude <= 2^63 while the q96 modulus is above 2^95. */
    const acc96_u128 reduced = magnitude;
    value = value >= reduced ? value - reduced : modulus - (reduced - value);
  } else {
    value += (uint64_t)centered_sample;
    if (value >= modulus) value -= modulus;
  }
  acc96_coeff_from_u128(out, value);
}

static void coefficient_add_signed_product(acc96_coeff_t *out,
                                           const acc96_coeff_t *sample,
                                           int64_t scalar) {
  static const uint64_t moduli[ACC96_LIMBS] = {ACC96_P0, ACC96_P1, ACC96_P2};
  const uint64_t magnitude =
      scalar < 0 ? (uint64_t)(-(scalar + 1)) + 1 : (uint64_t)scalar;
  if (magnitude == 0) return;
  for (size_t lane = 0; lane < ACC96_LIMBS; lane++) {
    /* The generated truncated-gadget basis has |scalar| <= 13. */
    const uint64_t product = (sample->limb[lane] * magnitude) % moduli[lane];
    if (scalar < 0)
      out->limb[lane] = out->limb[lane] >= product
                            ? out->limb[lane] - product
                            : moduli[lane] - (product - out->limb[lane]);
    else {
      out->limb[lane] += product;
      if (out->limb[lane] >= moduli[lane]) out->limb[lane] -= moduli[lane];
    }
  }
}

int acc96_klein_sample_exact(
    acc96_poly_t y[ACC_TRUNC_PARAM_K], const acc96_poly_t *w,
    const acc96_klein_exact_table_t *table,
    acc96_residual_gaussian_fn residual_sampler, void *opaque,
    const acc96_ctx_t *ctx) {
  acc96_i128 values[ACC_TRUNC_PARAM_K];
  acc96_coeff_t integer_part, sample;
  if (y == NULL || w == NULL || table == NULL || residual_sampler == NULL ||
      ctx == NULL || !ctx->initialized)
    return 0;
  for (size_t row = 0; row < ACC_TRUNC_PARAM_K; row++)
    acc96_poly_zero(&y[row]);
  for (size_t reverse = ACC96_N * ACC_TRUNC_PARAM_K; reverse > 0; reverse--) {
    const size_t column = (reverse - 1) / ACC96_N;
    const size_t coefficient = (reverse - 1) % ACC96_N;
    for (size_t row = 0; row < ACC_TRUNC_PARAM_K; row++)
      values[row] = acc96_coeff_to_centered(&y[row].coeff[coefficient], ctx);
    const acc96_i128 syndrome =
        acc96_coeff_to_centered(&w->coeff[coefficient], ctx);
    double fractional_center;
    if (!acc96_klein_center_split(table, column, syndrome, values,
                                  &integer_part, &fractional_center))
      return 0;
    const int64_t residual = residual_sampler(
        fractional_center, acc_trunc_klein_widths[column], opaque);
    acc96_klein_center_recombine(&sample, &integer_part, residual, ctx);
    for (size_t row = 0; row < ACC_TRUNC_PARAM_K; row++)
      coefficient_add_signed_product(
          &y[row].coeff[coefficient], &sample,
          acc_trunc_gadget_basis[row][column]);
  }
  acc96_poly_add(&y[0], &y[0], w);
  return 1;
}

int acc96_klein_sample_exact_rational(
    acc96_poly_t y[ACC_TRUNC_PARAM_K], const acc96_poly_t *w,
    const acc96_klein_exact_table_t *table,
    acc96_residual_gaussian_rational_fn residual_sampler, void *opaque,
    const acc96_ctx_t *ctx) {
  acc96_i128 values[ACC_TRUNC_PARAM_K];
  acc96_coeff_t integer_part, sample;
  klein_rational_scratch_t scratch;
  int ok = 0;

  if (y == NULL || w == NULL || table == NULL || residual_sampler == NULL ||
      ctx == NULL || !ctx->initialized)
    return 0;
  klein_rational_scratch_init(&scratch);
  for (size_t row = 0; row < ACC_TRUNC_PARAM_K; row++)
    acc96_poly_zero(&y[row]);
  for (size_t reverse = ACC96_N * ACC_TRUNC_PARAM_K; reverse > 0; reverse--) {
    const size_t column = (reverse - 1) / ACC96_N;
    const size_t coefficient = (reverse - 1) % ACC96_N;
    int64_t residual;

    for (size_t row = 0; row < ACC_TRUNC_PARAM_K; row++)
      values[row] = acc96_coeff_to_centered(&y[row].coeff[coefficient], ctx);
    if (!klein_center_split_rational_scratch(
            table, column,
            acc96_coeff_to_centered(&w->coeff[coefficient], ctx), values,
            &integer_part, scratch.residual_numerator,
            scratch.residual_denominator, &scratch) ||
        !residual_sampler(
            &residual, scratch.residual_numerator,
            scratch.residual_denominator,
            acc_trunc_klein_widths_high_precision[column], opaque)) {
      goto cleanup;
    }
    acc96_klein_center_recombine(&sample, &integer_part, residual, ctx);
    for (size_t row = 0; row < ACC_TRUNC_PARAM_K; row++)
      coefficient_add_signed_product(
          &y[row].coeff[coefficient], &sample,
          acc_trunc_gadget_basis[row][column]);
  }
  acc96_poly_add(&y[0], &y[0], w);
  ok = 1;
cleanup:
  klein_rational_scratch_clear(&scratch);
  return ok;
}

int acc96_klein_sample_exact_rational_parallel(
    acc96_poly_t y[ACC_TRUNC_PARAM_K], const acc96_poly_t *w,
    const acc96_klein_exact_table_t *table,
    acc96_residual_gaussian_rational_fn residual_sampler,
    void *const thread_opaque[], size_t thread_count,
    const acc96_ctx_t *ctx) {
  int all_ok = 1;
  klein_rational_scratch_t *scratch = NULL;

  if (y == NULL || w == NULL || table == NULL || residual_sampler == NULL ||
      thread_opaque == NULL || thread_count == 0 || ctx == NULL ||
      !ctx->initialized)
    return 0;
  scratch = calloc(thread_count, sizeof(*scratch));
  if (scratch == NULL) return 0;
  for (size_t i = 0; i < thread_count; i++)
    klein_rational_scratch_init(&scratch[i]);
  for (size_t row = 0; row < ACC_TRUNC_PARAM_K; row++)
    acc96_poly_zero(&y[row]);

#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(static) reduction(&:all_ok)
#endif
  for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
    acc96_i128 values[ACC_TRUNC_PARAM_K] = {0};
    acc96_coeff_t integer_part, sample;
    int chain_ok = 1;
    size_t worker = 0;
#ifdef _OPENMP
    worker = (size_t)omp_get_thread_num();
#endif
    klein_rational_scratch_t *local = &scratch[worker];
    for (size_t reverse = ACC_TRUNC_PARAM_K; reverse > 0; reverse--) {
      const size_t column = reverse - 1;
      int64_t residual;
      if (!klein_center_split_rational_scratch(
              table, column,
              acc96_coeff_to_centered(&w->coeff[coefficient], ctx), values,
              &integer_part, local->residual_numerator,
              local->residual_denominator, local) ||
          !residual_sampler(
              &residual, local->residual_numerator,
              local->residual_denominator,
              acc_trunc_klein_widths_high_precision[column],
              thread_opaque[worker])) {
        chain_ok = 0;
        break;
      }
      acc96_klein_center_recombine(&sample, &integer_part, residual, ctx);
      for (size_t row = 0; row < ACC_TRUNC_PARAM_K; row++) {
        const int64_t scalar = acc_trunc_gadget_basis[row][column];
        if (scalar != 0) {
          coefficient_add_signed_product(
              &y[row].coeff[coefficient], &sample, scalar);
          values[row] =
              acc96_coeff_to_centered(&y[row].coeff[coefficient], ctx);
        }
      }
    }
    all_ok &= chain_ok;
  }
  for (size_t i = 0; i < thread_count; i++)
    klein_rational_scratch_clear(&scratch[i]);
  free(scratch);
  if (!all_ok) return 0;
  acc96_poly_add(&y[0], &y[0], w);
  return 1;
}
