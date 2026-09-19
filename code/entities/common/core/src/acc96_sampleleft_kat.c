#include "acc96_sampleleft_kat.h"

#include <stdlib.h>
#include <string.h>

static int decompose_base13(int32_t digits[ACC96_TRUNC_K],
                            const acc96_coeff_t *value,
                            const acc96_ctx_t *ctx) {
  acc96_i128 remaining = acc96_coeff_to_centered(value, ctx);
  for (size_t i = 0; i < ACC96_TRUNC_K; i++) {
    acc96_i128 digit = remaining % (acc96_i128)ACC96_TRUNC_BASE;
    if (digit < 0) digit += ACC96_TRUNC_BASE;
    if (digit > (acc96_i128)(ACC96_TRUNC_BASE / 2))
      digit -= ACC96_TRUNC_BASE;
    digits[i] = (int32_t)digit;
    remaining = (remaining - digit) / (acc96_i128)ACC96_TRUNC_BASE;
  }
  return remaining == 0;
}

static void set_signed(acc96_coeff_t *out, int32_t value) {
  if (value < 0)
    acc96_coeff_from_u128(out, acc96_modulus() - (acc96_u128)(-value));
  else
    acc96_coeff_from_u128(out, (acc96_u128)value);
}

static void add_product(acc96_poly_t *acc, const acc96_poly_t *lhs,
                        const acc96_poly_t *rhs,
                        const acc96_ctx_t *ctx) {
  acc96_poly_t *product = malloc(sizeof(*product));
  acc96_poly_t *next = malloc(sizeof(*next));
  if (product == NULL || next == NULL) abort();
  acc96_poly_mul(product, lhs, rhs, ctx);
  acc96_poly_add(next, acc, product);
  acc96_poly_set(acc, next);
  free(next);
  free(product);
}

int acc96_sampleleft_algebraic_kat(
    acc96_relation_witness_t *out, const acc96_relation_pp_t *pp,
    const acc96_trapdoor_t *trapdoor, const acc96_poly_t *accumulator,
    uint32_t handle, void *opaque, const acc96_ctx_t *ctx) {
  acc96_poly_t *y = NULL;
  acc96_poly_t *lower = NULL;
  acc96_poly_t *scalar = NULL;
  (void)pp;
  (void)handle;
  (void)opaque;
  if (out == NULL || trapdoor == NULL || accumulator == NULL || ctx == NULL)
    return 0;
  y = calloc(ACC96_TRUNC_K, sizeof(*y));
  lower = calloc(1, sizeof(*lower));
  scalar = calloc(1, sizeof(*scalar));
  if (y == NULL || lower == NULL || scalar == NULL) goto fail;
  memset(out, 0, sizeof(*out));
  for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
    int32_t digits[ACC96_TRUNC_K];
    if (!decompose_base13(digits, &accumulator->coeff[coefficient], ctx))
      goto fail;
    for (size_t digit = 0; digit < ACC96_TRUNC_K; digit++)
      set_signed(&y[digit].coeff[coefficient], digits[digit]);
  }

  acc96_poly_set(&out->left, &y[0]);
  acc96_poly_zero(scalar);
  acc96_coeff_from_u128(&scalar->coeff[0], ACC96_TRUNC_BASE);
  acc96_poly_mul(lower, &y[1], scalar, ctx);
  acc96_poly_add(&out->left, &out->left, lower);
  acc96_poly_zero(&out->v12);
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++) {
    const acc96_poly_t *tail = &y[ACC96_TRUNC_L + i];
    acc96_poly_set(&out->v2[i], tail);
    add_product(&out->left, &trapdoor->r[0][i], tail, ctx);
    add_product(&out->v12, &trapdoor->r[1][i], tail, ctx);
  }
  acc96_poly_zero(&out->tag);
  acc96_coeff_from_u128(&out->tag.coeff[0], 1);
  free(scalar);
  free(lower);
  free(y);
  return 1;

fail:
  free(scalar);
  free(lower);
  free(y);
  if (out != NULL) memset(out, 0, sizeof(*out));
  return 0;
}
