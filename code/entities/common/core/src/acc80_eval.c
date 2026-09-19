#include "acc80_eval.h"

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
  acc80_poly_t block[ACC80_GADGET_DIGITS][ACC80_GADGET_DIGITS];
} acc80_ginv_t;

static int handle_bit(uint32_t value, size_t bit) {
  return (int)((value >> bit) & 1U);
}

static acc80_u128 gadget_power(size_t digit) {
  acc80_u128 value = 1;
  while (digit-- != 0) value *= ACC80_GADGET_BASE;
  return value;
}

void acc80_gadget_identity(acc80_gmat_t *out) {
  memset(out, 0, sizeof(*out));
  for (size_t digit = 0; digit < ACC80_GADGET_DIGITS; digit++) {
    acc80_coeff_from_u128(&out->block[digit].coeff[0],
                          gadget_power(digit));
  }
}

static void di(acc80_gmat_t *out, const acc80_ba_t *ba, size_t bit,
               int selected) {
  if (selected) {
    memcpy(out, &ba->bit[bit], sizeof(*out));
    return;
  }
  acc80_gmat_t identity;
  acc80_gadget_identity(&identity);
  for (size_t digit = 0; digit < ACC80_GADGET_DIGITS; digit++) {
    acc80_poly_sub(&out->block[digit], &identity.block[digit],
                   &ba->bit[bit].block[digit]);
  }
}

static void ginv(acc80_ginv_t *out, const acc80_gmat_t *in,
                 const acc80_ctx_t *ctx) {
  /* Matrix path: every T_i coefficient uses the congruence-class minimum. */
  memset(out, 0, sizeof(*out));
  for (size_t input = 0; input < ACC80_GADGET_DIGITS; input++) {
    for (size_t coeff = 0; coeff < ACC80_N; coeff++) {
      int32_t digits[ACC80_GADGET_DIGITS];
      if (!acc80_balanced_decompose(digits, &in->block[input].coeff[coeff],
                                    ctx)) {
        abort();
      }
      for (size_t output = 0; output < ACC80_GADGET_DIGITS; output++) {
        acc80_i128 value = digits[output];
        if (value < 0) value += (acc80_i128)acc80_modulus();
        acc80_coeff_from_u128(&out->block[output][input].coeff[coeff],
                              (acc80_u128)value);
      }
    }
  }
}

int acc80_ginv_reconstructs(const acc80_gmat_t *value,
                            const acc80_ctx_t *ctx) {
  acc80_ginv_t *inverse = malloc(sizeof(*inverse));
  if (inverse == NULL) return 0;
  ginv(inverse, value, ctx);
  for (size_t input = 0; input < ACC80_GADGET_DIGITS; input++) {
    for (size_t coeff = 0; coeff < ACC80_N; coeff++) {
      acc80_i128 reconstructed = 0;
      acc80_i128 power = 1;
      for (size_t digit = 0; digit < ACC80_GADGET_DIGITS; digit++) {
        reconstructed += acc80_coeff_to_centered(
            &inverse->block[digit][input].coeff[coeff], ctx) * power;
        power *= ACC80_GADGET_BASE;
      }
      const acc80_i128 modulus = (acc80_i128)acc80_modulus();
      reconstructed %= modulus;
      if (reconstructed < 0) reconstructed += modulus;
      if ((acc80_u128)reconstructed != acc80_coeff_to_u128(
                                               &value->block[input].coeff[coeff],
                                               ctx)) {
        free(inverse);
        return 0;
      }
    }
  }
  free(inverse);
  return 1;
}

static void mul_column(acc80_gmat_t *out, const acc80_gmat_t *lhs,
                       const acc80_ginv_t *rhs, size_t column,
                       const acc80_ctx_t *ctx) {
  acc80_poly_t *product = malloc(sizeof(*product));
  acc80_poly_t *sum = malloc(sizeof(*sum));
  acc80_poly_t *next = malloc(sizeof(*next));
  if (product == NULL || sum == NULL || next == NULL) abort();
  acc80_poly_zero(sum);
  for (size_t inner = 0; inner < ACC80_GADGET_DIGITS; inner++) {
    acc80_poly_mul(product, &lhs->block[inner], &rhs->block[inner][column],
                   ctx);
    acc80_poly_add(next, sum, product);
    acc80_poly_set(sum, next);
  }
  acc80_poly_set(&out->block[column], sum);
  free(next);
  free(sum);
  free(product);
}

static void compose(acc80_gmat_t *out, const acc80_gmat_t *lhs,
                    const acc80_ginv_t *rhs, const acc80_ctx_t *ctx) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(ACC80_GADGET_DIGITS)
#endif
  for (int column = 0; column < (int)ACC80_GADGET_DIGITS; column++) {
    mul_column(out, lhs, rhs, (size_t)column, ctx);
  }
}

static void ginv_vec(acc80_poly_t out[ACC80_GADGET_DIGITS],
                     const acc80_poly_t *in, const acc80_ctx_t *ctx) {
  /* Vector path: u_a uses the same minimum-norm G^-1 independently. */
  memset(out, 0, ACC80_GADGET_DIGITS * sizeof(*out));
  for (size_t coeff = 0; coeff < ACC80_N; coeff++) {
    int32_t digits[ACC80_GADGET_DIGITS];
    if (!acc80_balanced_decompose(digits, &in->coeff[coeff], ctx)) abort();
    for (size_t digit = 0; digit < ACC80_GADGET_DIGITS; digit++) {
      acc80_i128 value = digits[digit];
      if (value < 0) value += (acc80_i128)acc80_modulus();
      acc80_coeff_from_u128(&out[digit].coeff[coeff], (acc80_u128)value);
    }
  }
}

static void ginv_apply(acc80_poly_t out[ACC80_GADGET_DIGITS],
                       const acc80_ginv_t *matrix,
                       const acc80_poly_t in[ACC80_GADGET_DIGITS],
                       const acc80_ctx_t *ctx) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(ACC80_GADGET_DIGITS)
#endif
  for (int row = 0; row < (int)ACC80_GADGET_DIGITS; row++) {
    acc80_poly_t *product = malloc(sizeof(*product));
    acc80_poly_t *next = malloc(sizeof(*next));
    if (product == NULL || next == NULL) abort();
    acc80_poly_zero(&out[row]);
    for (size_t column = 0; column < ACC80_GADGET_DIGITS; column++) {
      acc80_poly_mul(product, &matrix->block[row][column], &in[column], ctx);
      acc80_poly_add(next, &out[row], product);
      acc80_poly_set(&out[row], next);
    }
    free(next);
    free(product);
  }
}

void acc80_evalf_apply(acc80_poly_t *out, const acc80_gmat_t *matrix,
                       const acc80_poly_t *value,
                       const acc80_ctx_t *ctx) {
  acc80_poly_t *digits = calloc(ACC80_GADGET_DIGITS, sizeof(*digits));
  acc80_poly_t *product = malloc(sizeof(*product));
  acc80_poly_t *next = malloc(sizeof(*next));
  if (digits == NULL || product == NULL || next == NULL) abort();
  ginv_vec(digits, value, ctx);
  acc80_poly_zero(out);
  for (size_t digit = 0; digit < ACC80_GADGET_DIGITS; digit++) {
    acc80_poly_mul(product, &matrix->block[digit], &digits[digit], ctx);
    acc80_poly_add(next, out, product);
    acc80_poly_set(out, next);
  }
  free(next);
  free(product);
  free(digits);
}

void acc80_evalf(acc80_gmat_t *out, const acc80_ba_t *ba,
                 uint32_t selector, const acc80_ctx_t *ctx) {
  acc80_gmat_t *current = malloc(sizeof(*current));
  acc80_gmat_t *choice = malloc(sizeof(*choice));
  acc80_gmat_t *next = malloc(sizeof(*next));
  acc80_ginv_t *inverse = malloc(sizeof(*inverse));
  if (current == NULL || choice == NULL || next == NULL || inverse == NULL) {
    abort();
  }
  di(current, ba, ACC80_HANDLE_BITS - 1,
     handle_bit(selector, ACC80_HANDLE_BITS - 1));
  for (size_t offset = 1; offset < ACC80_HANDLE_BITS; offset++) {
    const size_t bit = ACC80_HANDLE_BITS - 1 - offset;
    ginv(inverse, current, ctx);
    di(choice, ba, bit, handle_bit(selector, bit));
    compose(next, choice, inverse, ctx);
    acc80_gmat_t *swap = current;
    current = next;
    next = swap;
  }
  memcpy(out, current, sizeof(*out));
  free(inverse);
  free(next);
  free(choice);
  free(current);
}

int acc80_evalf_with_tail(acc80_gmat_t *out, acc80_gmat_t *tail_checkpoint,
                          size_t tail_bit, const acc80_ba_t *ba,
                          uint32_t selector, const acc80_ctx_t *ctx) {
  acc80_gmat_t *current = malloc(sizeof(*current));
  acc80_gmat_t *choice = malloc(sizeof(*choice));
  acc80_gmat_t *next = malloc(sizeof(*next));
  acc80_ginv_t *inverse = malloc(sizeof(*inverse));
  int captured = 0;
  if (out == NULL || tail_checkpoint == NULL || ba == NULL || ctx == NULL ||
      tail_bit == 0 || tail_bit >= ACC80_HANDLE_BITS || current == NULL ||
      choice == NULL || next == NULL || inverse == NULL)
    goto cleanup;
  di(current, ba, ACC80_HANDLE_BITS - 1,
     handle_bit(selector, ACC80_HANDLE_BITS - 1));
  if (tail_bit == ACC80_HANDLE_BITS - 1) {
    memcpy(tail_checkpoint, current, sizeof(*current));
    captured = 1;
  }
  for (size_t bit = ACC80_HANDLE_BITS - 1; bit > 0;) {
    bit--;
    ginv(inverse, current, ctx);
    di(choice, ba, bit, handle_bit(selector, bit));
    compose(next, choice, inverse, ctx);
    {
      acc80_gmat_t *swap = current;
      current = next;
      next = swap;
    }
    if (bit == tail_bit) {
      memcpy(tail_checkpoint, current, sizeof(*current));
      captured = 1;
    }
  }
  if (captured) memcpy(out, current, sizeof(*out));
cleanup:
  free(inverse);
  free(next);
  free(choice);
  free(current);
  return captured;
}

static int evalfx_prefix_impl(
    acc80_gmat_t *indicator,
    acc80_poly_t out[][ACC80_GADGET_DIGITS], size_t layers,
    const acc80_ba_t *ba, uint32_t selector, const acc80_poly_t *value,
    const acc80_ctx_t *ctx) {
  acc80_gmat_t *tails = calloc(ACC80_HANDLE_BITS, sizeof(*tails));
  acc80_gmat_t *choice = malloc(sizeof(*choice));
  acc80_gmat_t *next = malloc(sizeof(*next));
  acc80_ginv_t *inverse = malloc(sizeof(*inverse));
  acc80_poly_t *input = calloc(ACC80_GADGET_DIGITS, sizeof(*input));
  acc80_poly_t *applied = calloc(ACC80_GADGET_DIGITS, sizeof(*applied));
  if (out == NULL || layers == 0 || layers > ACC80_HANDLE_BITS || ba == NULL ||
      value == NULL || ctx == NULL) {
    free(applied);
    free(input);
    free(inverse);
    free(next);
    free(choice);
    free(tails);
    return 0;
  }
  if (tails == NULL || choice == NULL || next == NULL || inverse == NULL ||
      input == NULL || applied == NULL) abort();
  memset(out, 0, layers * ACC80_GADGET_DIGITS * sizeof(**out));

  di(&tails[ACC80_HANDLE_BITS - 1], ba, ACC80_HANDLE_BITS - 1,
     handle_bit(selector, ACC80_HANDLE_BITS - 1));
  for (size_t offset = 1; offset < ACC80_HANDLE_BITS - 1; offset++) {
    const size_t bit = ACC80_HANDLE_BITS - 1 - offset;
    ginv(inverse, &tails[bit + 1], ctx);
    di(choice, ba, bit, handle_bit(selector, bit));
    compose(next, choice, inverse, ctx);
    memcpy(&tails[bit], next, sizeof(*next));
  }
  if (indicator != NULL) {
    ginv(inverse, &tails[1], ctx);
    di(choice, ba, 0, handle_bit(selector, 0));
    compose(indicator, choice, inverse, ctx);
  }
  ginv_vec(input, value, ctx);
  for (size_t bit = 0; bit < layers; bit++) {
    if (bit + 1 < ACC80_HANDLE_BITS) {
      ginv(inverse, &tails[bit + 1], ctx);
      ginv_apply(applied, inverse, input, ctx);
    } else {
      memcpy(applied, input, ACC80_GADGET_DIGITS * sizeof(*input));
    }
    for (size_t digit = 0; digit < ACC80_GADGET_DIGITS; digit++) {
      if (handle_bit(selector, bit)) {
        acc80_poly_set(&out[bit][digit], &applied[digit]);
      } else {
        acc80_poly_neg(&out[bit][digit], &applied[digit]);
      }
    }
  }
  free(applied);
  free(input);
  free(inverse);
  free(next);
  free(choice);
  free(tails);
  return 1;
}

int acc80_evalfx_prefix(
    acc80_poly_t out[][ACC80_GADGET_DIGITS], size_t layers,
    const acc80_ba_t *ba, uint32_t selector, const acc80_poly_t *value,
    const acc80_ctx_t *ctx) {
  return evalfx_prefix_impl(NULL, out, layers, ba, selector, value, ctx);
}

int acc80_evalf_evalfx_prefix(
    acc80_gmat_t *indicator,
    acc80_poly_t out[][ACC80_GADGET_DIGITS], size_t layers,
    const acc80_ba_t *ba, uint32_t selector, const acc80_poly_t *value,
    const acc80_ctx_t *ctx) {
  if (indicator == NULL) return 0;
  return evalfx_prefix_impl(indicator, out, layers, ba, selector, value, ctx);
}

int acc80_evalf_tail_checkpoint(acc80_gmat_t *out,
                                const acc80_ba_t *ba,
                                uint32_t selector, size_t tail_bit,
                                const acc80_ctx_t *ctx) {
  acc80_gmat_t *current = malloc(sizeof(*current));
  acc80_gmat_t *choice = malloc(sizeof(*choice));
  acc80_gmat_t *next = malloc(sizeof(*next));
  acc80_ginv_t *inverse = malloc(sizeof(*inverse));
  int ok = 0;
  if (out == NULL || ba == NULL || ctx == NULL || tail_bit == 0 ||
      tail_bit >= ACC80_HANDLE_BITS || current == NULL || choice == NULL ||
      next == NULL || inverse == NULL)
    goto cleanup;
  di(current, ba, ACC80_HANDLE_BITS - 1,
     handle_bit(selector, ACC80_HANDLE_BITS - 1));
  for (size_t bit = ACC80_HANDLE_BITS - 1; bit > tail_bit;) {
    bit--;
    ginv(inverse, current, ctx);
    di(choice, ba, bit, handle_bit(selector, bit));
    compose(next, choice, inverse, ctx);
    {
      acc80_gmat_t *swap = current;
      current = next;
      next = swap;
    }
  }
  memcpy(out, current, sizeof(*out));
  ok = 1;
cleanup:
  free(inverse);
  free(next);
  free(choice);
  free(current);
  return ok;
}

int acc80_evalfx_from_tail(acc80_evalfx_delta_t *out,
                           const acc80_ba_t *ba,
                           uint32_t selector, uint32_t member,
                           const acc80_poly_t *value, size_t tail_bit,
                           const acc80_gmat_t *tail_checkpoint,
                           const acc80_ctx_t *ctx) {
  acc80_gmat_t *tails = calloc(tail_bit + 1, sizeof(*tails));
  acc80_gmat_t *choice = malloc(sizeof(*choice));
  acc80_gmat_t *next = malloc(sizeof(*next));
  acc80_ginv_t *inverse = malloc(sizeof(*inverse));
  acc80_poly_t *input = calloc(ACC80_GADGET_DIGITS, sizeof(*input));
  acc80_poly_t *applied = calloc(ACC80_GADGET_DIGITS, sizeof(*applied));
  uint32_t difference = selector ^ member;
  size_t first_difference = 0;
  int ok = 0;
  if (out == NULL || ba == NULL || value == NULL || tail_checkpoint == NULL ||
      ctx == NULL || tail_bit == 0 || tail_bit >= ACC80_HANDLE_BITS ||
      difference == 0 || tails == NULL || choice == NULL || next == NULL ||
      inverse == NULL || input == NULL || applied == NULL)
    goto cleanup;
  while ((difference & 1U) == 0U) {
    difference >>= 1;
    first_difference++;
  }
  if (first_difference >= tail_bit) goto cleanup;

  memset(out, 0, sizeof(*out));
  memcpy(&tails[tail_bit], tail_checkpoint, sizeof(*tail_checkpoint));
  for (size_t bit = tail_bit; bit > 1;) {
    bit--;
    ginv(inverse, &tails[bit + 1], ctx);
    di(choice, ba, bit, handle_bit(selector, bit));
    compose(next, choice, inverse, ctx);
    memcpy(&tails[bit], next, sizeof(*next));
  }
  ginv_vec(input, value, ctx);
  for (size_t bit = 0; bit <= first_difference; bit++) {
    ginv(inverse, &tails[bit + 1], ctx);
    ginv_apply(applied, inverse, input, ctx);
    for (size_t digit = 0; digit < ACC80_GADGET_DIGITS; digit++) {
      if (handle_bit(selector, bit))
        acc80_poly_set(&out->right[bit][digit], &applied[digit]);
      else
        acc80_poly_neg(&out->right[bit][digit], &applied[digit]);
    }
  }
  ok = 1;
cleanup:
  free(applied);
  free(input);
  free(inverse);
  free(next);
  free(choice);
  free(tails);
  return ok;
}

void acc80_evalfx(acc80_evalfx_delta_t *out, const acc80_ba_t *ba,
                  uint32_t selector, uint32_t member,
                  const acc80_poly_t *value, const acc80_ctx_t *ctx) {
  uint32_t difference = selector ^ member;
  size_t layers = ACC80_HANDLE_BITS;

  memset(out, 0, sizeof(*out));
  if (difference != 0) {
    layers = 1;
    while ((difference & 1U) == 0U) {
      difference >>= 1;
      layers++;
    }
  }
  if (!acc80_evalfx_prefix(out->right, layers, ba, selector, value, ctx)) {
    abort();
  }
}

void acc80_member_matrix_apply(acc80_poly_t *out, const acc80_ba_t *ba,
                               uint32_t member,
                               const acc80_evalfx_delta_t *delta,
                               const acc80_ctx_t *ctx) {
  acc80_gmat_t *identity = malloc(sizeof(*identity));
  acc80_poly_t *matrix = malloc(sizeof(*matrix));
  acc80_poly_t *product = malloc(sizeof(*product));
  acc80_poly_t *next = malloc(sizeof(*next));
  if (identity == NULL || matrix == NULL || product == NULL || next == NULL)
    abort();
  acc80_gadget_identity(identity);
  acc80_poly_zero(out);
  for (size_t bit = 0; bit < ACC80_HANDLE_BITS; bit++) {
    for (size_t digit = 0; digit < ACC80_GADGET_DIGITS; digit++) {
      if (handle_bit(member, bit)) {
        acc80_poly_sub(matrix, &ba->bit[bit].block[digit],
                       &identity->block[digit]);
      } else {
        acc80_poly_set(matrix, &ba->bit[bit].block[digit]);
      }
      acc80_poly_mul(product, matrix, &delta->right[bit][digit], ctx);
      acc80_poly_add(next, out, product);
      acc80_poly_set(out, next);
    }
  }
  free(next);
  free(product);
  free(matrix);
  free(identity);
}
