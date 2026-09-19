#include "acc96_relation.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "fips202.h"

static acc96_u128 trunc_power(size_t exponent) {
  acc96_u128 value = 1;
  while (exponent-- != 0) value *= ACC96_TRUNC_BASE;
  return value;
}

static void absorb_u32_le(keccak_state *state, uint32_t value) {
  const uint8_t encoded[4] = {
      (uint8_t)value, (uint8_t)(value >> 8),
      (uint8_t)(value >> 16), (uint8_t)(value >> 24)};
  shake256_absorb(state, encoded, sizeof(encoded));
}

static void expand_poly(acc96_poly_t *out,
                        const uint8_t seed[ACC96_SEED_BYTES],
                        uint32_t domain, uint32_t first, uint32_t second) {
  keccak_state state;
  uint8_t encoded[ACC96_COEFF_PACKED_BYTES];

  shake256_init(&state);
  shake256_absorb(&state, seed, ACC96_SEED_BYTES);
  absorb_u32_le(&state, domain);
  absorb_u32_le(&state, first);
  absorb_u32_le(&state, second);
  shake256_finalize(&state);
  for (size_t coeff = 0; coeff < ACC96_N; coeff++) {
    acc96_u128 value;
    do {
      value = 0;
      shake256_squeeze(encoded, sizeof(encoded), &state);
      for (size_t byte = 0; byte < sizeof(encoded); byte++)
        value |= (acc96_u128)encoded[byte] << (8 * byte);
    } while (value >= acc96_modulus());
    acc96_coeff_from_u128(&out->coeff[coeff], value);
  }
}

void acc96_relation_uniform_poly(
    acc96_poly_t *out, const uint8_t seed[ACC96_SEED_BYTES],
    uint32_t domain, uint32_t first, uint32_t second) {
  expand_poly(out, seed, domain, first, second);
}

static void sample_ternary_poly(acc96_poly_t *out,
                                const uint8_t seed[ACC96_SEED_BYTES],
                                uint32_t domain, uint32_t first,
                                uint32_t second) {
  keccak_state state;
  uint8_t byte = 0;

  acc96_poly_zero(out);
  shake256_init(&state);
  shake256_absorb(&state, seed, ACC96_SEED_BYTES);
  absorb_u32_le(&state, domain);
  absorb_u32_le(&state, first);
  absorb_u32_le(&state, second);
  shake256_finalize(&state);
  for (size_t coeff = 0; coeff < ACC96_N; coeff++) {
    do {
      shake256_squeeze(&byte, 1, &state);
      byte &= 3U;
    } while (byte == 3U);
    if (byte == 2U)
      acc96_coeff_from_u128(&out->coeff[coeff], acc96_modulus() - 1);
    else
      acc96_coeff_from_u128(&out->coeff[coeff], byte);
  }
}

int acc96_relation_pp_expand(acc96_relation_pp_t *pp,
                             const uint8_t seed[ACC96_SEED_BYTES],
                             const acc96_ctx_t *ctx) {
  if (pp == NULL || seed == NULL || ctx == NULL || !ctx->initialized)
    return 0;
  expand_poly(&pp->atrap, seed, 0x41393631U, 0, 0);
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++)
    expand_poly(&pp->trunc_b[i], seed, 0x42393631U, (uint32_t)i, 0);
  for (size_t bit = 0; bit < ACC96_HANDLE_BITS; bit++) {
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++) {
      expand_poly(&pp->ba.bit[bit].block[digit], seed, 0x44393631U,
                  (uint32_t)bit, (uint32_t)digit);
    }
  }
  expand_poly(&pp->ua, seed, 0x55393631U, 0, 0);
  return 1;
}

int acc96_relation_pp_reconstruct(
    acc96_relation_pp_t *pp,
    const uint8_t seed[ACC96_SEED_BYTES],
    const acc96_poly_t trunc_b[ACC96_TRUNC_RIGHT],
    const acc96_ctx_t *ctx) {
  if (pp == NULL || seed == NULL || trunc_b == NULL ||
      !acc96_relation_pp_expand(pp, seed, ctx))
    return 0;
  memcpy(pp->trunc_b, trunc_b, sizeof(pp->trunc_b));
  return 1;
}

static void add_product(acc96_poly_t *acc, const acc96_poly_t *lhs,
                        const acc96_poly_t *rhs, int subtract,
                        const acc96_ctx_t *ctx) {
  acc96_poly_t *product = malloc(sizeof(*product));
  acc96_poly_t *next = malloc(sizeof(*next));
  if (product == NULL || next == NULL) abort();
  acc96_poly_mul(product, lhs, rhs, ctx);
  if (subtract)
    acc96_poly_sub(next, acc, product);
  else
    acc96_poly_add(next, acc, product);
  acc96_poly_set(acc, next);
  free(next);
  free(product);
}

static void add_member_right(acc96_poly_t *out,
                             const acc96_relation_pp_t *pp,
                             uint32_t handle,
                             const acc96_evalfx_delta_t *right,
                             const acc96_ctx_t *ctx) {
  acc96_poly_t *contribution = malloc(sizeof(*contribution));
  if (contribution == NULL) abort();
  acc96_member_matrix_apply(contribution, &pp->ba, handle, right, ctx);
  acc96_poly_add(out, out, contribution);
  free(contribution);
}

void acc96_relation_evaluate(acc96_poly_t *out,
                             const acc96_relation_pp_t *pp,
                             uint32_t handle,
                             const acc96_relation_witness_t *witness,
                             const acc96_ctx_t *ctx) {
  acc96_poly_t *matrix = malloc(sizeof(*matrix));
  acc96_poly_t *scalar = malloc(sizeof(*scalar));
  acc96_poly_t *tag_term = malloc(sizeof(*tag_term));
  if (matrix == NULL || scalar == NULL || tag_term == NULL) abort();

  acc96_poly_set(out, &witness->left);
  add_product(out, &pp->atrap, &witness->v12, 0, ctx);
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++) {
    acc96_poly_zero(scalar);
    acc96_coeff_from_u128(&scalar->coeff[0], trunc_power(ACC96_TRUNC_L + i));
    acc96_poly_mul(tag_term, &witness->tag, scalar, ctx);
    acc96_poly_sub(matrix, &pp->trunc_b[i], tag_term);
    add_product(out, matrix, &witness->v2[i], 1, ctx);
  }
  add_member_right(out, pp, handle, &witness->right, ctx);
  free(tag_term);
  free(scalar);
  free(matrix);
}

static long double poly_norm2(const acc96_poly_t *poly,
                              const acc96_ctx_t *ctx) {
  long double result = 0;
  for (size_t i = 0; i < ACC96_N; i++) {
    const acc96_i128 value = acc96_coeff_to_centered(&poly->coeff[i], ctx);
    const long double converted = (long double)value;
    result += converted * converted;
  }
  return result;
}

long double acc96_relation_witness_norm2(
    const acc96_relation_witness_t *witness, const acc96_ctx_t *ctx) {
  long double result = poly_norm2(&witness->left, ctx) +
                       poly_norm2(&witness->v12, ctx) +
                       poly_norm2(&witness->tag, ctx);
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++)
    result += poly_norm2(&witness->v2[i], ctx);
  for (size_t bit = 0; bit < ACC96_HANDLE_BITS; bit++)
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++)
      result += poly_norm2(&witness->right.right[bit][digit], ctx);
  return result;
}

int acc96_relation_memver(const acc96_relation_pp_t *pp,
                          const acc96_poly_t *accumulator,
                          uint32_t handle,
                          const acc96_relation_witness_t *witness,
                          long double max_norm_sq,
                          const acc96_ctx_t *ctx) {
  acc96_poly_t *evaluated = malloc(sizeof(*evaluated));
  int ok;
  if (evaluated == NULL) return 0;
  acc96_relation_evaluate(evaluated, pp, handle, witness, ctx);
  ok = acc96_poly_equal(evaluated, accumulator);
  if (ok && max_norm_sq > 0)
    ok = acc96_relation_witness_norm2(witness, ctx) <= max_norm_sq;
  free(evaluated);
  return ok;
}

void acc96_relation_delete(acc96_poly_t *accumulator,
                           acc96_relation_update_t *update,
                           const acc96_relation_pp_t *pp,
                           uint32_t revoked_handle,
                           uint32_t next_sequence,
                           const acc96_ctx_t *ctx) {
  acc96_gmat_t *indicator = malloc(sizeof(*indicator));
  acc96_poly_t *delta = malloc(sizeof(*delta));
  acc96_poly_t *next = malloc(sizeof(*next));
  if (indicator == NULL || delta == NULL || next == NULL) abort();
  acc96_evalf(indicator, &pp->ba, revoked_handle, ctx);
  acc96_evalf_apply(delta, indicator, &pp->ua, ctx);
  acc96_poly_add(next, accumulator, delta);
  acc96_poly_set(accumulator, next);
  update->sequence = next_sequence;
  update->revoked_handle = revoked_handle;
  free(next);
  free(delta);
  free(indicator);
}

static int relation_delete_delta_prefix(
    acc96_poly_t *accumulator, acc96_relation_update_t *update,
    acc96_poly_t delta[][ACC96_GADGET_DIGITS], size_t layers,
    const acc96_relation_pp_t *pp,
    uint32_t revoked_handle, uint32_t next_sequence,
    const acc96_ctx_t *ctx) {
  acc96_gmat_t *indicator = malloc(sizeof(*indicator));
  acc96_poly_t *acc_delta = malloc(sizeof(*acc_delta));
  acc96_poly_t *next = malloc(sizeof(*next));
  int ok = 0;
  if (accumulator == NULL || update == NULL || delta == NULL || pp == NULL ||
      layers == 0 || layers > ACC96_HANDLE_BITS ||
      ctx == NULL || indicator == NULL || acc_delta == NULL || next == NULL)
    goto cleanup;
  if (!acc96_evalf_evalfx_prefix(indicator, delta, layers, &pp->ba,
                                 revoked_handle, &pp->ua, ctx))
    goto cleanup;
  acc96_evalf_apply(acc_delta, indicator, &pp->ua, ctx);
  acc96_poly_add(next, accumulator, acc_delta);
  acc96_poly_set(accumulator, next);
  update->sequence = next_sequence;
  update->revoked_handle = revoked_handle;
  ok = 1;
cleanup:
  free(next);
  free(acc_delta);
  free(indicator);
  return ok;
}

int acc96_relation_delete_delta4(
    acc96_poly_t *accumulator, acc96_relation_update_t *update,
    acc96_delta4_t *delta, const acc96_relation_pp_t *pp,
    uint32_t revoked_handle, uint32_t next_sequence,
    const acc96_ctx_t *ctx) {
  if (delta == NULL) return 0;
  return relation_delete_delta_prefix(accumulator, update, delta->right,
                                      ACC96_DELTA4_LAYERS, pp,
                                      revoked_handle, next_sequence, ctx);
}

int acc96_relation_delete_delta6(
    acc96_poly_t *accumulator, acc96_relation_update_t *update,
    acc96_delta6_t *delta, const acc96_relation_pp_t *pp,
    uint32_t revoked_handle, uint32_t next_sequence,
    const acc96_ctx_t *ctx) {
  if (delta == NULL) return 0;
  return relation_delete_delta_prefix(accumulator, update, delta->right,
                                      ACC96_DELTA6_LAYERS, pp,
                                      revoked_handle, next_sequence, ctx);
}

int acc96_relation_delete_tail4(
    acc96_poly_t *accumulator, acc96_relation_update_t *update,
    acc96_gmat_t *tail4, const acc96_relation_pp_t *pp,
    uint32_t revoked_handle, uint32_t next_sequence,
    const acc96_ctx_t *ctx) {
  acc96_gmat_t *indicator = malloc(sizeof(*indicator));
  acc96_poly_t *acc_delta = malloc(sizeof(*acc_delta));
  acc96_poly_t *next = malloc(sizeof(*next));
  int ok = 0;
  if (accumulator == NULL || update == NULL || tail4 == NULL || pp == NULL ||
      ctx == NULL || indicator == NULL || acc_delta == NULL || next == NULL)
    goto cleanup;
  if (!acc96_evalf_with_tail(indicator, tail4, ACC96_TAIL4_BIT, &pp->ba,
                             revoked_handle, ctx))
    goto cleanup;
  acc96_evalf_apply(acc_delta, indicator, &pp->ua, ctx);
  acc96_poly_add(next, accumulator, acc_delta);
  acc96_poly_set(accumulator, next);
  update->sequence = next_sequence;
  update->revoked_handle = revoked_handle;
  ok = 1;
cleanup:
  free(next);
  free(acc_delta);
  free(indicator);
  return ok;
}

int acc96_relation_memwitup(acc96_relation_witness_t *witness,
                            const acc96_relation_pp_t *pp,
                            uint32_t member_handle,
                            const acc96_relation_update_t *update,
                            const acc96_ctx_t *ctx) {
  acc96_evalfx_delta_t *delta;
  if (witness == NULL || pp == NULL || update == NULL ||
      update->revoked_handle == member_handle)
    return 0;
  delta = malloc(sizeof(*delta));
  if (delta == NULL) return 0;
  acc96_evalfx(delta, &pp->ba, update->revoked_handle, member_handle,
               &pp->ua, ctx);
  for (size_t bit = 0; bit < ACC96_HANDLE_BITS; bit++) {
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++) {
      acc96_poly_add(&witness->right.right[bit][digit],
                     &witness->right.right[bit][digit],
                     &delta->right[bit][digit]);
    }
  }
  free(delta);
  return 1;
}

int acc96_relation_delta4(acc96_delta4_t *delta,
                          const acc96_relation_pp_t *pp,
                          uint32_t revoked_handle,
                          const acc96_ctx_t *ctx) {
  if (delta == NULL || pp == NULL || ctx == NULL) return 0;
  return acc96_evalfx_prefix(delta->right, ACC96_DELTA4_LAYERS, &pp->ba,
                             revoked_handle, &pp->ua, ctx);
}

int acc96_relation_delta6(acc96_delta6_t *delta,
                          const acc96_relation_pp_t *pp,
                          uint32_t revoked_handle,
                          const acc96_ctx_t *ctx) {
  if (delta == NULL || pp == NULL || ctx == NULL) return 0;
  return acc96_evalfx_prefix(delta->right, ACC96_DELTA6_LAYERS, &pp->ba,
                             revoked_handle, &pp->ua, ctx);
}

static int relation_memwitup_delta_prefix(
    acc96_relation_witness_t *witness, uint32_t member_handle,
    const acc96_relation_update_t *update,
    const acc96_poly_t delta[][ACC96_GADGET_DIGITS], size_t available_layers,
    size_t *applied_layers) {
  uint32_t difference;
  size_t first_difference = 0;

  if (applied_layers != NULL) *applied_layers = 0;
  if (witness == NULL || update == NULL || delta == NULL ||
      available_layers == 0 || available_layers > ACC96_HANDLE_BITS ||
      update->revoked_handle == member_handle) {
    return 0;
  }
  difference = update->revoked_handle ^ member_handle;
  while ((difference & 1U) == 0U) {
    difference >>= 1;
    first_difference++;
  }
  if (first_difference >= available_layers) return 0;
  for (size_t bit = 0; bit <= first_difference; bit++) {
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++) {
      acc96_poly_add(&witness->right.right[bit][digit],
                     &witness->right.right[bit][digit],
                     &delta[bit][digit]);
    }
  }
  if (applied_layers != NULL) *applied_layers = first_difference + 1;
  return 1;
}

int acc96_relation_memwitup_delta4(
    acc96_relation_witness_t *witness, uint32_t member_handle,
    const acc96_relation_update_t *update, const acc96_delta4_t *delta,
    size_t *applied_layers) {
  if (delta == NULL) return 0;
  return relation_memwitup_delta_prefix(
      witness, member_handle, update, delta->right, ACC96_DELTA4_LAYERS,
      applied_layers);
}

int acc96_relation_memwitup_delta6(
    acc96_relation_witness_t *witness, uint32_t member_handle,
    const acc96_relation_update_t *update, const acc96_delta6_t *delta,
    size_t available_layers, size_t *applied_layers) {
  if (delta == NULL || available_layers > ACC96_DELTA6_LAYERS) return 0;
  return relation_memwitup_delta_prefix(
      witness, member_handle, update, delta->right, available_layers,
      applied_layers);
}

int acc96_relation_memwitup_tail4(
    acc96_relation_witness_t *witness, const acc96_relation_pp_t *pp,
    uint32_t member_handle, const acc96_relation_update_t *update,
    const acc96_gmat_t *tail4, const acc96_ctx_t *ctx) {
  acc96_evalfx_delta_t *delta;
  uint32_t difference;
  size_t first_difference = 0;
  if (witness == NULL || pp == NULL || update == NULL || tail4 == NULL ||
      ctx == NULL || update->revoked_handle == member_handle)
    return 0;
  difference = update->revoked_handle ^ member_handle;
  while ((difference & 1U) == 0U) {
    difference >>= 1;
    first_difference++;
  }
  if (first_difference >= ACC96_TAIL4_BIT) return 0;
  delta = malloc(sizeof(*delta));
  if (delta == NULL) return 0;
  if (!acc96_evalfx_from_tail(delta, &pp->ba, update->revoked_handle,
                              member_handle, &pp->ua, ACC96_TAIL4_BIT,
                              tail4, ctx)) {
    free(delta);
    return 0;
  }
  for (size_t bit = 0; bit <= first_difference; bit++)
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++)
      acc96_poly_add(&witness->right.right[bit][digit],
                     &witness->right.right[bit][digit],
                     &delta->right[bit][digit]);
  free(delta);
  return 1;
}

void acc96_relation_make_kat_instance(
    acc96_poly_t *accumulator, acc96_relation_witness_t *witness,
    const acc96_relation_pp_t *pp, uint32_t handle,
    const uint8_t seed[ACC96_SEED_BYTES], const acc96_ctx_t *ctx) {
  memset(witness, 0, sizeof(*witness));
  sample_ternary_poly(&witness->left, seed, 0x4c393631U, 0, 0);
  sample_ternary_poly(&witness->v12, seed, 0x56393631U, 0, 0);
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++)
    sample_ternary_poly(&witness->v2[i], seed, 0x76393631U, (uint32_t)i, 0);
  acc96_poly_zero(&witness->tag);
  acc96_coeff_from_u128(&witness->tag.coeff[0], 1);
  for (size_t bit = 0; bit < ACC96_HANDLE_BITS; bit++) {
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++) {
      sample_ternary_poly(&witness->right.right[bit][digit], seed,
                          0x52393631U, (uint32_t)bit, (uint32_t)digit);
    }
  }
  acc96_relation_evaluate(accumulator, pp, handle, witness, ctx);
}
