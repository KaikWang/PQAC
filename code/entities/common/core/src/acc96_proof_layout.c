#include "acc96_proof_layout.h"

size_t acc96_proof_right_slot(size_t bit, size_t digit) {
  if (bit >= ACC96_HANDLE_BITS || digit >= ACC96_GADGET_DIGITS)
    return ACC96_PROOF_WITNESS_POLYS;
  return ACC96_PROOF_SLOT_RIGHT_FIRST +
         bit * ACC96_GADGET_DIGITS + digit;
}

const acc96_poly_t *acc96_proof_witness_poly_at(
    const acc96_relation_witness_t *witness, size_t slot) {
  if (witness == NULL || slot >= ACC96_PROOF_WITNESS_POLYS) return NULL;
  if (slot == ACC96_PROOF_SLOT_LEFT) return &witness->left;
  if (slot == ACC96_PROOF_SLOT_V12) return &witness->v12;
  if (slot >= ACC96_PROOF_SLOT_V2_FIRST && slot < ACC96_PROOF_SLOT_TAG)
    return &witness->v2[slot - ACC96_PROOF_SLOT_V2_FIRST];
  if (slot == ACC96_PROOF_SLOT_TAG) return &witness->tag;
  slot -= ACC96_PROOF_SLOT_RIGHT_FIRST;
  return &witness->right.right[slot / ACC96_GADGET_DIGITS]
                                [slot % ACC96_GADGET_DIGITS];
}

acc96_poly_t *acc96_proof_witness_poly_at_mut(
    acc96_relation_witness_t *witness, size_t slot) {
  if (witness == NULL || slot >= ACC96_PROOF_WITNESS_POLYS) return NULL;
  if (slot == ACC96_PROOF_SLOT_LEFT) return &witness->left;
  if (slot == ACC96_PROOF_SLOT_V12) return &witness->v12;
  if (slot >= ACC96_PROOF_SLOT_V2_FIRST && slot < ACC96_PROOF_SLOT_TAG)
    return &witness->v2[slot - ACC96_PROOF_SLOT_V2_FIRST];
  if (slot == ACC96_PROOF_SLOT_TAG) return &witness->tag;
  slot -= ACC96_PROOF_SLOT_RIGHT_FIRST;
  return &witness->right.right[slot / ACC96_GADGET_DIGITS]
                                [slot % ACC96_GADGET_DIGITS];
}

int acc96_proof_slot_is_tag(size_t slot) {
  return slot == ACC96_PROOF_SLOT_TAG;
}

int acc96_lift_proof_limb_decompose(
    int32_t limbs[ACC96_LIFT_PROOF_LIMBS], acc96_i128 coefficient) {
  const acc96_i128 radix = ACC96_LIFT_PROOF_LIMB_RADIX;
  acc96_i128 remaining = coefficient;
  const acc96_u128 magnitude = coefficient < 0
      ? (acc96_u128)(-(coefficient + 1)) + 1U
      : (acc96_u128)coefficient;
  if (limbs == NULL || magnitude > ACC96_LIFT_COEFF_BOUND) return 0;
  for (size_t limb = 0; limb < ACC96_LIFT_PROOF_LIMBS; limb++) {
    acc96_i128 digit = remaining % radix;
    if (digit < 0) digit += radix;
    if (digit > ACC96_LIFT_PROOF_LIMB_BOUND) digit -= radix;
    limbs[limb] = (int32_t)digit;
    remaining = (remaining - digit) / radix;
  }
  return remaining == 0;
}

acc96_i128 acc96_lift_proof_limb_reconstruct(
    const int32_t limbs[ACC96_LIFT_PROOF_LIMBS]) {
  acc96_i128 result = 0;
  acc96_i128 power = 1;
  for (size_t limb = 0; limb < ACC96_LIFT_PROOF_LIMBS; limb++) {
    result += (acc96_i128)limbs[limb] * power;
    power *= ACC96_LIFT_PROOF_LIMB_RADIX;
  }
  return result;
}

int acc96_proof_limb_decompose(
    int32_t limbs[ACC96_PROOF_LIMBS],
    const acc96_coeff_t *coefficient,
    const acc96_ctx_t *ctx) {
  acc96_i128 remaining;
  const acc96_i128 radix = (acc96_i128)ACC96_PROOF_LIMB_RADIX;

  if (limbs == NULL || coefficient == NULL) return 0;
  remaining = acc96_coeff_to_centered(coefficient, ctx);
  for (size_t i = 0; i < ACC96_PROOF_LIMBS; i++) {
    acc96_i128 digit = remaining % radix;
    if (digit < 0) digit += radix;
    if (digit > ACC96_PROOF_LIMB_BOUND) digit -= radix;
    limbs[i] = (int32_t)digit;
    remaining = (remaining - digit) / radix;
  }
  return remaining == 0;
}

acc96_i128 acc96_proof_limb_reconstruct(
    const int32_t limbs[ACC96_PROOF_LIMBS]) {
  acc96_i128 result = 0;
  acc96_i128 power = 1;
  for (size_t i = 0; i < ACC96_PROOF_LIMBS; i++) {
    result += (acc96_i128)limbs[i] * power;
    power *= (acc96_i128)ACC96_PROOF_LIMB_RADIX;
  }
  return result;
}

int acc96_proof_poly_decompose(
    int32_t limbs[ACC96_PROOF_LIMBS][ACC96_N],
    const acc96_poly_t *polynomial, const acc96_ctx_t *ctx) {
  if (limbs == NULL || polynomial == NULL) return 0;
  for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
    int32_t digits[ACC96_PROOF_LIMBS];
    if (!acc96_proof_limb_decompose(digits, &polynomial->coeff[coefficient],
                                    ctx))
      return 0;
    for (size_t limb = 0; limb < ACC96_PROOF_LIMBS; limb++)
      limbs[limb][coefficient] = digits[limb];
  }
  return 1;
}

uint64_t acc96_proof_poly_max_chunk_norm2(
    const int32_t limbs[ACC96_PROOF_LIMBS][ACC96_N]) {
  uint64_t maximum = 0;
  for (size_t limb = 0; limb < ACC96_PROOF_LIMBS; limb++) {
    for (size_t phase = 0; phase < ACC96_LABRADOR_PHASES; phase++) {
      uint64_t norm = 0;
      for (size_t i = 0; i < ACC96_LABRADOR_BASE_DEGREE; i++) {
        const int64_t value = limbs[limb]
            [acc96_proof_polyphase_source_index(phase, i)];
        norm += (uint64_t)(value * value);
      }
      if (norm > maximum) maximum = norm;
    }
  }
  return maximum;
}

int acc96_window_proof_limb_decompose(
    int32_t limbs[ACC96_WINDOW_PROOF_LIMBS],
    const acc96_coeff_t *coefficient, const acc96_ctx_t *ctx) {
  acc96_i128 remaining;
  const acc96_i128 radix = ACC96_WINDOW_PROOF_LIMB_RADIX;
  if (limbs == NULL || coefficient == NULL || ctx == NULL) return 0;
  remaining = acc96_coeff_to_centered(coefficient, ctx);
  if (remaining < -ACC96_WINDOW_COEFF_BOUND ||
      remaining > ACC96_WINDOW_COEFF_BOUND)
    return 0;
  for (size_t limb = 0; limb < ACC96_WINDOW_PROOF_LIMBS; limb++) {
    acc96_i128 digit = remaining % radix;
    if (digit < 0) digit += radix;
    if (digit > ACC96_WINDOW_PROOF_LIMB_BOUND) digit -= radix;
    limbs[limb] = (int32_t)digit;
    remaining = (remaining - digit) / radix;
  }
  return remaining == 0;
}

int64_t acc96_window_proof_limb_reconstruct(
    const int32_t limbs[ACC96_WINDOW_PROOF_LIMBS]) {
  int64_t result = 0;
  int64_t power = 1;
  for (size_t limb = 0; limb < ACC96_WINDOW_PROOF_LIMBS; limb++) {
    result += (int64_t)limbs[limb] * power;
    power *= ACC96_WINDOW_PROOF_LIMB_RADIX;
  }
  return result;
}

int acc96_window_proof_poly_decompose(
    int32_t limbs[ACC96_WINDOW_PROOF_LIMBS][ACC96_N],
    const acc96_poly_t *polynomial, const acc96_ctx_t *ctx) {
  if (limbs == NULL || polynomial == NULL || ctx == NULL) return 0;
  for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
    int32_t digits[ACC96_WINDOW_PROOF_LIMBS];
    if (!acc96_window_proof_limb_decompose(
            digits, &polynomial->coeff[coefficient], ctx))
      return 0;
    for (size_t limb = 0; limb < ACC96_WINDOW_PROOF_LIMBS; limb++)
      limbs[limb][coefficient] = digits[limb];
  }
  return 1;
}

uint64_t acc96_window_proof_poly_max_chunk_norm2(
    const int32_t limbs[ACC96_WINDOW_PROOF_LIMBS][ACC96_N]) {
  uint64_t maximum = 0;
  for (size_t limb = 0; limb < ACC96_WINDOW_PROOF_LIMBS; limb++) {
    for (size_t phase = 0; phase < ACC96_LABRADOR_PHASES; phase++) {
      uint64_t norm = 0;
      for (size_t i = 0; i < ACC96_LABRADOR_BASE_DEGREE; i++) {
        const int64_t value = limbs[limb]
            [acc96_proof_polyphase_source_index(phase, i)];
        norm += (uint64_t)(value * value);
      }
      if (norm > maximum) maximum = norm;
    }
  }
  return maximum;
}

size_t acc96_proof_polyphase_source_index(size_t phase,
                                          size_t coefficient) {
  if (phase >= ACC96_LABRADOR_PHASES ||
      coefficient >= ACC96_LABRADOR_BASE_DEGREE)
    return ACC96_N;
  return ACC96_LABRADOR_PHASES * coefficient + phase;
}

static int window_check_poly(const acc96_poly_t *polynomial,
                             const acc96_ctx_t *ctx,
                             uint64_t *maximum) {
  int valid = 1;
  for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
    const acc96_i128 centered =
        acc96_coeff_to_centered(&polynomial->coeff[coefficient], ctx);
    const acc96_u128 magnitude = centered < 0
        ? (acc96_u128)(-centered) : (acc96_u128)centered;
    if (magnitude > UINT64_MAX) {
      *maximum = UINT64_MAX;
      valid = 0;
    } else if ((uint64_t)magnitude > *maximum) {
      *maximum = (uint64_t)magnitude;
    }
    if (magnitude > (acc96_u128)ACC96_WINDOW_COEFF_BOUND) valid = 0;
  }
  return valid;
}

int acc96_window_proof_witness_check(
    const acc96_relation_witness_t *witness,
    const acc96_ctx_t *ctx,
    uint64_t *max_abs) {
  uint64_t maximum = 0;

  if (witness == NULL || ctx == NULL || !ctx->initialized) return 0;
  int valid = window_check_poly(&witness->left, ctx, &maximum);
  valid &= window_check_poly(&witness->v12, ctx, &maximum);
  valid &= window_check_poly(&witness->tag, ctx, &maximum);
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++)
    valid &= window_check_poly(&witness->v2[i], ctx, &maximum);
  for (size_t bit = 0; bit < ACC96_HANDLE_BITS; bit++)
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++)
      valid &= window_check_poly(&witness->right.right[bit][digit], ctx,
                                 &maximum);
  if (max_abs != NULL) *max_abs = maximum;
  return valid;
}
