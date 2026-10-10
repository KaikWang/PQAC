#include "revac_q96_link_wire.h"

#include <string.h>

static void put_u32(uint8_t out[4], uint32_t value) {
  out[0] = (uint8_t)value;
  out[1] = (uint8_t)(value >> 8);
  out[2] = (uint8_t)(value >> 16);
  out[3] = (uint8_t)(value >> 24);
}

static void put_i32(uint8_t out[4], int32_t value) {
  put_u32(out, (uint32_t)value);
}

static int decompose_i64(int32_t out[REVAC_Q96_LINK_DECOMP_LIMBS],
                         int64_t value) {
  int64_t remaining = value;
  for (size_t limb = 0; limb < REVAC_Q96_LINK_DECOMP_LIMBS; limb++) {
    int64_t digit = remaining % REVAC_Q96_LINK_DECOMP_BASE;
    if (digit > REVAC_Q96_LINK_DECOMP_BASE / 2) {
      digit -= REVAC_Q96_LINK_DECOMP_BASE;
    } else if (digit < -(REVAC_Q96_LINK_DECOMP_BASE / 2)) {
      digit += REVAC_Q96_LINK_DECOMP_BASE;
    }
    out[limb] = (int32_t)digit;
    remaining = (remaining - digit) / REVAC_Q96_LINK_DECOMP_BASE;
  }
  return remaining == 0;
}

int revac_q96_link_wire_encode(
    uint8_t *out, size_t out_len,
    const revac_public_state_t *public_state, uint32_t handle,
    const poly_q_vec_d randomness[2], const poly_q_vec_d commitment) {
  int64_t lift[PARAM_D][PARAM_N];
  uint8_t accumulator_digest[REVAC_CONTEXT_HASH_BYTES];
  size_t offset = 0;

  if (out == NULL || public_state == NULL || randomness == NULL ||
      handle >= (UINT32_C(1) << ACC96_HANDLE_BITS) ||
      commitment == NULL || out_len != REVAC_Q96_LINK_WIRE_BYTES ||
      public_state->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      !revac_public_accumulator_digest(accumulator_digest, public_state) ||
      !revac_same_handle_exact_lift(lift, public_state, handle, randomness,
                                    commitment)) {
    return 0;
  }
  memcpy(out + offset, REVAC_Q96_LINK_MAGIC, REVAC_Q96_LINK_MAGIC_BYTES);
  offset += REVAC_Q96_LINK_MAGIC_BYTES;
  put_u32(out + offset, public_state->window_id);
  offset += 4;
  memcpy(out + offset, public_state->context_pp_hash,
         REVAC_CONTEXT_HASH_BYTES);
  offset += REVAC_CONTEXT_HASH_BYTES;
  memcpy(out + offset, public_state->parameter_digest,
         REVAC_CONTEXT_HASH_BYTES);
  offset += REVAC_CONTEXT_HASH_BYTES;
  memcpy(out + offset, accumulator_digest, REVAC_CONTEXT_HASH_BYTES);
  offset += REVAC_CONTEXT_HASH_BYTES;

  for (size_t row = 0; row < PARAM_D; row++) {
    for (size_t coefficient = 0; coefficient < PARAM_N; coefficient++) {
      put_u32(out + offset, (uint32_t)poly_q_get_coeff(
                                commitment->entries[row], coefficient));
      offset += 4;
    }
  }
  for (size_t part = 0; part < 2; part++) {
    for (size_t row = 0; row < PARAM_D; row++) {
      for (size_t coefficient = 0; coefficient < PARAM_N; coefficient++) {
        put_i32(out + offset, (int32_t)poly_q_get_coeff_centered(
                                  randomness[part]->entries[row],
                                  coefficient));
        offset += 4;
      }
    }
  }
  for (size_t row = 0; row < PARAM_D; row++) {
    const size_t row_offset = offset;
    for (size_t coefficient = 0; coefficient < PARAM_N; coefficient++) {
      int32_t digits[REVAC_Q96_LINK_DECOMP_LIMBS];
      if (!decompose_i64(digits, lift[row][coefficient])) return 0;
      for (size_t limb = 0; limb < REVAC_Q96_LINK_DECOMP_LIMBS; limb++) {
        put_i32(out + row_offset +
                    (limb * PARAM_N + coefficient) * 4U,
                digits[limb]);
      }
    }
    offset += REVAC_Q96_LINK_DECOMP_LIMBS * REVAC_Q96_LINK_POLY_BYTES;
  }
  return offset == out_len;
}

int revac_q96_instance_wire_encode(
    uint8_t *out, size_t out_len, const revac_public_state_t *public_state,
    const poly_q_vec_d commitment) {
  size_t offset = 0;
  uint8_t accumulator_digest[REVAC_CONTEXT_HASH_BYTES];
  if (out == NULL || public_state == NULL || commitment == NULL ||
      public_state->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      public_state->acc96_acc == NULL ||
      !revac_public_accumulator_digest(accumulator_digest, public_state) ||
      out_len != REVAC_Q96_INSTANCE_WIRE_BYTES)
    return 0;
  memcpy(out + offset, REVAC_Q96_INSTANCE_MAGIC,
         REVAC_Q96_LINK_MAGIC_BYTES);
  offset += REVAC_Q96_LINK_MAGIC_BYTES;
  put_u32(out + offset, public_state->acc96_seq);
  offset += 4;
  put_u32(out + offset, public_state->window_id);
  offset += 4;
  offset += acc96_poly_pack(out + offset, public_state->acc96_acc);
  memcpy(out + offset, public_state->context_pp_hash,
         REVAC_CONTEXT_HASH_BYTES);
  offset += REVAC_CONTEXT_HASH_BYTES;
  memcpy(out + offset, public_state->parameter_digest,
         REVAC_CONTEXT_HASH_BYTES);
  offset += REVAC_CONTEXT_HASH_BYTES;
  memcpy(out + offset, accumulator_digest, REVAC_CONTEXT_HASH_BYTES);
  offset += REVAC_CONTEXT_HASH_BYTES;
  for (size_t row = 0; row < PARAM_D; row++) {
    for (size_t coefficient = 0; coefficient < PARAM_N; coefficient++) {
      put_u32(out + offset, (uint32_t)poly_q_get_coeff(
                                commitment->entries[row], coefficient));
      offset += 4;
    }
  }
  return offset == out_len;
}
