#include "acc96_labrador_wire.h"

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

static int tag_is_unit(const acc96_poly_t *tag) {
  for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
    const uint64_t expected = coefficient == 0 ? 1U : 0U;
    for (size_t lane = 0; lane < ACC96_LIMBS; lane++)
      if (tag->coeff[coefficient].limb[lane] != expected) return 0;
  }
  return 1;
}

int acc96_labrador_public_encode(
    uint8_t *out, size_t out_len,
    const uint8_t public_seed[ACC96_SEED_BYTES],
    const acc96_relation_pp_t *pp) {
  size_t offset = 0;
  if (out == NULL || public_seed == NULL || pp == NULL ||
      out_len != ACC96_LABRADOR_PUBLIC_BYTES)
    return 0;
  memcpy(out + offset, ACC96_LABRADOR_PUBLIC_MAGIC,
         ACC96_LABRADOR_MAGIC_BYTES);
  offset += ACC96_LABRADOR_MAGIC_BYTES;
  memcpy(out + offset, public_seed, ACC96_SEED_BYTES);
  offset += ACC96_SEED_BYTES;
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++) {
    offset += acc96_poly_pack(out + offset, &pp->trunc_b[i]);
  }
  return offset == out_len;
}

int acc96_labrador_public_decode(
    acc96_relation_pp_t *pp,
    uint8_t public_seed[ACC96_SEED_BYTES],
    const uint8_t *encoded, size_t encoded_len,
    const acc96_ctx_t *ctx) {
  acc96_poly_t trunc_b[ACC96_TRUNC_RIGHT];
  size_t offset = ACC96_LABRADOR_MAGIC_BYTES;
  if (pp == NULL || public_seed == NULL || encoded == NULL || ctx == NULL ||
      encoded_len != ACC96_LABRADOR_PUBLIC_BYTES ||
      memcmp(encoded, ACC96_LABRADOR_PUBLIC_MAGIC,
             ACC96_LABRADOR_MAGIC_BYTES) != 0)
    return 0;
  memcpy(public_seed, encoded + offset, ACC96_SEED_BYTES);
  offset += ACC96_SEED_BYTES;
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++) {
    if (!acc96_poly_unpack(&trunc_b[i], encoded + offset,
                           ACC96_POLY_PACKED_BYTES))
      return 0;
    offset += ACC96_POLY_PACKED_BYTES;
  }
  return offset == encoded_len &&
         acc96_relation_pp_reconstruct(pp, public_seed, trunc_b, ctx);
}

int acc96_labrador_prover_encode(
    uint8_t *out, size_t out_len, uint32_t sequence, uint32_t handle,
    const acc96_poly_t *accumulator,
    const acc96_relation_witness_t *witness,
    const acc96_exact_lift_t *lift,
    const acc96_ctx_t *ctx) {
  size_t offset = 0;
  if (out == NULL || accumulator == NULL || witness == NULL || lift == NULL ||
      handle >= (UINT32_C(1) << ACC96_HANDLE_BITS) ||
      ctx == NULL || out_len != ACC96_LABRADOR_PROVER_BYTES ||
      !tag_is_unit(&witness->tag) ||
      !acc96_window_proof_witness_check(witness, ctx, NULL))
    return 0;
  memcpy(out + offset, ACC96_LABRADOR_PROVER_MAGIC,
         ACC96_LABRADOR_MAGIC_BYTES);
  offset += ACC96_LABRADOR_MAGIC_BYTES;
  put_u32(out + offset, sequence);
  offset += 4;
  put_u32(out + offset, handle);
  offset += 4;
  offset += acc96_poly_pack(out + offset, accumulator);
  for (size_t bit = 0; bit < ACC96_HANDLE_BITS; bit++)
    out[offset++] = (uint8_t)((handle >> bit) & 1U);

  {
    const size_t witness_offset = offset;
    int encode_ok = 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(&:encode_ok)
#endif
    for (size_t slot = 0; slot < ACC96_PROOF_WITNESS_POLYS; slot++) {
      const acc96_poly_t *poly;
      size_t ordinal;
      size_t slot_offset;
      if (acc96_proof_slot_is_tag(slot)) continue;
      ordinal = slot - (slot > ACC96_PROOF_SLOT_TAG ? 1U : 0U);
      slot_offset = witness_offset +
          ordinal * ACC96_WINDOW_PROOF_LIMBS * ACC96_N * 4U;
      poly = acc96_proof_witness_poly_at(witness, slot);
      if (poly == NULL) {
        encode_ok = 0;
        continue;
      }
    for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
      int32_t digits[ACC96_WINDOW_PROOF_LIMBS];
      if (!acc96_window_proof_limb_decompose(
              digits, &poly->coeff[coefficient], ctx)) {
        encode_ok = 0;
        continue;
      }
      const size_t phase = coefficient % ACC96_LABRADOR_PHASES;
      const size_t phase_coefficient = coefficient / ACC96_LABRADOR_PHASES;
      const size_t phased =
          phase * ACC96_LABRADOR_BASE_DEGREE + phase_coefficient;
      for (size_t limb = 0; limb < ACC96_WINDOW_PROOF_LIMBS; limb++)
        put_i32(out + slot_offset +
                (limb * ACC96_N + phased) * 4U, digits[limb]);
    }
    }
    if (!encode_ok) return 0;
    offset += ACC96_LABRADOR_WITNESS_LIMB_BYTES;
  }
  {
    const size_t lift_offset = offset;
    int encode_ok = 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(&:encode_ok)
#endif
    for (size_t lane = 0; lane < ACC96_LIMBS; lane++) {
    const size_t lane_offset = lift_offset +
        lane * ACC96_LIFT_PROOF_LIMBS * ACC96_N * 4U;
    for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
      int32_t digits[ACC96_LIFT_PROOF_LIMBS];
      if (!acc96_lift_proof_limb_decompose(
              digits, lift->coeff[lane][coefficient])) {
        encode_ok = 0;
        continue;
      }
      const size_t phase = coefficient % ACC96_LABRADOR_PHASES;
      const size_t phase_coefficient = coefficient / ACC96_LABRADOR_PHASES;
      const size_t phased =
          phase * ACC96_LABRADOR_BASE_DEGREE + phase_coefficient;
      for (size_t limb = 0; limb < ACC96_LIFT_PROOF_LIMBS; limb++)
        put_i32(out + lane_offset +
                (limb * ACC96_N + phased) * 4U, digits[limb]);
    }
    }
    if (!encode_ok) return 0;
    offset += ACC96_LABRADOR_LIFT_LIMB_BYTES;
  }
  return offset == out_len;
}
