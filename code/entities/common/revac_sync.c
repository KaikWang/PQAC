#include "revac_sync.h"

#include <stdlib.h>
#include <string.h>

#include "fips202.h"
#include "randombytes.h"

#define REVAC_SYNC_MAGIC "RVSYN1"
#define REVAC_SYNC_MAGIC_BYTES 6
#define REVAC_SYNC_DOMAIN "REVAC-MEMWITSYNC-V1"
#define REVAC_Q96_SYNC_MAGIC "RVSY9601"
#define REVAC_Q96_SYNC_MAGIC_BYTES 8
#define REVAC_Q96_SYNC_DOMAIN "REVAC-Q96-MEMWITSYNC-V1"
#define REVAC_Q96_SYNC_WIRE_VERSION 1

static void store_u32(uint8_t out[4], uint32_t value) {
  out[0] = (uint8_t)value;
  out[1] = (uint8_t)(value >> 8);
  out[2] = (uint8_t)(value >> 16);
  out[3] = (uint8_t)(value >> 24);
}

static uint32_t load_u32(const uint8_t in[4]) {
  return ((uint32_t)in[0]) |
         ((uint32_t)in[1] << 8) |
         ((uint32_t)in[2] << 16) |
         ((uint32_t)in[3] << 24);
}

static void store_i32(uint8_t out[4], int32_t value) {
  store_u32(out, (uint32_t)value);
}

static int32_t load_i32(const uint8_t in[4]) {
  return (int32_t)load_u32(in);
}

static void absorb_u32(keccak_state *state, uint32_t value) {
  uint8_t bytes[4];
  store_u32(bytes, value);
  shake256_absorb(state, bytes, sizeof(bytes));
}

static void write_bits(uint8_t *out, size_t *bitpos,
                       uint64_t value, unsigned int nbits) {
  for (unsigned int i = 0; i < nbits; i++) {
    size_t pos = *bitpos + i;
    unsigned int shift = nbits - 1 - i;
    out[pos / 8] |= (uint8_t)(((value >> shift) & 1U) << (7 - pos % 8));
  }
  *bitpos += nbits;
}

static int read_bits(uint64_t *value, const uint8_t *in, size_t in_len,
                     size_t *bitpos, unsigned int nbits) {
  uint64_t result = 0;

  if (*bitpos > in_len * 8 || nbits > in_len * 8 - *bitpos) {
    return 0;
  }
  for (unsigned int i = 0; i < nbits; i++) {
    size_t pos = *bitpos + i;
    result = (result << 1) |
             ((uint64_t)(in[pos / 8] >> (7 - pos % 8)) & 1U);
  }
  *bitpos += nbits;
  *value = result;
  return 1;
}

static void pack_poly_q(uint8_t *out, size_t *bitpos, const poly_q poly) {
  for (size_t i = 0; i < PARAM_N; i++) {
    write_bits(out, bitpos, (uint64_t)poly_q_get_coeff(poly, i),
               PARAM_Q_BITLEN);
  }
}

static int unpack_poly_q(poly_q poly, const uint8_t *in, size_t in_len,
                         size_t *bitpos) {
  uint64_t value;

  for (size_t i = 0; i < PARAM_N; i++) {
    if (!read_bits(&value, in, in_len, bitpos, PARAM_Q_BITLEN) ||
        value >= PARAM_Q) {
      return 0;
    }
    poly_q_set_coeff(poly, i, (coeff_q)value);
  }
  return 1;
}

static void pack_vec_q(uint8_t *out, size_t *bitpos,
                       const poly_q_vec_d vec) {
  for (size_t i = 0; i < PARAM_D; i++) {
    pack_poly_q(out, bitpos, vec->entries[i]);
  }
}

static int unpack_vec_q(poly_q_vec_d vec, const uint8_t *in, size_t in_len,
                        size_t *bitpos) {
  for (size_t i = 0; i < PARAM_D; i++) {
    if (!unpack_poly_q(vec->entries[i], in, in_len, bitpos)) {
      return 0;
    }
  }
  return 1;
}

static int pack_witness(uint8_t out[REVAC_SYNC_WITNESS_BYTES],
                        const acc_genpp_wit_t *witness) {
  size_t bitpos = 0;

  memset(out, 0, REVAC_SYNC_WITNESS_BYTES);
  pack_vec_q(out, &bitpos, witness->left);
  pack_vec_q(out, &bitpos, witness->v12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    pack_vec_q(out, &bitpos, witness->v2[i]);
  }
  for (size_t i = 0; i < PARAM_N; i++) {
    coeff_q value = poly_q_get_coeff(witness->tag, i);
    if (value != 0 && value != 1) {
      return 0;
    }
    write_bits(out, &bitpos, (uint64_t)value, 1);
  }
  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      pack_vec_q(out, &bitpos, witness->right[bit][j]);
    }
  }
  if (ACC_GENPP_RANDOMIZED_WITNESS_REFRESH) {
    for (size_t bit = 0; bit < ACC_GENPP_REFRESH_BITS_MAX; bit++) {
      write_bits(out, &bitpos,
                 (witness->refresh_r[bit >> 3] >> (bit & 7u)) & 1u, 1);
    }
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      pack_vec_q(out, &bitpos, witness->refresh_z[j]);
    }
  }
  return bitpos == REVAC_SYNC_WITNESS_BYTES * 8;
}

static int unpack_witness(acc_genpp_wit_t *witness,
                          const uint8_t in[REVAC_SYNC_WITNESS_BYTES]) {
  size_t bitpos = 0;
  uint64_t value;

  if (!unpack_vec_q(witness->left, in, REVAC_SYNC_WITNESS_BYTES, &bitpos) ||
      !unpack_vec_q(witness->v12, in, REVAC_SYNC_WITNESS_BYTES, &bitpos)) {
    return 0;
  }
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    if (!unpack_vec_q(witness->v2[i], in, REVAC_SYNC_WITNESS_BYTES,
                      &bitpos)) {
      return 0;
    }
  }
  poly_q_zero(witness->tag);
  for (size_t i = 0; i < PARAM_N; i++) {
    if (!read_bits(&value, in, REVAC_SYNC_WITNESS_BYTES, &bitpos, 1)) {
      return 0;
    }
    if (value != 0) {
      poly_q_set_coeff(witness->tag, i, 1);
    }
  }
  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      if (!unpack_vec_q(witness->right[bit][j], in,
                        REVAC_SYNC_WITNESS_BYTES, &bitpos)) {
        return 0;
      }
    }
  }
  memset(witness->refresh_r, 0, sizeof(witness->refresh_r));
  if (ACC_GENPP_RANDOMIZED_WITNESS_REFRESH) {
    for (size_t bit = 0; bit < ACC_GENPP_REFRESH_BITS_MAX; bit++) {
      if (!read_bits(&value, in, REVAC_SYNC_WITNESS_BYTES, &bitpos, 1)) {
        return 0;
      }
      witness->refresh_r[bit >> 3] |= (uint8_t)(value << (bit & 7u));
    }
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      if (!unpack_vec_q(witness->refresh_z[j], in,
                        REVAC_SYNC_WITNESS_BYTES, &bitpos)) {
        return 0;
      }
    }
  } else {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_vec_d_zero(witness->refresh_z[j]);
    }
  }
  return bitpos == REVAC_SYNC_WITNESS_BYTES * 8;
}

static void pack_signature(uint8_t out[REVAC_CONTEXT_SIGNATURE_BYTES],
                           const sep_sig_t *signature) {
  size_t bitpos = 0;

  memset(out, 0, REVAC_CONTEXT_SIGNATURE_BYTES);
  pack_poly_q(out, &bitpos, signature->tag);
  pack_vec_q(out, &bitpos, signature->v12);
  for (size_t i = 0; i < PARAM_KH; i++) {
    pack_vec_q(out, &bitpos, signature->v2[i]);
  }
  pack_poly_q(out, &bitpos, signature->v3);
}

static int unpack_signature(
    sep_sig_t *signature,
    const uint8_t in[REVAC_CONTEXT_SIGNATURE_BYTES]) {
  size_t bitpos = 0;

  if (!unpack_poly_q(signature->tag, in, REVAC_CONTEXT_SIGNATURE_BYTES,
                     &bitpos) ||
      !unpack_vec_q(signature->v12, in, REVAC_CONTEXT_SIGNATURE_BYTES,
                    &bitpos)) {
    return 0;
  }
  for (size_t i = 0; i < PARAM_KH; i++) {
    if (!unpack_vec_q(signature->v2[i], in,
                      REVAC_CONTEXT_SIGNATURE_BYTES, &bitpos)) {
      return 0;
    }
  }
  if (!unpack_poly_q(signature->v3, in, REVAC_CONTEXT_SIGNATURE_BYTES,
                     &bitpos)) {
    return 0;
  }
  return bitpos == REVAC_CONTEXT_SIGNATURE_BYTES * 8;
}

static void accumulator_hash(uint8_t out[REVAC_CONTEXT_HASH_BYTES],
                             const poly_q_vec_d acc) {
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, (const uint8_t *)"REVAC-ACC-STATE-V1", 18);
  for (size_t row = 0; row < PARAM_D; row++) {
    for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
      absorb_u32(&state,
                 (uint32_t)poly_q_get_coeff(acc->entries[row], coeff));
    }
  }
  shake256_finalize(&state);
  shake256_squeeze(out, REVAC_CONTEXT_HASH_BYTES, &state);
}

static void sync_message(
    uint8_t out[PARAM_M * PARAM_N / 8],
    uint32_t seq,
    uint32_t handle,
    const uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t nonce[REVAC_SYNC_NONCE_BYTES],
    const uint8_t acc_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t witness_hash[REVAC_CONTEXT_HASH_BYTES]) {
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, (const uint8_t *)REVAC_SYNC_DOMAIN,
                  sizeof(REVAC_SYNC_DOMAIN) - 1);
  absorb_u32(&state, REVAC_SYNC_WIRE_VERSION);
  absorb_u32(&state, ACC_GENPP_GADGET_BASE);
  absorb_u32(&state, ACC_GENPP_GADGET_DIM);
  absorb_u32(&state, ACC_GENPP_HANDLE_BITS);
  absorb_u32(&state, PARAM_Q);
  absorb_u32(&state, PARAM_Q_BITLEN);
  absorb_u32(&state, seq);
  absorb_u32(&state, handle);
  absorb_u32(&state, REVAC_SYNC_WITNESS_BYTES);
  absorb_u32(&state, REVAC_CONTEXT_SIGNATURE_BYTES);
  shake256_absorb(&state, pp_hash, REVAC_CONTEXT_HASH_BYTES);
  shake256_absorb(&state, nonce, REVAC_SYNC_NONCE_BYTES);
  shake256_absorb(&state, acc_hash, REVAC_CONTEXT_HASH_BYTES);
  shake256_absorb(&state, witness_hash, REVAC_CONTEXT_HASH_BYTES);
  shake256_finalize(&state);
  shake256_squeeze(out, PARAM_M * PARAM_N / 8, &state);
}

int revac_sync_request_create(revac_sync_request_t *request,
                              uint32_t handle,
                              uint32_t last_seq,
                              const revac_public_state_t *public_state) {
  if (request == NULL || public_state == NULL ||
      public_state->ta_pk == NULL || last_seq > public_state->acc_seq ||
      (public_state->acc_profile == REVAC_ACC_PROFILE_Q19 &&
       (public_state->acc_pp == NULL || public_state->acc == NULL)) ||
      (public_state->acc_profile == REVAC_ACC_PROFILE_Q96 &&
       (public_state->acc96_ctx == NULL || public_state->acc96_pp == NULL ||
        public_state->acc96_acc == NULL))) {
    return 0;
  }
  request->handle = handle;
  request->last_seq = last_seq;
  memcpy(request->pp_hash, public_state->context_pp_hash,
         REVAC_CONTEXT_HASH_BYTES);
  randombytes(request->nonce, sizeof(request->nonce));
  return 1;
}

static int q96_pack_poly_short(uint8_t *out, const acc96_poly_t *poly,
                               const acc96_ctx_t *ctx) {
  for (size_t i = 0; i < ACC96_N; i++) {
    const acc96_i128 centered = acc96_coeff_to_centered(&poly->coeff[i], ctx);
    if (centered < INT32_MIN || centered > INT32_MAX) return 0;
    store_i32(out + 4U * i, (int32_t)centered);
  }
  return 1;
}

static void q96_unpack_poly_short(acc96_poly_t *poly, const uint8_t *in) {
  const acc96_u128 modulus = acc96_modulus();
  for (size_t i = 0; i < ACC96_N; i++) {
    const int32_t encoded = load_i32(in + 4U * i);
    const acc96_u128 canonical = encoded < 0
        ? modulus - (acc96_u128)(-(int64_t)encoded)
        : (acc96_u128)encoded;
    acc96_coeff_from_u128(&poly->coeff[i], canonical);
  }
}

static int q96_pack_witness(uint8_t *out,
                            const acc96_relation_witness_t *witness,
                            const acc96_ctx_t *ctx) {
  size_t offset = 0;
  if (out == NULL || witness == NULL || ctx == NULL) return 0;
#define PACK_Q96_POLY(poly_ptr) do {                                  \
    if (!q96_pack_poly_short(out + offset, (poly_ptr), ctx)) return 0;\
    offset += ACC96_N * 4U;                                          \
  } while (0)
  PACK_Q96_POLY(&witness->left);
  PACK_Q96_POLY(&witness->v12);
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++)
    PACK_Q96_POLY(&witness->v2[i]);
  PACK_Q96_POLY(&witness->tag);
  for (size_t bit = 0; bit < ACC96_HANDLE_BITS; bit++)
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++)
      PACK_Q96_POLY(&witness->right.right[bit][digit]);
  memcpy(out + offset, witness->refresh_r, ACC96_REFRESH_BYTES);
  offset += ACC96_REFRESH_BYTES;
  for (size_t block = 0; block < ACC96_REFRESH_BLOCKS; block++)
    PACK_Q96_POLY(&witness->refresh_z[block]);
#undef PACK_Q96_POLY
  return offset == REVAC_Q96_SYNC_WITNESS_BYTES;
}

static int q96_unpack_witness(acc96_relation_witness_t *witness,
                              const uint8_t *in) {
  size_t offset = 0;
  const uint8_t *encoded = in;
  if (witness == NULL || in == NULL) return 0;
#define UNPACK_Q96_POLY(poly_ptr) do {                                \
    q96_unpack_poly_short((poly_ptr), encoded + offset);              \
    offset += ACC96_N * 4U;                                          \
  } while (0)
  UNPACK_Q96_POLY(&witness->left);
  UNPACK_Q96_POLY(&witness->v12);
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++)
    UNPACK_Q96_POLY(&witness->v2[i]);
  UNPACK_Q96_POLY(&witness->tag);
  for (size_t bit = 0; bit < ACC96_HANDLE_BITS; bit++)
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++)
      UNPACK_Q96_POLY(&witness->right.right[bit][digit]);
  memcpy(witness->refresh_r, encoded + offset, ACC96_REFRESH_BYTES);
  offset += ACC96_REFRESH_BYTES;
  for (size_t block = 0; block < ACC96_REFRESH_BLOCKS; block++)
    UNPACK_Q96_POLY(&witness->refresh_z[block]);
#undef UNPACK_Q96_POLY
  return offset == REVAC_Q96_SYNC_WITNESS_BYTES;
}

static void q96_accumulator_hash(
    uint8_t out[REVAC_CONTEXT_HASH_BYTES], const acc96_poly_t *accumulator) {
  uint8_t encoded[ACC96_POLY_PACKED_BYTES];
  (void)acc96_poly_pack(encoded, accumulator);
  sha3_256(out, encoded, sizeof(encoded));
}

static void q96_sync_message(
    uint8_t out[PARAM_M * PARAM_N / 8], uint32_t seq, uint32_t handle,
    const uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t nonce[REVAC_SYNC_NONCE_BYTES],
    const uint8_t acc_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t witness_hash[REVAC_CONTEXT_HASH_BYTES]) {
  keccak_state state;
  shake256_init(&state);
  shake256_absorb(&state, (const uint8_t *)REVAC_Q96_SYNC_DOMAIN,
                  sizeof(REVAC_Q96_SYNC_DOMAIN) - 1);
  absorb_u32(&state, REVAC_Q96_SYNC_WIRE_VERSION);
  absorb_u32(&state, REVAC_ACC_PROFILE_Q96);
  absorb_u32(&state, ACC96_N);
  absorb_u32(&state, ACC96_MODULUS_BITS);
  absorb_u32(&state, ACC96_GADGET_BASE);
  absorb_u32(&state, ACC96_GADGET_DIGITS);
  absorb_u32(&state, seq);
  absorb_u32(&state, handle);
  absorb_u32(&state, (uint32_t)REVAC_Q96_SYNC_WITNESS_BYTES);
  absorb_u32(&state, REVAC_CONTEXT_SIGNATURE_BYTES);
  shake256_absorb(&state, pp_hash, REVAC_CONTEXT_HASH_BYTES);
  shake256_absorb(&state, nonce, REVAC_SYNC_NONCE_BYTES);
  shake256_absorb(&state, acc_hash, REVAC_CONTEXT_HASH_BYTES);
  shake256_absorb(&state, witness_hash, REVAC_CONTEXT_HASH_BYTES);
  shake256_finalize(&state);
  shake256_squeeze(out, PARAM_M * PARAM_N / 8, &state);
}

int revac_q96_sync_response_pack(
    uint8_t **wire, size_t *wire_len,
    const revac_sync_request_t *request, uint32_t acc_seq,
    const acc96_poly_t *accumulator,
    const acc96_relation_witness_t *witness,
    const acc96_ctx_t *ctx,
    const sep_sk_t *ta_sk, const sep_pk_t *ta_pk) {
  uint8_t message[PARAM_M * PARAM_N / 8];
  uint8_t acc_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t witness_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t state[STATE_BYTES];
  sep_sig_t signature;
  const size_t total = REVAC_Q96_SYNC_WIRE_HEADER_BYTES +
                       REVAC_CONTEXT_SIGNATURE_BYTES +
                       REVAC_Q96_SYNC_WITNESS_BYTES;
  uint8_t *out;
  uint8_t *witness_bytes;

  if (wire == NULL || wire_len == NULL) return 0;
  *wire = NULL;
  *wire_len = 0;
  if (request == NULL || accumulator == NULL || witness == NULL || ctx == NULL ||
      ta_sk == NULL || ta_pk == NULL || acc_seq < request->last_seq ||
      REVAC_Q96_SYNC_WITNESS_BYTES > UINT32_MAX ||
      !acc96_relation_fresh_witness_check(witness, ctx)) return 0;
  out = calloc(total, 1);
  if (out == NULL) return 0;
  witness_bytes = out + REVAC_Q96_SYNC_WIRE_HEADER_BYTES +
                  REVAC_CONTEXT_SIGNATURE_BYTES;
  if (!q96_pack_witness(witness_bytes, witness, ctx)) {
    free(out);
    return 0;
  }
  q96_accumulator_hash(acc_hash, accumulator);
  sha3_256(witness_hash, witness_bytes, REVAC_Q96_SYNC_WITNESS_BYTES);
  q96_sync_message(message, acc_seq, request->handle, request->pp_hash,
                   request->nonce, acc_hash, witness_hash);
  sep_sig_init(&signature);
  randombytes(state, sizeof(state));
  sep_sign(&signature, state, ta_sk, ta_pk, message);

  memcpy(out, REVAC_Q96_SYNC_MAGIC, REVAC_Q96_SYNC_MAGIC_BYTES);
  out[8] = REVAC_Q96_SYNC_WIRE_VERSION;
  out[9] = REVAC_ACC_PROFILE_Q96;
  store_u32(out + 12, acc_seq);
  store_u32(out + 16, request->handle);
  store_u32(out + 20, (uint32_t)REVAC_Q96_SYNC_WITNESS_BYTES);
  store_u32(out + 24, (uint32_t)REVAC_CONTEXT_SIGNATURE_BYTES);
  memcpy(out + 28, request->pp_hash, REVAC_CONTEXT_HASH_BYTES);
  memcpy(out + 60, request->nonce, REVAC_SYNC_NONCE_BYTES);
  memcpy(out + 92, acc_hash, REVAC_CONTEXT_HASH_BYTES);
  memcpy(out + 124, witness_hash, REVAC_CONTEXT_HASH_BYTES);
  pack_signature(out + REVAC_Q96_SYNC_WIRE_HEADER_BYTES, &signature);
  sep_sig_clear(&signature);
  *wire = out;
  *wire_len = total;
  return 1;
}

int revac_q96_sync_response_unpack(
    acc96_relation_witness_t *witness, uint32_t *response_seq,
    const uint8_t *wire, size_t wire_len,
    const revac_sync_request_t *request,
    const revac_public_state_t *public_state) {
  uint8_t message[PARAM_M * PARAM_N / 8];
  uint8_t expected_acc_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t expected_witness_hash[REVAC_CONTEXT_HASH_BYTES];
  const size_t expected_len = REVAC_Q96_SYNC_WIRE_HEADER_BYTES +
                              REVAC_CONTEXT_SIGNATURE_BYTES +
                              REVAC_Q96_SYNC_WITNESS_BYTES;
  const uint8_t *witness_bytes;
  acc96_relation_witness_t *decoded = NULL;
  sep_sig_t signature;
  uint32_t seq;
  int ok = 0;

  if (witness == NULL || response_seq == NULL || wire == NULL ||
      request == NULL || public_state == NULL ||
      public_state->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      public_state->ta_pk == NULL || public_state->acc96_ctx == NULL ||
      public_state->acc96_pp == NULL || public_state->acc96_acc == NULL ||
      wire_len != expected_len ||
      memcmp(wire, REVAC_Q96_SYNC_MAGIC, REVAC_Q96_SYNC_MAGIC_BYTES) != 0 ||
      wire[8] != REVAC_Q96_SYNC_WIRE_VERSION ||
      wire[9] != REVAC_ACC_PROFILE_Q96 || wire[10] != 0 || wire[11] != 0 ||
      load_u32(wire + 16) != request->handle ||
      load_u32(wire + 20) != REVAC_Q96_SYNC_WITNESS_BYTES ||
      load_u32(wire + 24) != REVAC_CONTEXT_SIGNATURE_BYTES ||
      memcmp(wire + 28, request->pp_hash, REVAC_CONTEXT_HASH_BYTES) != 0 ||
      memcmp(wire + 28, public_state->context_pp_hash,
             REVAC_CONTEXT_HASH_BYTES) != 0 ||
      memcmp(wire + 60, request->nonce, REVAC_SYNC_NONCE_BYTES) != 0)
    return 0;
  seq = load_u32(wire + 12);
  if (seq != public_state->acc96_seq || seq < request->last_seq) return 0;
  q96_accumulator_hash(expected_acc_hash, public_state->acc96_acc);
  if (memcmp(wire + 92, expected_acc_hash, REVAC_CONTEXT_HASH_BYTES) != 0)
    return 0;
  witness_bytes = wire + REVAC_Q96_SYNC_WIRE_HEADER_BYTES +
                  REVAC_CONTEXT_SIGNATURE_BYTES;
  sha3_256(expected_witness_hash, witness_bytes,
           REVAC_Q96_SYNC_WITNESS_BYTES);
  if (memcmp(wire + 124, expected_witness_hash,
             REVAC_CONTEXT_HASH_BYTES) != 0) return 0;
  q96_sync_message(message, seq, request->handle, wire + 28, wire + 60,
                   wire + 92, wire + 124);
  sep_sig_init(&signature);
  if (!unpack_signature(&signature,
                        wire + REVAC_Q96_SYNC_WIRE_HEADER_BYTES)) {
    sep_sig_clear(&signature);
    return 0;
  }
  ok = sep_verify(&signature, message, public_state->ta_pk);
  sep_sig_clear(&signature);
  if (!ok) return 0;
  decoded = calloc(1, sizeof(*decoded));
  if (decoded == NULL) return 0;
  ok = q96_unpack_witness(decoded, witness_bytes) &&
       acc96_relation_fresh_witness_check(decoded, public_state->acc96_ctx) &&
       acc96_relation_memver(public_state->acc96_pp,
                             public_state->acc96_acc, request->handle,
                             decoded,
                             public_state->acc96_witness_norm_sq_bound,
                             public_state->acc96_ctx);
  if (ok) {
    memcpy(witness, decoded, sizeof(*witness));
    *response_seq = seq;
  }
  free(decoded);
  return ok;
}

int revac_sync_response_pack(uint8_t **wire,
                             size_t *wire_len,
                             const revac_sync_request_t *request,
                             uint32_t acc_seq,
                             const poly_q_vec_d acc,
                             const acc_genpp_wit_t *witness,
                             const sep_sk_t *ta_sk,
                             const sep_pk_t *ta_pk) {
  uint8_t message[PARAM_M * PARAM_N / 8];
  uint8_t acc_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t witness_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t state[STATE_BYTES];
  sep_sig_t signature;
  size_t total;
  uint8_t *out;
  uint8_t *witness_bytes;

  if (wire == NULL || wire_len == NULL) {
    return 0;
  }
  *wire = NULL;
  *wire_len = 0;
  if (request == NULL || acc == NULL || witness == NULL || ta_sk == NULL ||
      ta_pk == NULL || acc_seq < request->last_seq ||
      REVAC_SYNC_WITNESS_BYTES != acc_genpp_wit_full_size_bytes()) {
    return 0;
  }
  total = REVAC_SYNC_WIRE_HEADER_BYTES + REVAC_CONTEXT_SIGNATURE_BYTES +
          REVAC_SYNC_WITNESS_BYTES;
  out = calloc(total, 1);
  if (out == NULL) {
    return 0;
  }
  witness_bytes = out + REVAC_SYNC_WIRE_HEADER_BYTES +
                  REVAC_CONTEXT_SIGNATURE_BYTES;
  if (!pack_witness(witness_bytes, witness)) {
    free(out);
    return 0;
  }
  accumulator_hash(acc_hash, acc);
  sha3_256(witness_hash, witness_bytes, REVAC_SYNC_WITNESS_BYTES);
  sync_message(message, acc_seq, request->handle, request->pp_hash,
               request->nonce, acc_hash, witness_hash);
  sep_sig_init(&signature);
  randombytes(state, sizeof(state));
  sep_sign(&signature, state, ta_sk, ta_pk, message);

  memcpy(out, REVAC_SYNC_MAGIC, REVAC_SYNC_MAGIC_BYTES);
  out[6] = REVAC_SYNC_WIRE_VERSION;
  out[7] = ACC_GENPP_GADGET_BASE;
  out[8] = ACC_GENPP_GADGET_DIM;
  out[9] = ACC_GENPP_HANDLE_BITS;
  out[10] = PARAM_Q_BITLEN;
  store_u32(out + 12, acc_seq);
  store_u32(out + 16, request->handle);
  store_u32(out + 20, (uint32_t)REVAC_SYNC_WITNESS_BYTES);
  store_u32(out + 24, (uint32_t)REVAC_CONTEXT_SIGNATURE_BYTES);
  memcpy(out + 28, request->pp_hash, REVAC_CONTEXT_HASH_BYTES);
  memcpy(out + 60, request->nonce, REVAC_SYNC_NONCE_BYTES);
  memcpy(out + 92, acc_hash, REVAC_CONTEXT_HASH_BYTES);
  memcpy(out + 124, witness_hash, REVAC_CONTEXT_HASH_BYTES);
  pack_signature(out + REVAC_SYNC_WIRE_HEADER_BYTES, &signature);
  sep_sig_clear(&signature);
  *wire = out;
  *wire_len = total;
  return 1;
}

int revac_sync_response_unpack(acc_genpp_wit_t *witness,
                               uint32_t *response_seq,
                               const uint8_t *wire,
                               size_t wire_len,
                               const revac_sync_request_t *request,
                               const revac_public_state_t *public_state) {
  uint8_t message[PARAM_M * PARAM_N / 8];
  uint8_t expected_acc_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t expected_witness_hash[REVAC_CONTEXT_HASH_BYTES];
  sep_sig_t signature;
  acc_genpp_wit_t decoded;
  uint32_t seq;
  size_t expected_len = REVAC_SYNC_WIRE_HEADER_BYTES +
                        REVAC_CONTEXT_SIGNATURE_BYTES +
                        REVAC_SYNC_WITNESS_BYTES;
  const uint8_t *witness_bytes;
  int signature_ok;

  if (witness == NULL || response_seq == NULL || wire == NULL ||
      request == NULL || public_state == NULL ||
      public_state->ta_pk == NULL || public_state->acc_pp == NULL ||
      public_state->acc == NULL || wire_len != expected_len ||
      memcmp(wire, REVAC_SYNC_MAGIC, REVAC_SYNC_MAGIC_BYTES) != 0 ||
      wire[6] != REVAC_SYNC_WIRE_VERSION ||
      wire[7] != ACC_GENPP_GADGET_BASE ||
      wire[8] != ACC_GENPP_GADGET_DIM ||
      wire[9] != ACC_GENPP_HANDLE_BITS || wire[10] != PARAM_Q_BITLEN ||
      wire[11] != 0 ||
      load_u32(wire + 16) != request->handle ||
      load_u32(wire + 20) != REVAC_SYNC_WITNESS_BYTES ||
      load_u32(wire + 24) != REVAC_CONTEXT_SIGNATURE_BYTES ||
      memcmp(wire + 28, request->pp_hash, REVAC_CONTEXT_HASH_BYTES) != 0 ||
      memcmp(wire + 28, public_state->context_pp_hash,
             REVAC_CONTEXT_HASH_BYTES) != 0 ||
      memcmp(wire + 60, request->nonce, REVAC_SYNC_NONCE_BYTES) != 0) {
    return 0;
  }
  seq = load_u32(wire + 12);
  if (seq != public_state->acc_seq || seq < request->last_seq) {
    return 0;
  }
  accumulator_hash(expected_acc_hash, public_state->acc);
  if (memcmp(wire + 92, expected_acc_hash, REVAC_CONTEXT_HASH_BYTES) != 0) {
    return 0;
  }
  witness_bytes = wire + REVAC_SYNC_WIRE_HEADER_BYTES +
                  REVAC_CONTEXT_SIGNATURE_BYTES;
  sha3_256(expected_witness_hash, witness_bytes, REVAC_SYNC_WITNESS_BYTES);
  if (memcmp(wire + 124, expected_witness_hash,
             REVAC_CONTEXT_HASH_BYTES) != 0) {
    return 0;
  }
  sync_message(message, seq, request->handle, wire + 28, wire + 60,
               wire + 92, wire + 124);
  sep_sig_init(&signature);
  if (!unpack_signature(&signature,
                        wire + REVAC_SYNC_WIRE_HEADER_BYTES)) {
    sep_sig_clear(&signature);
    return 0;
  }
  signature_ok = sep_verify(&signature, message, public_state->ta_pk);
  sep_sig_clear(&signature);
  if (!signature_ok) {
    return 0;
  }

  /* Decode and verify into scratch storage so every rejected response is
     non-mutating even for callers below the UAV wrapper layer. */
  acc_genpp_wit_init(&decoded);
  if (!unpack_witness(&decoded, witness_bytes) ||
      !acc_genpp_memver(public_state->acc_pp, public_state->acc,
                        request->handle, &decoded)) {
    acc_genpp_wit_clear(&decoded);
    return 0;
  }
  revac_sync_witness_set(witness, &decoded);
  acc_genpp_wit_clear(&decoded);
  *response_seq = seq;
  return 1;
}

void revac_sync_witness_set(acc_genpp_wit_t *dst,
                            const acc_genpp_wit_t *src) {
  poly_q_vec_d_set(dst->left, src->left);
  poly_q_vec_d_set(dst->v12, src->v12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_vec_d_set(dst->v2[i], src->v2[i]);
  }
  poly_q_set(dst->tag, src->tag);
  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_vec_d_set(dst->right[bit][j], src->right[bit][j]);
    }
  }
  memcpy(dst->refresh_r, src->refresh_r, sizeof(dst->refresh_r));
  for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
    poly_q_vec_d_set(dst->refresh_z[j], src->refresh_z[j]);
  }
}
