#include "revac_context.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "fips202.h"
#include "randombytes.h"

#define REVAC_CONTEXT_MAGIC "RVDLT4"
#define REVAC_CONTEXT_MAGIC_BYTES 6
#define REVAC_CONTEXT_DOMAIN "REVAC-DELTA4-CONTEXT-V3"
#define REVAC_Q96_UPDATE_MAGIC "RVQ96U"
#define REVAC_Q96_UPDATE_MAGIC_BYTES 6
#define REVAC_Q96_UPDATE_DOMAIN "REVAC-Q96-SIGNED-TAIL4-V5"
#define REVAC_Q96_PREFETCH_MAGIC "RVQ96P"
#define REVAC_Q96_PREFETCH_MAGIC_BYTES 6
#define REVAC_Q96_PREFETCH_DOMAIN "REVAC-Q96-TAIL4-PREFETCH-V1"
#define REVAC_Q96_DELTA6_MANIFEST_MAGIC "RVQ88M"
#define REVAC_Q96_DELTA6_DOMAIN "REVAC-Q88-SIGNED-DELTA6-MANIFEST-V1"

static uint8_t encode_update_sign(int sign) {
  return sign > 0 ? UINT8_C(1) : sign < 0 ? UINT8_C(2) : UINT8_C(0);
}

static int decode_update_sign(uint8_t encoded) {
  return encoded == UINT8_C(1) ? 1 : encoded == UINT8_C(2) ? -1 : 0;
}

static void store_u32(uint8_t out[4], uint32_t x) {
  out[0] = (uint8_t)x;
  out[1] = (uint8_t)(x >> 8);
  out[2] = (uint8_t)(x >> 16);
  out[3] = (uint8_t)(x >> 24);
}

static uint32_t load_u32(const uint8_t in[4]) {
  return ((uint32_t)in[0]) | ((uint32_t)in[1] << 8) |
         ((uint32_t)in[2] << 16) | ((uint32_t)in[3] << 24);
}

static void absorb_u32(keccak_state *state, uint32_t x) {
  uint8_t buf[4];
  store_u32(buf, x);
  shake256_absorb(state, buf, sizeof(buf));
}

static void write_bits(uint8_t *out, size_t *bitpos,
                       uint64_t value, unsigned int nbits) {
  for (unsigned int i = 0; i < nbits; i++) {
    const size_t pos = *bitpos + i;
    const unsigned int shift = nbits - 1U - i;
    out[pos / 8] |=
        (uint8_t)(((value >> shift) & UINT64_C(1)) << (7U - pos % 8));
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
    const size_t pos = *bitpos + i;
    result = (result << 1) |
             ((uint64_t)(in[pos / 8] >> (7U - pos % 8)) & UINT64_C(1));
  }
  *bitpos += nbits;
  *value = result;
  return 1;
}

static void pack_poly(uint8_t *out, size_t *bitpos, const poly_q poly) {
  for (size_t i = 0; i < PARAM_N; i++) {
    write_bits(out, bitpos, (uint64_t)poly_q_get_coeff(poly, i),
               PARAM_Q_BITLEN);
  }
}

static int unpack_poly(poly_q poly, const uint8_t *in, size_t in_len,
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

static void pack_signature(uint8_t out[REVAC_CONTEXT_SIGNATURE_BYTES],
                           const sep_sig_t *signature) {
  size_t bitpos = 0;

  memset(out, 0, REVAC_CONTEXT_SIGNATURE_BYTES);
  pack_poly(out, &bitpos, signature->tag);
  for (size_t i = 0; i < PARAM_D; i++) {
    pack_poly(out, &bitpos, signature->v12->entries[i]);
  }
  for (size_t k = 0; k < PARAM_KH; k++) {
    for (size_t i = 0; i < PARAM_D; i++) {
      pack_poly(out, &bitpos, signature->v2[k]->entries[i]);
    }
  }
  pack_poly(out, &bitpos, signature->v3);
}

static int unpack_signature(sep_sig_t *signature, const uint8_t *in) {
  size_t bitpos = 0;

  if (!unpack_poly(signature->tag, in, REVAC_CONTEXT_SIGNATURE_BYTES,
                   &bitpos)) {
    return 0;
  }
  for (size_t i = 0; i < PARAM_D; i++) {
    if (!unpack_poly(signature->v12->entries[i], in,
                     REVAC_CONTEXT_SIGNATURE_BYTES, &bitpos)) {
      return 0;
    }
  }
  for (size_t k = 0; k < PARAM_KH; k++) {
    for (size_t i = 0; i < PARAM_D; i++) {
      if (!unpack_poly(signature->v2[k]->entries[i], in,
                       REVAC_CONTEXT_SIGNATURE_BYTES, &bitpos)) {
        return 0;
      }
    }
  }
  return unpack_poly(signature->v3, in, REVAC_CONTEXT_SIGNATURE_BYTES,
                     &bitpos) &&
         bitpos == REVAC_CONTEXT_SIGNATURE_BYTES * 8;
}

static uint32_t zigzag_encode(coeff_q value) {
  const int64_t wide = value;
  return value >= 0 ? (uint32_t)(2 * wide)
                    : (uint32_t)(-2 * wide - 1);
}

static coeff_q zigzag_decode(uint32_t value) {
  return (value & 1U) != 0U ? -(coeff_q)((value + 1U) / 2U)
                            : (coeff_q)(value / 2U);
}

static int delta_pack(
    uint8_t **out, size_t *out_len,
    const poly_q_vec_d delta[][ACC_GENPP_GADGET_DIM], size_t layers) {
  size_t escape_count = 0;
  size_t bitpos = 0;
  size_t bytes;
  uint8_t *encoded;

  if (out == NULL || out_len == NULL || delta == NULL || layers == 0 ||
      layers > ACC_GENPP_HANDLE_BITS) {
    return 0;
  }
  *out = NULL;
  *out_len = 0;
  for (size_t layer = 0; layer < layers; layer++) {
    for (size_t gadget = 0; gadget < ACC_GENPP_GADGET_DIM; gadget++) {
      for (size_t row = 0; row < PARAM_D; row++) {
        for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
          const coeff_q value = poly_q_get_coeff_centered(
              delta[layer][gadget]->entries[row], coeff);
          if (zigzag_encode(value) >= REVAC_CONTEXT_DELTA_ESCAPE) {
            escape_count++;
          }
        }
      }
    }
  }
  bytes = ((size_t)layers * ACC_GENPP_GADGET_DIM * PARAM_D * PARAM_N *
               REVAC_CONTEXT_DELTA_FAST_BITS +
           escape_count * PARAM_Q_BITLEN + 7) /
          8;
  encoded = calloc(bytes, 1);
  if (encoded == NULL) {
    return 0;
  }

  for (size_t layer = 0; layer < layers; layer++) {
    for (size_t gadget = 0; gadget < ACC_GENPP_GADGET_DIM; gadget++) {
      for (size_t row = 0; row < PARAM_D; row++) {
        for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
          const coeff_q value = poly_q_get_coeff_centered(
              delta[layer][gadget]->entries[row], coeff);
          const uint32_t zigzag = zigzag_encode(value);

          if (zigzag < REVAC_CONTEXT_DELTA_ESCAPE) {
            write_bits(encoded, &bitpos, zigzag,
                       REVAC_CONTEXT_DELTA_FAST_BITS);
          } else {
            write_bits(encoded, &bitpos, REVAC_CONTEXT_DELTA_ESCAPE,
                       REVAC_CONTEXT_DELTA_FAST_BITS);
            write_bits(encoded, &bitpos,
                       (uint32_t)poly_q_get_coeff(
                           delta[layer][gadget]->entries[row], coeff),
                       PARAM_Q_BITLEN);
          }
        }
      }
    }
  }
  *out = encoded;
  *out_len = bytes;
  return 1;
}

static int delta_unpack(
    poly_q_vec_d delta[][ACC_GENPP_GADGET_DIM], size_t layers,
    const uint8_t *encoded, size_t encoded_len) {
  size_t bitpos = 0;
  uint64_t word;

  if (delta == NULL || layers == 0 || layers > ACC_GENPP_HANDLE_BITS ||
      encoded == NULL || encoded_len == 0) {
    return 0;
  }
  for (size_t layer = 0; layer < layers; layer++) {
    for (size_t gadget = 0; gadget < ACC_GENPP_GADGET_DIM; gadget++) {
      for (size_t row = 0; row < PARAM_D; row++) {
        for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
          coeff_q value;

          if (!read_bits(&word, encoded, encoded_len, &bitpos,
                         REVAC_CONTEXT_DELTA_FAST_BITS)) {
            return 0;
          }
          if (word == REVAC_CONTEXT_DELTA_ESCAPE) {
            if (!read_bits(&word, encoded, encoded_len, &bitpos,
                           PARAM_Q_BITLEN) ||
                word >= PARAM_Q) {
              return 0;
            }
            value = (coeff_q)word;
            {
              const coeff_q centered =
                  value > PARAM_Q / 2 ? value - PARAM_Q : value;
              if (zigzag_encode(centered) < REVAC_CONTEXT_DELTA_ESCAPE) {
                return 0;
              }
            }
          } else {
            value = zigzag_decode((uint32_t)word);
          }
          poly_q_set_coeff(delta[layer][gadget]->entries[row], coeff,
                           value);
        }
      }
    }
  }
  if ((bitpos + 7) / 8 != encoded_len) {
    return 0;
  }
  while (bitpos < encoded_len * 8) {
    if (!read_bits(&word, encoded, encoded_len, &bitpos, 1) || word != 0) {
      return 0;
    }
  }
  return 1;
}

static void manifest_message(
    uint8_t out[PARAM_M * PARAM_N / 8], uint32_t seq, uint32_t handle,
    int sign,
    uint32_t payload_len, const uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t delta_hash[REVAC_CONTEXT_HASH_BYTES]) {
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, (const uint8_t *)REVAC_CONTEXT_DOMAIN,
                  sizeof(REVAC_CONTEXT_DOMAIN) - 1);
  absorb_u32(&state, REVAC_CONTEXT_WIRE_VERSION);
  absorb_u32(&state, REVAC_CONTEXT_PROFILE_DELTA4);
  absorb_u32(&state, REVAC_CONTEXT_DELTA_CODEC_ZIGZAG9_ESCAPE19);
  absorb_u32(&state, ACC_GENPP_GADGET_BASE);
  absorb_u32(&state, ACC_GENPP_GADGET_DIM);
  absorb_u32(&state, ACC_GENPP_HANDLE_BITS);
  absorb_u32(&state, PARAM_Q);
  absorb_u32(&state, PARAM_Q_BITLEN);
  absorb_u32(&state, seq);
  absorb_u32(&state, handle);
  absorb_u32(&state, encode_update_sign(sign));
  absorb_u32(&state, ACC_GENPP_DELTA4_LAYERS);
  absorb_u32(&state, (uint32_t)REVAC_CONTEXT_DELTA_COEFFS);
  absorb_u32(&state, payload_len);
  shake256_absorb(&state, pp_hash, REVAC_CONTEXT_HASH_BYTES);
  shake256_absorb(&state, delta_hash, REVAC_CONTEXT_HASH_BYTES);
  shake256_finalize(&state);
  shake256_squeeze(out, PARAM_M * PARAM_N / 8, &state);
}

static void q96_update_manifest_message(
    uint8_t out[PARAM_M * PARAM_N / 8], uint32_t sequence,
    uint32_t revoked_handle,
    const uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t accumulator_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t delta_hash[REVAC_CONTEXT_HASH_BYTES],
    uint32_t delta_len) {
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, (const uint8_t *)REVAC_Q96_UPDATE_DOMAIN,
                  sizeof(REVAC_Q96_UPDATE_DOMAIN) - 1);
  absorb_u32(&state, REVAC_Q96_UPDATE_WIRE_VERSION);
  absorb_u32(&state, REVAC_ACC_PROFILE_Q96);
  absorb_u32(&state, ACC96_MODULUS_BITS);
  absorb_u32(&state, ACC96_N);
  absorb_u32(&state, ACC96_GADGET_BASE);
  absorb_u32(&state, ACC96_GADGET_DIGITS);
  absorb_u32(&state, ACC96_TAIL4_BIT);
  absorb_u32(&state, REVAC_Q96_UPDATE_CODEC_TAIL4_CRT96);
  absorb_u32(&state, delta_len);
  absorb_u32(&state, sequence);
  absorb_u32(&state, revoked_handle);
  shake256_absorb(&state, pp_hash, REVAC_CONTEXT_HASH_BYTES);
  shake256_absorb(&state, accumulator_hash, REVAC_CONTEXT_HASH_BYTES);
  shake256_absorb(&state, delta_hash, REVAC_CONTEXT_HASH_BYTES);
  shake256_finalize(&state);
  shake256_squeeze(out, PARAM_M * PARAM_N / 8, &state);
}

static void q96_prefetch_message(
    uint8_t out[PARAM_M * PARAM_N / 8], uint32_t candidate_handle,
    const uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t tail_hash[REVAC_CONTEXT_HASH_BYTES]) {
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, (const uint8_t *)REVAC_Q96_PREFETCH_DOMAIN,
                  sizeof(REVAC_Q96_PREFETCH_DOMAIN) - 1);
  absorb_u32(&state, REVAC_Q96_PREFETCH_WIRE_VERSION);
  absorb_u32(&state, REVAC_ACC_PROFILE_Q96);
  absorb_u32(&state, ACC96_MODULUS_BITS);
  absorb_u32(&state, ACC96_N);
  absorb_u32(&state, ACC96_GADGET_BASE);
  absorb_u32(&state, ACC96_GADGET_DIGITS);
  absorb_u32(&state, ACC96_TAIL4_BIT);
  absorb_u32(&state, REVAC_Q96_UPDATE_CODEC_TAIL4_CRT96);
  absorb_u32(&state, REVAC_Q96_TAIL4_BYTES);
  absorb_u32(&state, candidate_handle);
  shake256_absorb(&state, pp_hash, REVAC_CONTEXT_HASH_BYTES);
  shake256_absorb(&state, tail_hash, REVAC_CONTEXT_HASH_BYTES);
  shake256_finalize(&state);
  shake256_squeeze(out, PARAM_M * PARAM_N / 8, &state);
}

static void q96_delta6_manifest_message(
    uint8_t out[PARAM_M * PARAM_N / 8],
    const revac_q96_delta6_entry_t *entry) {
  keccak_state state;
  shake256_init(&state);
  shake256_absorb(&state, (const uint8_t *)REVAC_Q96_DELTA6_DOMAIN,
                  sizeof(REVAC_Q96_DELTA6_DOMAIN) - 1);
  absorb_u32(&state, REVAC_Q96_DELTA6_MANIFEST_WIRE_VERSION);
  absorb_u32(&state, REVAC_ACC_PROFILE_Q96);
  absorb_u32(&state, ACC96_MODULUS_BITS);
  absorb_u32(&state, ACC96_N);
  absorb_u32(&state, ACC96_GADGET_BASE);
  absorb_u32(&state, ACC96_GADGET_DIGITS);
  absorb_u32(&state, REVAC_Q96_DELTA_LAYERS);
  absorb_u32(&state, entry->update.sequence);
  absorb_u32(&state, entry->update.revoked_handle);
  shake256_absorb(&state, entry->pp_hash, REVAC_CONTEXT_HASH_BYTES);
  shake256_absorb(&state, entry->accumulator_hash, REVAC_CONTEXT_HASH_BYTES);
  for (size_t i = 0; i < REVAC_Q96_DELTA_LAYERS; i++) {
    absorb_u32(&state, entry->layer_len[i]);
    shake256_absorb(&state, entry->layer_hash[i], REVAC_CONTEXT_HASH_BYTES);
  }
  shake256_finalize(&state);
  shake256_squeeze(out, PARAM_M * PARAM_N / 8, &state);
}

const char *revac_context_profile_name(revac_context_profile_t profile) {
  switch (profile) {
    case REVAC_CONTEXT_PROFILE_LOCAL:
      return "LOCAL";
    case REVAC_CONTEXT_PROFILE_DELTA4:
      return "DELTA4";
    default:
      return "invalid";
  }
}

size_t revac_context_profile_payload_bytes(revac_context_profile_t profile) {
  return profile == REVAC_CONTEXT_PROFILE_DELTA4
             ? REVAC_CONTEXT_DELTA_FAST_BYTES
             : 0;
}

size_t revac_context_profile_wire_bytes(revac_context_profile_t profile) {
  return profile == REVAC_CONTEXT_PROFILE_DELTA4
             ? REVAC_CONTEXT_WIRE_HEADER_BYTES +
                   REVAC_CONTEXT_SIGNATURE_BYTES +
                   REVAC_CONTEXT_DELTA_FAST_BYTES
             : 0;
}

double revac_context_profile_expected_composites(
    revac_context_profile_t profile) {
  switch (profile) {
    case REVAC_CONTEXT_PROFILE_LOCAL:
      return 30.0;
    case REVAC_CONTEXT_PROFILE_DELTA4:
      return 1.875;
    default:
      return INFINITY;
  }
}

revac_context_profile_t revac_context_profile_select(double goodput_mbps,
                                                      double recurrence_ms) {
  double local_seconds;
  double delta_seconds;

  if (!(goodput_mbps > 0.0) || !isfinite(goodput_mbps)) {
    return REVAC_CONTEXT_PROFILE_LOCAL;
  }
  if (!(recurrence_ms > 0.0) || !isfinite(recurrence_ms)) {
    recurrence_ms = 40.0;
  }
  local_seconds = 30.0 * recurrence_ms / 1000.0;
  delta_seconds =
      (double)revac_context_profile_wire_bytes(REVAC_CONTEXT_PROFILE_DELTA4) *
          8.0 / (goodput_mbps * 1000000.0) +
      revac_context_profile_expected_composites(
          REVAC_CONTEXT_PROFILE_DELTA4) *
          recurrence_ms / 1000.0;
  return delta_seconds < local_seconds ? REVAC_CONTEXT_PROFILE_DELTA4
                                       : REVAC_CONTEXT_PROFILE_LOCAL;
}

int revac_context_packet_peek_profile(const uint8_t *wire, size_t wire_len,
                                      revac_context_profile_t *profile) {
  if (wire == NULL || profile == NULL || wire_len < 8 ||
      memcmp(wire, REVAC_CONTEXT_MAGIC, REVAC_CONTEXT_MAGIC_BYTES) != 0 ||
      wire[6] != REVAC_CONTEXT_WIRE_VERSION ||
      wire[7] != REVAC_CONTEXT_PROFILE_DELTA4) {
    return 0;
  }
  *profile = REVAC_CONTEXT_PROFILE_DELTA4;
  return 1;
}

void revac_delta_context_init(revac_delta_context_t *ctx) {
  memset(ctx, 0, sizeof(*ctx));
  ctx->profile = REVAC_CONTEXT_PROFILE_LOCAL;
  for (size_t layer = 0; layer < ACC_GENPP_DELTA4_LAYERS; layer++) {
    for (size_t gadget = 0; gadget < ACC_GENPP_GADGET_DIM; gadget++) {
      poly_q_vec_d_init(ctx->delta[layer][gadget]);
    }
  }
}

void revac_delta_context_clear(revac_delta_context_t *ctx) {
  for (size_t layer = 0; layer < ACC_GENPP_DELTA4_LAYERS; layer++) {
    for (size_t gadget = 0; gadget < ACC_GENPP_GADGET_DIM; gadget++) {
      poly_q_vec_d_clear(ctx->delta[layer][gadget]);
    }
  }
  memset(ctx, 0, sizeof(*ctx));
}

void revac_context_compute_pp_hash(
    uint8_t out[REVAC_CONTEXT_HASH_BYTES], const acc_genpp_pp_t *pp) {
  static const uint8_t domain[] = "REVAC-ACC-PP-V2";
  static const uint8_t profile_id[REVAC_PARAM_PROFILE_ID_BYTES] =
      REVAC_PARAM_PROFILE_ID_BYTES_INIT;
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, domain, sizeof(domain) - 1);
  shake256_absorb(&state, profile_id, sizeof(profile_id));
  absorb_u32(&state, ACC_GENPP_GADGET_BASE);
  absorb_u32(&state, ACC_GENPP_GADGET_DIM);
  absorb_u32(&state, ACC_GENPP_HANDLE_BITS);
  absorb_u32(&state, PARAM_Q);
  absorb_u32(&state, PARAM_Q_BITLEN);
  shake256_absorb(&state, pp->A.seed, sizeof(pp->A.seed));
  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
    for (size_t gadget = 0; gadget < ACC_GENPP_GADGET_DIM; gadget++) {
      for (size_t row = 0; row < PARAM_D; row++) {
        for (size_t col = 0; col < PARAM_D; col++) {
          for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
            absorb_u32(&state, (uint32_t)poly_q_get_coeff(
                                   pp->Ba.bit[bit].blocks[gadget]
                                       ->rows[row]
                                       ->entries[col],
                                   coeff));
          }
        }
      }
    }
  }
  for (size_t row = 0; row < PARAM_D; row++) {
    for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
      absorb_u32(&state,
                 (uint32_t)poly_q_get_coeff(pp->ua->entries[row], coeff));
    }
  }
  shake256_finalize(&state);
  shake256_squeeze(out, REVAC_CONTEXT_HASH_BYTES, &state);
}

void revac_bml_entry_clear(revac_bml_entry_t *entry) {
  if (entry == NULL) {
    return;
  }
  free(entry->signature);
  free(entry->delta_data);
  memset(entry, 0, sizeof(*entry));
}

void revac_sync_delta_cache_entry_clear(
    revac_sync_delta_cache_entry_t *entry) {
  if (entry == NULL) {
    return;
  }
  free(entry->delta_data);
  memset(entry, 0, sizeof(*entry));
}

void revac_q96_bml_entry_clear(revac_q96_bml_entry_t *entry) {
  if (entry == NULL) return;
  free(entry->signature);
  free(entry->delta_data);
  memset(entry, 0, sizeof(*entry));
}

void revac_q96_delta6_entry_clear(revac_q96_delta6_entry_t *entry) {
  if (entry == NULL) return;
  for (size_t i = 0; i < REVAC_Q96_DELTA_LAYERS; i++) {
    free(entry->layer_data[i]);
    entry->layer_data[i] = NULL;
    entry->layer_len[i] = 0;
    memset(entry->layer_hash[i], 0, REVAC_CONTEXT_HASH_BYTES);
  }
  free(entry->signature);
  memset(entry, 0, sizeof(*entry));
}

int revac_q96_delta6_entry_build(
    revac_q96_delta6_entry_t *entry,
    const acc96_relation_update_t *update,
    const acc96_delta6_t *delta,
    const acc96_ctx_t *ctx,
    const uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t accumulator_hash[REVAC_CONTEXT_HASH_BYTES],
    const sep_sk_t *ta_sk, const sep_pk_t *ta_pk) {
  uint8_t message[PARAM_M * PARAM_N / 8];
  uint8_t state[STATE_BYTES];
  sep_sig_t signature;
  if (entry == NULL || update == NULL || update->sequence == 0 ||
      delta == NULL || ctx == NULL || pp_hash == NULL ||
      accumulator_hash == NULL || ta_sk == NULL || ta_pk == NULL)
    return 0;
  revac_q96_delta6_entry_clear(entry);
  entry->update = *update;
  memcpy(entry->pp_hash, pp_hash, REVAC_CONTEXT_HASH_BYTES);
  memcpy(entry->accumulator_hash, accumulator_hash, REVAC_CONTEXT_HASH_BYTES);
  for (size_t layer = 0; layer < REVAC_Q96_DELTA_LAYERS; layer++) {
    size_t packed_len = 0;
    if (!acc96_delta6_pack_layer_alloc(&entry->layer_data[layer],
                                       &packed_len,
                                       delta, layer, ctx)) {
      revac_q96_delta6_entry_clear(entry);
      return 0;
    }
    if (packed_len > UINT32_MAX) {
      revac_q96_delta6_entry_clear(entry);
      return 0;
    }
    entry->layer_len[layer] = (uint32_t)packed_len;
    sha3_256(entry->layer_hash[layer], entry->layer_data[layer],
             entry->layer_len[layer]);
  }
  entry->signature = malloc(REVAC_CONTEXT_SIGNATURE_BYTES);
  if (entry->signature == NULL) {
    revac_q96_delta6_entry_clear(entry);
    return 0;
  }
  q96_delta6_manifest_message(message, entry);
  sep_sig_init(&signature);
  randombytes(state, sizeof(state));
  sep_sign(&signature, state, ta_sk, ta_pk, message);
  pack_signature(entry->signature, &signature);
  sep_sig_clear(&signature);
  entry->signature_len = REVAC_CONTEXT_SIGNATURE_BYTES;
  entry->used = 1;
  return 1;
}

int revac_q96_delta6_manifest_pack(
    uint8_t **out, size_t *out_len,
    const revac_q96_delta6_entry_t *entry) {
  const size_t total = REVAC_Q96_DELTA6_MANIFEST_HEADER_BYTES +
                       REVAC_CONTEXT_SIGNATURE_BYTES;
  uint8_t *wire;
  if (out == NULL || out_len == NULL) return 0;
  *out = NULL; *out_len = 0;
  if (entry == NULL || !entry->used || entry->signature == NULL ||
      entry->signature_len != REVAC_CONTEXT_SIGNATURE_BYTES)
    return 0;
  wire = calloc(total, 1);
  if (wire == NULL) return 0;
  memcpy(wire, REVAC_Q96_DELTA6_MANIFEST_MAGIC, 6);
  wire[6] = REVAC_Q96_DELTA6_MANIFEST_WIRE_VERSION;
  wire[7] = REVAC_ACC_PROFILE_Q96;
  store_u32(wire + 8, entry->update.sequence);
  store_u32(wire + 12, entry->update.revoked_handle);
  store_u32(wire + 16, (uint32_t)REVAC_CONTEXT_SIGNATURE_BYTES);
  store_u32(wire + 20, REVAC_Q96_DELTA6_MANIFEST_HEADER_BYTES);
  memcpy(wire + 24, entry->pp_hash, REVAC_CONTEXT_HASH_BYTES);
  memcpy(wire + 56, entry->accumulator_hash, REVAC_CONTEXT_HASH_BYTES);
  store_u32(wire + 88, REVAC_Q96_DELTA_LAYERS);
  store_u32(wire + 92, ACC96_DELTA4_CODEC_VERSION);
  for (size_t i = 0; i < REVAC_Q96_DELTA_LAYERS; i++) {
    const size_t off = 104 + i * 36;
    store_u32(wire + off, entry->layer_len[i]);
    memcpy(wire + off + 4, entry->layer_hash[i], REVAC_CONTEXT_HASH_BYTES);
  }
  memcpy(wire + REVAC_Q96_DELTA6_MANIFEST_HEADER_BYTES, entry->signature,
         REVAC_CONTEXT_SIGNATURE_BYTES);
  *out = wire; *out_len = total;
  return 1;
}

int revac_q96_delta6_layer_pack(
    uint8_t **out, size_t *out_len,
    const revac_q96_delta6_entry_t *entry, size_t layer) {
  size_t total;
  uint8_t *wire;
  if (out == NULL || out_len == NULL) return 0;
  *out = NULL; *out_len = 0;
  if (entry == NULL || !entry->used || layer >= REVAC_Q96_DELTA_LAYERS ||
      entry->layer_data[layer] == NULL)
    return 0;
  total = REVAC_Q96_DELTA6_LAYER_HEADER_BYTES + entry->layer_len[layer];
  wire = calloc(total, 1);
  if (wire == NULL) return 0;
  memcpy(wire, REVAC_Q96_DELTA6_LAYER_MAGIC,
         REVAC_Q96_DELTA6_LAYER_MAGIC_BYTES);
  wire[6] = REVAC_Q96_DELTA6_MANIFEST_WIRE_VERSION;
  wire[7] = REVAC_ACC_PROFILE_Q96;
  store_u32(wire + 8, entry->update.sequence);
  store_u32(wire + 12, (uint32_t)layer);
  store_u32(wire + 16, entry->layer_len[layer]);
  memcpy(wire + REVAC_Q96_DELTA6_LAYER_HEADER_BYTES, entry->layer_data[layer],
         entry->layer_len[layer]);
  *out = wire; *out_len = total;
  return 1;
}

int revac_q96_delta6_layer_unpack(
    uint32_t *sequence, size_t *layer,
    const uint8_t **payload, size_t *payload_len,
    const uint8_t *wire, size_t wire_len) {
  uint32_t n, index;
  if (sequence == NULL || layer == NULL || payload == NULL ||
      payload_len == NULL || wire == NULL ||
      wire_len < REVAC_Q96_DELTA6_LAYER_HEADER_BYTES ||
      memcmp(wire, REVAC_Q96_DELTA6_LAYER_MAGIC, 6) != 0 ||
      wire[6] != REVAC_Q96_DELTA6_MANIFEST_WIRE_VERSION ||
      wire[7] != REVAC_ACC_PROFILE_Q96)
    return 0;
  index = load_u32(wire + 12);
  n = load_u32(wire + 16);
  if (index >= REVAC_Q96_DELTA_LAYERS ||
      n != wire_len - REVAC_Q96_DELTA6_LAYER_HEADER_BYTES)
    return 0;
  *layer = index;
  *sequence = load_u32(wire + 8);
  *payload = wire + REVAC_Q96_DELTA6_LAYER_HEADER_BYTES;
  *payload_len = n;
  return 1;
}

int revac_q96_delta6_manifest_unpack(
    acc96_relation_update_t *update,
    uint8_t layer_hash[REVAC_Q96_DELTA_LAYERS][REVAC_CONTEXT_HASH_BYTES],
    uint32_t layer_len[REVAC_Q96_DELTA_LAYERS],
    const uint8_t *wire, size_t wire_len,
    const revac_public_state_t *post_delete_state) {
  uint8_t message[PARAM_M * PARAM_N / 8];
  sep_sig_t signature;
  uint32_t sequence, handle, signature_len, header_len;
  int ok;
  if (update == NULL || layer_hash == NULL || layer_len == NULL ||
      wire == NULL || post_delete_state == NULL || post_delete_state->ta_pk == NULL ||
      wire_len != REVAC_Q96_DELTA6_MANIFEST_HEADER_BYTES +
                      REVAC_CONTEXT_SIGNATURE_BYTES ||
      memcmp(wire, REVAC_Q96_DELTA6_MANIFEST_MAGIC, 6) != 0 ||
      wire[6] != REVAC_Q96_DELTA6_MANIFEST_WIRE_VERSION ||
      wire[7] != REVAC_ACC_PROFILE_Q96)
    return 0;
  sequence = load_u32(wire + 8); handle = load_u32(wire + 12);
  signature_len = load_u32(wire + 16); header_len = load_u32(wire + 20);
  if (post_delete_state->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      sequence == 0 || sequence != post_delete_state->acc96_seq ||
      signature_len != REVAC_CONTEXT_SIGNATURE_BYTES ||
      header_len != REVAC_Q96_DELTA6_MANIFEST_HEADER_BYTES ||
      load_u32(wire + 88) != REVAC_Q96_DELTA_LAYERS ||
      load_u32(wire + 92) != ACC96_DELTA4_CODEC_VERSION ||
      load_u32(wire + 96) != 0 || load_u32(wire + 100) != 0 ||
      memcmp(wire + 24, post_delete_state->context_pp_hash,
             REVAC_CONTEXT_HASH_BYTES) != 0)
    return 0;
  for (size_t i = 0; i < REVAC_Q96_DELTA_LAYERS; i++) {
    const size_t off = 104 + i * 36;
    layer_len[i] = load_u32(wire + off);
    memcpy(layer_hash[i], wire + off + 4, REVAC_CONTEXT_HASH_BYTES);
  }
  {
    uint8_t accumulator_hash[REVAC_CONTEXT_HASH_BYTES];
    if (!revac_public_accumulator_digest(accumulator_hash, post_delete_state) ||
        memcmp(wire + 56, accumulator_hash, REVAC_CONTEXT_HASH_BYTES) != 0)
      return 0;
  }
  {
    revac_q96_delta6_entry_t view;
    memset(&view, 0, sizeof(view));
    view.update.sequence = sequence; view.update.revoked_handle = handle;
    memcpy(view.pp_hash, wire + 24, REVAC_CONTEXT_HASH_BYTES);
    memcpy(view.accumulator_hash, wire + 56, REVAC_CONTEXT_HASH_BYTES);
    for (size_t i = 0; i < REVAC_Q96_DELTA_LAYERS; i++) {
      view.layer_len[i] = layer_len[i];
      memcpy(view.layer_hash[i], layer_hash[i], REVAC_CONTEXT_HASH_BYTES);
    }
    q96_delta6_manifest_message(message, &view);
  }
  sep_sig_init(&signature);
  ok = unpack_signature(&signature,
                        wire + REVAC_Q96_DELTA6_MANIFEST_HEADER_BYTES) &&
       sep_verify(&signature, message, post_delete_state->ta_pk);
  sep_sig_clear(&signature);
  if (!ok) return 0;
  update->sequence = sequence; update->revoked_handle = handle;
  return 1;
}

static int q96_tail4_pack(uint8_t out[REVAC_Q96_TAIL4_BYTES],
                          const acc96_gmat_t *tail4) {
  size_t offset = 0;
  if (out == NULL || tail4 == NULL) return 0;
  for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++)
    offset += acc96_poly_pack(out + offset, &tail4->block[digit]);
  return offset == REVAC_Q96_TAIL4_BYTES;
}

static int q96_tail4_unpack(acc96_gmat_t *tail4, const uint8_t *in,
                            size_t in_len) {
  size_t offset = 0;
  if (tail4 == NULL || in == NULL || in_len != REVAC_Q96_TAIL4_BYTES)
    return 0;
  for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++) {
    if (!acc96_poly_unpack(&tail4->block[digit], in + offset,
                           ACC96_POLY_PACKED_BYTES))
      return 0;
    offset += ACC96_POLY_PACKED_BYTES;
  }
  return offset == in_len;
}

int revac_q96_tail4_prefetch_pack(
    uint8_t **out, size_t *out_len, uint32_t candidate_handle,
    const acc96_gmat_t *tail4,
    const uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const sep_sk_t *ta_sk, const sep_pk_t *ta_pk) {
  const size_t total = REVAC_Q96_PREFETCH_WIRE_HEADER_BYTES +
                       REVAC_CONTEXT_SIGNATURE_BYTES +
                       REVAC_Q96_TAIL4_BYTES;
  uint8_t message[PARAM_M * PARAM_N / 8];
  uint8_t state[STATE_BYTES];
  uint8_t *wire;
  uint8_t *payload;
  sep_sig_t signature;

  if (out == NULL || out_len == NULL || tail4 == NULL || pp_hash == NULL ||
      ta_sk == NULL || ta_pk == NULL)
    return 0;
  *out = NULL;
  *out_len = 0;
  wire = calloc(total, 1);
  if (wire == NULL) return 0;
  payload = wire + REVAC_Q96_PREFETCH_WIRE_HEADER_BYTES +
            REVAC_CONTEXT_SIGNATURE_BYTES;
  if (!q96_tail4_pack(payload, tail4)) {
    free(wire);
    return 0;
  }
  memcpy(wire, REVAC_Q96_PREFETCH_MAGIC, REVAC_Q96_PREFETCH_MAGIC_BYTES);
  wire[6] = REVAC_Q96_PREFETCH_WIRE_VERSION;
  wire[7] = REVAC_ACC_PROFILE_Q96;
  store_u32(wire + 8, candidate_handle);
  store_u32(wire + 12, (uint32_t)REVAC_CONTEXT_SIGNATURE_BYTES);
  store_u32(wire + 16, (uint32_t)REVAC_Q96_PREFETCH_WIRE_HEADER_BYTES);
  store_u32(wire + 20, (uint32_t)REVAC_Q96_TAIL4_BYTES);
  wire[24] = ACC96_TAIL4_BIT;
  wire[25] = REVAC_Q96_UPDATE_CODEC_TAIL4_CRT96;
  memcpy(wire + 32, pp_hash, REVAC_CONTEXT_HASH_BYTES);
  sha3_256(wire + 64, payload, REVAC_Q96_TAIL4_BYTES);
  q96_prefetch_message(message, candidate_handle, wire + 32, wire + 64);
  sep_sig_init(&signature);
  randombytes(state, sizeof(state));
  sep_sign(&signature, state, ta_sk, ta_pk, message);
  pack_signature(wire + REVAC_Q96_PREFETCH_WIRE_HEADER_BYTES, &signature);
  sep_sig_clear(&signature);
  *out = wire;
  *out_len = total;
  return 1;
}

int revac_q96_tail4_prefetch_unpack(
    uint32_t *candidate_handle, acc96_gmat_t *tail4,
    uint8_t tail_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t *wire, size_t wire_len,
    const revac_public_state_t *public_state) {
  const size_t expected = REVAC_Q96_PREFETCH_WIRE_HEADER_BYTES +
                          REVAC_CONTEXT_SIGNATURE_BYTES +
                          REVAC_Q96_TAIL4_BYTES;
  const uint8_t *payload;
  uint8_t computed_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t message[PARAM_M * PARAM_N / 8];
  sep_sig_t signature;
  uint32_t handle;
  int ok;

  if (candidate_handle == NULL || tail4 == NULL || tail_hash == NULL ||
      wire == NULL || public_state == NULL ||
      public_state->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      public_state->ta_pk == NULL || wire_len != expected ||
      memcmp(wire, REVAC_Q96_PREFETCH_MAGIC,
             REVAC_Q96_PREFETCH_MAGIC_BYTES) != 0 ||
      wire[6] != REVAC_Q96_PREFETCH_WIRE_VERSION ||
      wire[7] != REVAC_ACC_PROFILE_Q96 ||
      load_u32(wire + 12) != REVAC_CONTEXT_SIGNATURE_BYTES ||
      load_u32(wire + 16) != REVAC_Q96_PREFETCH_WIRE_HEADER_BYTES ||
      load_u32(wire + 20) != REVAC_Q96_TAIL4_BYTES ||
      wire[24] != ACC96_TAIL4_BIT ||
      wire[25] != REVAC_Q96_UPDATE_CODEC_TAIL4_CRT96 || wire[26] != 0 ||
      wire[27] != 0 || wire[28] != 0 || wire[29] != 0 || wire[30] != 0 ||
      wire[31] != 0 ||
      memcmp(wire + 32, public_state->context_pp_hash,
             REVAC_CONTEXT_HASH_BYTES) != 0)
    return 0;
  payload = wire + REVAC_Q96_PREFETCH_WIRE_HEADER_BYTES +
            REVAC_CONTEXT_SIGNATURE_BYTES;
  sha3_256(computed_hash, payload, REVAC_Q96_TAIL4_BYTES);
  if (memcmp(computed_hash, wire + 64, REVAC_CONTEXT_HASH_BYTES) != 0)
    return 0;
  handle = load_u32(wire + 8);
  q96_prefetch_message(message, handle, wire + 32, wire + 64);
  sep_sig_init(&signature);
  ok = unpack_signature(
           &signature, wire + REVAC_Q96_PREFETCH_WIRE_HEADER_BYTES) &&
       sep_verify(&signature, message, public_state->ta_pk);
  sep_sig_clear(&signature);
  if (!ok || !q96_tail4_unpack(tail4, payload, REVAC_Q96_TAIL4_BYTES))
    return 0;
  *candidate_handle = handle;
  memcpy(tail_hash, wire + 64, REVAC_CONTEXT_HASH_BYTES);
  return 1;
}

int revac_q96_bml_entry_build(
    revac_q96_bml_entry_t *entry,
    const acc96_relation_update_t *update,
    const acc96_gmat_t *tail4,
    const acc96_ctx_t *ctx,
    const uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t accumulator_hash[REVAC_CONTEXT_HASH_BYTES],
    const sep_sk_t *ta_sk, const sep_pk_t *ta_pk) {
  uint8_t message[PARAM_M * PARAM_N / 8];
  uint8_t state[STATE_BYTES];
  sep_sig_t signature;

  if (entry == NULL || update == NULL || update->sequence == 0 ||
      tail4 == NULL || ctx == NULL ||
      pp_hash == NULL || accumulator_hash == NULL || ta_sk == NULL ||
      ta_pk == NULL) {
    return 0;
  }
  revac_q96_bml_entry_clear(entry);
  entry->delta_data_len = REVAC_Q96_TAIL4_BYTES;
  entry->delta_data = malloc(entry->delta_data_len);
  entry->signature = malloc(REVAC_CONTEXT_SIGNATURE_BYTES);
  if (entry->delta_data == NULL || entry->signature == NULL ||
      !q96_tail4_pack(entry->delta_data, tail4)) {
    revac_q96_bml_entry_clear(entry);
    return 0;
  }
  sha3_256(entry->delta_hash, entry->delta_data, entry->delta_data_len);
  entry->update = *update;
  memcpy(entry->pp_hash, pp_hash, REVAC_CONTEXT_HASH_BYTES);
  memcpy(entry->accumulator_hash, accumulator_hash,
         REVAC_CONTEXT_HASH_BYTES);
  q96_update_manifest_message(message, update->sequence,
                              update->revoked_handle, pp_hash,
                              accumulator_hash, entry->delta_hash,
                              (uint32_t)entry->delta_data_len);
  sep_sig_init(&signature);
  randombytes(state, sizeof(state));
  sep_sign(&signature, state, ta_sk, ta_pk, message);
  pack_signature(entry->signature, &signature);
  sep_sig_clear(&signature);
  entry->signature_len = REVAC_CONTEXT_SIGNATURE_BYTES;
  entry->used = 1;
  return 1;
}

int revac_q96_bml_entry_packet_pack(
    uint8_t **out, size_t *out_len, const revac_q96_bml_entry_t *entry) {
  uint8_t *wire;
  const size_t total = REVAC_Q96_UPDATE_WIRE_HEADER_BYTES +
                       REVAC_CONTEXT_SIGNATURE_BYTES + entry->delta_data_len;

  if (out == NULL || out_len == NULL) return 0;
  *out = NULL;
  *out_len = 0;
  if (entry == NULL || !entry->used || entry->update.sequence == 0 ||
      entry->signature == NULL ||
      entry->signature_len != REVAC_CONTEXT_SIGNATURE_BYTES ||
      entry->delta_data == NULL ||
      entry->delta_data_len != REVAC_Q96_TAIL4_BYTES) {
    return 0;
  }
  wire = calloc(total, 1);
  if (wire == NULL) return 0;
  memcpy(wire, REVAC_Q96_UPDATE_MAGIC, REVAC_Q96_UPDATE_MAGIC_BYTES);
  wire[6] = REVAC_Q96_UPDATE_WIRE_VERSION;
  wire[7] = REVAC_ACC_PROFILE_Q96;
  store_u32(wire + 8, entry->update.sequence);
  store_u32(wire + 12, entry->update.revoked_handle);
  store_u32(wire + 16, (uint32_t)REVAC_CONTEXT_SIGNATURE_BYTES);
  store_u32(wire + 20, (uint32_t)REVAC_Q96_UPDATE_WIRE_HEADER_BYTES);
  memcpy(wire + 24, entry->pp_hash, REVAC_CONTEXT_HASH_BYTES);
  memcpy(wire + 56, entry->accumulator_hash, REVAC_CONTEXT_HASH_BYTES);
  memcpy(wire + 88, entry->delta_hash, REVAC_CONTEXT_HASH_BYTES);
  store_u32(wire + 120, (uint32_t)entry->delta_data_len);
  wire[124] = ACC96_TAIL4_BIT;
  wire[125] = REVAC_Q96_UPDATE_CODEC_TAIL4_CRT96;
  memcpy(wire + REVAC_Q96_UPDATE_WIRE_HEADER_BYTES, entry->signature,
         REVAC_CONTEXT_SIGNATURE_BYTES);
  memcpy(wire + REVAC_Q96_UPDATE_WIRE_HEADER_BYTES +
                    REVAC_CONTEXT_SIGNATURE_BYTES,
         entry->delta_data, entry->delta_data_len);
  *out = wire;
  *out_len = total;
  return 1;
}

int revac_q96_bml_entry_manifest_pack(
    uint8_t **out, size_t *out_len, const revac_q96_bml_entry_t *entry) {
  uint8_t *wire;
  const size_t manifest_len = REVAC_Q96_UPDATE_WIRE_HEADER_BYTES +
                              REVAC_CONTEXT_SIGNATURE_BYTES;
  if (out == NULL || out_len == NULL) return 0;
  *out = NULL;
  *out_len = 0;
  if (entry == NULL || !entry->used ||
      entry->signature_len != REVAC_CONTEXT_SIGNATURE_BYTES ||
      entry->delta_data == NULL ||
      entry->delta_data_len != REVAC_Q96_TAIL4_BYTES) {
    return 0;
  }
  wire = calloc(manifest_len, 1);
  if (wire == NULL) return 0;
  memcpy(wire, REVAC_Q96_UPDATE_MAGIC, REVAC_Q96_UPDATE_MAGIC_BYTES);
  wire[6] = REVAC_Q96_UPDATE_WIRE_VERSION;
  wire[7] = REVAC_ACC_PROFILE_Q96;
  store_u32(wire + 8, entry->update.sequence);
  store_u32(wire + 12, entry->update.revoked_handle);
  store_u32(wire + 16, (uint32_t)REVAC_CONTEXT_SIGNATURE_BYTES);
  store_u32(wire + 20, (uint32_t)REVAC_Q96_UPDATE_WIRE_HEADER_BYTES);
  memcpy(wire + 24, entry->pp_hash, REVAC_CONTEXT_HASH_BYTES);
  memcpy(wire + 56, entry->accumulator_hash, REVAC_CONTEXT_HASH_BYTES);
  memcpy(wire + 88, entry->delta_hash, REVAC_CONTEXT_HASH_BYTES);
  store_u32(wire + 120, (uint32_t)entry->delta_data_len);
  wire[124] = ACC96_TAIL4_BIT;
  wire[125] = REVAC_Q96_UPDATE_CODEC_TAIL4_CRT96;
  memcpy(wire + REVAC_Q96_UPDATE_WIRE_HEADER_BYTES, entry->signature,
         REVAC_CONTEXT_SIGNATURE_BYTES);
  *out = wire;
  *out_len = manifest_len;
  return 1;
}

static int q96_update_manifest_verify(
    acc96_relation_update_t *update,
    uint8_t tail_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t *wire, size_t wire_len,
    const revac_public_state_t *post_delete_state) {
  uint8_t message[PARAM_M * PARAM_N / 8];
  uint8_t accumulator_hash[REVAC_CONTEXT_HASH_BYTES];
  sep_sig_t signature;
  uint32_t sequence, revoked_handle, signature_len, header_len, delta_len;
  int ok;
  if (update == NULL || tail_hash == NULL || wire == NULL ||
      post_delete_state == NULL ||
      post_delete_state->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      post_delete_state->ta_pk == NULL ||
      wire_len != REVAC_Q96_UPDATE_WIRE_HEADER_BYTES +
                      REVAC_CONTEXT_SIGNATURE_BYTES ||
      memcmp(wire, REVAC_Q96_UPDATE_MAGIC,
             REVAC_Q96_UPDATE_MAGIC_BYTES) != 0 ||
      wire[6] != REVAC_Q96_UPDATE_WIRE_VERSION ||
      wire[7] != REVAC_ACC_PROFILE_Q96)
    return 0;
  sequence = load_u32(wire + 8);
  revoked_handle = load_u32(wire + 12);
  signature_len = load_u32(wire + 16);
  header_len = load_u32(wire + 20);
  delta_len = load_u32(wire + 120);
  if (sequence == 0 || sequence != post_delete_state->acc96_seq ||
      signature_len != REVAC_CONTEXT_SIGNATURE_BYTES ||
      header_len != REVAC_Q96_UPDATE_WIRE_HEADER_BYTES ||
      delta_len != REVAC_Q96_TAIL4_BYTES ||
      wire[124] != ACC96_TAIL4_BIT ||
      wire[125] != REVAC_Q96_UPDATE_CODEC_TAIL4_CRT96 || wire[126] != 0 ||
      wire[127] != 0 ||
      memcmp(wire + 24, post_delete_state->context_pp_hash,
             REVAC_CONTEXT_HASH_BYTES) != 0 ||
      !revac_public_accumulator_digest(accumulator_hash,
                                       post_delete_state) ||
      memcmp(wire + 56, accumulator_hash, REVAC_CONTEXT_HASH_BYTES) != 0)
    return 0;
  q96_update_manifest_message(message, sequence, revoked_handle, wire + 24,
                              wire + 56, wire + 88, delta_len);
  sep_sig_init(&signature);
  ok = unpack_signature(
           &signature, wire + REVAC_Q96_UPDATE_WIRE_HEADER_BYTES) &&
       sep_verify(&signature, message, post_delete_state->ta_pk);
  sep_sig_clear(&signature);
  if (!ok) return 0;
  update->sequence = sequence;
  update->revoked_handle = revoked_handle;
  memcpy(tail_hash, wire + 88, REVAC_CONTEXT_HASH_BYTES);
  return 1;
}

int revac_q96_update_manifest_unpack(
    acc96_relation_update_t *update,
    uint8_t tail_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t *wire, size_t wire_len,
    const revac_public_state_t *post_delete_state) {
  return q96_update_manifest_verify(update, tail_hash, wire, wire_len,
                                    post_delete_state);
}

static int q96_update_packet_unpack_impl(
    acc96_relation_update_t *update, acc96_gmat_t *tail4,
    const uint8_t *wire, size_t wire_len,
    const revac_public_state_t *post_delete_state) {
  uint8_t message[PARAM_M * PARAM_N / 8];
  uint8_t accumulator_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t delta_hash[REVAC_CONTEXT_HASH_BYTES];
  sep_sig_t signature;
  uint32_t sequence, revoked_handle, signature_len, header_len, delta_len;
  const uint8_t *delta_data;
  int ok;

  if (update == NULL || wire == NULL || post_delete_state == NULL ||
      post_delete_state->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      post_delete_state->ta_pk == NULL ||
      wire_len < REVAC_Q96_UPDATE_WIRE_HEADER_BYTES +
                     REVAC_CONTEXT_SIGNATURE_BYTES ||
      memcmp(wire, REVAC_Q96_UPDATE_MAGIC,
             REVAC_Q96_UPDATE_MAGIC_BYTES) != 0 ||
      wire[6] != REVAC_Q96_UPDATE_WIRE_VERSION ||
      wire[7] != REVAC_ACC_PROFILE_Q96) {
    return 0;
  }
  sequence = load_u32(wire + 8);
  revoked_handle = load_u32(wire + 12);
  signature_len = load_u32(wire + 16);
  header_len = load_u32(wire + 20);
  delta_len = load_u32(wire + 120);
  delta_data = wire + REVAC_Q96_UPDATE_WIRE_HEADER_BYTES +
               REVAC_CONTEXT_SIGNATURE_BYTES;
  if (sequence == 0 || sequence != post_delete_state->acc96_seq ||
      signature_len != REVAC_CONTEXT_SIGNATURE_BYTES ||
      header_len != REVAC_Q96_UPDATE_WIRE_HEADER_BYTES ||
      delta_len != REVAC_Q96_TAIL4_BYTES ||
      wire_len != REVAC_Q96_UPDATE_WIRE_HEADER_BYTES +
                      REVAC_CONTEXT_SIGNATURE_BYTES + delta_len ||
      wire[124] != ACC96_TAIL4_BIT ||
      wire[125] != REVAC_Q96_UPDATE_CODEC_TAIL4_CRT96 || wire[126] != 0 ||
      wire[127] != 0 ||
      memcmp(wire + 24, post_delete_state->context_pp_hash,
             REVAC_CONTEXT_HASH_BYTES) != 0 ||
      !revac_public_accumulator_digest(accumulator_hash,
                                       post_delete_state) ||
      memcmp(wire + 56, accumulator_hash,
             REVAC_CONTEXT_HASH_BYTES) != 0) {
    return 0;
  }
  sha3_256(delta_hash, delta_data, delta_len);
  if (memcmp(wire + 88, delta_hash, REVAC_CONTEXT_HASH_BYTES) != 0)
    return 0;
  q96_update_manifest_message(message, sequence, revoked_handle, wire + 24,
                              wire + 56, wire + 88, delta_len);
  sep_sig_init(&signature);
  ok = unpack_signature(
           &signature, wire + REVAC_Q96_UPDATE_WIRE_HEADER_BYTES) &&
       sep_verify(&signature, message, post_delete_state->ta_pk);
  sep_sig_clear(&signature);
  if (!ok) return 0;
  if (tail4 != NULL && !q96_tail4_unpack(tail4, delta_data, delta_len))
    return 0;
  update->sequence = sequence;
  update->revoked_handle = revoked_handle;
  return 1;
}

int revac_q96_update_packet_unpack(
    acc96_relation_update_t *update, const uint8_t *wire, size_t wire_len,
    const revac_public_state_t *post_delete_state) {
  return q96_update_packet_unpack_impl(update, NULL, wire, wire_len,
                                       post_delete_state);
}

int revac_q96_update_packet_unpack_tail4(
    acc96_relation_update_t *update, acc96_gmat_t *tail4,
    const uint8_t *wire, size_t wire_len,
    const revac_public_state_t *post_delete_state) {
  return q96_update_packet_unpack_impl(update, tail4, wire, wire_len,
                                       post_delete_state);
}

int revac_bml_entry_and_sync_cache_build(
    revac_bml_entry_t *bml_entry,
    revac_sync_delta_cache_entry_t *sync_entry,
    uint32_t seq, uint32_t handle, int sign,
    const uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const acc_genpp_gmat_t tails[REVAC_TA_SYNC_DELTA_LAYERS],
    const poly_q_vec_d u, const sep_sk_t *ta_sk, const sep_pk_t *ta_pk) {
  uint8_t message[PARAM_M * PARAM_N / 8];
  uint8_t state[STATE_BYTES];
  poly_q_vec_d delta[REVAC_TA_SYNC_DELTA_LAYERS][ACC_GENPP_GADGET_DIM];
  sep_sig_t signature;
  int ok = 0;

  if (bml_entry == NULL || sync_entry == NULL || pp_hash == NULL ||
      tails == NULL || ta_sk == NULL || ta_pk == NULL) {
    return 0;
  }
  if (sign != 1 && sign != -1) {
    return 0;
  }
  revac_bml_entry_clear(bml_entry);
  revac_sync_delta_cache_entry_clear(sync_entry);
  for (size_t layer = 0; layer < REVAC_TA_SYNC_DELTA_LAYERS; layer++) {
    for (size_t gadget = 0; gadget < ACC_GENPP_GADGET_DIM; gadget++) {
      poly_q_vec_d_init(delta[layer][gadget]);
    }
  }

  acc_genpp_delta_from_tails(delta, REVAC_TA_SYNC_DELTA_LAYERS, tails,
                            handle, u);

  sync_entry->seq = seq;
  sync_entry->handle = handle;
  sync_entry->sign = (int8_t)sign;
  memcpy(sync_entry->pp_hash, pp_hash, REVAC_CONTEXT_HASH_BYTES);
  if (!delta_pack(&sync_entry->delta_data, &sync_entry->delta_data_len,
                  delta, REVAC_TA_SYNC_DELTA_LAYERS)) {
    goto cleanup;
  }
  sha3_256(sync_entry->delta_hash, sync_entry->delta_data,
           sync_entry->delta_data_len);
  sync_entry->used = 1;

  bml_entry->signature = malloc(REVAC_CONTEXT_SIGNATURE_BYTES);
  if (bml_entry->signature == NULL ||
      !delta_pack(&bml_entry->delta_data, &bml_entry->delta_data_len,
                  delta, ACC_GENPP_DELTA4_LAYERS)) {
    goto cleanup;
  }
  bml_entry->update.seq = seq;
  bml_entry->update.op = ACC_GENPP_OP_DEL;
  bml_entry->update.handle = handle;
  bml_entry->update.sign = (int8_t)sign;
  bml_entry->signature_len = REVAC_CONTEXT_SIGNATURE_BYTES;
  memcpy(bml_entry->pp_hash, pp_hash, REVAC_CONTEXT_HASH_BYTES);
  sha3_256(bml_entry->delta_hash, bml_entry->delta_data,
           bml_entry->delta_data_len);
  manifest_message(message, seq, handle, sign,
                   (uint32_t)bml_entry->delta_data_len,
                   bml_entry->pp_hash, bml_entry->delta_hash);
  sep_sig_init(&signature);
  randombytes(state, sizeof(state));
  sep_sign(&signature, state, ta_sk, ta_pk, message);
  pack_signature(bml_entry->signature, &signature);
  sep_sig_clear(&signature);
  bml_entry->used = 1;
  ok = 1;

cleanup:
  for (size_t layer = 0; layer < REVAC_TA_SYNC_DELTA_LAYERS; layer++) {
    for (size_t gadget = 0; gadget < ACC_GENPP_GADGET_DIM; gadget++) {
      poly_q_vec_d_clear(delta[layer][gadget]);
    }
  }
  if (!ok) {
    revac_bml_entry_clear(bml_entry);
    revac_sync_delta_cache_entry_clear(sync_entry);
  }
  return ok;
}

int revac_sync_delta_cache_entry_apply(
    acc_genpp_wit_t *witness,
    const revac_sync_delta_cache_entry_t *entry,
    const uint8_t expected_pp_hash[REVAC_CONTEXT_HASH_BYTES],
    uint32_t member_handle,
    const acc_genpp_update_t *update, size_t *local_composites) {
  poly_q_vec_d delta[REVAC_TA_SYNC_DELTA_LAYERS][ACC_GENPP_GADGET_DIM];
  uint8_t computed_hash[REVAC_CONTEXT_HASH_BYTES];
  int ok = 0;

  if (local_composites != NULL) {
    *local_composites = 0;
  }
  if (witness == NULL || entry == NULL || expected_pp_hash == NULL ||
      update == NULL ||
      !entry->used || entry->seq != update->seq ||
      entry->handle != update->handle || entry->sign != update->sign ||
      !acc_genpp_update_sign_is_explicit(update) ||
      entry->delta_data == NULL ||
      entry->delta_data_len < REVAC_SYNC_DELTA_FAST_BYTES ||
      memcmp(entry->pp_hash, expected_pp_hash,
             REVAC_CONTEXT_HASH_BYTES) != 0) {
    return 0;
  }
  sha3_256(computed_hash, entry->delta_data, entry->delta_data_len);
  if (memcmp(computed_hash, entry->delta_hash,
             REVAC_CONTEXT_HASH_BYTES) != 0) {
    return 0;
  }
  for (size_t layer = 0; layer < REVAC_TA_SYNC_DELTA_LAYERS; layer++) {
    for (size_t gadget = 0; gadget < ACC_GENPP_GADGET_DIM; gadget++) {
      poly_q_vec_d_init(delta[layer][gadget]);
    }
  }
  if (delta_unpack(delta, REVAC_TA_SYNC_DELTA_LAYERS, entry->delta_data,
                   entry->delta_data_len)) {
    ok = acc_genpp_memwitup_delta(witness, member_handle, update, delta,
                                  REVAC_TA_SYNC_DELTA_LAYERS, NULL);
  }

  for (size_t layer = 0; layer < REVAC_TA_SYNC_DELTA_LAYERS; layer++) {
    for (size_t gadget = 0; gadget < ACC_GENPP_GADGET_DIM; gadget++) {
      poly_q_vec_d_clear(delta[layer][gadget]);
    }
  }
  return ok;
}

int revac_bml_entry_packet_pack(uint8_t **out, size_t *out_len,
                                const revac_bml_entry_t *entry) {
  size_t total;
  size_t off;
  uint8_t *wire;

  if (out == NULL || out_len == NULL) {
    return 0;
  }
  *out = NULL;
  *out_len = 0;
  if (entry == NULL || !entry->used || entry->delta_data == NULL ||
      entry->delta_data_len > UINT32_MAX ||
      entry->signature_len != REVAC_CONTEXT_SIGNATURE_BYTES) {
    return 0;
  }
  total = REVAC_CONTEXT_WIRE_HEADER_BYTES + REVAC_CONTEXT_SIGNATURE_BYTES +
          entry->delta_data_len;
  wire = calloc(total, 1);
  if (wire == NULL) {
    return 0;
  }
  memcpy(wire, REVAC_CONTEXT_MAGIC, REVAC_CONTEXT_MAGIC_BYTES);
  wire[6] = REVAC_CONTEXT_WIRE_VERSION;
  wire[7] = REVAC_CONTEXT_PROFILE_DELTA4;
  wire[8] = ACC_GENPP_GADGET_BASE;
  wire[9] = ACC_GENPP_GADGET_DIM;
  wire[10] = ACC_GENPP_HANDLE_BITS;
  wire[11] = PARAM_Q_BITLEN;
  wire[12] = REVAC_CONTEXT_DELTA_CODEC_ZIGZAG9_ESCAPE19;
  wire[13] = ACC_GENPP_DELTA4_LAYERS;
  wire[14] = encode_update_sign(entry->update.sign);
  store_u32(wire + 16, entry->update.seq);
  store_u32(wire + 20, entry->update.handle);
  store_u32(wire + 24, (uint32_t)REVAC_CONTEXT_DELTA_COEFFS);
  store_u32(wire + 28, (uint32_t)entry->delta_data_len);
  store_u32(wire + 32, (uint32_t)REVAC_CONTEXT_SIGNATURE_BYTES);
  store_u32(wire + 36, (uint32_t)REVAC_CONTEXT_WIRE_HEADER_BYTES);
  memcpy(wire + 40, entry->pp_hash, REVAC_CONTEXT_HASH_BYTES);
  memcpy(wire + 72, entry->delta_hash, REVAC_CONTEXT_HASH_BYTES);
  off = REVAC_CONTEXT_WIRE_HEADER_BYTES;
  memcpy(wire + off, entry->signature, REVAC_CONTEXT_SIGNATURE_BYTES);
  off += REVAC_CONTEXT_SIGNATURE_BYTES;
  memcpy(wire + off, entry->delta_data, entry->delta_data_len);
  *out = wire;
  *out_len = total;
  return 1;
}

int revac_context_packet_unpack(
    revac_delta_context_t *ctx, const uint8_t *wire, size_t wire_len,
    const acc_genpp_update_t *expected_update,
    const uint8_t expected_pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const sep_pk_t *ta_pk) {
  uint8_t message[PARAM_M * PARAM_N / 8];
  uint8_t computed_hash[REVAC_CONTEXT_HASH_BYTES];
  sep_sig_t signature;
  uint32_t seq, handle, coeffs, payload_len, signature_len, header_len;
  size_t expected_len;
  const uint8_t *payload;
  int signature_ok;
  int sign;

  if (ctx == NULL || wire == NULL || expected_update == NULL ||
      expected_pp_hash == NULL || ta_pk == NULL ||
      wire_len < REVAC_CONTEXT_WIRE_HEADER_BYTES ||
      memcmp(wire, REVAC_CONTEXT_MAGIC, REVAC_CONTEXT_MAGIC_BYTES) != 0 ||
      wire[6] != REVAC_CONTEXT_WIRE_VERSION ||
      wire[7] != REVAC_CONTEXT_PROFILE_DELTA4 ||
      wire[8] != ACC_GENPP_GADGET_BASE ||
      wire[9] != ACC_GENPP_GADGET_DIM ||
      wire[10] != ACC_GENPP_HANDLE_BITS || wire[11] != PARAM_Q_BITLEN ||
      wire[12] != REVAC_CONTEXT_DELTA_CODEC_ZIGZAG9_ESCAPE19 ||
      wire[13] != ACC_GENPP_DELTA4_LAYERS || wire[15] != 0) {
    return 0;
  }
  sign = decode_update_sign(wire[14]);
  seq = load_u32(wire + 16);
  handle = load_u32(wire + 20);
  coeffs = load_u32(wire + 24);
  payload_len = load_u32(wire + 28);
  signature_len = load_u32(wire + 32);
  header_len = load_u32(wire + 36);
  if (expected_update->op != ACC_GENPP_OP_DEL ||
      seq != expected_update->seq || handle != expected_update->handle ||
      sign == 0 || sign != acc_genpp_update_sign(expected_update) ||
      !acc_genpp_update_sign_is_explicit(expected_update) ||
      coeffs != REVAC_CONTEXT_DELTA_COEFFS || payload_len == 0 ||
      signature_len != REVAC_CONTEXT_SIGNATURE_BYTES ||
      header_len != REVAC_CONTEXT_WIRE_HEADER_BYTES ||
      memcmp(wire + 40, expected_pp_hash, REVAC_CONTEXT_HASH_BYTES) != 0) {
    return 0;
  }
  expected_len = REVAC_CONTEXT_WIRE_HEADER_BYTES + signature_len +
                 (size_t)payload_len;
  if (wire_len != expected_len) {
    return 0;
  }
  payload = wire + REVAC_CONTEXT_WIRE_HEADER_BYTES + signature_len;
  sha3_256(computed_hash, payload, payload_len);
  if (memcmp(computed_hash, wire + 72, REVAC_CONTEXT_HASH_BYTES) != 0) {
    return 0;
  }
  manifest_message(message, seq, handle, sign, payload_len, wire + 40,
                   wire + 72);
  sep_sig_init(&signature);
  if (!unpack_signature(&signature,
                        wire + REVAC_CONTEXT_WIRE_HEADER_BYTES)) {
    sep_sig_clear(&signature);
    return 0;
  }
  signature_ok = sep_verify(&signature, message, ta_pk);
  sep_sig_clear(&signature);
  if (!signature_ok ||
      !delta_unpack(ctx->delta, ACC_GENPP_DELTA4_LAYERS,
                    payload, payload_len)) {
    return 0;
  }
  ctx->profile = REVAC_CONTEXT_PROFILE_DELTA4;
  ctx->seq = seq;
  ctx->handle = handle;
  ctx->sign = (int8_t)sign;
  return 1;
}

int revac_context_packet_verify(
    const uint8_t *wire, size_t wire_len,
    const acc_genpp_update_t *expected_update,
    const uint8_t expected_pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const sep_pk_t *ta_pk) {
  revac_delta_context_t ctx;
  int ok;

  revac_delta_context_init(&ctx);
  ok = revac_context_packet_unpack(&ctx, wire, wire_len, expected_update,
                                   expected_pp_hash, ta_pk);
  revac_delta_context_clear(&ctx);
  return ok;
}

static void variance_absorb_u32(keccak_state *state, uint32_t value) {
  uint8_t encoded[4];
  store_u32(encoded, value);
  shake256_absorb(state, encoded, sizeof(encoded));
}

static void variance_absorb_u64(keccak_state *state, uint64_t value) {
  uint8_t encoded[8];
  for (size_t i = 0; i < sizeof(encoded); i++) {
    encoded[i] = (uint8_t)(value >> (8 * i));
  }
  shake256_absorb(state, encoded, sizeof(encoded));
}

void revac_window_variance_ledger_init(
    revac_window_variance_ledger_t *ledger,
    const revac_public_state_t *public_state,
    uint32_t member_handle, uint32_t window_id,
    uint32_t target_deletions) {
  static const uint8_t domain[] = "REVAC-WINDOW-VARIANCE-LEDGER-V1";
  keccak_state state;

  if (ledger == NULL) {
    return;
  }
  memset(ledger, 0, sizeof(*ledger));
  if (public_state == NULL || public_state->ta_pk == NULL ||
      public_state->acc_pp == NULL || target_deletions == 0) {
    return;
  }
  ledger->version = REVAC_WINDOW_VARIANCE_LEDGER_VERSION;
  ledger->window_id = window_id;
  ledger->member_handle = member_handle;
  ledger->target_deletions = target_deletions;
  memcpy(ledger->pp_hash, public_state->context_pp_hash,
         REVAC_CONTEXT_HASH_BYTES);
  shake256_init(&state);
  shake256_absorb(&state, domain, sizeof(domain) - 1);
  shake256_absorb(&state, ledger->pp_hash, REVAC_CONTEXT_HASH_BYTES);
  variance_absorb_u32(&state, ledger->version);
  variance_absorb_u32(&state, window_id);
  variance_absorb_u32(&state, member_handle);
  variance_absorb_u32(&state, target_deletions);
  shake256_finalize(&state);
  shake256_squeeze(ledger->transcript_hash, REVAC_CONTEXT_HASH_BYTES, &state);
  ledger->initialized = 1;
}

static int witness_right_distance_norm2(
    uint64_t *out, const acc_genpp_wit_t *before,
    const acc_genpp_wit_t *after) {
  poly_q_vec_d difference;
  uint64_t total = 0;
  int ok = 1;

  if (out == NULL || before == NULL || after == NULL) {
    return 0;
  }
  poly_q_vec_d_init(difference);
  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS && ok; bit++) {
    for (size_t digit = 0; digit < ACC_GENPP_GADGET_DIM; digit++) {
      uint64_t term;
      poly_q_vec_d_sub(difference, after->right[bit][digit],
                       before->right[bit][digit]);
      term = poly_q_vec_d_norm2(difference);
      if (UINT64_MAX - total < term) {
        ok = 0;
        break;
      }
      total += term;
    }
  }
  poly_q_vec_d_clear(difference);
  if (ok) {
    *out = total;
  }
  return ok;
}

int revac_window_variance_ledger_record(
    revac_window_variance_ledger_t *ledger,
    const revac_public_state_t *public_state,
    const acc_genpp_update_t *update,
    const uint8_t *wire, size_t wire_len,
    const acc_genpp_wit_t *before, const acc_genpp_wit_t *after) {
  static const uint8_t domain[] = "REVAC-WINDOW-VARIANCE-STEP-V1";
  uint8_t next_hash[REVAC_CONTEXT_HASH_BYTES];
  uint64_t increment_norm2;
  keccak_state state;

  if (ledger == NULL || !ledger->initialized || public_state == NULL ||
      update == NULL || update->op != ACC_GENPP_OP_DEL ||
      !acc_genpp_update_sign_is_explicit(update) ||
      ledger->deletion_count >= ledger->target_deletions ||
      memcmp(ledger->pp_hash, public_state->context_pp_hash,
             REVAC_CONTEXT_HASH_BYTES) != 0 ||
      update->handle == ledger->member_handle ||
      (ledger->deletion_count != 0 && update->seq != ledger->last_seq + 1) ||
      !revac_context_packet_verify(wire, wire_len, update,
                                   public_state->context_pp_hash,
                                   public_state->ta_pk) ||
      !witness_right_distance_norm2(&increment_norm2, before, after) ||
      UINT64_MAX - ledger->quadratic_variation < increment_norm2) {
    return 0;
  }

  shake256_init(&state);
  shake256_absorb(&state, domain, sizeof(domain) - 1);
  shake256_absorb(&state, ledger->transcript_hash,
                  REVAC_CONTEXT_HASH_BYTES);
  variance_absorb_u32(&state, update->seq);
  variance_absorb_u32(&state, update->handle);
  variance_absorb_u32(&state, encode_update_sign(update->sign));
  variance_absorb_u64(&state, increment_norm2);
  shake256_finalize(&state);
  shake256_squeeze(next_hash, sizeof(next_hash), &state);

  if (ledger->deletion_count == 0) {
    ledger->first_seq = update->seq;
  }
  ledger->last_seq = update->seq;
  ledger->deletion_count++;
  ledger->quadratic_variation += increment_norm2;
  if (increment_norm2 > ledger->max_increment_norm2) {
    ledger->max_increment_norm2 = increment_norm2;
  }
  memcpy(ledger->transcript_hash, next_hash, sizeof(next_hash));
  return 1;
}
