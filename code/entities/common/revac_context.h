#ifndef REVAC_CONTEXT_H
#define REVAC_CONTEXT_H

#include <stddef.h>
#include <stdint.h>

#include "revocable_ac.h"
#include "acc96_delta4_codec.h"

#define REVAC_CONTEXT_WIRE_VERSION REVAC_PARAM_CONTEXT_WIRE_VERSION
#define REVAC_CONTEXT_DELTA_CODEC_ZIGZAG9_ESCAPE19 1
#define REVAC_CONTEXT_DELTA_FAST_BITS 9
#define REVAC_CONTEXT_DELTA_ESCAPE ((1U << REVAC_CONTEXT_DELTA_FAST_BITS) - 1U)
#define REVAC_CONTEXT_SIGNATURE_COEFFS                                 \
  ((size_t)PARAM_N * (2 + PARAM_D * (1 + PARAM_KH)))
#define REVAC_CONTEXT_SIGNATURE_BYTES                                  \
  ((REVAC_CONTEXT_SIGNATURE_COEFFS * PARAM_Q_BITLEN + 7) / 8)
#define REVAC_CONTEXT_DELTA_COEFFS                                     \
  ((size_t)ACC_GENPP_DELTA4_LAYERS * ACC_GENPP_GADGET_DIM * PARAM_D * \
   PARAM_N)
#define REVAC_CONTEXT_DELTA_FAST_BYTES                                 \
  ((REVAC_CONTEXT_DELTA_COEFFS * REVAC_CONTEXT_DELTA_FAST_BITS + 7) / 8)
#define REVAC_SYNC_DELTA_COEFFS                                       \
  ((size_t)REVAC_TA_SYNC_DELTA_LAYERS * ACC_GENPP_GADGET_DIM *        \
   PARAM_D * PARAM_N)
#define REVAC_SYNC_DELTA_FAST_BYTES                                   \
  ((REVAC_SYNC_DELTA_COEFFS * REVAC_CONTEXT_DELTA_FAST_BITS + 7) / 8)
#define REVAC_CONTEXT_WIRE_HEADER_BYTES 104
#define REVAC_WINDOW_VARIANCE_LEDGER_VERSION 1
#define REVAC_Q96_UPDATE_WIRE_VERSION 5
#define REVAC_Q96_UPDATE_CODEC_TAIL4_CRT96 4
#define REVAC_Q96_TAIL4_BYTES                                        \
  ((size_t)ACC96_GADGET_DIGITS * ACC96_POLY_PACKED_BYTES)
#define REVAC_Q96_UPDATE_WIRE_HEADER_BYTES 128
#define REVAC_Q96_PREFETCH_WIRE_VERSION 1
#define REVAC_Q96_PREFETCH_WIRE_HEADER_BYTES 96
#define REVAC_Q96_DELTA6_MANIFEST_WIRE_VERSION 1
#define REVAC_Q96_DELTA6_MANIFEST_HEADER_BYTES 320
#define REVAC_Q96_DELTA6_LAYER_MAGIC "RVQ88L"
#define REVAC_Q96_DELTA6_LAYER_MAGIC_BYTES 6
#define REVAC_Q96_DELTA6_LAYER_HEADER_BYTES 20

/* Auditable transcript statistics for one survivor in one refresh window.
   Every recorded step must first pass the TA-signed DELTA4 packet check.
   These are observed values used to check/tune a pre-declared V_max/L_max;
   they are deliberately not presented as an a-priori security bound. */
typedef struct {
  uint32_t version;
  uint32_t window_id;
  uint32_t member_handle;
  uint32_t target_deletions;
  uint32_t first_seq;
  uint32_t last_seq;
  uint32_t deletion_count;
  uint64_t quadratic_variation;
  uint64_t max_increment_norm2;
  uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t transcript_hash[REVAC_CONTEXT_HASH_BYTES];
  int initialized;
} revac_window_variance_ledger_t;

const char *revac_context_profile_name(revac_context_profile_t profile);
size_t revac_context_profile_payload_bytes(revac_context_profile_t profile);
size_t revac_context_profile_wire_bytes(revac_context_profile_t profile);
double revac_context_profile_expected_composites(
    revac_context_profile_t profile);
revac_context_profile_t revac_context_profile_select(
    double goodput_mbps,
    double recurrence_ms);
int revac_context_packet_peek_profile(
    const uint8_t *wire,
    size_t wire_len,
    revac_context_profile_t *profile);

void revac_delta_context_init(revac_delta_context_t *ctx);
void revac_delta_context_clear(revac_delta_context_t *ctx);

void revac_context_compute_pp_hash(
    uint8_t out[REVAC_CONTEXT_HASH_BYTES],
    const acc_genpp_pp_t *pp);

void revac_bml_entry_clear(revac_bml_entry_t *entry);
void revac_sync_delta_cache_entry_clear(
    revac_sync_delta_cache_entry_t *entry);

void revac_q96_bml_entry_clear(revac_q96_bml_entry_t *entry);
void revac_q96_delta6_entry_clear(revac_q96_delta6_entry_t *entry);
int revac_q96_delta6_entry_build(
    revac_q96_delta6_entry_t *entry,
    const acc96_relation_update_t *update,
    const acc96_delta6_t *delta,
    const acc96_ctx_t *ctx,
    const uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t accumulator_hash[REVAC_CONTEXT_HASH_BYTES],
    const sep_sk_t *ta_sk, const sep_pk_t *ta_pk);
int revac_q96_delta6_manifest_pack(
    uint8_t **out, size_t *out_len,
    const revac_q96_delta6_entry_t *entry);
int revac_q96_delta6_layer_pack(
    uint8_t **out, size_t *out_len,
    const revac_q96_delta6_entry_t *entry, size_t layer);
int revac_q96_delta6_manifest_unpack(
    acc96_relation_update_t *update,
    uint8_t layer_hash[REVAC_Q96_DELTA_LAYERS][REVAC_CONTEXT_HASH_BYTES],
    uint32_t layer_len[REVAC_Q96_DELTA_LAYERS],
    const uint8_t *wire, size_t wire_len,
    const revac_public_state_t *post_delete_state);
int revac_q96_delta6_layer_unpack(
    uint32_t *sequence, size_t *layer,
    const uint8_t **payload, size_t *payload_len,
    const uint8_t *wire, size_t wire_len);
int revac_q96_bml_entry_build(
    revac_q96_bml_entry_t *entry,
    const acc96_relation_update_t *update,
    const acc96_gmat_t *tail4,
    const acc96_ctx_t *ctx,
    const uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t accumulator_hash[REVAC_CONTEXT_HASH_BYTES],
    const sep_sk_t *ta_sk, const sep_pk_t *ta_pk);
int revac_q96_bml_entry_packet_pack(
    uint8_t **out, size_t *out_len, const revac_q96_bml_entry_t *entry);
int revac_q96_bml_entry_manifest_pack(
    uint8_t **out, size_t *out_len, const revac_q96_bml_entry_t *entry);
/* Verify the TA signature and its binding to the supplied post-deletion
   public state, then return the authenticated update. */
int revac_q96_update_packet_unpack(
    acc96_relation_update_t *update, const uint8_t *wire, size_t wire_len,
    const revac_public_state_t *post_delete_state);
int revac_q96_update_packet_unpack_tail4(
    acc96_relation_update_t *update, acc96_gmat_t *tail4,
    const uint8_t *wire, size_t wire_len,
    const revac_public_state_t *post_delete_state);
int revac_q96_update_manifest_unpack(
    acc96_relation_update_t *update,
    uint8_t tail_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t *wire, size_t wire_len,
    const revac_public_state_t *post_delete_state);
/* An idle-time packet authenticated independently of any deletion sequence or
   accumulator state.  The later signed update manifest binds its tail hash to
   the actual post-delete state. */
int revac_q96_tail4_prefetch_pack(
    uint8_t **out, size_t *out_len, uint32_t candidate_handle,
    const acc96_gmat_t *tail4,
    const uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const sep_sk_t *ta_sk, const sep_pk_t *ta_pk);
int revac_q96_tail4_prefetch_unpack(
    uint32_t *candidate_handle, acc96_gmat_t *tail4,
    uint8_t tail_hash[REVAC_CONTEXT_HASH_BYTES],
    const uint8_t *wire, size_t wire_len,
    const revac_public_state_t *public_state);

/* Derive DELTA6 once, retain all six layers in TA's private sync cache,
   and optionally create the signed public DELTA4 prefix. */
int revac_bml_entry_and_sync_cache_build(
    revac_bml_entry_t *bml_entry,
    revac_sync_delta_cache_entry_t *sync_entry,
    uint32_t seq,
    uint32_t handle,
    int sign,
    const uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const acc_genpp_gmat_t tails[REVAC_TA_SYNC_DELTA_LAYERS],
    const poly_q_vec_d u,
    const sep_sk_t *ta_sk,
    const sep_pk_t *ta_pk);

/* Apply the TA-private DELTA6 cache without exposing it on the wire. */
int revac_sync_delta_cache_entry_apply(
    acc_genpp_wit_t *witness,
    const revac_sync_delta_cache_entry_t *entry,
    const uint8_t expected_pp_hash[REVAC_CONTEXT_HASH_BYTES],
    uint32_t member_handle,
    const acc_genpp_update_t *update,
    size_t *local_composites);

int revac_bml_entry_packet_pack(
    uint8_t **out,
    size_t *out_len,
    const revac_bml_entry_t *entry);

/* Verify a mandatory BML DELTA4 packet without applying it to a witness. */
int revac_context_packet_verify(
    const uint8_t *wire,
    size_t wire_len,
    const acc_genpp_update_t *expected_update,
    const uint8_t expected_pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const sep_pk_t *ta_pk);

int revac_context_packet_unpack(
    revac_delta_context_t *ctx,
    const uint8_t *wire,
    size_t wire_len,
    const acc_genpp_update_t *expected_update,
    const uint8_t expected_pp_hash[REVAC_CONTEXT_HASH_BYTES],
    const sep_pk_t *ta_pk);

void revac_window_variance_ledger_init(
    revac_window_variance_ledger_t *ledger,
    const revac_public_state_t *public_state,
    uint32_t member_handle,
    uint32_t window_id,
    uint32_t target_deletions);

/* Authenticate wire/update, then append the exact centered L2 norm squared
   of after.right-before.right.  Sequence gaps, overflow, a revoked member,
   or a forged sign/packet are rejected without changing the ledger. */
int revac_window_variance_ledger_record(
    revac_window_variance_ledger_t *ledger,
    const revac_public_state_t *public_state,
    const acc_genpp_update_t *update,
    const uint8_t *wire,
    size_t wire_len,
    const acc_genpp_wit_t *before,
    const acc_genpp_wit_t *after);

#endif
