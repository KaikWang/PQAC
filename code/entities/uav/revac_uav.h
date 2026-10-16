#ifndef REVAC_UAV_H
#define REVAC_UAV_H

#include "../common/revocable_ac.h"

/* The callbacks return malloc-compatible packets owned by the UAV caller.
   A sync transport must authenticate the peer as request->handle and protect
   the returned per-member witness in transit (for example, mTLS or QUIC). */
typedef int (*revac_context_fetch_fn)(void *opaque,
                                      uint32_t seq,
                                      revac_context_profile_t profile,
                                      uint8_t **wire,
                                      size_t *wire_len);
typedef int (*revac_sync_fetch_fn)(void *opaque,
                                   const revac_sync_request_t *request,
                                   uint8_t **wire,
                                   size_t *wire_len);

void revac_uav_init(revac_uav_t *uav);
void revac_uav_clear(revac_uav_t *uav);
void revac_uav_keygen(revac_uav_t *uav, uint64_t id, const uint8_t seed[SEED_BYTES]);

int revac_uav_apply_update_q96(revac_uav_t *uav,
                               const revac_public_state_t *public_state,
                               const acc96_relation_update_t *update);
int revac_uav_apply_update_q96_delta4(
    revac_uav_t *uav, const revac_public_state_t *public_state,
    const acc96_relation_update_t *update, const acc96_delta4_t *delta,
    size_t *applied_layers);
int revac_uav_apply_update_q96_wire(
    revac_uav_t *uav, const revac_public_state_t *public_state,
    const uint8_t *wire, size_t wire_len);
int revac_uav_cache_q96_tail4(revac_uav_t *uav,
                              const revac_public_state_t *public_state,
                              const uint8_t *full_wire, size_t full_wire_len);
int revac_uav_cache_q96_tail4_prefetch(
    revac_uav_t *uav, const revac_public_state_t *public_state,
    const uint8_t *prefetch_wire, size_t prefetch_wire_len);
int revac_uav_apply_update_q96_manifest(
    revac_uav_t *uav, const revac_public_state_t *public_state,
    const uint8_t *manifest_wire, size_t manifest_wire_len);
/* Selective DELTA6 path.  The manifest is always verified first.  For
 * first-difference j<6, layer_packets[0..j] are required; for j>=6 the
 * function ignores packets and performs the authenticated LOCAL recurrence. */
int revac_uav_apply_update_q96_delta6(
    revac_uav_t *uav, const revac_public_state_t *public_state,
    const uint8_t *manifest_wire, size_t manifest_wire_len,
    const uint8_t *layer_packets[REVAC_Q96_DELTA_LAYERS],
    const size_t layer_packet_lens[REVAC_Q96_DELTA_LAYERS],
    size_t supplied_layers, size_t *applied_layers, int *used_local);
int revac_uav_memver_q96(const revac_uav_t *uav,
                         const revac_public_state_t *public_state);
int revac_uav_q88_cached_lift_copy(acc96_exact_lift_t *out, unsigned *maximum_abs_bits,
    const revac_uav_t *uav, const revac_public_state_t *public_state);
int revac_uav_memver_cache_prewarm(revac_uav_t *uav,
                                   const revac_public_state_t *public_state);

void revac_uav_apply_update(revac_uav_t *uav,
                            const revac_public_state_t *public_state,
                            const acc_genpp_update_t *upmsg);

/* LOCAL computation path.  It first authenticates the BML-carried DELTA4,
   then performs the local recurrence; sequence gaps are rejected. */
int revac_uav_apply_update_local(revac_uav_t *uav,
                                 const revac_public_state_t *public_state,
                                 const acc_genpp_update_t *upmsg,
                                 size_t *local_composites);

/* Verify and consume a TA DELTA4 packet; optionally fall back to LOCAL when
   the two handles share all four low bits. */
int revac_uav_apply_update_context(revac_uav_t *uav,
                                   const revac_public_state_t *public_state,
                                   const acc_genpp_update_t *upmsg,
                                   const uint8_t *wire,
                                   size_t wire_len,
                                   int allow_local_fallback,
                                   revac_context_profile_t *used_profile,
                                   size_t *local_composites);

/* Every update first receives and verifies its BML-carried DELTA4.  The
   selector then chooses whether the witness itself uses LOCAL or DELTA4.
   fetch_context/fetch_opaque remain ABI-compatible but are no longer used. */
int revac_uav_apply_update_auto(revac_uav_t *uav,
                                const revac_public_state_t *public_state,
                                const acc_genpp_update_t *upmsg,
                                double goodput_mbps,
                                double recurrence_ms,
                                revac_context_fetch_fn fetch_context,
                                void *fetch_opaque,
                                revac_context_profile_t *used_profile,
                                size_t *local_composites,
                                size_t *download_bytes);

int revac_uav_recover_updates(revac_uav_t *uav,
                              const revac_public_state_t *public_state,
                              revac_sync_fetch_fn fetch_sync,
                              void *sync_opaque);

int revac_uav_recover_updates_auto(revac_uav_t *uav,
                                   const revac_public_state_t *public_state,
                                   double goodput_mbps,
                                   double recurrence_ms,
                                   revac_context_fetch_fn fetch_context,
                                   void *fetch_opaque,
                                   revac_sync_fetch_fn fetch_sync,
                                   void *sync_opaque,
                                   size_t *download_bytes,
                                   size_t *local_composites);

int revac_uav_show_prove(revac_show_proof_t *proof,
                         const revac_public_state_t *public_state,
                         revac_uav_t *uav,
                         const uint8_t msg[PARAM_M * PARAM_N / 8],
                         const uint8_t nonce[REVAC_SHOW_NONCE_BYTES]);

int revac_uav_show_prove_wire(uint8_t **wire,
                              size_t *wire_len,
                              const revac_public_state_t *public_state,
                              revac_uav_t *uav,
                              const uint8_t msg[PARAM_M * PARAM_N / 8],
                              const uint8_t nonce[REVAC_SHOW_NONCE_BYTES]);

#endif
