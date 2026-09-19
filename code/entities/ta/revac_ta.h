#ifndef REVAC_TA_H
#define REVAC_TA_H

#include "../common/revocable_ac.h"

void revac_ta_init(revac_ta_t *ta);
void revac_ta_clear(revac_ta_t *ta);
/* Production default: certified q96 accumulator.  Setup can fail while
   constructing/certifying the trapdoor sampler, so callers must check it. */
int revac_ta_setup(revac_ta_t *ta);

/* Explicit compatibility/ablation profile; never selected implicitly. */
void revac_ta_setup_q19(revac_ta_t *ta);

/* Explicit spelling retained for profile-aware deployments. */
int revac_ta_setup_q96(revac_ta_t *ta);
int revac_ta_configure_q96_sampler(revac_ta_t *ta,
                                   acc96_sampleleft_fn sampler,
                                   void *sampler_opaque,
                                   long double witness_norm_sq_bound);
void revac_ta_public_state(revac_public_state_t *public_state,
                           const revac_ta_t *ta);

int revac_ta_register_uav(revac_ta_t *ta, revac_uav_t *uav);
int revac_ta_register_gcs(revac_ta_t *ta, revac_gcs_t *gcs);

int revac_ta_issue_uav(revac_ta_t *ta,
                       revac_uav_t *uav,
                       const uint8_t msg[PARAM_M * PARAM_N / 8]);

int revac_ta_issue_uav_q96(revac_ta_t *ta,
                            revac_uav_t *uav,
                            const uint8_t msg[PARAM_M * PARAM_N / 8]);

int revac_ta_revoke_uav_q96(revac_ta_t *ta,
                             uint32_t handle,
                             acc96_relation_update_t *update);
int revac_ta_export_q96_update(const revac_ta_t *ta, uint32_t sequence,
                               uint8_t **wire, size_t *wire_len);
int revac_ta_export_q96_update_manifest(const revac_ta_t *ta,
                                        uint32_t sequence,
                                        uint8_t **wire, size_t *wire_len);
int revac_ta_export_q96_delta6_manifest(const revac_ta_t *ta,
                                        uint32_t sequence,
                                        uint8_t **wire, size_t *wire_len);
int revac_ta_export_q96_delta6_layer(const revac_ta_t *ta,
                                     uint32_t sequence, size_t layer,
                                     uint8_t **wire, size_t *wire_len);
/* Build an authenticated T4 cache object before the candidate is revoked.
   This is intended for registration/idle-time distribution. */
int revac_ta_export_q96_tail4_prefetch(const revac_ta_t *ta,
                                       uint32_t candidate_handle,
                                       uint8_t **wire, size_t *wire_len);

int revac_ta_revoke_uav(revac_ta_t *ta,
                        uint32_t handle,
                        acc_genpp_update_t *upmsg);

/* Compatibility entry point: every revocation now retains TA-private DELTA6
   and inserts a mandatory authenticated DELTA4 record into the public BML. */
int revac_ta_revoke_uav_context(revac_ta_t *ta,
                                uint32_t handle,
                                acc_genpp_update_t *upmsg);

/* Export the mandatory DELTA4 carried by one public BML record. */
int revac_ta_export_context(const revac_ta_t *ta,
                            uint32_t seq,
                            revac_context_profile_t profile,
                            uint8_t **wire,
                            size_t *wire_len);

/* Reconstruct a current member witness while reusing private DELTA6 entries.
   Missing, invalid, collected, or unsupported entries transparently fall
   back to the original local recurrence.
   REVAC_TA_DISABLE_SYNC_DELTA_CACHE=1 forces that fallback for ablation
   tests; the former REVAC_TA_DISABLE_TAIL_CACHE name remains accepted. */
int revac_ta_memwitsync_cached(acc_genpp_wit_t *witness,
                               const revac_ta_t *ta,
                               uint32_t handle,
                               size_t *cached_updates,
                               size_t *local_composites);

/* Drop only private delta entries at or below the minimum issuance/checkpoint
   baseline of all active, non-revoked members. */
size_t revac_ta_sync_delta_cache_gc(revac_ta_t *ta);

int revac_ta_export_memwitsync(const revac_ta_t *ta,
                               /* Identity established by the secure session. */
                               uint32_t authenticated_handle,
                               const revac_sync_request_t *request,
                               uint8_t **wire,
                               size_t *wire_len);

#endif
