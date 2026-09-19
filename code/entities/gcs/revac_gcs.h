#ifndef REVAC_GCS_H
#define REVAC_GCS_H

#include <pthread.h>

#include "../common/revocable_ac.h"

#define REVAC_GCS_REPLAY_CACHE_ENTRIES 64u
#define REVAC_GCS_MAX_SHOW_WIRE_BYTES (1024u * 1024u)

/* Local admission state for the opt-in parallel verification endpoint.
   It is deliberately not part of revac_gcs_t or any wire object. */
typedef struct {
  pthread_mutex_t lock;
  double tokens;
  double capacity;
  double refill_per_second;
  double last_refill;
  uint8_t replay_nonce[REVAC_GCS_REPLAY_CACHE_ENTRIES]
                      [REVAC_SHOW_NONCE_BYTES];
  size_t replay_count;
  size_t replay_next;
  /* Opaque, process-local immutable signature public-embedding cache. */
  void *signature_cache;
  int initialized;
} revac_gcs_verify_guard_t;

void revac_gcs_init(revac_gcs_t *gcs);
void revac_gcs_clear(revac_gcs_t *gcs);
void revac_gcs_keygen(revac_gcs_t *gcs, uint64_t id, const uint8_t seed[SEED_BYTES], uint64_t timestamp);

int revac_gcs_verify_uav(const revac_public_state_t *public_state,
                         const revac_gcs_t *gcs,
                         const revac_uav_t *uav,
                         const uint8_t msg[PARAM_M * PARAM_N / 8]);

int revac_gcs_verify_show_signature(
                                    const revac_public_state_t *public_state,
                                    const revac_gcs_t *gcs,
                                    const revac_show_context_t *ctx,
                                    const revac_show_proof_t *proof);

int revac_gcs_verify_show(const revac_public_state_t *public_state,
                          const revac_gcs_t *gcs,
                          const revac_show_context_t *ctx,
                          const revac_show_proof_t *proof);

/* Component-only gate for profiling/pipelines that have already accepted
   revac_gcs_verify_show_signature on the same immutable proof.  Applications
   that have not done so must call revac_gcs_verify_show instead. */
int revac_gcs_verify_show_accumulator(
    const revac_public_state_t *public_state,
    const revac_gcs_t *gcs,
    const revac_show_context_t *ctx,
    const revac_show_proof_t *proof);

int revac_gcs_verify_show_signature_wire(
                                         const revac_public_state_t *public_state,
                                         const revac_gcs_t *gcs,
                                         const revac_show_context_t *ctx,
                                         const uint8_t *wire,
                                         size_t wire_len);

int revac_gcs_verify_show_wire(const revac_public_state_t *public_state,
                               const revac_gcs_t *gcs,
                               const revac_show_context_t *ctx,
                               const uint8_t *wire,
                               size_t wire_len);

/* Trusted GCS endpoint: performs cheap structural/state checks, token-bucket
   admission and nonce replay rejection before verifying the independent
   signature and accumulator proof components concurrently.  The guard is
   local operational state; enabling this API does not change proof bytes. */
int revac_gcs_verify_guard_init(revac_gcs_verify_guard_t *guard,
                                double capacity,
                                double refill_per_second);
void revac_gcs_verify_guard_clear(revac_gcs_verify_guard_t *guard);
int revac_gcs_verify_show_wire_guarded_parallel(
    const revac_public_state_t *public_state,
    const revac_gcs_t *gcs,
    const revac_show_context_t *ctx,
    revac_gcs_verify_guard_t *guard,
    const uint8_t *wire,
    size_t wire_len);

int revac_gcs_verify_show_accumulator_wire(
    const revac_public_state_t *public_state,
    const revac_gcs_t *gcs,
    const revac_show_context_t *ctx,
    const uint8_t *wire,
    size_t wire_len);

#endif
