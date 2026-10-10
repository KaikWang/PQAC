#include "revac_gcs.h"

#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "revac_labrador_bridge.h"
#include "revac_proof_codec.h"
#include "revac_uav.h"
#include "poly_q_sampling.h"

typedef struct {
  const revac_public_state_t *public_state;
  const revac_gcs_t *gcs;
  const revac_show_context_t *ctx;
  const revac_show_proof_t *proof;
  const void *signature_cache;
  int result;
} revac_gcs_verify_task_t;

typedef struct {
  const revac_public_state_t *public_state;
  const sep_pk_t *ta_pk;
  poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D];
  poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K];
  poly_qshow_mat_k_k A3_embed[PARAM_D];
  poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D];
  poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M];
  poly_qshow_mat_k_k Dx_embed[PARAM_D];
  poly_qshow_vec_k u_embed[PARAM_D];
  poly_qshow_mat_k_k A_link_embed[PARAM_D][PARAM_D];
  poly_qshow_mat_k_k D_link_embed[PARAM_D];
} revac_gcs_signature_cache_t;

static double revac_gcs_monotonic_seconds(void) {
  struct timespec now;

  if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) {
    return 0.0;
  }
  return (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
}

static unsigned int revac_gcs_parallel_signature_threads(void) {
  const char *value = getenv("REVAC_SHOW_PARALLEL_VERIFY_THREADS");

  if (value != NULL && value[0] != '\0') {
    char *end = NULL;
    long requested = strtol(value, &end, 10);
    if (end != value && *end == '\0' && requested >= 1 && requested <= 64) {
      return (unsigned int)requested;
    }
  }
  return 4;
}

static int revac_gcs_show_prefilter(const revac_public_state_t *public_state,
                                    const revac_gcs_t *gcs,
                                    const revac_show_context_t *ctx,
                                    const revac_show_proof_t *proof) {
  return public_state != NULL &&
         (public_state->acc_profile == REVAC_ACC_PROFILE_Q19 ||
          public_state->acc_profile == REVAC_ACC_PROFILE_Q96) &&
         gcs != NULL && ctx != NULL && proof != NULL &&
         gcs->registered && gcs->cert_valid == REVAC_GCS_CERT_VALID &&
         public_state->ta_pk != NULL &&
         ((public_state->acc_profile == REVAC_ACC_PROFILE_Q19 &&
           public_state->acc != NULL) ||
          (public_state->acc_profile == REVAC_ACC_PROFILE_Q96 &&
           public_state->acc96_acc != NULL)) &&
         revac_show_context_matches_public(ctx, public_state) &&
         (public_state->acc_profile == REVAC_ACC_PROFILE_Q96 ||
          poly_q_vec_d_equal(ctx->acc, public_state->acc)) &&
         ((public_state->acc_profile == REVAC_ACC_PROFILE_Q19 &&
           (proof->acc_zk_kind == REVAC_ACC_ZK_LABRADOR_MEMVER ||
            proof->acc_zk_kind == REVAC_ACC_ZK_LABRADOR_MEMVER_DEFLATE)) ||
          (public_state->acc_profile == REVAC_ACC_PROFILE_Q96 &&
           (proof->acc_zk_kind == REVAC_ACC_ZK_Q96_LABRADOR_MEMVER ||
            proof->acc_zk_kind == REVAC_ACC_ZK_Q96_LABRADOR_MEMVER_DEFLATE))) &&
         proof->acc_zk_proof != NULL && proof->acc_zk_proof_len != 0 &&
         proof->acc_zk_proof_len <= REVAC_GCS_MAX_SHOW_WIRE_BYTES;
}

static int revac_gcs_guard_admit(revac_gcs_verify_guard_t *guard,
                                 const uint8_t nonce[REVAC_SHOW_NONCE_BYTES]) {
  const double now = revac_gcs_monotonic_seconds();
  int admitted = 0;

  if (guard == NULL || !guard->initialized ||
      pthread_mutex_lock(&guard->lock) != 0) {
    return 0;
  }
  if (now > guard->last_refill) {
    guard->tokens += (now - guard->last_refill) * guard->refill_per_second;
    if (guard->tokens > guard->capacity) {
      guard->tokens = guard->capacity;
    }
    guard->last_refill = now;
  }
  for (size_t i = 0; i < guard->replay_count; i++) {
    if (memcmp(guard->replay_nonce[i], nonce, REVAC_SHOW_NONCE_BYTES) == 0) {
      goto out;
    }
  }
  if (guard->tokens < 1.0) {
    goto out;
  }
  guard->tokens -= 1.0;
  memcpy(guard->replay_nonce[guard->replay_next], nonce,
         REVAC_SHOW_NONCE_BYTES);
  guard->replay_next =
      (guard->replay_next + 1u) % REVAC_GCS_REPLAY_CACHE_ENTRIES;
  if (guard->replay_count < REVAC_GCS_REPLAY_CACHE_ENTRIES) {
    guard->replay_count++;
  }
  admitted = 1;
out:
  (void)pthread_mutex_unlock(&guard->lock);
  return admitted;
}

static void revac_gcs_signature_cache_free(
    revac_gcs_signature_cache_t *cache) {
  if (cache == NULL) {
    return;
  }
  for (size_t row = 0; row < PARAM_D; row++) {
    for (size_t col = 0; col < PARAM_D; col++) {
      poly_qshow_mat_k_k_clear(cache->A_link_embed[row][col]);
    }
    poly_qshow_mat_k_k_clear(cache->D_link_embed[row]);
  }
  revac_show_lang_clear(cache->A_embed, cache->B_embed, cache->A3_embed,
                        cache->Ds_embed, cache->D_embed, cache->Dx_embed,
                        cache->u_embed);
  free(cache);
}

static revac_gcs_signature_cache_t *revac_gcs_signature_cache_build(
    const revac_public_state_t *public_state,
    const revac_show_context_t *ctx) {
  revac_gcs_signature_cache_t *cache;
  poly_q_mat_d_d A_link;
  poly_q_mat_d_d D_link;

  cache = calloc(1, sizeof(*cache));
  if (cache == NULL) {
    return NULL;
  }
  cache->public_state = public_state;
  cache->ta_pk = public_state->ta_pk;
  revac_show_lang_init(cache->A_embed, cache->B_embed, cache->A3_embed,
                       cache->Ds_embed, cache->D_embed, cache->Dx_embed,
                       cache->u_embed);
  poly_q_mat_d_d_init(A_link);
  poly_q_mat_d_d_init(D_link);
  revac_same_handle_public_matrices(A_link, D_link, public_state);
  for (size_t row = 0; row < PARAM_D; row++) {
    poly_qshow_mat_k_k_init(cache->D_link_embed[row]);
    poly_qshow_subring_embed_mat_k_k(cache->D_link_embed[row],
                                    D_link->rows[row]->entries[0],
                                    PARAM_Q1_SHOW);
    for (size_t col = 0; col < PARAM_D; col++) {
      poly_qshow_mat_k_k_init(cache->A_link_embed[row][col]);
      poly_qshow_subring_embed_mat_k_k(cache->A_link_embed[row][col],
                                      A_link->rows[row]->entries[col],
                                      PARAM_Q1_SHOW);
    }
  }
  revac_show_public_embed(cache->A_embed, cache->B_embed, cache->A3_embed,
                          cache->Ds_embed, cache->D_embed, cache->Dx_embed,
                          cache->u_embed, public_state->ta_pk, ctx->user_seed);
  poly_q_mat_d_d_clear(D_link);
  poly_q_mat_d_d_clear(A_link);
  return cache;
}

static const revac_gcs_signature_cache_t *revac_gcs_signature_cache_get(
    revac_gcs_verify_guard_t *guard,
    const revac_public_state_t *public_state,
    const revac_show_context_t *ctx) {
  revac_gcs_signature_cache_t *cache;

  if (pthread_mutex_lock(&guard->lock) != 0) {
    return NULL;
  }
  cache = (revac_gcs_signature_cache_t *)guard->signature_cache;
  if (cache == NULL) {
    cache = revac_gcs_signature_cache_build(public_state, ctx);
    guard->signature_cache = cache;
  } else if (cache->public_state != public_state ||
             cache->ta_pk != public_state->ta_pk) {
    cache = NULL;
  }
  (void)pthread_mutex_unlock(&guard->lock);
  return cache;
}

static int revac_gcs_verify_show_signature_cached(
    const revac_public_state_t *public_state,
    const revac_show_context_t *ctx,
    const revac_show_proof_t *proof,
    const revac_gcs_signature_cache_t *cache) {
  poly_qshow_vec_k C_link_embed[PARAM_D];
  poly_qshow_mat_k_k Dx_embed[PARAM_D];
  poly_q_mat_d_d Dx;
  uint8_t expected_crs[CRS_SEED_BYTES];
  int ok;

  if (cache == NULL) {
    return 0;
  }
  revac_show_derive_crs_seed(expected_crs, public_state->ta_pk, ctx->user_seed,
                             proof->nonce, ctx->acc_profile, ctx->acc_seq,
                             ctx->window_id, ctx->parameter_digest,
                             ctx->acc_digest);
  /* Dx is the only public signature matrix keyed by the per-user seed.  Keep
     the much larger TA/static embedding cache shareable across all UAVs and
     expand only these four matrices per request. */
  poly_q_mat_d_d_init(Dx);
  poly_q_mat_d_d_uniform(Dx, ctx->user_seed, DOMAIN_SEPARATOR_DX, 0);
  for (size_t row = 0; row < PARAM_D; row++) {
    poly_qshow_vec_k_init(C_link_embed[row]);
    poly_qshow_mat_k_k_init(Dx_embed[row]);
    poly_qshow_subring_embed_vec_k(C_link_embed[row],
                                  proof->handle_commitment->entries[row],
                                  PARAM_Q1_SHOW);
    poly_qshow_subring_embed_mat_k_k(Dx_embed[row],
                                    Dx->rows[row]->entries[0],
                                    PARAM_Q1_SHOW);
  }
  ok = show_verify_handle(&proof->sig_show, cache->A_embed, cache->B_embed,
                          cache->A3_embed, cache->Ds_embed, cache->D_embed,
                          Dx_embed, cache->A_link_embed,
                          cache->D_link_embed, C_link_embed, cache->u_embed,
                          expected_crs, ctx->user_seed);
  for (size_t row = 0; row < PARAM_D; row++) {
    poly_qshow_mat_k_k_clear(Dx_embed[row]);
    poly_qshow_vec_k_clear(C_link_embed[row]);
  }
  poly_q_mat_d_d_clear(Dx);
  return ok;
}

static void *revac_gcs_signature_verify_thread(void *opaque) {
  revac_gcs_verify_task_t *task = (revac_gcs_verify_task_t *)opaque;

  show_verify_set_thread_cap(revac_gcs_parallel_signature_threads());
  if (task->signature_cache != NULL) {
    task->result = revac_gcs_verify_show_signature_cached(
        task->public_state, task->ctx, task->proof,
        (const revac_gcs_signature_cache_t *)task->signature_cache);
  } else {
    task->result = revac_gcs_verify_show_signature(
        task->public_state, task->gcs, task->ctx, task->proof);
  }
  show_verify_set_thread_cap(0);
  return NULL;
}

void revac_gcs_init(revac_gcs_t *gcs) {
  memset(gcs, 0, sizeof(*gcs));
  user_keys_init(&gcs->pk, &gcs->sk);
}

void revac_gcs_clear(revac_gcs_t *gcs) {
  user_keys_clear(&gcs->pk, &gcs->sk);
}

void revac_gcs_keygen(revac_gcs_t *gcs, uint64_t id, const uint8_t seed[SEED_BYTES], uint64_t timestamp) {
  gcs->id = id;
  gcs->timestamp = timestamp;
  osig_user_keygen(&gcs->pk, &gcs->sk, seed);
}

int revac_gcs_verify_uav(const revac_public_state_t *public_state,
                         const revac_gcs_t *gcs,
                         const revac_uav_t *uav,
                         const uint8_t msg[PARAM_M * PARAM_N / 8]) {
  if (!gcs->registered || gcs->cert_valid != REVAC_GCS_CERT_VALID || !uav->issued) {
    return 0;
  }
  if (public_state == NULL || public_state->ta_pk == NULL ||
      !osig_user_verify_handle(&uav->sig, public_state->ta_pk, &uav->pk,
                               msg, uav->handle)) {
    return 0;
  }
  if (public_state->acc_profile == REVAC_ACC_PROFILE_Q96) {
    return revac_uav_memver_q96(uav, public_state);
  }
  if (public_state->acc_profile != REVAC_ACC_PROFILE_Q19 ||
      public_state->acc_pp == NULL || public_state->acc == NULL) {
    return 0;
  }
  return acc_genpp_memver(public_state->acc_pp, public_state->acc,
                          uav->handle, &uav->witness);
}

int revac_gcs_verify_show_signature(
                                    const revac_public_state_t *public_state,
                                    const revac_gcs_t *gcs,
                                    const revac_show_context_t *ctx,
                                    const revac_show_proof_t *proof) {
  poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D];
  poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K];
  poly_qshow_mat_k_k A3_embed[PARAM_D];
  poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D];
  poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M];
  poly_qshow_mat_k_k Dx_embed[PARAM_D];
  poly_qshow_mat_k_k A_link_embed[PARAM_D][PARAM_D];
  poly_qshow_mat_k_k D_link_embed[PARAM_D];
  poly_qshow_vec_k C_link_embed[PARAM_D];
  poly_qshow_vec_k u_embed[PARAM_D];
  poly_q_mat_d_d A_link, D_link;
  uint8_t expected_crs[CRS_SEED_BYTES];
  int ok;

  if (!gcs->registered || gcs->cert_valid != REVAC_GCS_CERT_VALID ||
      public_state == NULL ||
      (public_state->acc_profile != REVAC_ACC_PROFILE_Q19 &&
       public_state->acc_profile != REVAC_ACC_PROFILE_Q96)) {
    return 0;
  }
  if (public_state->ta_pk == NULL ||
      (public_state->acc_profile == REVAC_ACC_PROFILE_Q19 &&
       public_state->acc == NULL) ||
      (public_state->acc_profile == REVAC_ACC_PROFILE_Q96 &&
       public_state->acc96_acc == NULL) ||
      !revac_show_context_matches_public(ctx, public_state) ||
      (public_state->acc_profile == REVAC_ACC_PROFILE_Q19 &&
       !poly_q_vec_d_equal(ctx->acc, public_state->acc))) {
    return 0;
  }
  revac_show_derive_crs_seed(expected_crs, public_state->ta_pk, ctx->user_seed,
                             proof->nonce, ctx->acc_profile, ctx->acc_seq,
                             ctx->window_id, ctx->parameter_digest,
                             ctx->acc_digest);

  revac_show_lang_init(A_embed, B_embed, A3_embed, Ds_embed, D_embed, Dx_embed, u_embed);
  poly_q_mat_d_d_init(A_link);
  poly_q_mat_d_d_init(D_link);
  for (size_t i = 0; i < PARAM_D; i++) {
    poly_qshow_mat_k_k_init(D_link_embed[i]);
    poly_qshow_vec_k_init(C_link_embed[i]);
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_init(A_link_embed[i][j]);
    }
  }
  revac_same_handle_public_matrices(A_link, D_link, public_state);
  show_link_public_embed(A_link_embed, D_link_embed, C_link_embed,
                         A_link, D_link, proof->handle_commitment);
  revac_show_public_embed(A_embed, B_embed, A3_embed, Ds_embed, D_embed,
                          Dx_embed, u_embed, public_state->ta_pk,
                          ctx->user_seed);
  ok = show_verify_handle(&proof->sig_show, A_embed, B_embed, A3_embed,
                          Ds_embed, D_embed, Dx_embed, A_link_embed,
                          D_link_embed, C_link_embed, u_embed,
                          expected_crs, ctx->user_seed);
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_clear(A_link_embed[i][j]);
    }
    poly_qshow_vec_k_clear(C_link_embed[i]);
    poly_qshow_mat_k_k_clear(D_link_embed[i]);
  }
  poly_q_mat_d_d_clear(D_link);
  poly_q_mat_d_d_clear(A_link);
  revac_show_lang_clear(A_embed, B_embed, A3_embed, Ds_embed, D_embed, Dx_embed, u_embed);
  return ok;
}

int revac_gcs_verify_show(const revac_public_state_t *public_state,
                          const revac_gcs_t *gcs,
                          const revac_show_context_t *ctx,
                          const revac_show_proof_t *proof) {
  if (!revac_gcs_verify_show_signature(public_state, gcs, ctx, proof)) {
    return 0;
  }
  return revac_gcs_verify_show_accumulator(public_state, gcs, ctx, proof);
}

int revac_gcs_verify_show_accumulator(
    const revac_public_state_t *public_state,
    const revac_gcs_t *gcs,
    const revac_show_context_t *ctx,
    const revac_show_proof_t *proof) {
  if (public_state == NULL || gcs == NULL || ctx == NULL || proof == NULL ||
      (public_state->acc_profile != REVAC_ACC_PROFILE_Q19 &&
       public_state->acc_profile != REVAC_ACC_PROFILE_Q96) ||
      !gcs->registered || gcs->cert_valid != REVAC_GCS_CERT_VALID ||
      (public_state->acc_profile == REVAC_ACC_PROFILE_Q19 &&
       public_state->acc == NULL) ||
      (public_state->acc_profile == REVAC_ACC_PROFILE_Q96 &&
       public_state->acc96_acc == NULL) ||
      !revac_show_context_matches_public(ctx, public_state) ||
      (public_state->acc_profile == REVAC_ACC_PROFILE_Q19 &&
       !poly_q_vec_d_equal(ctx->acc, public_state->acc))) {
    return 0;
  }
  if (public_state->acc_profile == REVAC_ACC_PROFILE_Q96) {
    if ((proof->acc_zk_kind != REVAC_ACC_ZK_Q96_LABRADOR_MEMVER &&
         proof->acc_zk_kind != REVAC_ACC_ZK_Q96_LABRADOR_MEMVER_DEFLATE) ||
        proof->acc_zk_proof == NULL || proof->acc_zk_proof_len == 0)
      return 0;
    return revac_labrador_q96_memver_verify(public_state, ctx, proof);
  }
  if ((proof->acc_zk_kind != REVAC_ACC_ZK_LABRADOR_MEMVER &&
       proof->acc_zk_kind != REVAC_ACC_ZK_LABRADOR_MEMVER_DEFLATE) ||
      proof->acc_zk_proof == NULL ||
      proof->acc_zk_proof_len == 0) {
    return 0;
  }
  return revac_labrador_memver_verify(public_state, ctx, proof);
}

int revac_gcs_verify_show_signature_wire(
                                         const revac_public_state_t *public_state,
                                         const revac_gcs_t *gcs,
                                         const revac_show_context_t *ctx,
                                         const uint8_t *wire,
                                         size_t wire_len) {
  revac_show_proof_t proof;
  int ok = 0;

  if (wire == NULL || wire_len == 0) {
    return 0;
  }
  revac_show_proof_init(&proof);
  if (revac_show_proof_unpack_compact(&proof, wire, wire_len)) {
    ok = revac_gcs_verify_show_signature(public_state, gcs, ctx, &proof);
  }
  revac_show_proof_clear(&proof);
  return ok;
}

int revac_gcs_verify_show_wire(const revac_public_state_t *public_state,
                               const revac_gcs_t *gcs,
                               const revac_show_context_t *ctx,
                               const uint8_t *wire,
                               size_t wire_len) {
  revac_show_proof_t proof;
  int ok = 0;

  if (wire == NULL || wire_len == 0) {
    return 0;
  }
  revac_show_proof_init(&proof);
  if (revac_show_proof_unpack_compact(&proof, wire, wire_len)) {
    ok = revac_gcs_verify_show(public_state, gcs, ctx, &proof);
  }
  revac_show_proof_clear(&proof);
  return ok;
}

int revac_gcs_verify_guard_init(revac_gcs_verify_guard_t *guard,
                                double capacity,
                                double refill_per_second) {
  if (guard == NULL || capacity < 1.0 || refill_per_second <= 0.0) {
    return 0;
  }
  memset(guard, 0, sizeof(*guard));
  if (pthread_mutex_init(&guard->lock, NULL) != 0) {
    return 0;
  }
  guard->tokens = capacity;
  guard->capacity = capacity;
  guard->refill_per_second = refill_per_second;
  guard->last_refill = revac_gcs_monotonic_seconds();
  guard->initialized = 1;
  return 1;
}

void revac_gcs_verify_guard_clear(revac_gcs_verify_guard_t *guard) {
  if (guard == NULL) {
    return;
  }
  if (guard->initialized) {
    revac_gcs_signature_cache_free(
        (revac_gcs_signature_cache_t *)guard->signature_cache);
    guard->signature_cache = NULL;
    (void)pthread_mutex_destroy(&guard->lock);
  }
  memset(guard, 0, sizeof(*guard));
}

int revac_gcs_verify_show_wire_guarded_parallel(
    const revac_public_state_t *public_state,
    const revac_gcs_t *gcs,
    const revac_show_context_t *ctx,
    revac_gcs_verify_guard_t *guard,
    const uint8_t *wire,
    size_t wire_len) {
  revac_show_proof_t proof;
  revac_gcs_verify_task_t signature_task;
  pthread_t signature_thread;
  int accumulator_ok;
  int ok = 0;

  /* Reject oversized/untrusted input before allocation or worker dispatch. */
  if (wire == NULL || wire_len == 0 ||
      wire_len > REVAC_GCS_MAX_SHOW_WIRE_BYTES) {
    return 0;
  }
  revac_show_proof_init(&proof);
  if (!revac_show_proof_unpack_compact(&proof, wire, wire_len) ||
      !revac_gcs_show_prefilter(public_state, gcs, ctx, &proof) ||
      !revac_gcs_guard_admit(guard, proof.nonce)) {
    goto out;
  }
  signature_task.public_state = public_state;
  signature_task.gcs = gcs;
  signature_task.ctx = ctx;
  signature_task.proof = &proof;
  signature_task.signature_cache =
      revac_gcs_signature_cache_get(guard, public_state, ctx);
  signature_task.result = 0;
  if (pthread_create(&signature_thread, NULL,
                     revac_gcs_signature_verify_thread,
                     &signature_task) == 0) {
    /* Admission already performed the GCS/state/type checks.  Calling the
       bridge directly avoids repeating accumulator equality and proof-shape
       scans on the accepted hot path. */
    if (public_state->acc_profile == REVAC_ACC_PROFILE_Q96) {
      accumulator_ok =
          revac_labrador_q96_memver_verify(public_state, ctx, &proof);
    } else {
      accumulator_ok = revac_labrador_memver_verify(public_state, ctx, &proof);
    }
    (void)pthread_join(signature_thread, NULL);
    ok = accumulator_ok && signature_task.result;
  } else {
    /* Resource exhaustion remains correct: use the established serial path. */
    ok = revac_gcs_verify_show(public_state, gcs, ctx, &proof);
  }
out:
  revac_show_proof_clear(&proof);
  return ok;
}

int revac_gcs_verify_show_accumulator_wire(
    const revac_public_state_t *public_state,
    const revac_gcs_t *gcs,
    const revac_show_context_t *ctx,
    const uint8_t *wire,
    size_t wire_len) {
  revac_show_proof_t proof;
  int ok = 0;

  if (wire == NULL || wire_len == 0) {
    return 0;
  }
  revac_show_proof_init(&proof);
  if (revac_show_proof_unpack_compact(&proof, wire, wire_len)) {
    ok = revac_gcs_verify_show_accumulator(public_state, gcs, ctx, &proof);
  }
  revac_show_proof_clear(&proof);
  return ok;
}
