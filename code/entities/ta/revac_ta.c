#include "revac_ta.h"

#include <stdlib.h>
#include <string.h>

#include "randombytes.h"
#include "fips202.h"
#include "acc96_perturbation_cert.h"
#include "acc96_sampleleft.h"
#include "acc96_trapdoor_cert.h"
#include "acc96_truncated_params_generated.h"
#include "revac_context.h"
#include "revac_sync.h"

static void revac_copy_user_pk(user_pk_t *dst, const user_pk_t *src) {
  poly_q_vec_d_set(dst->t, src->t);
  memcpy(dst->seed, src->seed, sizeof(dst->seed));
}

static void revac_ta_acc96_secret_free(acc96_accumulator_t **state) {
  volatile unsigned char *bytes;

  if (state == NULL || *state == NULL) {
    return;
  }
  bytes = (volatile unsigned char *)(void *)*state;
  for (size_t i = 0; i < sizeof(**state); i++) {
    bytes[i] = 0;
  }
  free(*state);
  *state = NULL;
}

static int revac_find_uav_by_handle(const revac_ta_t *ta, uint32_t handle) {
  for (size_t i = 0; i < REVAC_MAX_UAVS; i++) {
    if (ta->uavs[i].used && ta->uavs[i].handle == handle) {
      return (int)i;
    }
  }
  return -1;
}

static const revac_bml_entry_t *revac_find_bml_entry_by_seq(
    const revac_ta_t *ta,
    uint32_t seq) {
  for (size_t i = 0; i < ACC_GENPP_BML_CAPACITY; i++) {
    const revac_bml_entry_t *entry = &ta->bml[i];

    if (entry->used && entry->update.seq == seq) {
      return entry;
    }
  }
  return NULL;
}

static const revac_sync_delta_cache_entry_t *revac_find_sync_delta_by_seq(
    const revac_ta_t *ta, uint32_t seq) {
  const revac_sync_delta_cache_entry_t *entry;

  if (seq == 0 || seq > ACC_GENPP_MAX_DELETIONS) {
    return NULL;
  }
  entry = &ta->sync_delta_cache[seq - 1];
  return entry->used && entry->seq == seq ? entry : NULL;
}

typedef struct {
  const revac_ta_t *ta;
  size_t cached_updates;
  size_t local_composites;
} revac_ta_sync_cache_t;

static int revac_ta_apply_cached_update(
    void *opaque,
    acc_genpp_wit_t *witness,
    const acc_genpp_pp_t *pp,
    uint32_t member_handle,
    const acc_genpp_update_t *update) {
  revac_ta_sync_cache_t *sync = opaque;
  const char *disable_cache = getenv("REVAC_TA_DISABLE_SYNC_DELTA_CACHE");
  const char *disable_legacy = getenv("REVAC_TA_DISABLE_TAIL_CACHE");
  const int cache_enabled =
      (disable_cache == NULL || disable_cache[0] == '\0' ||
       strcmp(disable_cache, "0") == 0) &&
      (disable_legacy == NULL || disable_legacy[0] == '\0' ||
       strcmp(disable_legacy, "0") == 0);
  const revac_sync_delta_cache_entry_t *entry =
      cache_enabled
          ? revac_find_sync_delta_by_seq(sync->ta, update->seq)
          : NULL;
  size_t composites = 0;
  int ok = 0;

  if (entry != NULL && entry->handle == update->handle) {
    ok = revac_sync_delta_cache_entry_apply(
        witness, entry, sync->ta->context_pp_hash, member_handle,
        update, &composites);
    if (ok) {
      sync->cached_updates++;
    }
  }

  if (!ok) {
    composites = 0;
    ok = acc_genpp_memwitup_tails(witness, pp, member_handle, update,
                                  NULL, NULL, 0, &composites);
  }
  sync->local_composites += composites;
  return ok;
}

size_t revac_ta_sync_delta_cache_gc(revac_ta_t *ta) {
  uint32_t watermark;
  int have_active_member = 0;
  size_t removed = 0;

  if (ta == NULL) {
    return 0;
  }
  watermark = ta->acc_st.ctr;
  for (size_t i = 0; i < REVAC_MAX_UAVS; i++) {
    const revac_uav_record_t *record = &ta->uavs[i];

    if (!record->used || !record->issued || record->revoked) {
      continue;
    }
    if (!have_active_member || record->baseline_seq < watermark) {
      watermark = record->baseline_seq;
    }
    have_active_member = 1;
  }
  if (!have_active_member) {
    watermark = ta->acc_st.ctr;
  }
  for (size_t i = 0; i < ACC_GENPP_MAX_DELETIONS; i++) {
    revac_sync_delta_cache_entry_t *entry = &ta->sync_delta_cache[i];

    if (entry->used && entry->seq <= watermark) {
      revac_sync_delta_cache_entry_clear(entry);
      removed++;
    }
  }
  return removed;
}

int revac_ta_memwitsync_cached(acc_genpp_wit_t *witness,
                               const revac_ta_t *ta,
                               uint32_t handle,
                               size_t *cached_updates,
                               size_t *local_composites) {
  revac_ta_sync_cache_t sync;
  int ok;

  if (witness == NULL || ta == NULL ||
      ta->acc_profile != REVAC_ACC_PROFILE_Q19) {
    return 0;
  }
  sync.ta = ta;
  sync.cached_updates = 0;
  sync.local_composites = 0;
  ok = acc_genpp_memwitsync_apply(
      witness, &ta->acc_pp, &ta->acc_td, &ta->acc_st, ta->acc, handle,
      revac_ta_apply_cached_update, &sync);
  if (cached_updates != NULL) {
    *cached_updates = sync.cached_updates;
  }
  if (local_composites != NULL) {
    *local_composites = sync.local_composites;
  }
  return ok;
}

void revac_ta_init(revac_ta_t *ta) {
  sep_keys_init(&ta->ta_pk, &ta->ta_sk);
  acc_genpp_pp_init(&ta->acc_pp);
  acc_genpp_td_init(&ta->acc_td);
  acc_genpp_state_init(&ta->acc_st);
  poly_q_vec_d_init(ta->acc);
  memset(ta->context_pp_hash, 0, sizeof(ta->context_pp_hash));
  memset(ta->bml, 0, sizeof(ta->bml));
  memset(ta->q96_bml, 0, sizeof(ta->q96_bml));
  memset(ta->q96_delta6_bml, 0, sizeof(ta->q96_delta6_bml));
  memset(ta->sync_delta_cache, 0, sizeof(ta->sync_delta_cache));
  memset(ta->uavs, 0, sizeof(ta->uavs));
  ta->acc_profile = REVAC_ACC_PROFILE_Q19;
  memset(&ta->acc96_ctx, 0, sizeof(ta->acc96_ctx));
  memset(ta->acc96_public_seed, 0, sizeof(ta->acc96_public_seed));
  ta->acc96_window_id = 0;
  ta->acc96_state = NULL;
  ta->acc96_sampler = NULL;
  ta->acc96_sampler_opaque = NULL;
  ta->acc96_secure_sampler = NULL;
  ta->acc96_witness_norm_sq_bound = 0.0L;
  for (size_t i = 0; i < REVAC_MAX_UAVS; i++) {
    poly_q_vec_d_init(ta->uavs[i].pk.t);
    memset(ta->uavs[i].pk.seed, 0, sizeof(ta->uavs[i].pk.seed));
  }
}

void revac_ta_clear(revac_ta_t *ta) {
  acc96_sampleleft_ctx_free(ta->acc96_secure_sampler);
  ta->acc96_secure_sampler = NULL;
  revac_ta_acc96_secret_free(&ta->acc96_state);
  memset(ta->acc96_public_seed, 0, sizeof(ta->acc96_public_seed));
  for (size_t i = 0; i < ACC_GENPP_BML_CAPACITY; i++) {
    revac_bml_entry_clear(&ta->bml[i]);
    revac_q96_bml_entry_clear(&ta->q96_bml[i]);
    revac_q96_delta6_entry_clear(&ta->q96_delta6_bml[i]);
  }
  for (size_t i = 0; i < ACC_GENPP_MAX_DELETIONS; i++) {
    revac_sync_delta_cache_entry_clear(&ta->sync_delta_cache[i]);
  }
  for (size_t i = 0; i < REVAC_MAX_UAVS; i++) {
    poly_q_vec_d_clear(ta->uavs[i].pk.t);
  }
  poly_q_vec_d_clear(ta->acc);
  acc_genpp_state_clear(&ta->acc_st);
  acc_genpp_td_clear(&ta->acc_td);
  acc_genpp_pp_clear(&ta->acc_pp);
  sep_keys_clear(&ta->ta_pk, &ta->ta_sk);
}

void revac_ta_setup_q19(revac_ta_t *ta) {
  acc96_sampleleft_ctx_free(ta->acc96_secure_sampler);
  ta->acc96_secure_sampler = NULL;
  revac_ta_acc96_secret_free(&ta->acc96_state);
  ta->acc96_sampler = NULL;
  ta->acc96_sampler_opaque = NULL;
  ta->acc96_witness_norm_sq_bound = 0.0L;
  memset(ta->acc96_public_seed, 0, sizeof(ta->acc96_public_seed));
  ta->acc96_window_id = 0;
  for (size_t i = 0; i < ACC_GENPP_BML_CAPACITY; i++) {
    revac_bml_entry_clear(&ta->bml[i]);
    revac_q96_bml_entry_clear(&ta->q96_bml[i]);
    revac_q96_delta6_entry_clear(&ta->q96_delta6_bml[i]);
  }
  for (size_t i = 0; i < ACC_GENPP_MAX_DELETIONS; i++) {
    revac_sync_delta_cache_entry_clear(&ta->sync_delta_cache[i]);
  }
  sep_keygen(&ta->ta_pk, &ta->ta_sk);
  acc_genpp_gen(&ta->acc_pp, &ta->acc_td, &ta->acc_st, ta->acc);
  revac_context_compute_pp_hash(ta->context_pp_hash, &ta->acc_pp);
  ta->acc_profile = REVAC_ACC_PROFILE_Q19;
}

int revac_ta_configure_q96_sampler(revac_ta_t *ta,
                                   acc96_sampleleft_fn sampler,
                                   void *sampler_opaque,
                                   long double witness_norm_sq_bound) {
  if (ta == NULL || sampler == NULL || witness_norm_sq_bound <= 0.0L) {
    return 0;
  }
  acc96_sampleleft_ctx_free(ta->acc96_secure_sampler);
  ta->acc96_secure_sampler = NULL;
  ta->acc96_sampler = sampler;
  ta->acc96_sampler_opaque = sampler_opaque;
  ta->acc96_witness_norm_sq_bound = witness_norm_sq_bound;
  return 1;
}

int revac_ta_setup_q96(revac_ta_t *ta) {
  static const uint8_t domain[] = "REVAC-Q96-ENTITY-PUBLIC-V1";
  uint8_t public_seed[ACC96_SEED_BYTES];
  uint8_t trapdoor_seed[ACC96_SEED_BYTES];
  uint8_t accumulator_seed[ACC96_SEED_BYTES];
  acc96_accumulator_t *state;
  acc96_sampleleft_ctx_t *secure_sampler = NULL;
  keccak_state hash;
  double certified_upper = 0.0;
  acc96_perturbation_covariance_cert_t covariance_certificate;

  if (ta == NULL || !acc96_ctx_init(&ta->acc96_ctx)) {
    return 0;
  }
  /* A failed re-setup must not leave a stale seed paired with new/partial
     q96 state. */
  memset(ta->acc96_public_seed, 0, sizeof(ta->acc96_public_seed));
  /* Keep credential-key generation in the same setup order as the q19
     production profile.  The later q96 setup/certification is independent. */
  sep_keygen(&ta->ta_pk, &ta->ta_sk);
  state = calloc(1, sizeof(*state));
  if (state == NULL) {
    return 0;
  }
  randombytes(public_seed, sizeof(public_seed));
  randombytes(trapdoor_seed, sizeof(trapdoor_seed));
  randombytes(accumulator_seed, sizeof(accumulator_seed));
  if (!acc96_accgen(state, public_seed, trapdoor_seed, accumulator_seed,
                    &ta->acc96_ctx) ||
      !acc96_trapdoor_spectral_certify_acb(
          &state->trapdoor, &ta->acc96_ctx, 160, &certified_upper) ||
      !acc96_perturbation_covariance_certify_acb(
          &state->trapdoor, &ta->acc96_ctx, 192,
          &covariance_certificate)) {
    revac_ta_acc96_secret_free(&state);
    return 0;
  }
  secure_sampler = acc96_sampleleft_ctx_new(
      &state->trapdoor, &ta->acc96_ctx, 384, 0);
  if (secure_sampler == NULL) {
    revac_ta_acc96_secret_free(&state);
    return 0;
  }

  acc96_sampleleft_ctx_free(ta->acc96_secure_sampler);
  ta->acc96_secure_sampler = NULL;
  revac_ta_acc96_secret_free(&ta->acc96_state);
  ta->acc96_state = state;
  ta->acc96_secure_sampler = secure_sampler;
  ta->acc96_sampler = acc96_sampleleft_secure;
  ta->acc96_sampler_opaque = secure_sampler;
  ta->acc96_witness_norm_sq_bound = ACC96_PUBLIC_WITNESS_NORM_SQ_BOUND;
  memcpy(ta->acc96_public_seed, public_seed, sizeof(ta->acc96_public_seed));
  ta->acc96_window_id = 0;
  ta->acc_profile = REVAC_ACC_PROFILE_Q96;
  for (size_t i = 0; i < ACC_GENPP_BML_CAPACITY; i++) {
    revac_bml_entry_clear(&ta->bml[i]);
    revac_q96_bml_entry_clear(&ta->q96_bml[i]);
    revac_q96_delta6_entry_clear(&ta->q96_delta6_bml[i]);
  }
  for (size_t i = 0; i < ACC_GENPP_MAX_DELETIONS; i++) {
    revac_sync_delta_cache_entry_clear(&ta->sync_delta_cache[i]);
  }
  shake256_init(&hash);
  shake256_absorb(&hash, domain, sizeof(domain) - 1);
  shake256_absorb(&hash, public_seed, sizeof(public_seed));
  shake256_absorb(&hash, accumulator_seed, sizeof(accumulator_seed));
  shake256_finalize(&hash);
  shake256_squeeze(ta->context_pp_hash, sizeof(ta->context_pp_hash), &hash);
  (void)certified_upper;
  (void)covariance_certificate;
  return 1;
}

int revac_ta_setup(revac_ta_t *ta) {
  return revac_ta_setup_q96(ta);
}

void revac_ta_public_state(revac_public_state_t *public_state,
                           const revac_ta_t *ta) {
  memset(public_state, 0, sizeof(*public_state));
  public_state->acc_profile = ta->acc_profile;
  public_state->ta_pk = &ta->ta_pk;
  public_state->window_id = ta->acc_profile == REVAC_ACC_PROFILE_Q96
                                ? ta->acc96_window_id
                                : 0;
  memcpy(public_state->parameter_digest, ta->context_pp_hash,
         sizeof(public_state->parameter_digest));
  if (ta->acc_profile == REVAC_ACC_PROFILE_Q96) {
    if (ta->acc96_state == NULL) {
      return;
    }
    public_state->acc96_ctx = &ta->acc96_ctx;
    public_state->acc96_pp = &ta->acc96_state->pp;
    public_state->acc96_acc = &ta->acc96_state->accumulator;
    memcpy(public_state->acc96_public_seed, ta->acc96_public_seed,
           sizeof(public_state->acc96_public_seed));
    public_state->acc96_seq = ta->acc96_state->sequence;
    public_state->acc_seq = ta->acc96_state->sequence;
    public_state->acc96_witness_norm_sq_bound =
        ta->acc96_witness_norm_sq_bound;
    memcpy(public_state->context_pp_hash, ta->context_pp_hash,
           sizeof(public_state->context_pp_hash));
    {
      const uint32_t count = ta->acc96_state->sequence <
                                     ACC_GENPP_BML_CAPACITY
                                 ? ta->acc96_state->sequence
                                 : ACC_GENPP_BML_CAPACITY;
      const uint32_t first = ta->acc96_state->sequence - count + 1;
      for (uint32_t sequence = first; sequence <= ta->acc96_state->sequence;
           sequence++) {
        const revac_q96_bml_entry_t *entry =
            &ta->q96_bml[(sequence - 1) % ACC_GENPP_BML_CAPACITY];
        if (!entry->used || entry->update.sequence != sequence) {
          public_state->q96_bml_len = 0;
          return;
        }
        public_state->q96_bml[public_state->q96_bml_len++] = entry;
      }
    }
    return;
  }
  public_state->acc_pp = &ta->acc_pp;
  public_state->acc = ta->acc;
  public_state->acc_seq = ta->acc_st.ctr;
  memcpy(public_state->context_pp_hash, ta->context_pp_hash,
         sizeof(public_state->context_pp_hash));
  for (size_t i = 0; i < ta->acc_st.bml_len; i++) {
    size_t pos = (ta->acc_st.bml_start + i) % ACC_GENPP_BML_CAPACITY;
    const acc_genpp_update_t *update = &ta->acc_st.bml[pos];
    const revac_bml_entry_t *entry =
        &ta->bml[(update->seq - 1) % ACC_GENPP_BML_CAPACITY];

    if (!entry->used || entry->update.seq != update->seq ||
        entry->update.op != update->op ||
        entry->update.handle != update->handle ||
        entry->update.sign != update->sign) {
      public_state->bml_len = 0;
      return;
    }
    public_state->bml[i] = entry;
  }
  public_state->bml_len = ta->acc_st.bml_len;
}

int revac_ta_register_uav(revac_ta_t *ta, revac_uav_t *uav) {
  uint32_t handle = revac_handle_from_identity(&uav->pk, uav->id);

  if (revac_find_uav_by_handle(ta, handle) >= 0) {
    return 0;
  }
  for (size_t i = 0; i < REVAC_MAX_UAVS; i++) {
    if (!ta->uavs[i].used) {
      ta->uavs[i].used = 1;
      ta->uavs[i].id = uav->id;
      ta->uavs[i].handle = handle;
      revac_copy_user_pk(&ta->uavs[i].pk, &uav->pk);
      uav->handle = handle;
      uav->registered = 1;
      return 1;
    }
  }
  return 0;
}

int revac_ta_register_gcs(revac_ta_t *ta, revac_gcs_t *gcs) {
  (void)ta;
  gcs->registered = 1;
  gcs->cert_valid = REVAC_GCS_CERT_VALID;
  return 1;
}

int revac_ta_issue_uav(revac_ta_t *ta,
                       revac_uav_t *uav,
                       const uint8_t msg[PARAM_M * PARAM_N / 8]) {
  uint8_t state[STATE_BYTES];
  poly_q_vec_d r[2], cmt, next_acc;
  acc_genpp_update_t upmsg;
  int db_pos;
  int ok = 0;

  if (ta != NULL && ta->acc_profile == REVAC_ACC_PROFILE_Q96) {
    return revac_ta_issue_uav_q96(ta, uav, msg);
  }

  if (!uav->registered) {
    return 0;
  }
  db_pos = revac_find_uav_by_handle(ta, uav->handle);
  if (db_pos < 0 || ta->uavs[db_pos].revoked) {
    return 0;
  }

  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  poly_q_vec_d_init(next_acc);

  randombytes(state, sizeof(state));
  osig_user_commit_handle(r, cmt, msg, &uav->pk, uav->handle);
  osig_signer_sign_commitment(&uav->sig, state, &ta->ta_sk, &ta->ta_pk, cmt);
  osig_user_sig_complete(&uav->sig, r);
  if (!osig_user_verify_handle(&uav->sig, &ta->ta_pk, &uav->pk, msg, uav->handle)) {
    goto cleanup;
  }
  if (!acc_genpp_add(next_acc, &uav->witness, &upmsg, &ta->acc_st,
                     &ta->acc_pp, &ta->acc_td, ta->acc, uav->handle)) {
    goto cleanup;
  }

  poly_q_vec_d_set(ta->acc, next_acc);
  uav->last_seq = ta->acc_st.ctr;
  uav->memver_cache_valid = 0;
  uav->memver_lift_cache_valid = 0;
  uav->issued = 1;
  ta->uavs[db_pos].issued = 1;
  ta->uavs[db_pos].baseline_seq = ta->acc_st.ctr;
  (void)revac_ta_sync_delta_cache_gc(ta);
  ok = 1;

cleanup:
  poly_q_vec_d_clear(next_acc);
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  return ok;
}

int revac_ta_issue_uav_q96(revac_ta_t *ta,
                            revac_uav_t *uav,
                            const uint8_t msg[PARAM_M * PARAM_N / 8]) {
  uint8_t state[STATE_BYTES];
  poly_q_vec_d r[2], cmt;
  int db_pos;
  int ok = 0;

  if (ta == NULL || uav == NULL || msg == NULL ||
      ta->acc_profile != REVAC_ACC_PROFILE_Q96 || ta->acc96_state == NULL ||
      ta->acc96_sampler == NULL || ta->acc96_witness_norm_sq_bound <= 0.0L ||
      !uav->registered) {
    return 0;
  }
  db_pos = revac_find_uav_by_handle(ta, uav->handle);
  if (db_pos < 0 || ta->uavs[db_pos].revoked) {
    return 0;
  }
  if (uav->acc96_witness == NULL) {
    uav->acc96_witness = calloc(1, sizeof(*uav->acc96_witness));
    if (uav->acc96_witness == NULL) {
      return 0;
    }
  }
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  randombytes(state, sizeof(state));
  osig_user_commit_handle(r, cmt, msg, &uav->pk, uav->handle);
  osig_signer_sign_commitment(&uav->sig, state, &ta->ta_sk, &ta->ta_pk, cmt);
  osig_user_sig_complete(&uav->sig, r);
  if (!osig_user_verify_handle(&uav->sig, &ta->ta_pk, &uav->pk, msg,
                               uav->handle) ||
      !acc96_add(ta->acc96_state, uav->acc96_witness, uav->handle,
                 ta->acc96_witness_norm_sq_bound, ta->acc96_sampler,
                 ta->acc96_sampler_opaque, &ta->acc96_ctx)) {
    goto cleanup;
  }
  uav->last_seq = ta->acc96_state->sequence;
  uav->issued = 1;
  ta->uavs[db_pos].issued = 1;
  ta->uavs[db_pos].baseline_seq = ta->acc96_state->sequence;
  ok = 1;
cleanup:
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  return ok;
}

static int revac_ta_revoke_uav_internal(revac_ta_t *ta,
                                        uint32_t handle,
                                        acc_genpp_update_t *upmsg) {
  poly_q_vec_d next_acc;
  acc_genpp_gmat_t tails[REVAC_TA_SYNC_DELTA_LAYERS];
  acc_genpp_gmat_t *capture[ACC_GENPP_HANDLE_BITS] = {NULL};
  acc_genpp_state_t state_before;
  revac_bml_entry_t bml_entry = {0};
  revac_sync_delta_cache_entry_t sync_entry = {0};
  int db_pos;
  int update_sign;
  int ok = 0;

  if (ta == NULL || ta->acc_profile != REVAC_ACC_PROFILE_Q19) {
    return 0;
  }
  db_pos = revac_find_uav_by_handle(ta, handle);
  if (db_pos < 0 || !ta->uavs[db_pos].issued || ta->uavs[db_pos].revoked) {
    return 0;
  }

  poly_q_vec_d_init(next_acc);
  for (size_t i = 0; i < REVAC_TA_SYNC_DELTA_LAYERS; i++) {
    acc_genpp_gmat_init(&tails[i]);
    capture[i + 1] = &tails[i];
  }

  /* TA already owns the authenticated registration database and accumulator
     dictionary.  Do not reconstruct a temporary witness for the revoked
     handle merely to feed the generic Delete API's MemVer gate. */
  /* The reviewed PKC'26 path uses the construction's deterministic +1
     update.  An authenticated Rademacher sign remains available only to an
     explicitly selected research profile because it needs a separate
     martingale/reduction argument. */
#if REVAC_PARAM_ACC_RANDOM_UPDATE_SIGNS
  {
    uint8_t sign_byte;
    randombytes(&sign_byte, sizeof(sign_byte));
    update_sign = (sign_byte & UINT8_C(1)) != 0 ? 1 : -1;
  }
#else
  update_sign = 1;
#endif
  state_before = ta->acc_st;
  if (!acc_genpp_delete_trusted_capture_signed(
          next_acc, upmsg, &ta->acc_st, &ta->acc_pp, ta->acc, handle,
          update_sign, capture)) {
    goto cleanup;
  }
  {
    const size_t bml_slot = (upmsg->seq - 1) % ACC_GENPP_BML_CAPACITY;
    const size_t sync_slot = upmsg->seq - 1;

    /* Every BML event is a mandatory signed DELTA4 broadcast; the same
       revocation also receives its TA-private DELTA6 acceleration entry. */
    if (!revac_bml_entry_and_sync_cache_build(
        &bml_entry, &sync_entry, upmsg->seq,
        upmsg->handle, upmsg->sign, ta->context_pp_hash, tails, ta->acc_pp.ua,
        &ta->ta_sk, &ta->ta_pk)) {
      ta->acc_st = state_before;
      goto cleanup;
    }
    revac_bml_entry_clear(&ta->bml[bml_slot]);
    revac_sync_delta_cache_entry_clear(&ta->sync_delta_cache[sync_slot]);
    ta->bml[bml_slot] = bml_entry;
    ta->sync_delta_cache[sync_slot] = sync_entry;
    memset(&bml_entry, 0, sizeof(bml_entry));
    memset(&sync_entry, 0, sizeof(sync_entry));
  }
  poly_q_vec_d_set(ta->acc, next_acc);
  ta->uavs[db_pos].revoked = 1;
  (void)revac_ta_sync_delta_cache_gc(ta);
  ok = 1;

cleanup:
  revac_bml_entry_clear(&bml_entry);
  revac_sync_delta_cache_entry_clear(&sync_entry);
  for (size_t i = 0; i < REVAC_TA_SYNC_DELTA_LAYERS; i++) {
    acc_genpp_gmat_clear(&tails[i]);
  }
  poly_q_vec_d_clear(next_acc);
  return ok;
}

int revac_ta_revoke_uav_q96(revac_ta_t *ta,
                             uint32_t handle,
                             acc96_relation_update_t *update) {
  acc96_poly_t accumulator_before;
  acc96_member_entry_t member_before;
  acc96_gmat_t *tail4 = NULL;
  acc96_delta6_t *delta6 = NULL;
  revac_q96_bml_entry_t entry = {0};
  revac_q96_delta6_entry_t delta6_entry = {0};
  revac_public_state_t public_state;
  uint8_t accumulator_hash[REVAC_CONTEXT_HASH_BYTES];
  uint32_t sequence_before;
  size_t slot;
  int db_pos;

  if (ta == NULL || update == NULL ||
      ta->acc_profile != REVAC_ACC_PROFILE_Q96 || ta->acc96_state == NULL) {
    return 0;
  }
  db_pos = revac_find_uav_by_handle(ta, handle);
  if (db_pos < 0 || !ta->uavs[db_pos].issued || ta->uavs[db_pos].revoked) {
    return 0;
  }
  sequence_before = ta->acc96_state->sequence;
  tail4 = malloc(sizeof(*tail4));
  delta6 = calloc(1, sizeof(*delta6));
  if (tail4 == NULL || delta6 == NULL) {
    free(tail4);
    free(delta6);
    return 0;
  }
  acc96_poly_set(&accumulator_before, &ta->acc96_state->accumulator);
  {
    size_t member_pos;
    for (member_pos = 0; member_pos < ACC96_MAX_MEMBERS; member_pos++) {
      if (ta->acc96_state->members[member_pos].status != ACC96_MEMBER_EMPTY &&
          ta->acc96_state->members[member_pos].handle == handle) {
        break;
      }
    }
    if (member_pos == ACC96_MAX_MEMBERS) {
      free(tail4);
      free(delta6);
      return 0;
    }
    member_before = ta->acc96_state->members[member_pos];
    /* Fuse the authoritative Delete with DELTA6 construction.  The legacy
       T4 checkpoint is then derived only from the high suffix, avoiding a
       second complete EvalF while preserving the old export API. */
    if (!acc96_delete_delta6(ta->acc96_state, update, delta6, handle,
                             &ta->acc96_ctx) ||
        !acc96_evalf_tail_checkpoint(tail4, &ta->acc96_state->pp.ba, handle,
                                     ACC96_TAIL4_BIT, &ta->acc96_ctx)) {
      acc96_poly_set(&ta->acc96_state->accumulator, &accumulator_before);
      ta->acc96_state->sequence = sequence_before;
      ta->acc96_state->members[member_pos] = member_before;
      memset(update, 0, sizeof(*update));
      free(tail4);
      free(delta6);
      return 0;
    }
    revac_ta_public_state(&public_state, ta);
    if (!revac_public_accumulator_digest(accumulator_hash, &public_state) ||
        !revac_q96_bml_entry_build(
            &entry, update, tail4, &ta->acc96_ctx, ta->context_pp_hash,
            accumulator_hash,
            &ta->ta_sk, &ta->ta_pk) ||
        !revac_q96_delta6_entry_build(
            &delta6_entry, update, delta6, &ta->acc96_ctx,
            ta->context_pp_hash, accumulator_hash,
            &ta->ta_sk, &ta->ta_pk)) {
      acc96_poly_set(&ta->acc96_state->accumulator, &accumulator_before);
      ta->acc96_state->sequence = sequence_before;
      ta->acc96_state->members[member_pos] = member_before;
      revac_q96_bml_entry_clear(&entry);
      revac_q96_delta6_entry_clear(&delta6_entry);
      memset(update, 0, sizeof(*update));
      free(tail4);
      free(delta6);
      return 0;
    }
  }
  slot = (update->sequence - 1) % ACC_GENPP_BML_CAPACITY;
  revac_q96_bml_entry_clear(&ta->q96_bml[slot]);
  ta->q96_bml[slot] = entry;
  revac_q96_delta6_entry_clear(&ta->q96_delta6_bml[slot]);
  ta->q96_delta6_bml[slot] = delta6_entry;
  memset(&entry, 0, sizeof(entry));
  memset(&delta6_entry, 0, sizeof(delta6_entry));
  ta->uavs[db_pos].revoked = 1;
  free(tail4);
  free(delta6);
  return 1;
}

int revac_ta_export_q96_update(const revac_ta_t *ta, uint32_t sequence,
                               uint8_t **wire, size_t *wire_len) {
  const revac_q96_bml_entry_t *entry;

  if (wire == NULL || wire_len == NULL) return 0;
  *wire = NULL;
  *wire_len = 0;
  if (ta == NULL || ta->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      ta->acc96_state == NULL || sequence == 0 ||
      sequence > ta->acc96_state->sequence) {
    return 0;
  }
  entry = &ta->q96_bml[(sequence - 1) % ACC_GENPP_BML_CAPACITY];
  if (!entry->used || entry->update.sequence != sequence) return 0;
  return revac_q96_bml_entry_packet_pack(wire, wire_len, entry);
}

int revac_ta_export_q96_update_manifest(const revac_ta_t *ta,
                                        uint32_t sequence,
                                        uint8_t **wire, size_t *wire_len) {
  const revac_q96_bml_entry_t *entry;
  if (wire == NULL || wire_len == NULL) return 0;
  *wire = NULL;
  *wire_len = 0;
  if (ta == NULL || ta->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      ta->acc96_state == NULL || sequence == 0 ||
      sequence > ta->acc96_state->sequence)
    return 0;
  entry = &ta->q96_bml[(sequence - 1) % ACC_GENPP_BML_CAPACITY];
  if (!entry->used || entry->update.sequence != sequence) return 0;
  return revac_q96_bml_entry_manifest_pack(wire, wire_len, entry);
}

int revac_ta_export_q96_delta6_manifest(const revac_ta_t *ta,
                                        uint32_t sequence,
                                        uint8_t **wire, size_t *wire_len) {
  const revac_q96_delta6_entry_t *entry;
  if (wire == NULL || wire_len == NULL) return 0;
  *wire = NULL; *wire_len = 0;
  if (ta == NULL || ta->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      ta->acc96_state == NULL || sequence == 0 ||
      sequence > ta->acc96_state->sequence)
    return 0;
  entry = &ta->q96_delta6_bml[(sequence - 1) % ACC_GENPP_BML_CAPACITY];
  if (!entry->used || entry->update.sequence != sequence) return 0;
  return revac_q96_delta6_manifest_pack(wire, wire_len, entry);
}

int revac_ta_export_q96_delta6_layer(const revac_ta_t *ta,
                                     uint32_t sequence, size_t layer,
                                     uint8_t **wire, size_t *wire_len) {
  const revac_q96_delta6_entry_t *entry;
  if (wire == NULL || wire_len == NULL) return 0;
  *wire = NULL; *wire_len = 0;
  if (ta == NULL || ta->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      ta->acc96_state == NULL || sequence == 0 ||
      sequence > ta->acc96_state->sequence || layer >= REVAC_Q96_DELTA_LAYERS)
    return 0;
  entry = &ta->q96_delta6_bml[(sequence - 1) % ACC_GENPP_BML_CAPACITY];
  if (!entry->used || entry->update.sequence != sequence) return 0;
  return revac_q96_delta6_layer_pack(wire, wire_len, entry, layer);
}

int revac_ta_export_q96_tail4_prefetch(const revac_ta_t *ta,
                                       uint32_t candidate_handle,
                                       uint8_t **wire, size_t *wire_len) {
  acc96_gmat_t *tail4;
  int db_pos;
  int ok;
  if (wire == NULL || wire_len == NULL) return 0;
  *wire = NULL;
  *wire_len = 0;
  if (ta == NULL || ta->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      ta->acc96_state == NULL)
    return 0;
  db_pos = revac_find_uav_by_handle(ta, candidate_handle);
  if (db_pos < 0 || !ta->uavs[db_pos].issued || ta->uavs[db_pos].revoked)
    return 0;
  tail4 = malloc(sizeof(*tail4));
  if (tail4 == NULL) return 0;
  ok = acc96_evalf_tail_checkpoint(
           tail4, &ta->acc96_state->pp.ba, candidate_handle,
           ACC96_TAIL4_BIT, &ta->acc96_ctx) &&
       revac_q96_tail4_prefetch_pack(
           wire, wire_len, candidate_handle, tail4,
           ta->context_pp_hash, &ta->ta_sk, &ta->ta_pk);
  free(tail4);
  return ok;
}

int revac_ta_revoke_uav(revac_ta_t *ta,
                        uint32_t handle,
                        acc_genpp_update_t *upmsg) {
  return revac_ta_revoke_uav_internal(ta, handle, upmsg);
}

int revac_ta_revoke_uav_context(revac_ta_t *ta,
                                uint32_t handle,
                                acc_genpp_update_t *upmsg) {
  return revac_ta_revoke_uav_internal(ta, handle, upmsg);
}

int revac_ta_export_context(const revac_ta_t *ta,
                            uint32_t seq,
                            revac_context_profile_t profile,
                            uint8_t **wire,
                            size_t *wire_len) {
  if (wire == NULL || wire_len == NULL) {
    return 0;
  }
  *wire = NULL;
  *wire_len = 0;
  if (ta == NULL || ta->acc_profile != REVAC_ACC_PROFILE_Q19 ||
      profile != REVAC_CONTEXT_PROFILE_DELTA4) {
    return 0;
  }
  {
    const revac_bml_entry_t *entry = revac_find_bml_entry_by_seq(ta, seq);

    if (entry != NULL) {
      return revac_bml_entry_packet_pack(wire, wire_len, entry);
    }
  }
  return 0;
}

int revac_ta_export_memwitsync(const revac_ta_t *ta,
                               uint32_t authenticated_handle,
                               const revac_sync_request_t *request,
                               uint8_t **wire,
                               size_t *wire_len) {
  acc_genpp_wit_t witness;
  acc96_relation_witness_t *witness96 = NULL;
  int db_pos;
  int ok = 0;

  if (wire == NULL || wire_len == NULL) {
    return 0;
  }
  *wire = NULL;
  *wire_len = 0;
  if (ta == NULL || request == NULL ||
      request->handle != authenticated_handle ||
      memcmp(request->pp_hash, ta->context_pp_hash,
             REVAC_CONTEXT_HASH_BYTES) != 0 ||
      (ta->acc_profile == REVAC_ACC_PROFILE_Q19 &&
       request->last_seq > ta->acc_st.ctr) ||
      (ta->acc_profile == REVAC_ACC_PROFILE_Q96 &&
       (ta->acc96_state == NULL || ta->acc96_sampler == NULL ||
        request->last_seq > ta->acc96_state->sequence))) {
    return 0;
  }
  db_pos = revac_find_uav_by_handle(ta, request->handle);
  if (db_pos < 0 || !ta->uavs[db_pos].issued || ta->uavs[db_pos].revoked) {
    return 0;
  }

  if (ta->acc_profile == REVAC_ACC_PROFILE_Q96) {
    witness96 = calloc(1, sizeof(*witness96));
    if (witness96 == NULL) return 0;
    if (ta->acc96_sampler(witness96, &ta->acc96_state->pp,
                          &ta->acc96_state->trapdoor,
                          &ta->acc96_state->accumulator, request->handle,
                          ta->acc96_sampler_opaque, &ta->acc96_ctx) &&
        acc96_relation_fresh_witness_check(witness96, &ta->acc96_ctx) &&
        acc96_relation_memver(&ta->acc96_state->pp,
                              &ta->acc96_state->accumulator,
                              request->handle, witness96,
                              ta->acc96_witness_norm_sq_bound,
                              &ta->acc96_ctx)) {
      ok = revac_q96_sync_response_pack(
          wire, wire_len, request, ta->acc96_state->sequence,
          &ta->acc96_state->accumulator, witness96,
          &ta->acc96_ctx,
          &ta->ta_sk, &ta->ta_pk);
    }
    free(witness96);
    return ok;
  }
  if (ta->acc_profile != REVAC_ACC_PROFILE_Q19) return 0;
  acc_genpp_wit_init(&witness);
  if (revac_ta_memwitsync_cached(&witness, ta, request->handle,
                                 NULL, NULL))
    ok = revac_sync_response_pack(wire, wire_len, request,
                                  ta->acc_st.ctr, ta->acc,
                                  &witness, &ta->ta_sk, &ta->ta_pk);
  acc_genpp_wit_clear(&witness);
  return ok;
}
