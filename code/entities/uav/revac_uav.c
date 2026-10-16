#include "revac_uav.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "randombytes.h"
#include "fips202.h"
#include "revac_context.h"
#include "revac_labrador_bridge.h"
#include "revac_proof_codec.h"
#include "revac_sync.h"

static double revac_uav_monotonic_seconds(void) {
  struct timespec now;

  (void)clock_gettime(CLOCK_MONOTONIC, &now);
  return (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
}

static void revac_uav_print_proof_timing(const char *label, double started) {
  printf("[证明计时] %-42s %10.6f 秒\n",
         label, revac_uav_monotonic_seconds() - started);
}

static int revac_uav_parallel_show_enabled(void) {
  const char *value = getenv("REVAC_SHOW_PARALLEL_COMPONENTS");
  long cpus;

  if (value != NULL && value[0] != '\0') {
    return value[0] != '0';
  }
  cpus = sysconf(_SC_NPROCESSORS_ONLN);
  return cpus >= 8;
}

static unsigned int revac_uav_parallel_signature_threads(void) {
  const char *value = getenv("REVAC_SHOW_PARALLEL_SIGNATURE_THREADS");

  if (value != NULL && value[0] != '\0') {
    char *end = NULL;
    long requested = strtol(value, &end, 10);
    if (end != value && *end == '\0' && requested >= 0 && requested <= 64) {
      return (unsigned int)requested;
    }
  }
  return 4;
}

static void revac_uav_invalidate_memver_caches(revac_uav_t *uav) {
  uav->memver_cache_valid = 0;
  uav->memver_lift_cache_valid = 0;
  uav->memver_q88_lift_valid = 0;
}

static void revac_uav_memver_absorb_u32(keccak_state *state, uint32_t value) {
  uint8_t encoded[4] = {
      (uint8_t)value, (uint8_t)(value >> 8),
      (uint8_t)(value >> 16), (uint8_t)(value >> 24)};
  shake256_absorb(state, encoded, sizeof(encoded));
}

static void revac_uav_memver_absorb_poly(keccak_state *state,
                                         const poly_q poly) {
  for (size_t i = 0; i < PARAM_N; i++) {
    revac_uav_memver_absorb_u32(
        state, (uint32_t)poly_q_get_coeff(poly, i));
  }
}

static void revac_uav_memver_absorb_vec(keccak_state *state,
                                        const poly_q_vec_d vec) {
  for (size_t i = 0; i < PARAM_D; i++) {
    revac_uav_memver_absorb_poly(state, vec->entries[i]);
  }
}

static void revac_uav_memver_absorb_acc96_poly(
    keccak_state *state, const acc96_poly_t *poly) {
  for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
    for (size_t lane = 0; lane < ACC96_LIMBS; lane++) {
      revac_uav_memver_absorb_u32(
          state, (uint32_t)poly->coeff[coefficient].limb[lane]);
    }
  }
}

static void revac_uav_memver_q96_cache_key(
    uint8_t out[REVAC_CONTEXT_HASH_BYTES],
    const revac_public_state_t *public_state, const revac_uav_t *uav) {
  static const uint8_t domain[] = "REVAC-UAV-Q96-MEMVER-CACHE-V1";
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, domain, sizeof(domain) - 1);
  char norm_bound[96];
  const int norm_bytes = snprintf(norm_bound, sizeof(norm_bound), "%.36La",
                                  public_state->acc96_witness_norm_sq_bound);
  shake256_absorb(&state, (const uint8_t *)norm_bound, (size_t)norm_bytes);
  shake256_absorb(&state, public_state->context_pp_hash,
                  REVAC_CONTEXT_HASH_BYTES);
  revac_uav_memver_absorb_u32(&state, public_state->acc96_seq);
  revac_uav_memver_absorb_u32(&state, uav->last_seq);
  revac_uav_memver_absorb_u32(&state, uav->handle);
  revac_uav_memver_absorb_acc96_poly(&state, public_state->acc96_acc);
  revac_uav_memver_absorb_acc96_poly(&state, &uav->acc96_witness->left);
  revac_uav_memver_absorb_acc96_poly(&state, &uav->acc96_witness->v12);
  for (size_t index = 0; index < ACC96_TRUNC_RIGHT; index++) {
    revac_uav_memver_absorb_acc96_poly(
        &state, &uav->acc96_witness->v2[index]);
  }
  revac_uav_memver_absorb_acc96_poly(&state, &uav->acc96_witness->tag);
  for (size_t bit = 0; bit < ACC96_HANDLE_BITS; bit++) {
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++) {
      revac_uav_memver_absorb_acc96_poly(
          &state, &uav->acc96_witness->right.right[bit][digit]);
    }
  }
  shake256_absorb(&state, uav->acc96_witness->refresh_r,
                  sizeof(uav->acc96_witness->refresh_r));
  for (size_t block = 0; block < ACC96_REFRESH_BLOCKS; block++) {
    revac_uav_memver_absorb_acc96_poly(
        &state, &uav->acc96_witness->refresh_z[block]);
  }
  shake256_finalize(&state);
  shake256_squeeze(out, REVAC_CONTEXT_HASH_BYTES, &state);
}

/* Canonical cache key for the exact clear MemVer input.  Hashing the roughly
   1 MiB witness is much cheaper than repeating all public matrix products,
   while still detecting direct in-memory mutation by callers/tests. */
static void revac_uav_memver_cache_key(
    uint8_t out[REVAC_CONTEXT_HASH_BYTES],
    const revac_public_state_t *public_state, const revac_uav_t *uav) {
  static const uint8_t domain[] = "REVAC-UAV-MEMVER-CACHE-V1";
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, domain, sizeof(domain) - 1);
  shake256_absorb(&state, public_state->context_pp_hash,
                  REVAC_CONTEXT_HASH_BYTES);
  revac_uav_memver_absorb_u32(&state, public_state->acc_seq);
  revac_uav_memver_absorb_u32(&state, uav->last_seq);
  revac_uav_memver_absorb_u32(&state, uav->handle);
  revac_uav_memver_absorb_vec(&state, public_state->acc);
  revac_uav_memver_absorb_vec(&state, uav->witness.left);
  revac_uav_memver_absorb_vec(&state, uav->witness.v12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    revac_uav_memver_absorb_vec(&state, uav->witness.v2[i]);
  }
  revac_uav_memver_absorb_poly(&state, uav->witness.tag);
  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      revac_uav_memver_absorb_vec(&state, uav->witness.right[bit][j]);
    }
  }
  shake256_absorb(&state, uav->witness.refresh_r,
                  sizeof(uav->witness.refresh_r));
  for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
    revac_uav_memver_absorb_vec(&state, uav->witness.refresh_z[j]);
  }
  shake256_finalize(&state);
  shake256_squeeze(out, REVAC_CONTEXT_HASH_BYTES, &state);
}

void revac_uav_init(revac_uav_t *uav) {
  memset(uav, 0, sizeof(*uav));
  user_keys_init(&uav->pk, &uav->sk);
  sep_sig_init(&uav->sig);
  acc_genpp_wit_init(&uav->witness);
  uav->acc96_witness = NULL;
}

void revac_uav_clear(revac_uav_t *uav) {
  for (size_t i = 0; i < ACC96_MAX_MEMBERS; i++) {
    free(uav->q96_tail_cache[i].tail4_packed);
    uav->q96_tail_cache[i].tail4_packed = NULL;
    uav->q96_tail_cache[i].used = 0;
  }
  free(uav->memver_q88_lift_cache);
  uav->memver_q88_lift_cache = NULL;
  free(uav->acc96_witness);
  uav->acc96_witness = NULL;
  acc_genpp_wit_clear(&uav->witness);
  sep_sig_clear(&uav->sig);
  user_keys_clear(&uav->pk, &uav->sk);
}

void revac_uav_keygen(revac_uav_t *uav, uint64_t id, const uint8_t seed[SEED_BYTES]) {
  uav->id = id;
  osig_user_keygen(&uav->pk, &uav->sk, seed);
}

void revac_uav_apply_update(revac_uav_t *uav,
                            const revac_public_state_t *public_state,
                            const acc_genpp_update_t *upmsg) {
  /* Keep the legacy void API, but do not let it bypass the mandatory BML
     DELTA4 authentication step.  A zero goodput selects LOCAL afterwards. */
  (void)revac_uav_apply_update_auto(uav, public_state, upmsg, 0.0, 0.0,
                                    NULL, NULL, NULL, NULL, NULL);
}

int revac_uav_memver_q96(const revac_uav_t *uav,
                         const revac_public_state_t *public_state) {
  if (uav == NULL || public_state == NULL || !uav->issued ||
      public_state->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      public_state->acc96_ctx == NULL || public_state->acc96_pp == NULL ||
      public_state->acc96_acc == NULL || uav->acc96_witness == NULL ||
      public_state->acc96_witness_norm_sq_bound <= 0.0L ||
      uav->last_seq != public_state->acc96_seq) {
    return 0;
  }
  return acc96_relation_memver(public_state->acc96_pp,
                               public_state->acc96_acc, uav->handle,
                               uav->acc96_witness,
                               public_state->acc96_witness_norm_sq_bound,
                               public_state->acc96_ctx);
}

static int revac_uav_q88_prepare_checked_lift(
    revac_uav_t *uav, const revac_public_state_t *public_state) {
  const char *mode = getenv("REVAC_Q88_CHECKED_LIFT_CACHE");
  if (!revac_labrador_online_enabled() || (mode && mode[0] == '0'))
    return revac_uav_memver_q96(uav, public_state);
  uav->memver_q88_lift_valid = 0;
  if (!uav->issued || !public_state->acc96_ctx || !public_state->acc96_pp ||
      !public_state->acc96_acc || !uav->acc96_witness ||
      uav->last_seq != public_state->acc96_seq ||
      !(public_state->acc96_witness_norm_sq_bound > 0.0L) ||
      !(acc96_relation_witness_norm2(uav->acc96_witness, public_state->acc96_ctx)
        <= public_state->acc96_witness_norm_sq_bound)) return 0;
  if (!uav->memver_q88_lift_cache)
    uav->memver_q88_lift_cache = malloc(sizeof(*uav->memver_q88_lift_cache));
  if (!uav->memver_q88_lift_cache || !acc96_memver_exact_lift(
      uav->memver_q88_lift_cache, &uav->memver_q88_lift_bits, public_state->acc96_pp,
      public_state->acc96_acc, uav->handle, uav->acc96_witness,
      public_state->acc96_ctx)) return 0;
  /* Divisibility in all three primes proves the same modular MemVer
     equation. The original aggregate norm bound was checked above. */
  revac_uav_memver_q96_cache_key(uav->memver_q88_lift_key, public_state, uav);
  uav->memver_q88_lift_valid = 1;
  return 1;
}

int revac_uav_q88_cached_lift_copy(acc96_exact_lift_t *out, unsigned *maximum_abs_bits,
    const revac_uav_t *uav, const revac_public_state_t *public_state) {
  uint8_t key[REVAC_CONTEXT_HASH_BYTES];
  if (!out || !uav || !public_state || !uav->memver_q88_lift_valid ||
      !uav->memver_q88_lift_cache || !uav->acc96_witness ||
      !public_state->acc96_acc || uav->last_seq != public_state->acc96_seq) return 0;
  revac_uav_memver_q96_cache_key(key, public_state, uav);
  if (memcmp(key, uav->memver_q88_lift_key, sizeof(key))) return 0;
  memcpy(out, uav->memver_q88_lift_cache, sizeof(*out));
  if (maximum_abs_bits) *maximum_abs_bits = uav->memver_q88_lift_bits;
  return 1;
}

int revac_uav_memver_cache_prewarm(
    revac_uav_t *uav, const revac_public_state_t *public_state) {
  uint8_t cache_key[REVAC_CONTEXT_HASH_BYTES];

  if (uav == NULL || public_state == NULL ||
      public_state->acc_profile != REVAC_ACC_PROFILE_Q96) {
    return 0;
  }
  revac_uav_memver_q96_cache_key(cache_key, public_state, uav);
  if (uav->memver_cache_valid &&
      uav->memver_cache_seq == public_state->acc96_seq &&
      memcmp(uav->memver_cache_key, cache_key, sizeof(cache_key)) == 0) {
    return 1;
  }
  if (!revac_uav_q88_prepare_checked_lift(uav, public_state)) {
    uav->memver_cache_valid = 0;
    return 0;
  }
  memcpy(uav->memver_cache_key, cache_key, sizeof(cache_key));
  uav->memver_cache_seq = public_state->acc96_seq;
  uav->memver_cache_valid = 1;
  return 1;
}

int revac_uav_apply_update_q96(revac_uav_t *uav,
                               const revac_public_state_t *public_state,
                               const acc96_relation_update_t *update) {
  acc96_relation_witness_t *candidate;
  int ok;

  if (uav == NULL || public_state == NULL || update == NULL ||
      public_state->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      public_state->acc96_ctx == NULL || public_state->acc96_pp == NULL ||
      public_state->acc96_acc == NULL || uav->acc96_witness == NULL ||
      update->sequence != uav->last_seq + 1 ||
      update->sequence > public_state->acc96_seq ||
      update->revoked_handle == uav->handle) {
    return 0;
  }
  candidate = malloc(sizeof(*candidate));
  if (candidate == NULL) {
    return 0;
  }
  memcpy(candidate, uav->acc96_witness, sizeof(*candidate));
  ok = acc96_relation_memwitup(candidate, public_state->acc96_pp,
                               uav->handle, update,
                               public_state->acc96_ctx) &&
       acc96_relation_memver(public_state->acc96_pp,
                             public_state->acc96_acc, uav->handle,
                             candidate,
                             public_state->acc96_witness_norm_sq_bound,
                             public_state->acc96_ctx);
  if (ok) {
    memcpy(uav->acc96_witness, candidate, sizeof(*candidate));
    uav->last_seq = update->sequence;
    revac_uav_invalidate_memver_caches(uav);
  }
  free(candidate);
  return ok;
}

int revac_uav_apply_update_q96_delta4(
    revac_uav_t *uav, const revac_public_state_t *public_state,
    const acc96_relation_update_t *update, const acc96_delta4_t *delta,
    size_t *applied_layers) {
  acc96_relation_witness_t *candidate;
  size_t used = 0;
  int ok;

  if (applied_layers != NULL) *applied_layers = 0;
  if (uav == NULL || public_state == NULL || update == NULL || delta == NULL ||
      public_state->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      public_state->acc96_ctx == NULL || public_state->acc96_pp == NULL ||
      public_state->acc96_acc == NULL || uav->acc96_witness == NULL ||
      update->sequence != uav->last_seq + 1 ||
      update->sequence > public_state->acc96_seq ||
      update->revoked_handle == uav->handle) {
    return 0;
  }
  candidate = malloc(sizeof(*candidate));
  if (candidate == NULL) return 0;
  memcpy(candidate, uav->acc96_witness, sizeof(*candidate));
  ok = acc96_relation_memwitup_delta4(candidate, uav->handle, update, delta,
                                      &used) &&
       acc96_relation_memver(public_state->acc96_pp,
                             public_state->acc96_acc, uav->handle,
                             candidate,
                             public_state->acc96_witness_norm_sq_bound,
                             public_state->acc96_ctx);
  if (ok) {
    memcpy(uav->acc96_witness, candidate, sizeof(*candidate));
    uav->last_seq = update->sequence;
    revac_uav_invalidate_memver_caches(uav);
    if (applied_layers != NULL) *applied_layers = used;
  }
  free(candidate);
  return ok;
}

/* The packet verifier has already authenticated the exact T4 checkpoint, the
 * public parameters, sequence and post-delete accumulator.  Given a witness
 * valid at sequence-1, the EvalF/EvalFx identity proves the updated algebraic
 * relation, so repeating the roughly 185-product clear MemVer here is not a
 * security check.  Keep the public norm gate; the next anonymous show still
 * performs its normal cached MemVer/proof checks. */
static int revac_uav_apply_authenticated_update_q96_tail4(
    revac_uav_t *uav, const revac_public_state_t *public_state,
    const acc96_relation_update_t *update, const acc96_gmat_t *tail4) {
  acc96_relation_witness_t *candidate;
  int ok;

  if (uav == NULL || public_state == NULL || update == NULL || tail4 == NULL ||
      public_state->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      public_state->acc96_ctx == NULL ||
      public_state->acc96_witness_norm_sq_bound <= 0.0L ||
      uav->acc96_witness == NULL || update->sequence != uav->last_seq + 1 ||
      update->revoked_handle == uav->handle)
    return 0;
  candidate = malloc(sizeof(*candidate));
  if (candidate == NULL) return 0;
  memcpy(candidate, uav->acc96_witness, sizeof(*candidate));
  ok = acc96_relation_memwitup_tail4(candidate, public_state->acc96_pp,
                                     uav->handle, update, tail4,
                                     public_state->acc96_ctx) &&
       acc96_relation_witness_norm2(candidate, public_state->acc96_ctx) <=
           public_state->acc96_witness_norm_sq_bound;
  if (ok) {
    memcpy(uav->acc96_witness, candidate, sizeof(*candidate));
    uav->last_seq = update->sequence;
    revac_uav_invalidate_memver_caches(uav);
  }
  free(candidate);
  return ok;
}

int revac_uav_apply_update_q96_wire(
    revac_uav_t *uav, const revac_public_state_t *public_state,
    const uint8_t *wire, size_t wire_len) {
  acc96_relation_update_t update;
  acc96_gmat_t *tail4;
  int ok;

  tail4 = malloc(sizeof(*tail4));
  if (tail4 == NULL) return 0;
  if (!revac_q96_update_packet_unpack_tail4(
          &update, tail4, wire, wire_len, public_state)) {
    free(tail4);
    return 0;
  }
  if (((update.revoked_handle ^ uav->handle) & UINT32_C(0xf)) != 0) {
    ok = revac_uav_apply_authenticated_update_q96_tail4(
        uav, public_state, &update, tail4);
  } else {
    ok = revac_uav_apply_update_q96(uav, public_state, &update);
  }
  free(tail4);
  return ok;
}

int revac_uav_cache_q96_tail4(revac_uav_t *uav,
                              const revac_public_state_t *public_state,
                              const uint8_t *full_wire, size_t full_wire_len) {
  acc96_relation_update_t update;
  acc96_gmat_t *tail4;
  size_t slot = ACC96_MAX_MEMBERS;
  if (uav == NULL || public_state == NULL || full_wire == NULL) return 0;
  tail4 = malloc(sizeof(*tail4));
  if (tail4 == NULL) return 0;
  if (!revac_q96_update_packet_unpack_tail4(
          &update, tail4, full_wire, full_wire_len, public_state)) {
    free(tail4);
    return 0;
  }
  for (size_t i = 0; i < ACC96_MAX_MEMBERS; i++) {
    if (uav->q96_tail_cache[i].used &&
        uav->q96_tail_cache[i].revoked_handle == update.revoked_handle) {
      slot = i;
      break;
    }
    if (slot == ACC96_MAX_MEMBERS && !uav->q96_tail_cache[i].used) slot = i;
  }
  if (slot == ACC96_MAX_MEMBERS) {
    free(tail4);
    return 0;
  }
  {
    uint8_t *packed = malloc(REVAC_Q96_TAIL4_BYTES);
    if (packed == NULL) {
      free(tail4);
      return 0;
    }
    memcpy(packed,
           full_wire + REVAC_Q96_UPDATE_WIRE_HEADER_BYTES +
               REVAC_CONTEXT_SIGNATURE_BYTES,
           REVAC_Q96_TAIL4_BYTES);
    free(uav->q96_tail_cache[slot].tail4_packed);
    uav->q96_tail_cache[slot].tail4_packed = packed;
  }
  uav->q96_tail_cache[slot].revoked_handle = update.revoked_handle;
  memcpy(uav->q96_tail_cache[slot].pp_hash, full_wire + 24,
         REVAC_CONTEXT_HASH_BYTES);
  memcpy(uav->q96_tail_cache[slot].tail_hash, full_wire + 88,
         REVAC_CONTEXT_HASH_BYTES);
  uav->q96_tail_cache[slot].used = 1;
  free(tail4);
  return 1;
}

int revac_uav_cache_q96_tail4_prefetch(
    revac_uav_t *uav, const revac_public_state_t *public_state,
    const uint8_t *prefetch_wire, size_t prefetch_wire_len) {
  uint32_t candidate_handle;
  uint8_t tail_hash[REVAC_CONTEXT_HASH_BYTES];
  acc96_gmat_t *tail4;
  size_t slot = ACC96_MAX_MEMBERS;
  if (uav == NULL || public_state == NULL || prefetch_wire == NULL)
    return 0;
  tail4 = malloc(sizeof(*tail4));
  if (tail4 == NULL) return 0;
  if (!revac_q96_tail4_prefetch_unpack(
          &candidate_handle, tail4, tail_hash, prefetch_wire,
          prefetch_wire_len, public_state)) {
    free(tail4);
    return 0;
  }
  for (size_t i = 0; i < ACC96_MAX_MEMBERS; i++) {
    if (uav->q96_tail_cache[i].used &&
        uav->q96_tail_cache[i].revoked_handle == candidate_handle) {
      slot = i;
      break;
    }
    if (slot == ACC96_MAX_MEMBERS && !uav->q96_tail_cache[i].used) slot = i;
  }
  if (slot == ACC96_MAX_MEMBERS) {
    free(tail4);
    return 0;
  }
  {
    uint8_t *packed = malloc(REVAC_Q96_TAIL4_BYTES);
    if (packed == NULL) {
      free(tail4);
      return 0;
    }
    memcpy(packed,
           prefetch_wire + REVAC_Q96_PREFETCH_WIRE_HEADER_BYTES +
               REVAC_CONTEXT_SIGNATURE_BYTES,
           REVAC_Q96_TAIL4_BYTES);
    free(uav->q96_tail_cache[slot].tail4_packed);
    uav->q96_tail_cache[slot].tail4_packed = packed;
  }
  uav->q96_tail_cache[slot].revoked_handle = candidate_handle;
  memcpy(uav->q96_tail_cache[slot].pp_hash,
         public_state->context_pp_hash, REVAC_CONTEXT_HASH_BYTES);
  memcpy(uav->q96_tail_cache[slot].tail_hash, tail_hash,
         REVAC_CONTEXT_HASH_BYTES);
  uav->q96_tail_cache[slot].used = 1;
  free(tail4);
  return 1;
}

static int revac_uav_unpack_cached_q96_tail4(
    acc96_gmat_t *tail4, const uint8_t *packed) {
  size_t offset = 0;
  if (tail4 == NULL || packed == NULL) return 0;
  for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++) {
    if (!acc96_poly_unpack(&tail4->block[digit], packed + offset,
                           ACC96_POLY_PACKED_BYTES))
      return 0;
    offset += ACC96_POLY_PACKED_BYTES;
  }
  return offset == REVAC_Q96_TAIL4_BYTES;
}

int revac_uav_apply_update_q96_manifest(
    revac_uav_t *uav, const revac_public_state_t *public_state,
    const uint8_t *manifest_wire, size_t manifest_wire_len) {
  acc96_relation_update_t update;
  uint8_t tail_hash[REVAC_CONTEXT_HASH_BYTES];
  if (uav == NULL || public_state == NULL || manifest_wire == NULL ||
      !revac_q96_update_manifest_unpack(&update, tail_hash, manifest_wire,
                                        manifest_wire_len, public_state))
    return 0;
  if (((update.revoked_handle ^ uav->handle) & UINT32_C(0xf)) == 0)
    return revac_uav_apply_update_q96(uav, public_state, &update);
  for (size_t i = 0; i < ACC96_MAX_MEMBERS; i++) {
    revac_q96_tail_cache_entry_t *entry = &uav->q96_tail_cache[i];
    if (entry->used && entry->revoked_handle == update.revoked_handle &&
        entry->tail4_packed != NULL &&
        memcmp(entry->pp_hash, public_state->context_pp_hash,
               REVAC_CONTEXT_HASH_BYTES) == 0 &&
        memcmp(entry->tail_hash, tail_hash, REVAC_CONTEXT_HASH_BYTES) == 0) {
      acc96_gmat_t *tail4 = malloc(sizeof(*tail4));
      int ok;
      if (tail4 == NULL) return 0;
      ok = revac_uav_unpack_cached_q96_tail4(tail4,
                                             entry->tail4_packed) &&
           revac_uav_apply_authenticated_update_q96_tail4(
               uav, public_state, &update, tail4);
      free(tail4);
      return ok;
    }
  }
  /* The manifest already authenticates the revoked handle and resulting
     accumulator.  A cache is only an acceleration hint, so a miss can always
     use the exact LOCAL recurrence without another network transfer. */
  return revac_uav_apply_update_q96(uav, public_state, &update);
}

int revac_uav_apply_update_q96_delta6(
    revac_uav_t *uav, const revac_public_state_t *public_state,
    const uint8_t *manifest_wire, size_t manifest_wire_len,
    const uint8_t *layer_packets[REVAC_Q96_DELTA_LAYERS],
    const size_t layer_packet_lens[REVAC_Q96_DELTA_LAYERS],
    size_t supplied_layers, size_t *applied_layers, int *used_local) {
  acc96_relation_update_t update;
  uint8_t layer_hash[REVAC_Q96_DELTA_LAYERS][REVAC_CONTEXT_HASH_BYTES];
  uint32_t layer_len[REVAC_Q96_DELTA_LAYERS];
  uint32_t difference;
  size_t first_difference = 0;
  if (applied_layers != NULL) *applied_layers = 0;
  if (used_local != NULL) *used_local = 0;
  if (uav == NULL || public_state == NULL || manifest_wire == NULL ||
      supplied_layers > REVAC_Q96_DELTA_LAYERS ||
      !revac_q96_delta6_manifest_unpack(&update, layer_hash, layer_len,
                                        manifest_wire, manifest_wire_len,
                                        public_state) ||
      update.sequence != uav->last_seq + 1 ||
      update.revoked_handle == uav->handle ||
      uav->acc96_witness == NULL || public_state->acc96_ctx == NULL ||
      public_state->acc96_pp == NULL)
    return 0;
  difference = update.revoked_handle ^ uav->handle;
  while ((difference & 1U) == 0U) {
    difference >>= 1;
    first_difference++;
  }
  if (first_difference >= REVAC_Q96_DELTA_LAYERS) {
    int ok = revac_uav_apply_update_q96(uav, public_state, &update);
    if (ok && used_local != NULL) *used_local = 1;
    return ok;
  }
  if (supplied_layers < first_difference + 1 || layer_packets == NULL ||
      layer_packet_lens == NULL)
    return 0;
  {
    acc96_delta6_t *delta = calloc(1, sizeof(*delta));
    acc96_relation_witness_t *candidate = malloc(sizeof(*candidate));
    int ok = 0;
    if (delta == NULL || candidate == NULL) {
      free(delta); free(candidate); return 0;
    }
    for (size_t i = 0; i <= first_difference; i++) {
      size_t packet_layer, payload_len;
      uint32_t packet_sequence;
      const uint8_t *payload;
      uint8_t digest[REVAC_CONTEXT_HASH_BYTES];
      if (layer_packets[i] == NULL ||
          !revac_q96_delta6_layer_unpack(&packet_sequence, &packet_layer,
                                         &payload, &payload_len,
                                         layer_packets[i], layer_packet_lens[i]) ||
          packet_sequence != update.sequence || packet_layer != i ||
          payload_len != layer_len[i]) {
        free(delta); free(candidate); return 0;
      }
      sha3_256(digest, payload, payload_len);
      if (memcmp(digest, layer_hash[i], REVAC_CONTEXT_HASH_BYTES) != 0 ||
          !acc96_delta6_unpack_layer(delta, NULL, payload, payload_len,
                                     public_state->acc96_ctx)) {
        free(delta); free(candidate); return 0;
      }
    }
    memcpy(candidate, uav->acc96_witness, sizeof(*candidate));
    ok = acc96_relation_memwitup_delta6(
             candidate, uav->handle, &update, delta, first_difference + 1,
             applied_layers) &&
         acc96_relation_witness_norm2(candidate, public_state->acc96_ctx) <=
             public_state->acc96_witness_norm_sq_bound;
    if (ok) {
      memcpy(uav->acc96_witness, candidate, sizeof(*candidate));
      uav->last_seq = update.sequence;
      revac_uav_invalidate_memver_caches(uav);
    }
    free(delta); free(candidate);
    return ok;
  }
}

static int revac_uav_update_sequence_ok(const revac_uav_t *uav,
                                        const acc_genpp_update_t *upmsg,
                                        int *already_applied) {
  *already_applied = upmsg->seq <= uav->last_seq;
  if (*already_applied) {
    return 1;
  }
  return upmsg->seq == uav->last_seq + 1;
}

static const revac_bml_entry_t *revac_public_find_bml_entry(
    const revac_public_state_t *public_state,
    const acc_genpp_update_t *upmsg) {
  if (public_state == NULL || upmsg == NULL) {
    return NULL;
  }
  for (size_t i = 0; i < public_state->bml_len; i++) {
    const revac_bml_entry_t *entry = public_state->bml[i];

    if (entry != NULL && entry->used &&
        entry->update.seq == upmsg->seq &&
        entry->update.op == upmsg->op &&
        entry->update.handle == upmsg->handle &&
        entry->update.sign == upmsg->sign) {
      return entry;
    }
  }
  return NULL;
}

/* A LOCAL computation must authenticate exactly the same BML broadcast as a
   DELTA4 computation.  This helper models receipt of that record in the
   in-process demo and verifies the wire-equivalent packet. */
static int revac_uav_verify_bml_delta4(
    const revac_public_state_t *public_state,
    const acc_genpp_update_t *upmsg) {
  const revac_bml_entry_t *entry;
  uint8_t *wire = NULL;
  size_t wire_len = 0;
  int ok;

  if (public_state == NULL || public_state->ta_pk == NULL || upmsg == NULL) {
    return 0;
  }
  entry = revac_public_find_bml_entry(public_state, upmsg);
  if (entry == NULL || !revac_bml_entry_packet_pack(&wire, &wire_len, entry)) {
    return 0;
  }
  ok = revac_context_packet_verify(wire, wire_len, upmsg,
                                   public_state->context_pp_hash,
                                   public_state->ta_pk);
  free(wire);
  return ok;
}

static int revac_uav_apply_update_local_unchecked(
    revac_uav_t *uav, const revac_public_state_t *public_state,
    const acc_genpp_update_t *upmsg, size_t *local_composites) {
  if (!acc_genpp_memwitup_tails(&uav->witness, public_state->acc_pp,
                                uav->handle, upmsg, NULL, NULL, 0,
                                local_composites)) {
    return 0;
  }
  uav->last_seq = upmsg->seq;
  revac_uav_invalidate_memver_caches(uav);
  return 1;
}

int revac_uav_apply_update_local(revac_uav_t *uav,
                                 const revac_public_state_t *public_state,
                                 const acc_genpp_update_t *upmsg,
                                 size_t *local_composites) {
  int already_applied;

  if (uav == NULL || public_state == NULL ||
      public_state->acc_profile != REVAC_ACC_PROFILE_Q19 ||
      public_state->acc_pp == NULL || upmsg == NULL ||
      !revac_uav_update_sequence_ok(uav, upmsg, &already_applied)) {
    return 0;
  }
  if (local_composites != NULL) {
    *local_composites = 0;
  }
  if (already_applied) {
    return 1;
  }
  if (!revac_uav_verify_bml_delta4(public_state, upmsg)) {
    return 0;
  }
  return revac_uav_apply_update_local_unchecked(uav, public_state, upmsg,
                                                local_composites);
}

int revac_uav_apply_update_context(revac_uav_t *uav,
                                   const revac_public_state_t *public_state,
                                   const acc_genpp_update_t *upmsg,
                                   const uint8_t *wire,
                                   size_t wire_len,
                                   int allow_local_fallback,
                                   revac_context_profile_t *used_profile,
                                   size_t *local_composites) {
  revac_delta_context_t ctx;
  int already_applied;
  int ok;

  if (used_profile != NULL) {
    *used_profile = REVAC_CONTEXT_PROFILE_LOCAL;
  }
  if (local_composites != NULL) {
    *local_composites = 0;
  }
  if (uav == NULL || public_state == NULL ||
      public_state->ta_pk == NULL || public_state->acc_pp == NULL ||
      upmsg == NULL ||
      !revac_uav_update_sequence_ok(uav, upmsg, &already_applied)) {
    return 0;
  }
  if (already_applied) {
    return 1;
  }

  revac_delta_context_init(&ctx);
  ok = wire != NULL &&
       revac_context_packet_unpack(&ctx, wire, wire_len, upmsg,
                                   public_state->context_pp_hash,
                                   public_state->ta_pk);
  if (!ok) {
    revac_delta_context_clear(&ctx);
    if (!allow_local_fallback) {
      return 0;
    }
    return revac_uav_apply_update_local(uav, public_state, upmsg,
                                        local_composites);
  }
  ok = acc_genpp_memwitup_delta4(&uav->witness, uav->handle, upmsg,
                                 ctx.delta, NULL);
  if (ok) {
    uav->last_seq = upmsg->seq;
    revac_uav_invalidate_memver_caches(uav);
    if (used_profile != NULL) {
      *used_profile = ctx.profile;
    }
  } else if (allow_local_fallback) {
    revac_delta_context_clear(&ctx);
    return revac_uav_apply_update_local(uav, public_state, upmsg,
                                        local_composites);
  }
  revac_delta_context_clear(&ctx);
  return ok;
}

int revac_uav_apply_update_auto(revac_uav_t *uav,
                                const revac_public_state_t *public_state,
                                const acc_genpp_update_t *upmsg,
                                double goodput_mbps,
                                double recurrence_ms,
                                revac_context_fetch_fn fetch_context,
                                void *fetch_opaque,
                                revac_context_profile_t *used_profile,
                                size_t *local_composites,
                                size_t *download_bytes) {
  revac_context_profile_t profile;
  revac_context_profile_t actual_profile = REVAC_CONTEXT_PROFILE_LOCAL;
  const revac_bml_entry_t *bml_entry;
  uint8_t *wire = NULL;
  size_t wire_len = 0;
  int already_applied;
  int ok;

  (void)fetch_context;
  (void)fetch_opaque;
  if (download_bytes != NULL) {
    *download_bytes = 0;
  }
  if (uav == NULL || public_state == NULL || upmsg == NULL ||
      public_state->acc_profile != REVAC_ACC_PROFILE_Q19) {
    return 0;
  }
  if (!revac_uav_update_sequence_ok(uav, upmsg, &already_applied)) {
    return 0;
  }
  if (already_applied) {
    if (used_profile != NULL) {
      *used_profile = REVAC_CONTEXT_PROFILE_LOCAL;
    }
    if (local_composites != NULL) {
      *local_composites = 0;
    }
    return 1;
  }
  bml_entry = revac_public_find_bml_entry(public_state, upmsg);
  if (bml_entry == NULL ||
      !revac_bml_entry_packet_pack(&wire, &wire_len, bml_entry)) {
    return 0;
  }
  if (download_bytes != NULL) {
    *download_bytes = wire_len;
  }
  if (!revac_context_packet_verify(wire, wire_len, upmsg,
                                   public_state->context_pp_hash,
                                   public_state->ta_pk)) {
    free(wire);
    return 0;
  }
  profile = revac_context_profile_select(goodput_mbps, recurrence_ms);
  if (profile == REVAC_CONTEXT_PROFILE_DELTA4) {
    ok = revac_uav_apply_update_context(uav, public_state, upmsg,
                                        wire, wire_len, 1,
                                        &actual_profile, local_composites);
  } else {
    /* The BML packet above has already been received and authenticated. */
    ok = revac_uav_apply_update_local_unchecked(
        uav, public_state, upmsg, local_composites);
  }
  if (used_profile != NULL) {
    *used_profile = actual_profile;
  }
  if (download_bytes != NULL) {
    *download_bytes = wire_len;
  }
  free(wire);
  return ok;
}

static int revac_public_get_missing(
    const revac_public_state_t *public_state,
    uint32_t last_seq,
    acc_genpp_update_t *missing,
    size_t missing_cap,
    size_t *missing_len) {
  size_t count = 0;

  if (public_state == NULL || missing_len == NULL ||
      last_seq > public_state->acc_seq) {
    return 0;
  }
  *missing_len = 0;
  if (last_seq == public_state->acc_seq) {
    return 1;
  }
  if (public_state->bml_len == 0 || public_state->bml[0] == NULL ||
      last_seq + 1 < public_state->bml[0]->update.seq) {
    return 0;
  }
  for (size_t i = 0; i < public_state->bml_len; i++) {
    if (public_state->bml[i] == NULL) {
      return 0;
    }
    if (public_state->bml[i]->update.seq <= last_seq) {
      continue;
    }
    if (count >= missing_cap) {
      return 0;
    }
    missing[count++] = public_state->bml[i]->update;
  }
  *missing_len = count;
  return count != 0 && missing[count - 1].seq == public_state->acc_seq;
}

static int revac_uav_remote_memwitsync(
    revac_uav_t *uav,
    const revac_public_state_t *public_state,
    revac_sync_fetch_fn fetch_sync,
    void *sync_opaque,
    size_t *download_bytes) {
  revac_sync_request_t request;
  acc_genpp_wit_t witness;
  acc96_relation_witness_t *witness96 = NULL;
  uint8_t *wire = NULL;
  size_t wire_len = 0;
  uint32_t response_seq = 0;
  int ok = 0;

  if (download_bytes != NULL) {
    *download_bytes = 0;
  }
  if (fetch_sync == NULL ||
      !revac_sync_request_create(&request, uav->handle, uav->last_seq,
                                 public_state) ||
      !fetch_sync(sync_opaque, &request, &wire, &wire_len)) {
    free(wire);
    return 0;
  }
  if (download_bytes != NULL) {
    *download_bytes = wire_len;
  }
  if (public_state->acc_profile == REVAC_ACC_PROFILE_Q96) {
    witness96 = calloc(1, sizeof(*witness96));
    if (witness96 != NULL &&
        revac_q96_sync_response_unpack(witness96, &response_seq, wire,
                                       wire_len, &request, public_state)) {
      if (uav->acc96_witness == NULL)
        uav->acc96_witness = calloc(1, sizeof(*uav->acc96_witness));
      if (uav->acc96_witness != NULL) {
        memcpy(uav->acc96_witness, witness96, sizeof(*witness96));
        uav->last_seq = response_seq;
        revac_uav_invalidate_memver_caches(uav);
        ok = 1;
      }
    }
    free(witness96);
    free(wire);
    return ok;
  }
  acc_genpp_wit_init(&witness);
  if (revac_sync_response_unpack(&witness, &response_seq, wire, wire_len,
                                 &request, public_state)) {
    revac_sync_witness_set(&uav->witness, &witness);
    uav->last_seq = response_seq;
    revac_uav_invalidate_memver_caches(uav);
    ok = 1;
  }
  acc_genpp_wit_clear(&witness);
  free(wire);
  return ok;
}

int revac_uav_recover_updates(revac_uav_t *uav,
                              const revac_public_state_t *public_state,
                              revac_sync_fetch_fn fetch_sync,
                              void *sync_opaque) {
  acc_genpp_update_t missing[ACC_GENPP_BML_CAPACITY];
  size_t missing_len = 0;

  if (uav == NULL || public_state == NULL) {
    return 0;
  }
  if (public_state->acc_profile == REVAC_ACC_PROFILE_Q96)
    return revac_uav_remote_memwitsync(uav, public_state, fetch_sync,
                                       sync_opaque, NULL);
  if (revac_public_get_missing(public_state, uav->last_seq, missing,
                               ACC_GENPP_BML_CAPACITY, &missing_len)) {
    for (size_t i = 0; i < missing_len; i++) {
      /* The non-auto recovery API deliberately uses LOCAL after receiving
         and authenticating the BML-carried DELTA4 for every update. */
      if (!revac_uav_apply_update_auto(uav, public_state, &missing[i], 0.0,
                                       0.0, NULL, NULL, NULL, NULL, NULL)) {
        return 0;
      }
    }
    return 1;
  }
  return revac_uav_remote_memwitsync(uav, public_state, fetch_sync,
                                     sync_opaque, NULL);
}

int revac_uav_recover_updates_auto(revac_uav_t *uav,
                                   const revac_public_state_t *public_state,
                                   double goodput_mbps,
                                   double recurrence_ms,
                                   revac_context_fetch_fn fetch_context,
                                   void *fetch_opaque,
                                   revac_sync_fetch_fn fetch_sync,
                                   void *sync_opaque,
                                   size_t *download_bytes,
                                   size_t *local_composites) {
  acc_genpp_update_t missing[ACC_GENPP_BML_CAPACITY];
  size_t missing_len = 0;
  size_t total_bytes = 0;
  size_t total_composites = 0;

  if (download_bytes != NULL) {
    *download_bytes = 0;
  }
  if (local_composites != NULL) {
    *local_composites = 0;
  }
  if (uav == NULL || public_state == NULL) {
    return 0;
  }
  if (public_state->acc_profile == REVAC_ACC_PROFILE_Q96)
    return revac_uav_remote_memwitsync(uav, public_state, fetch_sync,
                                       sync_opaque, download_bytes);
  if (revac_public_get_missing(public_state, uav->last_seq, missing,
                               ACC_GENPP_BML_CAPACITY, &missing_len)) {
    for (size_t i = 0; i < missing_len; i++) {
      size_t bytes = 0;
      size_t composites = 0;
      if (!revac_uav_apply_update_auto(uav, public_state, &missing[i],
                                       goodput_mbps,
                                       recurrence_ms, fetch_context,
                                       fetch_opaque, NULL, &composites,
                                       &bytes)) {
        return 0;
      }
      total_bytes += bytes;
      total_composites += composites;
    }
    if (download_bytes != NULL) {
      *download_bytes = total_bytes;
    }
    if (local_composites != NULL) {
      *local_composites = total_composites;
    }
    return 1;
  }

  /* BML exhaustion crosses the device boundary: TA computes and signs a
     fresh witness; the UAV only verifies and installs the response. */
  return revac_uav_remote_memwitsync(uav, public_state, fetch_sync,
                                     sync_opaque, download_bytes);
}

int revac_uav_show_prove(revac_show_proof_t *proof,
                         const revac_public_state_t *public_state,
                         revac_uav_t *uav,
                         const uint8_t msg[PARAM_M * PARAM_N / 8],
                         const uint8_t nonce[REVAC_SHOW_NONCE_BYTES]) {
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
  poly_qshow_vec_m1 s1;
  poly_q_mat_d_d A_link, D_link;
  poly_q_vec_d link_randomness[2];
  uint8_t crs_seed[CRS_SEED_BYTES];
  double step_started;
  double total_started = revac_uav_monotonic_seconds();
  uint8_t memver_key[REVAC_CONTEXT_HASH_BYTES];
  revac_labrador_prove_async_t *labrador_job = NULL;
  int memver_cache_hit;
  int ok = 0;

  if (!uav->issued || public_state == NULL ||
      (public_state->acc_profile != REVAC_ACC_PROFILE_Q19 &&
       public_state->acc_profile != REVAC_ACC_PROFILE_Q96)) {
    return 0;
  }
  if (public_state->ta_pk == NULL ||
      (public_state->acc_profile == REVAC_ACC_PROFILE_Q19 &&
       (public_state->acc_pp == NULL || public_state->acc == NULL)) ||
      (public_state->acc_profile == REVAC_ACC_PROFILE_Q96 &&
       (public_state->acc96_ctx == NULL || public_state->acc96_pp == NULL ||
        public_state->acc96_acc == NULL || uav->acc96_witness == NULL))) {
    return 0;
  }
  step_started = revac_uav_monotonic_seconds();
  const char *async_precheck_mode = getenv("REVAC_Q88_ASYNC_PRECHECK");
  const int async_precheck = public_state->acc_profile == REVAC_ACC_PROFILE_Q96 &&
      revac_labrador_online_enabled() && revac_uav_parallel_show_enabled() &&
      (async_precheck_mode != NULL && async_precheck_mode[0] == '1');
  if (async_precheck) {
    /* The proof thread checks the complete state/witness before dispatch.
       Signature work only reads the witness and does not touch these caches.
       No proof is returned until both tasks have succeeded. */
    memver_cache_hit = 0;
  } else if (public_state->acc_profile == REVAC_ACC_PROFILE_Q96) {
    revac_uav_memver_q96_cache_key(memver_key, public_state, uav);
    memver_cache_hit =
        uav->memver_cache_valid &&
        uav->memver_cache_seq == public_state->acc96_seq &&
        memcmp(uav->memver_cache_key, memver_key, sizeof(memver_key)) == 0;
    if (!memver_cache_hit) {
      if (!revac_uav_q88_prepare_checked_lift(uav, public_state)) {
        uav->memver_cache_valid = 0;
        return 0;
      }
      memcpy(uav->memver_cache_key, memver_key, sizeof(memver_key));
      uav->memver_cache_seq = public_state->acc96_seq;
      uav->memver_cache_valid = 1;
    }
  } else {
    revac_uav_memver_cache_key(memver_key, public_state, uav);
    memver_cache_hit =
        uav->memver_cache_valid &&
        uav->memver_cache_seq == public_state->acc_seq &&
        memcmp(uav->memver_cache_key, memver_key, sizeof(memver_key)) == 0;
    if (!memver_cache_hit) {
      if (!acc_genpp_memver(public_state->acc_pp, public_state->acc,
                            uav->handle, &uav->witness)) {
        uav->memver_cache_valid = 0;
        return 0;
      }
      memcpy(uav->memver_cache_key, memver_key, sizeof(memver_key));
      uav->memver_cache_seq = public_state->acc_seq;
      uav->memver_cache_valid = 1;
    }
  }
  revac_uav_print_proof_timing(
      async_precheck ? "UAV MemVer 预检查委托证明线程" :
      memver_cache_hit ? "UAV 累加器 MemVer 摘要缓存命中"
                       : "UAV 累加器 MemVer 预检查（缓存未命中）",
      step_started);

  step_started = revac_uav_monotonic_seconds();
  memcpy(proof->nonce, nonce, REVAC_SHOW_NONCE_BYTES);
  if (!revac_show_proof_set_acc_zk(proof, REVAC_ACC_ZK_NONE, NULL, 0)) {
    return 0;
  }
  {
    uint8_t acc_digest[REVAC_CONTEXT_HASH_BYTES];
    if (!revac_public_accumulator_digest(acc_digest, public_state)) return 0;
    revac_show_derive_crs_seed(crs_seed, public_state->ta_pk, uav->pk.seed,
                               proof->nonce, public_state->acc_profile,
                               public_state->acc_seq, public_state->window_id,
                               public_state->parameter_digest, acc_digest);
  }

  poly_qshow_vec_m1_init(s1);
  poly_q_mat_d_d_init(A_link);
  poly_q_mat_d_d_init(D_link);
  for (size_t i = 0; i < 2; i++) {
    poly_q_vec_d_init(link_randomness[i]);
  }
  revac_show_lang_init(A_embed, B_embed, A3_embed, Ds_embed, D_embed, Dx_embed, u_embed);
  for (size_t i = 0; i < PARAM_D; i++) {
    poly_qshow_mat_k_k_init(D_link_embed[i]);
    poly_qshow_vec_k_init(C_link_embed[i]);
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_init(A_link_embed[i][j]);
    }
  }
  revac_same_handle_public_matrices(A_link, D_link, public_state);
  revac_same_handle_commitment_generate(proof->handle_commitment,
                                        link_randomness, public_state,
                                        uav->handle);
  show_link_public_embed(A_link_embed, D_link_embed, C_link_embed,
                         A_link, D_link, proof->handle_commitment);
  revac_uav_print_proof_timing("CRS 派生与 show 语言初始化", step_started);

  step_started = revac_uav_monotonic_seconds();
  show_user_embed_handle(A_embed, B_embed, A3_embed, Ds_embed, D_embed,
                         Dx_embed, u_embed, s1, &uav->pk, &uav->sk,
                         public_state->ta_pk, &uav->sig, msg, uav->handle);
  show_user_embed_link_randomness(s1, link_randomness);
  revac_uav_print_proof_timing("匿名签名 statement/witness 嵌入", step_started);

  if (revac_labrador_online_enabled() &&
      revac_uav_parallel_show_enabled()) {
    (void)revac_labrador_memver_prove_async_start(
        &labrador_job, proof, public_state, uav, link_randomness);
  }

  /* Failed thread creation uses the same checks on the serial fallback. */
  if (async_precheck && labrador_job == NULL &&
      !revac_uav_memver_cache_prewarm(uav, public_state)) goto cleanup;

  step_started = revac_uav_monotonic_seconds();
  if (labrador_job != NULL) {
    show_user_prove_set_thread_cap(revac_uav_parallel_signature_threads());
  }
  show_user_prove_handle(&proof->sig_show, A_embed, B_embed, A3_embed,
                         Ds_embed, D_embed, Dx_embed, A_link_embed,
                         D_link_embed, C_link_embed, s1,
                         crs_seed, uav->pk.seed);
  show_user_prove_set_thread_cap(0);
  revac_uav_print_proof_timing("匿名签名 show_user_prove_handle", step_started);

  if (public_state->acc_profile == REVAC_ACC_PROFILE_Q96 &&
      !revac_labrador_online_enabled()) {
    goto cleanup;
  }
  if (revac_labrador_online_enabled()) {
    step_started = revac_uav_monotonic_seconds();
    if (labrador_job != NULL &&
        revac_labrador_memver_prove_async_finish(labrador_job, proof)) {
      labrador_job = NULL;
    } else {
      labrador_job = NULL;
      if (async_precheck &&
          !revac_uav_memver_cache_prewarm(uav, public_state)) goto cleanup;
      const int proved = public_state->acc_profile == REVAC_ACC_PROFILE_Q96
          ? revac_labrador_q96_memver_prove(
                proof, public_state, uav, link_randomness)
          : revac_labrador_memver_prove(
                proof, public_state, uav, link_randomness);
      if (!proved) {
        goto cleanup;
      }
    }
    if (proof->acc_zk_proof == NULL || proof->acc_zk_proof_len == 0) {
      goto cleanup;
    }
    revac_uav_print_proof_timing("累加器 Labrador MemVer 证明", step_started);
  }
  ok = 1;

cleanup:
  step_started = revac_uav_monotonic_seconds();
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_clear(A_link_embed[i][j]);
    }
    poly_qshow_vec_k_clear(C_link_embed[i]);
    poly_qshow_mat_k_k_clear(D_link_embed[i]);
  }
  revac_show_lang_clear(A_embed, B_embed, A3_embed, Ds_embed, D_embed, Dx_embed, u_embed);
  for (size_t i = 0; i < 2; i++) {
    poly_q_vec_d_clear(link_randomness[i]);
  }
  poly_q_mat_d_d_clear(D_link);
  poly_q_mat_d_d_clear(A_link);
  poly_qshow_vec_m1_clear(s1);
  revac_uav_print_proof_timing("show 证明临时对象清理", step_started);
  revac_uav_print_proof_timing("UAV show proof 生成总计", total_started);
  return ok;
}

int revac_uav_show_prove_wire(uint8_t **wire,
                              size_t *wire_len,
                              const revac_public_state_t *public_state,
                              revac_uav_t *uav,
                              const uint8_t msg[PARAM_M * PARAM_N / 8],
                              const uint8_t nonce[REVAC_SHOW_NONCE_BYTES]) {
  revac_show_proof_t proof;
  double step_started;
  double total_started = revac_uav_monotonic_seconds();
  int ok = 0;

  *wire = NULL;
  *wire_len = 0;
  revac_show_proof_init(&proof);
  if (!revac_uav_show_prove(&proof, public_state, uav, msg, nonce)) {
    goto cleanup;
  }
  step_started = revac_uav_monotonic_seconds();
  ok = revac_show_proof_pack_compact(wire, wire_len, &proof, 1);
  revac_uav_print_proof_timing("show packet 压缩与序列化", step_started);

cleanup:
  revac_show_proof_clear(&proof);
  revac_uav_print_proof_timing("show packet 生成总计", total_started);
  return ok;
}
