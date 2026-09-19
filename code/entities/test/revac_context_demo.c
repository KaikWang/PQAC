#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "arith.h"
#include "random.h"
#include "randombytes.h"
#include "revac_context.h"
#include "revac_ta.h"
#include "revac_uav.h"

static double monotonic_seconds(void) {
  struct timespec now;
  (void)clock_gettime(CLOCK_MONOTONIC, &now);
  return (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
}

static int expect(int condition, const char *message) {
  if (!condition) {
    fprintf(stderr, "FAILED: %s\n", message);
    return 0;
  }
  return 1;
}

static int fetch_from_ta(void *opaque,
                         uint32_t seq,
                         revac_context_profile_t profile,
                         uint8_t **wire,
                         size_t *wire_len) {
  return revac_ta_export_context((const revac_ta_t *)opaque, seq, profile,
                                 wire, wire_len);
}

static void witness_set(acc_genpp_wit_t *dst, const acc_genpp_wit_t *src) {
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

static int witness_equal(const acc_genpp_wit_t *lhs,
                         const acc_genpp_wit_t *rhs) {
  if (!poly_q_vec_d_equal(lhs->left, rhs->left) ||
      !poly_q_vec_d_equal(lhs->v12, rhs->v12) ||
      !poly_q_equal(lhs->tag, rhs->tag)) {
    return 0;
  }
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    if (!poly_q_vec_d_equal(lhs->v2[i], rhs->v2[i])) {
      return 0;
    }
  }
  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      if (!poly_q_vec_d_equal(lhs->right[bit][j], rhs->right[bit][j])) {
        return 0;
      }
    }
  }
  return 1;
}

static void update_clone(revac_uav_t *dst, const revac_uav_t *src) {
  dst->id = src->id;
  dst->handle = src->handle;
  dst->last_seq = src->last_seq;
  dst->registered = src->registered;
  dst->issued = src->issued;
  witness_set(&dst->witness, &src->witness);
}

static size_t sync_delta_cache_count(const revac_ta_t *ta) {
  size_t count = 0;

  for (size_t i = 0; i < ACC_GENPP_MAX_DELETIONS; i++) {
    if (ta->sync_delta_cache[i].used) {
      count++;
    }
  }
  return count;
}

static int member_baseline_is(const revac_ta_t *ta, uint32_t handle,
                              uint32_t baseline_seq) {
  for (size_t i = 0; i < REVAC_MAX_UAVS; i++) {
    if (ta->uavs[i].used && ta->uavs[i].handle == handle) {
      return ta->uavs[i].baseline_seq == baseline_seq;
    }
  }
  return 0;
}

static int run_demo(void) {
  revac_ta_t ta;
  revac_public_state_t public_state;
  revac_uav_t revoked, member, reference, candidate, probe_local, probe_delta;
  acc_genpp_wit_t cached_sync, fallback_sync;
  acc_genpp_update_t update;
  revac_window_variance_ledger_t variance_ledger;
  uint8_t message0[PARAM_M * PARAM_N / 8];
  uint8_t message1[PARAM_M * PARAM_N / 8];
  uint8_t *wire = NULL;
  size_t wire_len = 0;
  size_t local_composites = 0;
  double started;
  int ok = 0;

  revac_ta_init(&ta);
  revac_uav_init(&revoked);
  revac_uav_init(&member);
  revac_uav_init(&reference);
  revac_uav_init(&candidate);
  revac_uav_init(&probe_local);
  revac_uav_init(&probe_delta);
  acc_genpp_wit_init(&cached_sync);
  acc_genpp_wit_init(&fallback_sync);

  revac_ta_setup_q19(&ta);
  revac_ta_public_state(&public_state, &ta);
  revac_uav_keygen(&revoked, UINT64_C(0x7100000000000001), ta.ta_pk.seed);
  revac_uav_keygen(&member, UINT64_C(0x7100000000000002), ta.ta_pk.seed);
  randombytes(message0, sizeof(message0));
  randombytes(message1, sizeof(message1));
  if (!expect(revac_ta_register_uav(&ta, &revoked),
              "revoked UAV registration") ||
      !expect(revac_ta_register_uav(&ta, &member),
              "member UAV registration") ||
      !expect(revac_ta_issue_uav(&ta, &revoked, message0),
              "revoked UAV issuance") ||
      !expect(revac_ta_issue_uav(&ta, &member, message1),
              "member UAV issuance")) {
    goto cleanup;
  }
  revac_window_variance_ledger_init(&variance_ledger, &public_state,
                                    member.handle, 1, 500);
  if (!expect(variance_ledger.initialized,
              "window variance ledger initialization")) {
    goto cleanup;
  }

  update_clone(&reference, &member);
  started = monotonic_seconds();
  if (!expect(revac_ta_revoke_uav_context(&ta, revoked.handle, &update),
              "TA context revocation") ||
      !expect(update.seq == 1 && update.handle == revoked.handle &&
                  acc_genpp_update_sign_is_explicit(&update),
              "context update metadata and explicit authenticated sign")) {
    goto cleanup;
  }
  printf("TA Delete + context capture/cache: %.3f s\n",
         monotonic_seconds() - started);
  printf("TA private DELTA6 cache: %.2f KiB for seq=%u\n",
         ta.sync_delta_cache[update.seq - 1].delta_data_len / 1024.0,
         update.seq);
  revac_ta_public_state(&public_state, &ta);
  if (!expect(public_state.bml_len == 1 && public_state.bml[0] != NULL &&
                  public_state.bml[0]->update.seq == update.seq &&
                  public_state.bml[0]->update.handle == update.handle &&
                  public_state.bml[0]->delta_data != NULL &&
                  public_state.bml[0]->delta_data_len > 0,
              "BML carries the revocation update and DELTA4")) {
    goto cleanup;
  }

  started = monotonic_seconds();
  if (!expect(revac_uav_apply_update_local(&reference, &public_state, &update,
                                           &local_composites),
              "local update") ||
      !expect(acc_genpp_memver(&ta.acc_pp, ta.acc, reference.handle,
                               &reference.witness),
              "local reference MemVer")) {
    goto cleanup;
  }
  printf("local: %zu composites, %.3f s\n", local_composites,
         monotonic_seconds() - started);

  {
    const size_t slot = update.seq - 1;
    const int saved_used = ta.sync_delta_cache[slot].used;
    size_t cached_updates = 0;
    size_t cached_composites = 0;
    size_t fallback_updates = 0;
    size_t fallback_composites = 0;
    int fallback_ok;

    started = monotonic_seconds();
    if (!expect(revac_ta_memwitsync_cached(
                    &cached_sync, &ta, member.handle, &cached_updates,
                    &cached_composites),
                "TA cached MemWitSync") ||
        !expect(cached_updates == 1, "TA cache hit count") ||
        !expect(acc_genpp_memver(&ta.acc_pp, ta.acc, member.handle,
                                 &cached_sync),
                "TA cached MemWitSync witness")) {
      goto cleanup;
    }
    printf("TA DELTA6 MemWitSync: %zu cached update, %zu composites, %.3f s\n",
           cached_updates, cached_composites,
           monotonic_seconds() - started);

    ta.sync_delta_cache[slot].used = 0;
    started = monotonic_seconds();
    fallback_ok = revac_ta_memwitsync_cached(
        &fallback_sync, &ta, member.handle, &fallback_updates,
        &fallback_composites);
    ta.sync_delta_cache[slot].used = saved_used;
    if (!expect(fallback_ok, "TA cache-miss MemWitSync fallback") ||
        !expect(fallback_updates == 0, "TA fallback cache hit count") ||
        !expect(fallback_composites > 0, "TA fallback local composites") ||
        !expect(acc_genpp_memver(&ta.acc_pp, ta.acc, member.handle,
                                 &fallback_sync),
                "TA fallback MemWitSync witness")) {
      goto cleanup;
    }
    printf("TA forced DELTA6 miss: %zu cached updates, %zu composites, %.3f s\n",
           fallback_updates, fallback_composites,
           monotonic_seconds() - started);

    {
      revac_sync_delta_cache_entry_t *entry = &ta.sync_delta_cache[slot];
      const size_t last = entry->delta_data_len - 1;
      const uint8_t saved_byte = entry->delta_data[last];

      entry->delta_data[last] ^= UINT8_C(1);
      cached_updates = 0;
      cached_composites = 0;
      fallback_ok = revac_ta_memwitsync_cached(
          &fallback_sync, &ta, member.handle, &cached_updates,
          &cached_composites);
      entry->delta_data[last] = saved_byte;
      if (!expect(fallback_ok, "corrupt DELTA6 local fallback") ||
          !expect(cached_updates == 0 && cached_composites > 0,
                  "corrupt DELTA6 was not consumed") ||
          !expect(acc_genpp_memver(&ta.acc_pp, ta.acc, member.handle,
                                   &fallback_sync),
                  "corrupt DELTA6 fallback witness")) {
        goto cleanup;
      }
    }
  }

  {
    revac_context_profile_t requested = REVAC_CONTEXT_PROFILE_DELTA4;
    revac_context_profile_t used = REVAC_CONTEXT_PROFILE_LOCAL;
    size_t composites = 0;
    const int delta4_supported =
        ((member.handle ^ update.handle) & UINT32_C(0x0f)) != 0;

    update_clone(&candidate, &member);
    if (!expect(revac_ta_export_context(&ta, update.seq, requested,
                                        &wire, &wire_len),
                "TA DELTA4 export") ||
        !expect(wire_len >= revac_context_profile_wire_bytes(requested),
                "DELTA4 wire size")) {
      goto cleanup;
    }

    /* The sign is security-critical: changing +1 to -1 (or vice versa) in
       the reserved header byte must invalidate the signed manifest. */
    {
      const uint8_t saved_sign = wire[14];
      wire[14] = saved_sign == UINT8_C(1) ? UINT8_C(2) : UINT8_C(1);
      if (!expect(!revac_context_packet_verify(
                      wire, wire_len, &update, public_state.context_pp_hash,
                      public_state.ta_pk),
                  "tampered update-sign rejection")) {
        wire[14] = saved_sign;
        goto cleanup;
      }
      wire[14] = saved_sign;
    }
    {
      acc_genpp_update_t forged_update = update;
      forged_update.sign = (int8_t)-forged_update.sign;
      if (!expect(!revac_window_variance_ledger_record(
                      &variance_ledger, &public_state, &forged_update,
                      wire, wire_len, &member.witness, &reference.witness) &&
                      variance_ledger.deletion_count == 0,
                  "forged sign cannot enter variance ledger")) {
        goto cleanup;
      }
    }
    if (!expect(revac_window_variance_ledger_record(
                    &variance_ledger, &public_state, &update, wire, wire_len,
                    &member.witness, &reference.witness),
                "authenticated window variance record") ||
        !expect(variance_ledger.deletion_count == 1 &&
                    variance_ledger.first_seq == update.seq &&
                    variance_ledger.last_seq == update.seq &&
                    variance_ledger.quadratic_variation > 0 &&
                    variance_ledger.max_increment_norm2 ==
                        variance_ledger.quadratic_variation,
                "window variance V/L accounting")) {
      goto cleanup;
    }
    started = monotonic_seconds();
    if (!expect(revac_uav_apply_update_context(
                    &candidate, &public_state, &update, wire, wire_len, 1, &used,
                    &composites),
                "UAV DELTA4/local update") ||
        !expect(used == (delta4_supported ? requested
                                          : REVAC_CONTEXT_PROFILE_LOCAL),
                "UAV selected DELTA4 or safe LOCAL fallback") ||
        !expect(witness_equal(&candidate.witness, &reference.witness),
                "DELTA4 witness equals LOCAL witness") ||
        !expect(acc_genpp_memver(&ta.acc_pp, ta.acc, candidate.handle,
                                 &candidate.witness),
                "DELTA4 witness MemVer")) {
      goto cleanup;
    }
    printf("%-8s %8.2f KiB, %2zu local composites, %.3f s\n",
           revac_context_profile_name(requested), wire_len / 1024.0,
           composites, monotonic_seconds() - started);

    /* Exercise every DELTA4 stopping layer independently of the random
       registered handles.  The base witness need not match the synthetic
       handle for this exact transition-equivalence check. */
    for (size_t first_difference = 0;
         first_difference < ACC_GENPP_DELTA4_LAYERS;
         first_difference++) {
      const uint32_t synthetic_handle =
          update.handle ^ (UINT32_C(1) << first_difference);

      update_clone(&probe_local, &member);
      update_clone(&probe_delta, &member);
      probe_local.handle = synthetic_handle;
      probe_delta.handle = synthetic_handle;
      if (!expect(revac_uav_apply_update_local(
                      &probe_local, &public_state, &update, NULL),
                  "synthetic LOCAL transition") ||
          !expect(revac_uav_apply_update_context(
                      &probe_delta, &public_state, &update, wire, wire_len, 0,
                      &used, NULL),
                  "synthetic DELTA4 transition") ||
          !expect(used == REVAC_CONTEXT_PROFILE_DELTA4,
                  "synthetic transition used DELTA4") ||
          !expect(witness_equal(&probe_local.witness, &probe_delta.witness),
                  "every DELTA4 layer equals LOCAL")) {
        goto cleanup;
      }
    }

    /* A first difference at bit 4 is deliberately outside DELTA4.  Strict
       mode must reject without mutation; fallback must then equal LOCAL. */
    update_clone(&probe_local, &member);
    update_clone(&probe_delta, &member);
    probe_local.handle = update.handle ^ UINT32_C(0x10);
    probe_delta.handle = probe_local.handle;
    if (!expect(!revac_uav_apply_update_context(
                    &probe_delta, &public_state, &update, wire, wire_len, 0,
                    NULL, NULL),
                "DELTA4 unsupported-prefix strict rejection") ||
        !expect(probe_delta.last_seq == member.last_seq &&
                    witness_equal(&probe_delta.witness, &member.witness),
                "DELTA4 rejection is non-mutating") ||
        !expect(revac_uav_apply_update_local(
                    &probe_local, &public_state, &update, NULL),
                "unsupported-prefix LOCAL reference") ||
        !expect(revac_uav_apply_update_context(
                    &probe_delta, &public_state, &update, wire, wire_len, 1,
                    NULL, NULL),
                "unsupported-prefix LOCAL fallback") ||
        !expect(witness_equal(&probe_delta.witness, &probe_local.witness),
                "DELTA4 fallback equals LOCAL")) {
      goto cleanup;
    }
    free(wire);
    wire = NULL;
    wire_len = 0;
  }

  /* A corrupt authenticated chunk must not mutate state. */
  update_clone(&candidate, &member);
  if (!expect(revac_ta_export_context(&ta, update.seq,
                                      REVAC_CONTEXT_PROFILE_DELTA4,
                                      &wire, &wire_len),
              "tamper profile export")) {
    goto cleanup;
  }
  wire[wire_len - 1] ^= UINT8_C(1);
  if (!expect(!revac_uav_apply_update_context(
                  &candidate, &public_state, &update, wire, wire_len, 0, NULL, NULL),
              "tampered context rejection") ||
      !expect(candidate.last_seq == member.last_seq &&
                  witness_equal(&candidate.witness, &member.witness),
              "tampered context is non-mutating") ||
      !expect(revac_uav_apply_update_context(
                  &candidate, &public_state, &update, wire, wire_len, 1, NULL,
                  &local_composites),
              "tampered context local fallback") ||
      !expect(witness_equal(&candidate.witness, &reference.witness),
              "fallback witness equals local witness")) {
    goto cleanup;
  }
  free(wire);
  wire = NULL;
  wire_len = 0;

  /* Auto mode and BML recovery both use the same authenticated cache. */
  update_clone(&candidate, &member);
  {
    revac_context_profile_t selected =
        revac_context_profile_select(100.0, 158.6093566);
    revac_context_profile_t used = REVAC_CONTEXT_PROFILE_LOCAL;
    size_t bytes = 0;
    size_t composites = 0;
    if (!expect(revac_uav_apply_update_auto(
                    &candidate, &public_state, &update, 100.0, 158.6093566,
                    fetch_from_ta, &ta, &used, &composites, &bytes),
                "auto update") ||
        !expect(used == selected || used == REVAC_CONTEXT_PROFILE_LOCAL,
                "auto selected DELTA4 or fallback LOCAL") ||
        !expect(bytes > 0, "auto DELTA4 download byte count") ||
        !expect(witness_equal(&candidate.witness, &reference.witness),
                "auto witness equals local witness")) {
      goto cleanup;
    }
    printf("auto@100Mbps: %s, %.2f KiB, %zu composites\n",
           revac_context_profile_name(used), bytes / 1024.0, composites);
  }

  /* Even an explicitly LOCAL receiver consumes and authenticates the BML
     DELTA4 broadcast before choosing its local witness recurrence. */
  update_clone(&candidate, &member);
  {
    revac_context_profile_t used = REVAC_CONTEXT_PROFILE_DELTA4;
    size_t bytes = 0;
    size_t composites = 0;

    if (!expect(revac_uav_apply_update_auto(
                    &candidate, &public_state, &update, 0.0, 158.6093566,
                    NULL, NULL, &used, &composites, &bytes),
                "mandatory BML DELTA4 LOCAL update") ||
        !expect(used == REVAC_CONTEXT_PROFILE_LOCAL,
                "zero-goodput selector chose LOCAL") ||
        !expect(bytes >= revac_context_profile_wire_bytes(
                             REVAC_CONTEXT_PROFILE_DELTA4),
                "LOCAL receiver downloaded BML DELTA4") ||
        !expect(composites > 0,
                "LOCAL receiver did not execute local recurrence") ||
        !expect(witness_equal(&candidate.witness, &reference.witness),
                "mandatory-broadcast LOCAL witness equals reference")) {
      goto cleanup;
    }
    printf("auto@LOCAL: mandatory %.2f KiB DELTA4, %zu composites\n",
           bytes / 1024.0, composites);
  }

  update_clone(&candidate, &member);
  {
    size_t bytes = 0;
    size_t composites = 0;
    if (!expect(revac_uav_recover_updates_auto(
                    &candidate, &public_state, 20.0, 158.6093566,
                    fetch_from_ta, &ta, NULL, NULL, &bytes, &composites),
                "BML auto recovery") ||
        !expect(candidate.last_seq == update.seq,
                "BML auto sequence") ||
        !expect(witness_equal(&candidate.witness, &reference.witness),
                "BML auto witness equals local witness")) {
      goto cleanup;
    }
    printf("BML auto@20Mbps: %.2f KiB, %zu composites\n",
           bytes / 1024.0, composites);
  }

  printf("profile selector (only LOCAL/DELTA4):\n");
  {
    static const double links[] = {0.5, 1.0, 2.0, 5.0, 10.0,
                                   20.0, 50.0, 100.0};
    for (size_t i = 0; i < sizeof(links) / sizeof(links[0]); i++) {
      revac_context_profile_t profile =
          revac_context_profile_select(links[i], 158.6093566);
      printf("  %6.1f Mbps -> %-27s %.2f KiB\n", links[i],
             revac_context_profile_name(profile),
             revac_context_profile_wire_bytes(profile) / 1024.0);
    }
  }
  ok = 1;

cleanup:
  free(wire);
  acc_genpp_wit_clear(&fallback_sync);
  acc_genpp_wit_clear(&cached_sync);
  revac_uav_clear(&candidate);
  revac_uav_clear(&probe_delta);
  revac_uav_clear(&probe_local);
  revac_uav_clear(&reference);
  revac_uav_clear(&member);
  revac_uav_clear(&revoked);
  revac_ta_clear(&ta);
  return ok;
}

static int run_sync_delta_gc_demo(void) {
  revac_ta_t ta;
  revac_uav_t source, early, late;
  acc_genpp_update_t update1, update2, update3;
  uint8_t message[PARAM_M * PARAM_N / 8];
  uint8_t *wire = NULL;
  size_t wire_len = 0;
  int ok = 0;

  revac_ta_init(&ta);
  revac_uav_init(&source);
  revac_uav_init(&early);
  revac_uav_init(&late);
  revac_ta_setup_q19(&ta);
  revac_uav_keygen(&source, UINT64_C(0x7200000000000001), ta.ta_pk.seed);
  revac_uav_keygen(&early, UINT64_C(0x7200000000000002), ta.ta_pk.seed);
  revac_uav_keygen(&late, UINT64_C(0x7200000000000003), ta.ta_pk.seed);
  randombytes(message, sizeof(message));

  if (!expect(revac_ta_register_uav(&ta, &source),
              "GC source registration") ||
      !expect(revac_ta_register_uav(&ta, &early),
              "GC early member registration") ||
      !expect(revac_ta_register_uav(&ta, &late),
              "GC late member registration") ||
      !expect(revac_ta_issue_uav(&ta, &source, message),
              "GC source issuance") ||
      !expect(revac_ta_issue_uav(&ta, &early, message),
              "GC early member issuance") ||
      !expect(member_baseline_is(&ta, early.handle, 0),
              "early member baseline")) {
    goto cleanup;
  }

  if (!expect(revac_ta_revoke_uav_context(&ta, source.handle, &update1),
              "GC first revocation") ||
      !expect(update1.seq == 1 && ta.sync_delta_cache[0].used,
              "first DELTA6 retained for early member") ||
      !expect(revac_ta_issue_uav(&ta, &late, message),
              "GC late member issuance") ||
      !expect(member_baseline_is(&ta, late.handle, 1),
              "late member baseline")) {
    goto cleanup;
  }

  /* Revoking the only baseline-0 survivor advances the safe watermark to
     late.baseline_seq=1.  DELTA6(seq=1) can be collected, while seq=2 must
     remain for the late member.  The public DELTA4/BML history is separate. */
  if (!expect(revac_ta_revoke_uav_context(&ta, early.handle, &update2),
              "GC second revocation") ||
      !expect(update2.seq == 2, "second revocation sequence") ||
      !expect(!ta.sync_delta_cache[0].used &&
                  ta.sync_delta_cache[1].used &&
                  sync_delta_cache_count(&ta) == 1,
              "baseline watermark collected only seq=1") ||
      !expect(ta.acc_st.bml_len == 2,
              "private GC did not alter public BML") ||
      !expect(revac_ta_export_context(&ta, update1.seq,
                                      REVAC_CONTEXT_PROFILE_DELTA4,
                                      &wire, &wire_len),
              "private GC did not remove public DELTA4")) {
    goto cleanup;
  }
  free(wire);
  wire = NULL;
  wire_len = 0;

  /* With no active member left, every private replay delta is unreachable
     and may be collected immediately. */
  if (!expect(revac_ta_revoke_uav_context(&ta, late.handle, &update3),
              "GC final revocation") ||
      !expect(update3.seq == 3 && sync_delta_cache_count(&ta) == 0,
              "no-active-member GC cleared all private deltas") ||
      !expect(ta.acc_st.bml_len == 3,
              "final private GC kept public BML intact")) {
    goto cleanup;
  }

  printf("DELTA6 safe GC: baseline watermark and no-active-member cleanup passed\n");
  ok = 1;

cleanup:
  free(wire);
  revac_uav_clear(&late);
  revac_uav_clear(&early);
  revac_uav_clear(&source);
  revac_ta_clear(&ta);
  return ok;
}

int main(void) {
  int ok;

  arith_setup();
  random_init();
  ok = run_demo() && run_sync_delta_gc_demo();
  arith_teardown();
  if (!ok) {
    return 1;
  }
  printf("revac_context_demo passed.\n");
  return 0;
}
