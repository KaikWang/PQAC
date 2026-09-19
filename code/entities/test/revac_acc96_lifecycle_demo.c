#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "acc96_accumulator.h"
#include "acc96_sampleleft_kat.h"

static double seconds(void) {
  struct timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);
  return now.tv_sec + now.tv_nsec / 1e9;
}

int main(void) {
  const uint8_t public_seed[ACC96_SEED_BYTES] = {0x96, 0x50, 0x01};
  const uint8_t trapdoor_seed[ACC96_SEED_BYTES] = {0x96, 0x54, 0x01};
  const uint8_t accumulator_seed[ACC96_SEED_BYTES] = {0x96, 0x41, 0x01};
  const uint32_t revoked = UINT32_C(0x00203040);
  const uint32_t survivor = UINT32_C(0x00203041);
  const long double kat_bound = 1.0e40L;
  acc96_ctx_t ctx;
  acc96_accumulator_t *state = calloc(1, sizeof(*state));
  acc96_relation_witness_t *revoked_witness =
      calloc(1, sizeof(*revoked_witness));
  acc96_relation_witness_t *survivor_witness =
      calloc(1, sizeof(*survivor_witness));
  acc96_relation_update_t update;
  acc96_relation_update_t blocked_update;
  acc96_relation_update_t blocked_update_before;
  acc96_delta4_t *blocked_delta = calloc(1, sizeof(*blocked_delta));
  acc96_gmat_t *blocked_tail = calloc(1, sizeof(*blocked_tail));
  acc96_poly_t accumulator_at_limit;
  acc96_member_entry_t survivor_at_limit;
  double started;
  double setup_time;
  double add_time;
  double delete_time;
  double update_time;
  if (state == NULL || revoked_witness == NULL || survivor_witness == NULL ||
      blocked_delta == NULL || blocked_tail == NULL ||
      !acc96_ctx_init(&ctx))
    return 1;
  started = seconds();
  if (!acc96_accgen(state, public_seed, trapdoor_seed, accumulator_seed, &ctx))
    return 1;
  setup_time = seconds() - started;

  /* Missing samplers and witnesses over a caller-supplied too-small norm bound
   * must fail before dictionary insertion. */
  if (acc96_add(state, revoked_witness, revoked, kat_bound, NULL, NULL, &ctx) ||
      acc96_add(state, revoked_witness, revoked, 1.0L,
                acc96_sampleleft_algebraic_kat, NULL, &ctx))
    return 1;
  started = seconds();
  if (!acc96_add(state, revoked_witness, revoked, kat_bound,
                 acc96_sampleleft_algebraic_kat, NULL, &ctx) ||
      !acc96_add(state, survivor_witness, survivor, kat_bound,
                 acc96_sampleleft_algebraic_kat, NULL, &ctx))
    return 1;
  add_time = (seconds() - started) / 2.0;
  if (acc96_add(state, revoked_witness, revoked, kat_bound,
                acc96_sampleleft_algebraic_kat, NULL, &ctx))
    return 1;
  if (!acc96_memver(state, revoked, revoked_witness, kat_bound, &ctx) ||
      !acc96_memver(state, survivor, survivor_witness, kat_bound, &ctx))
    return 1;

  /* Model the last admissible deletion in a window.  The sequence is only
   * authenticated metadata; the algebraic update itself is independent of
   * the number of preceding events. */
  state->sequence = ACC96_MAX_DELETIONS_PER_WINDOW - 1;
  started = seconds();
  if (!acc96_delete(state, &update, revoked, &ctx)) return 1;
  delete_time = seconds() - started;
  if (state->sequence != ACC96_MAX_DELETIONS_PER_WINDOW ||
      update.sequence != ACC96_MAX_DELETIONS_PER_WINDOW)
    return 1;
  started = seconds();
  if (!acc96_memwitup(survivor_witness, state, survivor, &update, &ctx))
    return 1;
  update_time = seconds() - started;
  if (acc96_delete(state, &update, revoked, &ctx) ||
      acc96_memver(state, revoked, revoked_witness, kat_bound, &ctx) ||
      !acc96_memver(state, survivor, survivor_witness, kat_bound, &ctx))
    return 1;

  /* All deletion entry points must reject event 10001 before changing the
   * accumulator, member dictionary, sequence, or caller-owned update. */
  acc96_poly_set(&accumulator_at_limit, &state->accumulator);
  survivor_at_limit = state->members[1];
  memset(&blocked_update, 0xa5, sizeof(blocked_update));
  blocked_update_before = blocked_update;
  if (acc96_delete(state, &blocked_update, survivor, &ctx) ||
      acc96_delete_delta4(state, &blocked_update, blocked_delta, survivor,
                          &ctx) ||
      acc96_delete_tail4(state, &blocked_update, blocked_tail, survivor,
                         &ctx) ||
      state->sequence != ACC96_MAX_DELETIONS_PER_WINDOW ||
      !acc96_poly_equal(&state->accumulator, &accumulator_at_limit) ||
      memcmp(&state->members[1], &survivor_at_limit,
             sizeof(survivor_at_limit)) != 0 ||
      memcmp(&blocked_update, &blocked_update_before,
             sizeof(blocked_update)) != 0)
    return 1;

  const long double updated_norm =
      acc96_relation_witness_norm2(survivor_witness, &ctx);

  acc96_coeff_from_u128(&survivor_witness->left.coeff[0],
      (acc96_coeff_to_u128(&survivor_witness->left.coeff[0], &ctx) + 1) %
          acc96_modulus());
  if (acc96_memver(state, survivor, survivor_witness, kat_bound, &ctx))
    return 1;

  printf("q88 lifecycle passed: AccGen/Add/Delete/MemWitUp/MemVer; "
         "missing sampler, low norm, duplicate add/delete, revoked member "
         "and tampering all rejected; deletion 10000 accepted and deletion "
         "10001 rejected without state change. Setup %.3f s, algebraic SampleLeft/Add "
         "%.3f s/member, Delete %.3f s, MemWitUp %.3f s, updated norm^2 "
         "%.0Lf. Algebraic SampleLeft was KAT-only.\n",
         setup_time, add_time, delete_time, update_time, updated_norm);
  free(survivor_witness);
  free(revoked_witness);
  free(blocked_tail);
  free(blocked_delta);
  free(state);
  return 0;
}
