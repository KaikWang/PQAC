#include "acc96_accumulator.h"

#include <string.h>

static int handle_is_canonical(uint32_t handle) {
  return handle < (UINT32_C(1) << ACC96_HANDLE_BITS);
}

static int deletion_allowed(const acc96_accumulator_t *state) {
  return state != NULL &&
         state->sequence < ACC96_MAX_DELETIONS_PER_WINDOW;
}

static int find_member(const acc96_accumulator_t *state, uint32_t handle) {
  for (size_t i = 0; i < ACC96_MAX_MEMBERS; i++)
    if (state->members[i].status != ACC96_MEMBER_EMPTY &&
        state->members[i].handle == handle)
      return (int)i;
  return -1;
}

static int insert_member(acc96_accumulator_t *state, uint32_t handle) {
  for (size_t i = 0; i < ACC96_MAX_MEMBERS; i++) {
    if (state->members[i].status == ACC96_MEMBER_EMPTY) {
      state->members[i].handle = handle;
      state->members[i].issued_sequence = state->sequence;
      state->members[i].revoked_sequence = 0;
      state->members[i].status = ACC96_MEMBER_ACTIVE;
      return 1;
    }
  }
  return 0;
}

int acc96_accgen(acc96_accumulator_t *state,
                 const uint8_t public_seed[ACC96_SEED_BYTES],
                 const uint8_t trapdoor_seed[ACC96_SEED_BYTES],
                 const uint8_t accumulator_seed[ACC96_SEED_BYTES],
                 const acc96_ctx_t *ctx) {
  if (state == NULL || public_seed == NULL || trapdoor_seed == NULL ||
      accumulator_seed == NULL || ctx == NULL || !ctx->initialized)
    return 0;
  memset(state, 0, sizeof(*state));
  if (!acc96_relation_pp_expand(&state->pp, public_seed, ctx) ||
      !acc96_trapdoor_generate(&state->trapdoor, &state->pp,
                               trapdoor_seed, ctx))
    return 0;
  acc96_relation_uniform_poly(&state->accumulator, accumulator_seed,
                              UINT32_C(0x41393630), 0, 0);
  return 1;
}

int acc96_add(acc96_accumulator_t *state,
              acc96_relation_witness_t *witness,
              uint32_t handle,
              long double public_norm_sq_bound,
              acc96_sampleleft_fn sampler,
              void *sampler_opaque,
              const acc96_ctx_t *ctx) {
  if (state == NULL || witness == NULL || sampler == NULL ||
      !handle_is_canonical(handle) || public_norm_sq_bound <= 0 ||
      find_member(state, handle) >= 0)
    return 0;
  memset(witness, 0, sizeof(*witness));
  if (!sampler(witness, &state->pp, &state->trapdoor,
               &state->accumulator, handle, sampler_opaque, ctx) ||
      !acc96_relation_memver(&state->pp, &state->accumulator, handle,
                             witness, public_norm_sq_bound, ctx)) {
    memset(witness, 0, sizeof(*witness));
    return 0;
  }
  if (!insert_member(state, handle)) {
    memset(witness, 0, sizeof(*witness));
    return 0;
  }
  return 1;
}

int acc96_delete(acc96_accumulator_t *state,
                 acc96_relation_update_t *update,
                 uint32_t handle,
                 const acc96_ctx_t *ctx) {
  const int position = deletion_allowed(state) && handle_is_canonical(handle)
                           ? find_member(state, handle) : -1;
  if (position < 0 || state->members[position].status != ACC96_MEMBER_ACTIVE ||
      update == NULL)
    return 0;
  state->sequence++;
  acc96_relation_delete(&state->accumulator, update, &state->pp, handle,
                        state->sequence, ctx);
  state->members[position].status = ACC96_MEMBER_REVOKED;
  state->members[position].revoked_sequence = state->sequence;
  return 1;
}

int acc96_delete_delta4(acc96_accumulator_t *state,
                        acc96_relation_update_t *update,
                        acc96_delta4_t *delta,
                        uint32_t handle,
                        const acc96_ctx_t *ctx) {
  const int position = deletion_allowed(state) && handle_is_canonical(handle)
                           ? find_member(state, handle) : -1;
  if (position < 0 || state->members[position].status != ACC96_MEMBER_ACTIVE ||
      update == NULL || delta == NULL)
    return 0;
  if (!acc96_relation_delete_delta4(
          &state->accumulator, update, delta, &state->pp, handle,
          state->sequence + 1, ctx))
    return 0;
  state->sequence++;
  state->members[position].status = ACC96_MEMBER_REVOKED;
  state->members[position].revoked_sequence = state->sequence;
  return 1;
}

int acc96_delete_delta6(acc96_accumulator_t *state,
                        acc96_relation_update_t *update,
                        acc96_delta6_t *delta,
                        uint32_t handle,
                        const acc96_ctx_t *ctx) {
  const int position = deletion_allowed(state) && handle_is_canonical(handle)
                           ? find_member(state, handle) : -1;
  if (position < 0 || state->members[position].status != ACC96_MEMBER_ACTIVE ||
      update == NULL || delta == NULL)
    return 0;
  if (!acc96_relation_delete_delta6(
          &state->accumulator, update, delta, &state->pp, handle,
          state->sequence + 1, ctx))
    return 0;
  state->sequence++;
  state->members[position].status = ACC96_MEMBER_REVOKED;
  state->members[position].revoked_sequence = state->sequence;
  return 1;
}

int acc96_delete_tail4(acc96_accumulator_t *state,
                       acc96_relation_update_t *update,
                       acc96_gmat_t *tail4, uint32_t handle,
                       const acc96_ctx_t *ctx) {
  const int position = deletion_allowed(state) && handle_is_canonical(handle)
                           ? find_member(state, handle) : -1;
  if (state == NULL || update == NULL || tail4 == NULL ||
      !deletion_allowed(state))
    return 0;
  if (position < 0 ||
      state->members[position].status != ACC96_MEMBER_ACTIVE)
    return 0;
  if (!acc96_relation_delete_tail4(
          &state->accumulator, update, tail4, &state->pp, handle,
          state->sequence + 1, ctx))
    return 0;
  state->sequence++;
  state->members[position].status = ACC96_MEMBER_REVOKED;
  state->members[position].revoked_sequence = state->sequence;
  return 1;
}

int acc96_memwitup(acc96_relation_witness_t *witness,
                   const acc96_accumulator_t *state,
                   uint32_t member_handle,
                   const acc96_relation_update_t *update,
                   const acc96_ctx_t *ctx) {
  const int position = state != NULL && handle_is_canonical(member_handle)
                           ? find_member(state, member_handle) : -1;
  if (position < 0 || state->members[position].status != ACC96_MEMBER_ACTIVE ||
      update == NULL || update->sequence > state->sequence)
    return 0;
  return acc96_relation_memwitup(witness, &state->pp, member_handle,
                                 update, ctx);
}

int acc96_memver(const acc96_accumulator_t *state,
                 uint32_t handle,
                 const acc96_relation_witness_t *witness,
                 long double public_norm_sq_bound,
                 const acc96_ctx_t *ctx) {
  const int position = state != NULL && handle_is_canonical(handle)
                           ? find_member(state, handle) : -1;
  return position >= 0 &&
         state->members[position].status == ACC96_MEMBER_ACTIVE &&
         acc96_relation_memver(&state->pp, &state->accumulator, handle,
                               witness, public_norm_sq_bound, ctx);
}
