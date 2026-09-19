#ifndef ACC96_ACCUMULATOR_H
#define ACC96_ACCUMULATOR_H

#include <stddef.h>
#include <stdint.h>

#include "acc96_trapdoor.h"

#define ACC96_MAX_MEMBERS 10000U

typedef enum {
  ACC96_MEMBER_EMPTY = 0,
  ACC96_MEMBER_ACTIVE = 1,
  ACC96_MEMBER_REVOKED = 2
} acc96_member_status_t;

typedef struct {
  uint32_t handle;
  uint32_t issued_sequence;
  uint32_t revoked_sequence;
  acc96_member_status_t status;
} acc96_member_entry_t;

typedef struct {
  acc96_relation_pp_t pp;
  acc96_trapdoor_t trapdoor;
  acc96_poly_t accumulator;
  acc96_member_entry_t members[ACC96_MAX_MEMBERS];
  uint32_t sequence;
} acc96_accumulator_t;

/* Real q88 SampleLeft implementations must sample a short witness for the
 * supplied fixed accumulator and handle. */
typedef int (*acc96_sampleleft_fn)(
    acc96_relation_witness_t *out,
    const acc96_relation_pp_t *pp,
    const acc96_trapdoor_t *trapdoor,
    const acc96_poly_t *accumulator,
    uint32_t handle,
    void *opaque,
    const acc96_ctx_t *ctx);

int acc96_accgen(acc96_accumulator_t *state,
                 const uint8_t public_seed[ACC96_SEED_BYTES],
                 const uint8_t trapdoor_seed[ACC96_SEED_BYTES],
                 const uint8_t accumulator_seed[ACC96_SEED_BYTES],
                 const acc96_ctx_t *ctx);

int acc96_add(acc96_accumulator_t *state,
              acc96_relation_witness_t *witness,
              uint32_t handle,
              long double public_norm_sq_bound,
              acc96_sampleleft_fn sampler,
              void *sampler_opaque,
              const acc96_ctx_t *ctx);

int acc96_delete(acc96_accumulator_t *state,
                 acc96_relation_update_t *update,
                 uint32_t handle,
                 const acc96_ctx_t *ctx);
int acc96_delete_delta4(acc96_accumulator_t *state,
                        acc96_relation_update_t *update,
                        acc96_delta4_t *delta,
                        uint32_t handle,
                        const acc96_ctx_t *ctx);
int acc96_delete_delta6(acc96_accumulator_t *state,
                        acc96_relation_update_t *update,
                        acc96_delta6_t *delta,
                        uint32_t handle,
                        const acc96_ctx_t *ctx);
int acc96_delete_tail4(acc96_accumulator_t *state,
                       acc96_relation_update_t *update,
                       acc96_gmat_t *tail4, uint32_t handle,
                       const acc96_ctx_t *ctx);

int acc96_memwitup(acc96_relation_witness_t *witness,
                   const acc96_accumulator_t *state,
                   uint32_t member_handle,
                   const acc96_relation_update_t *update,
                   const acc96_ctx_t *ctx);

int acc96_memver(const acc96_accumulator_t *state,
                 uint32_t handle,
                 const acc96_relation_witness_t *witness,
                 long double public_norm_sq_bound,
                 const acc96_ctx_t *ctx);

#endif
