#ifndef ACC96_RELATION_H
#define ACC96_RELATION_H

#include <stddef.h>
#include <stdint.h>

#include "acc96_eval.h"

/* q88 rank-one specialization of the accumulator MemVer relation.  This
 * module deliberately contains no issuing sampler: a witness must be supplied
 * by the q88 truncated-gadget SampleLeft implementation.  Keeping that
 * boundary explicit prevents a relation-only test witness from being used as
 * a production credential. */
#define ACC96_TRUNC_BASE UINT64_C(13)
#define ACC96_TRUNC_K 24U
#define ACC96_TRUNC_L 2U
#define ACC96_TRUNC_RIGHT (ACC96_TRUNC_K - ACC96_TRUNC_L)
#define ACC96_SEED_BYTES 32U
#define ACC96_DELTA4_LAYERS 4U
#define ACC96_DELTA6_LAYERS 6U
#define ACC96_PUBLIC_DELTA_LAYERS ACC96_DELTA6_LAYERS
#define ACC96_TAIL4_BIT 4U
#define ACC96_REFRESH_BITS 128U
#define ACC96_REFRESH_BYTES (ACC96_REFRESH_BITS / 8U)
#define ACC96_REFRESH_BLOCKS ACC96_GADGET_DIGITS
/* Accepted fresh samples: 24 head polynomials are bounded by F; right/z
 * remain ternary. This is separate from the larger post-update bound. */
#define ACC96_FRESH_COEFF_BOUND ((int64_t)((UINT64_C(1) << 30) - 1))
/* The proof parameters and witness norm bound are instantiated for one
 * replacement-free window containing at most this many deletions.  A fresh
 * window must be created before accepting another deletion. */
#define ACC96_MAX_DELETIONS_PER_WINDOW 50000U

typedef struct {
  acc96_poly_t right[ACC96_DELTA4_LAYERS][ACC96_GADGET_DIGITS];
} acc96_delta4_t;

/* Six independently transportable leading EvalFx layers.  A receiver whose
 * first differing handle bit is j consumes exactly layers 0..j. */
typedef struct {
  acc96_poly_t right[ACC96_DELTA6_LAYERS][ACC96_GADGET_DIGITS];
} acc96_delta6_t;

typedef struct {
  acc96_poly_t atrap;
  acc96_poly_t trunc_b[ACC96_TRUNC_RIGHT];
  acc96_ba_t ba;
  acc96_poly_t ua;
  uint8_t refresh_seed[ACC96_SEED_BYTES];
} acc96_relation_pp_t;

typedef struct {
  acc96_poly_t left;
  acc96_poly_t v12;
  acc96_poly_t v2[ACC96_TRUNC_RIGHT];
  acc96_poly_t tag;
  acc96_evalfx_delta_t right;
  uint8_t refresh_r[ACC96_REFRESH_BYTES];
  acc96_poly_t refresh_z[ACC96_REFRESH_BLOCKS];
} acc96_relation_witness_t;

typedef struct {
  uint32_t sequence;
  uint32_t revoked_handle;
} acc96_relation_update_t;

/* Deterministic SHAKE256 expansion of the initially uniform public matrices.
 * TrapGen subsequently replaces trunc_b by A*R_1+R_0, so a verifier that did
 * not participate in Setup must additionally receive that public block. */
int acc96_relation_pp_expand(acc96_relation_pp_t *pp,
                             const uint8_t seed[ACC96_SEED_BYTES],
                             const acc96_ctx_t *ctx);

/* Reconstruct the final public parameters from the 32-byte uniform seed and
 * the trapdoor-dependent public trunc_b block.  This avoids transmitting Ba,
 * A and ua while never deriving the secret trapdoor from public material. */
int acc96_relation_pp_reconstruct(
    acc96_relation_pp_t *pp,
    const uint8_t seed[ACC96_SEED_BYTES],
    const acc96_poly_t trunc_b[ACC96_TRUNC_RIGHT],
    const acc96_ctx_t *ctx);

void acc96_relation_uniform_poly(
    acc96_poly_t *out,
    const uint8_t seed[ACC96_SEED_BYTES],
    uint32_t domain,
    uint32_t first,
    uint32_t second);

/* Remark 4.3 / Construction E.1 branch A_r=A_0+sum_i r_i A_i.  The q88
 * rank-one specialization has one row of ACC96_REFRESH_BLOCKS ring elements
 * in each seeded A_i and one short ring element per block in z. */
void acc96_relation_refresh_matrix_poly(
    acc96_poly_t *out, const acc96_relation_pp_t *pp,
    size_t matrix_index, size_t block);
void acc96_relation_refresh_apply(
    acc96_poly_t *out, const acc96_relation_pp_t *pp,
    const acc96_relation_witness_t *witness, const acc96_ctx_t *ctx);
void acc96_relation_refresh_sample(
    acc96_relation_witness_t *witness,
    const uint8_t seed[ACC96_SEED_BYTES]);

/* Check newly sampled witnesses only: head coefficients in [-F,F], right/z
 * ternary, tag=1, canonical CRT coefficients. Does not check membership or
 * certify sampling probabilities. Never apply this to an updated witness. */
int acc96_relation_fresh_witness_check(
    const acc96_relation_witness_t *witness, const acc96_ctx_t *ctx);

/* Evaluate the exact MemVer left-hand side in R_q. */
void acc96_relation_evaluate(acc96_poly_t *out,
                             const acc96_relation_pp_t *pp,
                             uint32_t handle,
                             const acc96_relation_witness_t *witness,
                             const acc96_ctx_t *ctx);

/* Algebraic verification plus an optional Euclidean-norm-squared bound.
 * max_norm_sq <= 0 disables the norm gate and is intended only for relation
 * KATs.  Production callers must pass the theorem-derived public bound. */
int acc96_relation_memver(const acc96_relation_pp_t *pp,
                          const acc96_poly_t *accumulator,
                          uint32_t handle,
                          const acc96_relation_witness_t *witness,
                          long double max_norm_sq,
                          const acc96_ctx_t *ctx);

long double acc96_relation_witness_norm2(
    const acc96_relation_witness_t *witness,
    const acc96_ctx_t *ctx);

/* TA deletion and UAV witness update.  The two functions implement the
 * EvalF/EvalFx indicator identity over all three CRT lanes. */
void acc96_relation_delete(acc96_poly_t *accumulator,
                           acc96_relation_update_t *update,
                           const acc96_relation_pp_t *pp,
                           uint32_t revoked_handle,
                           uint32_t next_sequence,
                           const acc96_ctx_t *ctx);
int acc96_relation_delete_delta4(
    acc96_poly_t *accumulator, acc96_relation_update_t *update,
    acc96_delta4_t *delta, const acc96_relation_pp_t *pp,
    uint32_t revoked_handle, uint32_t next_sequence,
    const acc96_ctx_t *ctx);
int acc96_relation_delete_delta6(
    acc96_poly_t *accumulator, acc96_relation_update_t *update,
    acc96_delta6_t *delta, const acc96_relation_pp_t *pp,
    uint32_t revoked_handle, uint32_t next_sequence,
    const acc96_ctx_t *ctx);
int acc96_relation_delete_tail4(
    acc96_poly_t *accumulator, acc96_relation_update_t *update,
    acc96_gmat_t *tail4, const acc96_relation_pp_t *pp,
    uint32_t revoked_handle, uint32_t next_sequence,
    const acc96_ctx_t *ctx);

int acc96_relation_memwitup(acc96_relation_witness_t *witness,
                            const acc96_relation_pp_t *pp,
                            uint32_t member_handle,
                            const acc96_relation_update_t *update,
                            const acc96_ctx_t *ctx);

int acc96_relation_delta4(acc96_delta4_t *delta,
                          const acc96_relation_pp_t *pp,
                          uint32_t revoked_handle,
                          const acc96_ctx_t *ctx);
int acc96_relation_delta6(acc96_delta6_t *delta,
                          const acc96_relation_pp_t *pp,
                          uint32_t revoked_handle,
                          const acc96_ctx_t *ctx);

/* Returns zero without changing witness when the first differing low-order
 * handle bit is outside DELTA4. */
int acc96_relation_memwitup_delta4(
    acc96_relation_witness_t *witness, uint32_t member_handle,
    const acc96_relation_update_t *update, const acc96_delta4_t *delta,
    size_t *applied_layers);
int acc96_relation_memwitup_delta6(
    acc96_relation_witness_t *witness, uint32_t member_handle,
    const acc96_relation_update_t *update, const acc96_delta6_t *delta,
    size_t available_layers, size_t *applied_layers);
int acc96_relation_memwitup_tail4(
    acc96_relation_witness_t *witness, const acc96_relation_pp_t *pp,
    uint32_t member_handle, const acc96_relation_update_t *update,
    const acc96_gmat_t *tail4, const acc96_ctx_t *ctx);

/* Test-only constructor for known-answer and preservation tests.  It samples
 * small relation variables and defines accumulator := MemVer(witness).  It is
 * not SampleLeft and must never be called by Issue. */
void acc96_relation_make_kat_instance(
    acc96_poly_t *accumulator,
    acc96_relation_witness_t *witness,
    const acc96_relation_pp_t *pp,
    uint32_t handle,
    const uint8_t seed[ACC96_SEED_BYTES],
    const acc96_ctx_t *ctx);

#endif
