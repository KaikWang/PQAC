#ifndef ACC_GENPP_H
#define ACC_GENPP_H

#include <stdint.h>
#include <stddef.h>

#include "params.h"
#include "revac_profile_generated.h"
#include "arith.h"
#include "sep.h"

/*
 * Prototype implementation of the AccGenpp.pdf accumulator algorithms.
 *
 * The PDF uses n = 1 and its own matrix dimensions.  This implementation
 * follows the existing code base and works over R_q^PARAM_D.  The accumulator
 * uses a balanced ternary gadget independently of the binary handle encoding
 * and of PARAM_Q_BITLEN, which remains a global modulus/serialization
 * parameter.  Handles are fixed 32-bit strings as in the PDF.
 */
#define ACC_GENPP_HANDLE_BITS REVAC_PARAM_ACC_HANDLE_BITS
#define ACC_GENPP_GADGET_BASE REVAC_PARAM_ACC_GADGET_BASE
#define ACC_GENPP_GADGET_DIM REVAC_PARAM_ACC_GADGET_DIM
#define ACC_GENPP_GADGET_SIGNED REVAC_PARAM_ACC_GADGET_SIGNED
#define ACC_GENPP_DELTA4_LAYERS REVAC_PARAM_PUBLIC_DELTA_LAYERS
#define ACC_GENPP_GADGET_MAX_CENTERED REVAC_PARAM_ACC_GADGET_MAX_CENTERED
#define ACC_GENPP_MAX_DICT_ENTRIES REVAC_PARAM_ACC_MAX_DICT_ENTRIES
#define ACC_GENPP_MAX_DELETIONS REVAC_PARAM_ACC_MAX_DELETIONS_STORAGE
#define ACC_GENPP_BML_CAPACITY REVAC_PARAM_BML_CAPACITY
#define ACC_GENPP_SEED_BYTES REVAC_PARAM_SEED_BYTES
#define ACC_GENPP_TRUNC_K REVAC_PARAM_ACC_TRUNC_K
#define ACC_GENPP_TRUNC_L REVAC_PARAM_ACC_TRUNC_L
#define ACC_GENPP_TRUNC_K_L (ACC_GENPP_TRUNC_K - ACC_GENPP_TRUNC_L)
#define ACC_GENPP_TRUNC_B REVAC_PARAM_ACC_TRUNC_BASE
#define ACC_GENPP_TRUNC_B_POW_L REVAC_PARAM_ACC_TRUNC_BASE_POW_L
#define ACC_GENPP_WIT_BETA_SQ REVAC_PARAM_ACC_WIT_BETA_SQ
#define ACC_GENPP_RANDOMIZED_WITNESS_REFRESH \
  REVAC_PARAM_ACC_RANDOMIZED_WITNESS_REFRESH
#define ACC_GENPP_REFRESH_BITS_MIN 128
#define ACC_GENPP_REFRESH_BITS_MAX 192
#define ACC_GENPP_REFRESH_BYTES_MAX ((ACC_GENPP_REFRESH_BITS_MAX + 7) / 8)

typedef enum {
  ACC_GENPP_OP_EMPTY = 0,
  ACC_GENPP_OP_ADD = 1,
  ACC_GENPP_OP_DEL = 2
} acc_genpp_op_t;

typedef struct {
  poly_q_mat_d_d blocks[ACC_GENPP_GADGET_DIM];
} acc_genpp_gmat_t;

typedef struct {
  acc_genpp_gmat_t bit[ACC_GENPP_HANDLE_BITS];
} acc_genpp_Ba_t;

typedef struct {
  uint32_t handle;
  acc_genpp_op_t op;
  uint32_t index;
  int8_t update_sign;
  int used;
} acc_genpp_dict_entry_t;

typedef struct {
  uint32_t seq;
  acc_genpp_op_t op;
  uint32_t handle;
  /* +1/-1 for authenticated randomized deletion.  Zero is accepted only by
     legacy in-process callers and is interpreted as +1; new wire records
     always carry an explicit nonzero sign. */
  int8_t sign;
} acc_genpp_update_t;

int acc_genpp_update_sign(const acc_genpp_update_t *update);
int acc_genpp_update_sign_is_explicit(const acc_genpp_update_t *update);

typedef struct {
  uint8_t ka[ACC_GENPP_SEED_BYTES];
  uint32_t ctr;
  poly_q_vec_d acc0;
  acc_genpp_dict_entry_t dict[ACC_GENPP_MAX_DICT_ENTRIES];
  acc_genpp_update_t bml[ACC_GENPP_BML_CAPACITY];
  size_t bml_start;
  size_t bml_len;
} acc_genpp_state_t;

typedef struct {
  sep_pk_t A;
  poly_q_mat_d_d trunc_B[ACC_GENPP_TRUNC_K_L];
  acc_genpp_Ba_t Ba;
  poly_q_vec_d ua;
  uint8_t refresh_seed[ACC_GENPP_SEED_BYTES];
  uint16_t refresh_bits;
} acc_genpp_pp_t;

typedef struct {
  sep_sk_t RA;
  poly_q_mat_d_d trunc_R[2][ACC_GENPP_TRUNC_K_L];
  poly_real_mat_2d_2d trunc_S;
} acc_genpp_td_t;

typedef struct {
  poly_q_vec_d left;
  poly_q_vec_d v12;
  poly_q_vec_d v2[ACC_GENPP_TRUNC_K_L];
  poly_q tag;
  poly_q_vec_d right[ACC_GENPP_HANDLE_BITS][ACC_GENPP_GADGET_DIM];
  uint8_t refresh_r[ACC_GENPP_REFRESH_BYTES_MAX];
  poly_q_vec_d refresh_z[ACC_GENPP_GADGET_DIM];
} acc_genpp_wit_t;

/* Optional update hook used by trusted authorities that can supply cached
   EvalF tail checkpoints.  Implementations must apply exactly one update to
   sx and return nonzero on success. */
typedef int (*acc_genpp_update_apply_fn)(
    void *opaque,
    acc_genpp_wit_t *sx,
    const acc_genpp_pp_t *pp,
    uint32_t member_handle,
    const acc_genpp_update_t *update);

/* Exact integer lifting backend for the accumulator MemVer relation. */
typedef enum {
  ACC_GENPP_LIFT_AUTO = 0,
  ACC_GENPP_LIFT_FORCE_INT64 = 1,
  ACC_GENPP_LIFT_FORCE_FLINT = 2,
  ACC_GENPP_LIFT_CROSSCHECK = 3
} acc_genpp_lift_mode_t;

typedef enum {
  ACC_GENPP_LIFT_BACKEND_NONE = 0,
  ACC_GENPP_LIFT_BACKEND_INT64 = 1,
  ACC_GENPP_LIFT_BACKEND_FLINT = 2
} acc_genpp_lift_backend_t;

typedef struct {
  acc_genpp_lift_backend_t backend;
  uint64_t public_abs_bound;
  int int64_bound_ok;
  int fallback_used;
  int crosscheck_performed;
  int crosscheck_passed;
} acc_genpp_lift_info_t;

void acc_genpp_gmat_init(acc_genpp_gmat_t *m);
void acc_genpp_gmat_clear(acc_genpp_gmat_t *m);
void acc_genpp_Ba_init(acc_genpp_Ba_t *B);
void acc_genpp_Ba_clear(acc_genpp_Ba_t *B);
void acc_genpp_state_init(acc_genpp_state_t *st);
void acc_genpp_state_clear(acc_genpp_state_t *st);
void acc_genpp_pp_init(acc_genpp_pp_t *pp);
void acc_genpp_pp_clear(acc_genpp_pp_t *pp);
void acc_genpp_td_init(acc_genpp_td_t *td);
void acc_genpp_td_clear(acc_genpp_td_t *td);
void acc_genpp_wit_init(acc_genpp_wit_t *wit);
void acc_genpp_wit_clear(acc_genpp_wit_t *wit);
size_t acc_genpp_acc_size_bytes(void);
size_t acc_genpp_wit_compact_size_bytes(void);
size_t acc_genpp_wit_full_size_bytes(void);
uint64_t acc_genpp_wit_norm2(const acc_genpp_wit_t *wit);
uint64_t acc_genpp_gadget_power(size_t j);

/* Appendix-E-style randomized refresh matrix A_r=A_0+sum r_i A_i. */
void acc_genpp_refresh_matrix_expand(acc_genpp_gmat_t *out,
                                     const uint8_t master[ACC_GENPP_SEED_BYTES],
                                     size_t matrix_index);
void acc_genpp_refresh_matrix_compose(acc_genpp_gmat_t *out,
                                      const acc_genpp_pp_t *pp,
                                      const uint8_t r[ACC_GENPP_REFRESH_BYTES_MAX]);

/* Canonical centered balanced-base-3 decomposition used by G_3^{-1}.
   Returns nonzero exactly when value is covered by ACC_GENPP_GADGET_DIM
   signed digits; on success digits are in {-1,0,1}, least significant first. */
int acc_genpp_balanced_ternary_decompose(
    coeff_q digits[ACC_GENPP_GADGET_DIM], coeff_q value);

void acc_genpp_Fk(uint8_t r[ACC_GENPP_SEED_BYTES],
                  const uint8_t ka[ACC_GENPP_SEED_BYTES],
                  uint32_t x,
                  const poly_q_vec_d acc);

void acc_genpp_evalf(acc_genpp_gmat_t *out,
                     const acc_genpp_Ba_t *B,
                     uint32_t selector);

/* Build the EvalF suffix only through T_1 (30 matrix composites). */
void acc_genpp_evalf_tail1(acc_genpp_gmat_t *out,
                           const acc_genpp_Ba_t *B,
                           uint32_t selector);

/*
 * Compute EvalF through T_0 and copy any requested T_i into pre-initialized
 * capture[i] objects.  capture may be NULL; entries 1..31 may independently
 * be NULL.  T_32 is the public gadget identity and is never captured.
 */
void acc_genpp_evalf_capture(
    acc_genpp_gmat_t *out,
    const acc_genpp_Ba_t *B,
    uint32_t selector,
    acc_genpp_gmat_t *capture[ACC_GENPP_HANDLE_BITS]);

/* Number of worker threads used by the independent gadget-column products.
   REVAC_ACC_THREADS=1 forces the reference serial path; otherwise the value is
   capped at ACC_GENPP_GADGET_DIM and the OpenMP runtime's available workers. */
size_t acc_genpp_parallel_threads(void);

void acc_genpp_evalfx_apply(acc_genpp_wit_t *delta,
                            const acc_genpp_Ba_t *B,
                            uint32_t selector,
                            uint32_t point,
                            const poly_q_vec_d u);

int acc_genpp_evalfx_apply_tails(
    acc_genpp_wit_t *delta,
    const acc_genpp_Ba_t *B,
    uint32_t selector,
    uint32_t point,
    const poly_q_vec_d u,
    const uint8_t *tail_indices,
    const acc_genpp_gmat_t *tail_values,
    size_t tail_count,
    size_t *local_composites);

/* Convert TA-private T1..T4 checkpoints for selector into the exact public
   witness increments consumed by EvalFX.  The output is independent of the
   member handle; that handle only determines how many leading layers apply. */
void acc_genpp_delta4_from_tails(
    poly_q_vec_d delta[ACC_GENPP_DELTA4_LAYERS][ACC_GENPP_GADGET_DIM],
    const acc_genpp_gmat_t tails[ACC_GENPP_DELTA4_LAYERS],
    uint32_t selector,
    const poly_q_vec_d u);

/* Convert an arbitrary leading T_1,...,T_layers prefix into reusable exact
   witness increments.  This is the generic primitive behind public DELTA4
   and TA-private, longer MemWitSync delta caches. */
void acc_genpp_delta_from_tails(
    poly_q_vec_d delta[][ACC_GENPP_GADGET_DIM],
    size_t layers,
    const acc_genpp_gmat_t tails[],
    uint32_t selector,
    const poly_q_vec_d u);

/* Apply an authenticated DELTA4 broadcast.  Returns zero without mutating sx
   when the member and selector share all four low bits, so the caller can
   safely choose the complete LOCAL recurrence instead. */
int acc_genpp_memwitup_delta4(
    acc_genpp_wit_t *sx,
    uint32_t member_handle,
    const acc_genpp_update_t *update,
    const poly_q_vec_d delta[ACC_GENPP_DELTA4_LAYERS]
                            [ACC_GENPP_GADGET_DIM],
    size_t *applied_layers);

/* Apply a reusable leading delta prefix.  Returns zero without mutating sx
   when the first differing handle bit lies outside the cached prefix. */
int acc_genpp_memwitup_delta(
    acc_genpp_wit_t *sx,
    uint32_t member_handle,
    const acc_genpp_update_t *update,
    const poly_q_vec_d delta[][ACC_GENPP_GADGET_DIM],
    size_t layers,
    size_t *applied_layers);

void acc_genpp_sampleleft(acc_genpp_wit_t *e,
                          const acc_genpp_pp_t *pp,
                          const acc_genpp_td_t *td,
                          const poly_q_vec_d u,
                          uint32_t x,
                          const uint8_t ra[ACC_GENPP_SEED_BYTES]);

void acc_genpp_gen(acc_genpp_pp_t *pp,
                   acc_genpp_td_t *td,
                   acc_genpp_state_t *st,
                   poly_q_vec_d A0);

int acc_genpp_add(poly_q_vec_d A_out,
                  acc_genpp_wit_t *sx,
                  acc_genpp_update_t *upmsg,
                  acc_genpp_state_t *st,
                  const acc_genpp_pp_t *pp,
                  const acc_genpp_td_t *td,
                  const poly_q_vec_d A,
                  uint32_t x);

int acc_genpp_delete(poly_q_vec_d A_out,
                     acc_genpp_update_t *upmsg,
                     acc_genpp_state_t *st,
                     const acc_genpp_pp_t *pp,
                     const acc_genpp_td_t *td,
                     const poly_q_vec_d A,
                     uint32_t x,
                     const acc_genpp_wit_t *sx);

int acc_genpp_delete_capture(
    poly_q_vec_d A_out,
    acc_genpp_update_t *upmsg,
    acc_genpp_state_t *st,
    const acc_genpp_pp_t *pp,
    const acc_genpp_td_t *td,
    const poly_q_vec_d A,
    uint32_t x,
    const acc_genpp_wit_t *sx,
    acc_genpp_gmat_t *capture[ACC_GENPP_HANDLE_BITS]);

/* TA-internal deletion path.  The caller must already be the authoritative
   revocation service and must enforce its own registration/issuance/revoked
   policy.  It retains the accumulator dictionary check and all EvalF/tail
   computation, but intentionally does not reconstruct or verify a member
   witness merely to authorize the TA's own deletion. */
int acc_genpp_delete_trusted_capture(
    poly_q_vec_d A_out,
    acc_genpp_update_t *upmsg,
    acc_genpp_state_t *st,
    const acc_genpp_pp_t *pp,
    const poly_q_vec_d A,
    uint32_t x,
    acc_genpp_gmat_t *capture[ACC_GENPP_HANDLE_BITS]);

/* Authenticated-randomized deletion used by the TA after the target handle is
   fixed.  sign must be exactly +1 or -1. */
int acc_genpp_delete_trusted_capture_signed(
    poly_q_vec_d A_out,
    acc_genpp_update_t *upmsg,
    acc_genpp_state_t *st,
    const acc_genpp_pp_t *pp,
    const poly_q_vec_d A,
    uint32_t x,
    int sign,
    acc_genpp_gmat_t *capture[ACC_GENPP_HANDLE_BITS]);

void acc_genpp_memwitup(acc_genpp_wit_t *sx,
                        const acc_genpp_pp_t *pp,
                        uint32_t x,
                        const acc_genpp_update_t *upmsg);

int acc_genpp_memwitup_tails(
    acc_genpp_wit_t *sx,
    const acc_genpp_pp_t *pp,
    uint32_t x,
    const acc_genpp_update_t *upmsg,
    const uint8_t *tail_indices,
    const acc_genpp_gmat_t *tail_values,
    size_t tail_count,
    size_t *local_composites);

int acc_genpp_bml_get_missing(const acc_genpp_state_t *st,
                              uint32_t last_seq,
                              acc_genpp_update_t *out,
                              size_t out_cap,
                              size_t *out_len);

int acc_genpp_memver(const acc_genpp_pp_t *pp,
                     const poly_q_vec_d A,
                     uint32_t x,
                     const acc_genpp_wit_t *sx);

/*
 * Compute lift = (A - MemVerExpression) / PARAM_Q coefficient-wise over
 * centered integer representatives in Z[X]/(X^PARAM_N + 1).  AUTO uses the
 * fixed-size int64 backend only when a public, parameter-only absolute bound
 * proves every intermediate safe; otherwise it falls back to FLINT.  The
 * CROSSCHECK mode runs both implementations and rejects any mismatch.
 */
int acc_genpp_memver_exact_lift(
    int64_t lift[PARAM_D][PARAM_N],
    const acc_genpp_pp_t *pp,
    const poly_q_vec_d A,
    uint32_t x,
    const acc_genpp_wit_t *sx,
    acc_genpp_lift_mode_t mode,
    acc_genpp_lift_info_t *info);

int acc_genpp_memwitsync(acc_genpp_wit_t *sx,
                         const acc_genpp_pp_t *pp,
                         const acc_genpp_td_t *td,
                         const acc_genpp_state_t *st,
                         const poly_q_vec_d A,
                         uint32_t x);

int acc_genpp_memwitsync_apply(acc_genpp_wit_t *sx,
                               const acc_genpp_pp_t *pp,
                               const acc_genpp_td_t *td,
                               const acc_genpp_state_t *st,
                               const poly_q_vec_d A,
                               uint32_t x,
                               acc_genpp_update_apply_fn apply_update,
                               void *apply_opaque);

#endif
