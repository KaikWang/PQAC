#ifndef ACC80_EVAL_H
#define ACC80_EVAL_H

#include <stdint.h>

#include "acc80_crt.h"

#define ACC80_HANDLE_BITS 24U

typedef struct {
  acc80_poly_t block[ACC80_GADGET_DIGITS];
} acc80_gmat_t;

typedef struct {
  acc80_gmat_t bit[ACC80_HANDLE_BITS];
} acc80_ba_t;

typedef struct {
  acc80_poly_t right[ACC80_HANDLE_BITS][ACC80_GADGET_DIGITS];
} acc80_evalfx_delta_t;

/* Rank-one q_acc implementation of the PKC'26 indicator-function chain. */
void acc80_gadget_identity(acc80_gmat_t *out);
int acc80_ginv_reconstructs(const acc80_gmat_t *value,
                            const acc80_ctx_t *ctx);
void acc80_evalf(acc80_gmat_t *out, const acc80_ba_t *ba,
                 uint32_t selector, const acc80_ctx_t *ctx);
int acc80_evalf_with_tail(acc80_gmat_t *out, acc80_gmat_t *tail_checkpoint,
                          size_t tail_bit, const acc80_ba_t *ba,
                          uint32_t selector, const acc80_ctx_t *ctx);
void acc80_evalf_apply(acc80_poly_t *out, const acc80_gmat_t *matrix,
                       const acc80_poly_t *value,
                       const acc80_ctx_t *ctx);
void acc80_evalfx(acc80_evalfx_delta_t *out, const acc80_ba_t *ba,
                  uint32_t selector, uint32_t member,
                  const acc80_poly_t *value, const acc80_ctx_t *ctx);
/* Compute the member-independent leading EvalFx increments.  A receiver only
 * consumes entries through the first low-order bit at which its handle
 * differs from selector. */
int acc80_evalfx_prefix(
    acc80_poly_t out[][ACC80_GADGET_DIGITS], size_t layers,
    const acc80_ba_t *ba, uint32_t selector, const acc80_poly_t *value,
    const acc80_ctx_t *ctx);
int acc80_evalf_evalfx_prefix(
    acc80_gmat_t *indicator,
    acc80_poly_t out[][ACC80_GADGET_DIGITS], size_t layers,
    const acc80_ba_t *ba, uint32_t selector, const acc80_poly_t *value,
    const acc80_ctx_t *ctx);
/* Export one exact high-to-low EvalF tail T_tail_bit.  A receiver whose first
 * differing bit is below tail_bit can resume the recurrence from this single
 * checkpoint instead of rebuilding all 31 upper levels. */
int acc80_evalf_tail_checkpoint(acc80_gmat_t *out,
                                const acc80_ba_t *ba,
                                uint32_t selector, size_t tail_bit,
                                const acc80_ctx_t *ctx);
int acc80_evalfx_from_tail(acc80_evalfx_delta_t *out,
                           const acc80_ba_t *ba,
                           uint32_t selector, uint32_t member,
                           const acc80_poly_t *value, size_t tail_bit,
                           const acc80_gmat_t *tail_checkpoint,
                           const acc80_ctx_t *ctx);
void acc80_member_matrix_apply(acc80_poly_t *out, const acc80_ba_t *ba,
                               uint32_t member,
                               const acc80_evalfx_delta_t *delta,
                               const acc80_ctx_t *ctx);

#endif
