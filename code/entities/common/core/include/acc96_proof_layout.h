#ifndef ACC96_PROOF_LAYOUT_H
#define ACC96_PROOF_LAYOUT_H

#include <stdint.h>

#include "acc96_relation.h"

/* LaBRADOR stores each witness-vector squared-norm bound in uint64_t.  Four
 * balanced 28-bit limbs cover every centered q88 coefficient while keeping a
 * degree-256 limb chunk below 2^62 in the worst case. */
#define ACC96_PROOF_LIMB_BITS 28U
#define ACC96_PROOF_LIMBS 4U
#define ACC96_PROOF_LIMB_RADIX (UINT32_C(1) << ACC96_PROOF_LIMB_BITS)
#define ACC96_PROOF_LIMB_BOUND (INT32_C(1) << (ACC96_PROOF_LIMB_BITS - 1))
#define ACC96_LABRADOR_BASE_DEGREE 256U
#define ACC96_LABRADOR_PHASES \
  (ACC96_N / ACC96_LABRADOR_BASE_DEGREE)
#define ACC96_PROOF_LIMB_CHUNK_NORM2_BOUND \
  (UINT64_C(1) << 62)

/* The generic 4x28 layout above can encode every centered q88 residue, but it
 * is unnecessarily large for a windowed credential witness.  EvalFx emits
 * balanced gadget digits bounded by 13003, hence 50000 updates add at most
 * A coefficient of ginv_apply is a sum of at most 6*4096 products of two
 * such digits.  Consequently one update is bounded by
 * 6*4096*13003^2 < 2^42, and 50000 updates plus a 2^30-1 fresh-witness gate
 * fit a common five-limb signed representation shared by all three CRT
 * equations.  Although LaBRADOR can store int16_t, using the maximum radix
 * makes the conservative JL norm too large for secure Dachshund parameters.
 * Radix 8191 (digits in [-4095,4095]) reduces that squared-norm envelope by
 * roughly 64x while adding only one witness limb.  This is a proof encoding
 * gate, not a replacement for the theorem-derived L2 norm gate. */
#define ACC96_WINDOW_MAX_DELETIONS ACC96_MAX_DELETIONS_PER_WINDOW
#define ACC96_WINDOW_UPDATE_COEFF_BOUND \
  (INT64_C(6) * ACC96_N * INT64_C(13003) * INT64_C(13003))
#define ACC96_WINDOW_FRESH_COEFF_BOUND \
  ((int64_t)((UINT64_C(1) << 30) - 1))
#define ACC96_WINDOW_COEFF_BOUND \
  ((int64_t)ACC96_WINDOW_MAX_DELETIONS * \
       ACC96_WINDOW_UPDATE_COEFF_BOUND + \
   ACC96_WINDOW_FRESH_COEFF_BOUND)
#define ACC96_WINDOW_PROOF_LIMB_BITS 12U
#define ACC96_WINDOW_PROOF_LIMBS 5U
#define ACC96_WINDOW_PROOF_LIMB_RADIX INT64_C(4095)
#define ACC96_WINDOW_PROOF_LIMB_BOUND INT32_C(2047)
#define ACC96_WINDOW_PROOF_CHUNK_NORM2_BOUND (UINT64_C(1) << 33)

/* Conservative exact-lift gate.  With tag fixed to 1, every lane contains
 * 173 public-polynomial convolutions; 174*N*W is a simple strict envelope
 * after division by the ~30-bit lane modulus, where W is the common witness
 * coefficient gate.  Six balanced radix-8191 limbs cover this <2^75 bound
 * while keeping the LaBRADOR JL norm parameterizable. */
#define ACC96_LIFT_COEFF_BOUND \
  ((acc96_u128)174U * ACC96_N * (acc96_u128)ACC96_WINDOW_COEFF_BOUND)
#define ACC96_LIFT_PROOF_LIMB_BITS ACC96_WINDOW_PROOF_LIMB_BITS
#define ACC96_LIFT_PROOF_LIMBS 7U
#define ACC96_LIFT_PROOF_LIMB_RADIX ACC96_WINDOW_PROOF_LIMB_RADIX
#define ACC96_LIFT_PROOF_LIMB_BOUND ACC96_WINDOW_PROOF_LIMB_BOUND
#define ACC96_LIFT_PROOF_CHUNK_NORM2_BOUND (UINT64_C(1) << 33)

/* Stable witness-slot ABI shared by the C exporter and the LaBRADOR worker.
 * Do not reorder these slots without changing the q88 proof format version.
 * The tag occupies one slot but uses dedicated short/binary constraints
 * rather than the generic radix-65535 signed-limb path. */
#define ACC96_PROOF_SLOT_LEFT 0U
#define ACC96_PROOF_SLOT_V12 1U
#define ACC96_PROOF_SLOT_V2_FIRST 2U
#define ACC96_PROOF_SLOT_TAG \
  (ACC96_PROOF_SLOT_V2_FIRST + ACC96_TRUNC_RIGHT)
#define ACC96_PROOF_SLOT_RIGHT_FIRST (ACC96_PROOF_SLOT_TAG + 1U)
#define ACC96_PROOF_SLOT_REFRESH_Z_FIRST \
  (ACC96_PROOF_SLOT_RIGHT_FIRST + \
   ACC96_HANDLE_BITS * ACC96_GADGET_DIGITS)
#define ACC96_PROOF_WITNESS_POLYS \
  (ACC96_PROOF_SLOT_REFRESH_Z_FIRST + ACC96_REFRESH_BLOCKS)
#define ACC96_PROOF_NON_TAG_POLYS (ACC96_PROOF_WITNESS_POLYS - 1U)

size_t acc96_proof_right_slot(size_t bit, size_t digit);
const acc96_poly_t *acc96_proof_witness_poly_at(
    const acc96_relation_witness_t *witness, size_t slot);
acc96_poly_t *acc96_proof_witness_poly_at_mut(
    acc96_relation_witness_t *witness, size_t slot);
int acc96_proof_slot_is_tag(size_t slot);

int acc96_lift_proof_limb_decompose(
    int32_t limbs[ACC96_LIFT_PROOF_LIMBS], acc96_i128 coefficient);
acc96_i128 acc96_lift_proof_limb_reconstruct(
    const int32_t limbs[ACC96_LIFT_PROOF_LIMBS]);

int acc96_proof_limb_decompose(
    int32_t limbs[ACC96_PROOF_LIMBS],
    const acc96_coeff_t *coefficient,
    const acc96_ctx_t *ctx);

acc96_i128 acc96_proof_limb_reconstruct(
    const int32_t limbs[ACC96_PROOF_LIMBS]);

int acc96_proof_poly_decompose(
    int32_t limbs[ACC96_PROOF_LIMBS][ACC96_N],
    const acc96_poly_t *polynomial,
    const acc96_ctx_t *ctx);

/* Maximum squared norm among the 4*16 degree-256 backend groups. */
uint64_t acc96_proof_poly_max_chunk_norm2(
    const int32_t limbs[ACC96_PROOF_LIMBS][ACC96_N]);

/* Compact representation used by the 50000-deletion window relation.  It
 * deliberately rejects a coefficient outside the public window gate instead
 * of silently falling back to the generic four-limb encoding. */
int acc96_window_proof_limb_decompose(
    int32_t limbs[ACC96_WINDOW_PROOF_LIMBS],
    const acc96_coeff_t *coefficient,
    const acc96_ctx_t *ctx);

int64_t acc96_window_proof_limb_reconstruct(
    const int32_t limbs[ACC96_WINDOW_PROOF_LIMBS]);

int acc96_window_proof_poly_decompose(
    int32_t limbs[ACC96_WINDOW_PROOF_LIMBS][ACC96_N],
    const acc96_poly_t *polynomial,
    const acc96_ctx_t *ctx);

uint64_t acc96_window_proof_poly_max_chunk_norm2(
    const int32_t limbs[ACC96_WINDOW_PROOF_LIMBS][ACC96_N]);

/* Polyphase order used by the degree-256 LaBRADOR backend:
 *   f(X) = sum_{r=0}^{15} X^r f_r(X^16).
 * The coefficient k of phase r is therefore source coefficient 16*k+r.
 * This layout preserves the q88 negacyclic product using 16 phase equations
 * while keeping scalar handle bits in degree-0 backend slots. */
size_t acc96_proof_polyphase_source_index(size_t phase,
                                          size_t coefficient);

/* Verify that every coefficient of a real MemVer witness has the common
 * signed representative required by the radix-65535 window relation. max_abs is
 * optional and reports the largest non-negative absolute coefficient. */
int acc96_window_proof_witness_check(
    const acc96_relation_witness_t *witness,
    const acc96_ctx_t *ctx,
    uint64_t *max_abs);

#endif
