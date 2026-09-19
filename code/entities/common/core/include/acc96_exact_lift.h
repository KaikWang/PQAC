#ifndef ACC96_EXACT_LIFT_H
#define ACC96_EXACT_LIFT_H

#include <stddef.h>

#include "acc96_proof_layout.h"

/* Exact integer quotients for the three CRT-lane MemVer equations:
 *
 *   accumulator - MemVer_integer(witness) = p_lane * lift_lane.
 *
 * The common signed q96 witness is used in every lane.  This is essential:
 * independently choosing a witness modulo each CRT prime would not prove
 * knowledge of one bounded q96 witness. */
typedef struct {
  acc96_i128 coeff[ACC96_LIMBS][ACC96_N];
} acc96_exact_lift_t;

int acc96_memver_exact_lift(
    acc96_exact_lift_t *lift,
    unsigned *maximum_abs_bits,
    const acc96_relation_pp_t *pp,
    const acc96_poly_t *accumulator,
    uint32_t handle,
    const acc96_relation_witness_t *witness,
    const acc96_ctx_t *ctx);

#endif
