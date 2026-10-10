#ifndef ACC96_LABRADOR_WIRE_H
#define ACC96_LABRADOR_WIRE_H

#include <stddef.h>
#include <stdint.h>

#include "acc96_exact_lift.h"

#define ACC96_LABRADOR_PUBLIC_MAGIC "RV88PP03"
#define ACC96_LABRADOR_PROVER_MAGIC "RV88PV04"
#define ACC96_LABRADOR_MAGIC_BYTES 8U

/* Static verifier/prover artifact.  A, Ba and ua are expanded from seed;
 * trunc_b is trapdoor-dependent public data and therefore remains explicit. */
#define ACC96_LABRADOR_PUBLIC_BYTES \
  (ACC96_LABRADOR_MAGIC_BYTES + ACC96_SEED_BYTES + \
   ACC96_TRUNC_RIGHT * ACC96_POLY_PACKED_BYTES)

/* Local-only prover payload.  The tag is omitted because the relation fixes
 * it to the constant polynomial one.  All signed witness/lift coefficients
 * are emitted as little-endian int16 limbs ready for native LaBRADOR slots.
 * The encoder checks the decomposition bound before narrowing.
 * Each limb is ordered phase-major under X^i -> X^(16k+r), not as 16
 * contiguous source chunks. */
#define ACC96_LABRADOR_WITNESS_LIMB_BYTES \
  (ACC96_PROOF_NON_TAG_POLYS * ACC96_WINDOW_PROOF_LIMBS * ACC96_N * 2U)
#define ACC96_LABRADOR_LIFT_LIMB_BYTES \
  (ACC96_LIMBS * ACC96_LIFT_PROOF_LIMBS * ACC96_N * 2U)
#define ACC96_LABRADOR_PROVER_BYTES \
  (ACC96_LABRADOR_MAGIC_BYTES + 2U * 4U + ACC96_POLY_PACKED_BYTES + \
   ACC96_HANDLE_BITS + ACC96_REFRESH_BYTES + \
   ACC96_LABRADOR_WITNESS_LIMB_BYTES + \
   ACC96_LABRADOR_LIFT_LIMB_BYTES)


int acc96_labrador_public_encode(
    uint8_t *out, size_t out_len,
    const uint8_t public_seed[ACC96_SEED_BYTES],
    const acc96_relation_pp_t *pp);

int acc96_labrador_public_decode(
    acc96_relation_pp_t *pp,
    uint8_t public_seed[ACC96_SEED_BYTES],
    const uint8_t *encoded, size_t encoded_len,
    const acc96_ctx_t *ctx);

int acc96_labrador_prover_encode(
    uint8_t *out, size_t out_len, uint32_t sequence, uint32_t handle,
    const acc96_poly_t *accumulator,
    const acc96_relation_witness_t *witness,
    const acc96_exact_lift_t *lift,
    const acc96_ctx_t *ctx);

#endif
