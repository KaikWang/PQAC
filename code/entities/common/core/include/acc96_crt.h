#ifndef REVAC_ACC96_CRT_H
#define REVAC_ACC96_CRT_H

/* Public name for the q96 accumulator arithmetic.  The implementation keeps
 * its original acc80_* symbols temporarily so older experiment code remains
 * source-compatible while the candidate is integrated. */
#include "acc80_crt.h"

typedef acc80_u128 acc96_u128;
typedef acc80_i128 acc96_i128;
typedef acc80_coeff_t acc96_coeff_t;
typedef acc80_poly_t acc96_poly_t;
typedef acc80_ctx_t acc96_ctx_t;

#define ACC96_N ACC80_N
#define ACC96_LIMBS ACC80_LIMBS
#define ACC96_P0 ACC80_P0
#define ACC96_P1 ACC80_P1
#define ACC96_P2 ACC80_P2
#define ACC96_MODULUS_BITS ACC80_MODULUS_BITS
#define ACC96_GADGET_BASE ACC80_GADGET_BASE
#define ACC96_GADGET_DIGITS ACC80_GADGET_DIGITS
#define ACC96_GADGET_BOUND ACC80_GADGET_BOUND
#define ACC96_COEFF_PACKED_BYTES ACC80_COEFF_PACKED_BYTES
#define ACC96_POLY_PACKED_BYTES ACC80_POLY_PACKED_BYTES

#define acc96_ctx_init acc80_ctx_init
#define acc96_modulus acc80_modulus
#define acc96_poly_zero acc80_poly_zero
#define acc96_poly_set acc80_poly_set
#define acc96_poly_add acc80_poly_add
#define acc96_poly_sub acc80_poly_sub
#define acc96_poly_neg acc80_poly_neg
#define acc96_poly_equal acc80_poly_equal
#define acc96_poly_mul acc80_poly_mul
#define acc96_poly_mul_workspace acc80_poly_mul_workspace
#define acc96_coeff_from_u128 acc80_coeff_from_u128
#define acc96_coeff_to_u128 acc80_coeff_to_u128
#define acc96_coeff_to_centered acc80_coeff_to_centered
#define acc96_balanced_decompose acc80_balanced_decompose
#define acc96_poly_pack acc80_poly_pack
#define acc96_poly_unpack acc80_poly_unpack

#endif
