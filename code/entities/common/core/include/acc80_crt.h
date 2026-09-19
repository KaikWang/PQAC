#ifndef ACC80_CRT_H
#define ACC80_CRT_H

#include <stddef.h>
#include <stdint.h>

/* Independent CRT accumulator arithmetic.  The historical acc80 symbol
 * prefix is retained for source compatibility; the candidate below uses an
 * 88-bit modulus.  Three <=31-bit lanes keep the accumulator relation within
 * the coefficient range supported by the existing LaBRADOR bridge.  Signature
 * arithmetic remains in R_405769. */
#define ACC80_N 4096U
#define ACC80_LOG_N 12U
#define ACC80_LIMBS 3U
#define ACC80_P0 UINT64_C(1073692673)
#define ACC80_P1 UINT64_C(536813569)
#define ACC80_P2 UINT64_C(536690689)
#define ACC80_GARNER_INV1 UINT64_C(536780804)
#define ACC80_GARNER_INV2 UINT64_C(317361870)
#define ACC80_MODULUS_BITS 88U
#define ACC80_GADGET_BASE UINT64_C(26007)
#define ACC80_GADGET_DIGITS 6U
#define ACC80_GADGET_BOUND INT32_C(13003)
#define ACC80_COEFF_PACKED_BYTES 11U
#define ACC80_POLY_PACKED_BYTES (ACC80_N * ACC80_COEFF_PACKED_BYTES)

typedef unsigned __int128 acc80_u128;
typedef __int128 acc80_i128;

typedef struct {
  uint64_t limb[ACC80_LIMBS];
} acc80_coeff_t;

typedef struct {
  acc80_coeff_t coeff[ACC80_N];
} acc80_poly_t;

typedef struct {
  uint64_t psi[ACC80_LIMBS];
  uint64_t psi_inv[ACC80_LIMBS];
  uint64_t omega[ACC80_LIMBS];
  uint64_t omega_inv[ACC80_LIMBS];
  uint64_t n_inv[ACC80_LIMBS];
  uint64_t mul_reciprocal[ACC80_LIMBS];
  uint64_t stage_root[ACC80_LIMBS][ACC80_LOG_N + 1U];
  uint64_t stage_root_inv[ACC80_LIMBS][ACC80_LOG_N + 1U];
  uint64_t garner_inv[ACC80_LIMBS];
  int initialized;
} acc80_ctx_t;

int acc80_ctx_init(acc80_ctx_t *ctx);
acc80_u128 acc80_modulus(void);

void acc80_poly_zero(acc80_poly_t *out);
void acc80_poly_set(acc80_poly_t *out, const acc80_poly_t *in);
void acc80_poly_add(acc80_poly_t *out, const acc80_poly_t *lhs,
                    const acc80_poly_t *rhs);
void acc80_poly_sub(acc80_poly_t *out, const acc80_poly_t *lhs,
                    const acc80_poly_t *rhs);
void acc80_poly_neg(acc80_poly_t *out, const acc80_poly_t *in);
int acc80_poly_equal(const acc80_poly_t *lhs, const acc80_poly_t *rhs);
void acc80_poly_mul(acc80_poly_t *out, const acc80_poly_t *lhs,
                    const acc80_poly_t *rhs, const acc80_ctx_t *ctx);
/* Workspace variant used by high-volume SampleLeft arithmetic.  The two
 * buffers must each contain ACC80_N uint64_t words and are private to the
 * caller/thread. */
void acc80_poly_mul_workspace(acc80_poly_t *out, const acc80_poly_t *lhs,
                              const acc80_poly_t *rhs,
                              const acc80_ctx_t *ctx, uint64_t *left,
                              uint64_t *right);

void acc80_coeff_from_u128(acc80_coeff_t *out, acc80_u128 value);
acc80_u128 acc80_coeff_to_u128(const acc80_coeff_t *value,
                               const acc80_ctx_t *ctx);
acc80_i128 acc80_coeff_to_centered(const acc80_coeff_t *value,
                                   const acc80_ctx_t *ctx);
/* Decompose the centered integer lift.  This entry point is retained for
 * regression tests and for measuring the effect of the congruence-class
 * optimization. */
int acc80_balanced_decompose_canonical(
    int32_t out[ACC80_GADGET_DIGITS], const acc80_coeff_t *value,
    const acc80_ctx_t *ctx);

/* Among every representable integer congruent to value modulo q, select a
 * balanced-base preimage with minimum squared Euclidean norm.  Ties are
 * resolved deterministically by l1 norm, Hamming weight, representative
 * magnitude, and finally lexicographic digit order. */
int acc80_balanced_decompose_min_norm(
    int32_t out[ACC80_GADGET_DIGITS], const acc80_coeff_t *value,
    const acc80_ctx_t *ctx);

/* Production G^-1 uses the minimum-norm congruence-class representative. */
int acc80_balanced_decompose(int32_t out[ACC80_GADGET_DIGITS],
                             const acc80_coeff_t *value,
                             const acc80_ctx_t *ctx);
size_t acc80_poly_pack(uint8_t out[ACC80_POLY_PACKED_BYTES],
                       const acc80_poly_t *poly);
int acc80_poly_unpack(acc80_poly_t *poly, const uint8_t *in, size_t in_len);

#endif
