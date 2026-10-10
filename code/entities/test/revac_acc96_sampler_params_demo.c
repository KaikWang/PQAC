#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "acc80_crt.h"
#include "acc80_eval.h"
#include "acc96_proof_layout.h"
#include "acc96_truncated_params_generated.h"

static acc80_u128 generated_modulus(void) {
  return ((acc80_u128)ACC_TRUNC_PARAM_Q_HI << 64) |
         (acc80_u128)ACC_TRUNC_PARAM_Q_LO;
}

static int check_kernel_basis(void) {
  const acc80_u128 q = generated_modulus();
  for (size_t column = 0; column < ACC_TRUNC_PARAM_K; column++) {
    __int128 sum = 0;
    __int128 power = 1;
    for (size_t row = 0; row < ACC_TRUNC_PARAM_K; row++) {
      sum += power * acc_trunc_gadget_basis[row][column];
      power *= ACC_TRUNC_PARAM_B;
    }
    __int128 reduced = sum % (__int128)q;
    if (reduced < 0) reduced += (__int128)q;
    if (reduced != 0) return 0;
  }
  return 1;
}

int main(void) {
  int ok = 1;
  ok &= ACC_TRUNC_PARAM_Q_BITS == ACC80_MODULUS_BITS;
  ok &= ACC_TRUNC_PARAM_N == ACC80_N;
  ok &= ACC_TRUNC_PARAM_D == 1;
  ok &= ACC_TRUNC_PARAM_B == 13;
  ok &= ACC_TRUNC_PARAM_K == 24;
  ok &= ACC_TRUNC_PARAM_L == 2;
  ok &= ACC_TRUNC_PARAM_STATISTICAL_BITS == 128;
  ok &= generated_modulus() == acc80_modulus();
  ok &= ACC80_HANDLE_BITS == 24;
  ok &= ACC80_GADGET_BASE == UINT64_C(26007);
  ok &= ACC80_GADGET_DIGITS == 6;
  ok &= ACC80_GADGET_BOUND == 13003;
  ok &= ACC80_COEFF_PACKED_BYTES == 11;
  ok &= ACC80_POLY_PACKED_BYTES == 44U * 1024U;
  ok &= ACC96_MAX_DELETIONS_PER_WINDOW == 50000U;
  ok &= ACC96_REFRESH_BITS == 128U;
  ok &= ACC96_PROOF_NON_TAG_POLYS == 174U;
  ok &= ACC96_PROOF_NON_TAG_POLYS * ACC80_POLY_PACKED_BYTES ==
        7656U * 1024U;
  ok &= ACC80_HANDLE_BITS * ACC80_GADGET_DIGITS *
            ACC80_POLY_PACKED_BYTES == 6336U * 1024U;
  ok &= ACC_TRUNC_S1SQ > ACC_TRUNC_S3SQ;
  ok &= ACC_TRUNC_S2SQ > 0.0;
  ok &= ACC_TRUNC_R_MAX_SQ_SPECTRAL_NORM > 0.0;
  for (size_t i = 0; i < ACC_TRUNC_PARAM_K; i++) {
    ok &= isfinite(acc_trunc_klein_widths[i]);
    ok &= acc_trunc_klein_widths[i] > 0.0;
  }
  ok &= check_kernel_basis();

  if (!ok) {
    fprintf(stderr, "q88 truncated-sampler parameter KAT failed.\n");
    return 1;
  }
  printf("q88 truncated-sampler parameters passed: q=%u bit, N=%u, "
         "base=%u, k=%u, l=%u, statistical=%u bit.\n",
         ACC_TRUNC_PARAM_Q_BITS, ACC_TRUNC_PARAM_N, ACC_TRUNC_PARAM_B,
         ACC_TRUNC_PARAM_K, ACC_TRUNC_PARAM_L,
         ACC_TRUNC_PARAM_STATISTICAL_BITS);
  return 0;
}
