/*************************************************
 * TSampler covariance computation
 * Ported from tag-friendly_sampler covariance.c
 *************************************************/

#include "covariance.h"
#include "arith.h"

#include <flint/arb.h>
#include <flint/acb_dft.h>

#define SPECTRAL_NORM_PRECISION 128

static acb_dft_pre_t pre;

void fft_precomp_setup(void) {
  acb_dft_precomp_init(pre, 2*PARAM_N, SPECTRAL_NORM_PRECISION);
}

void fft_precomp_teardown(void) {
  acb_dft_precomp_clear(pre);
}

void compute_covariance(poly_real_mat_2d_2d S, const poly_q_mat_d_d RRstar[2][2]) {
  size_t i,j,k;
  poly_real finv, tmppoly;

  poly_real_init(finv);
  poly_real_init(tmppoly);

  /* convert and scale: -s_H²*s_4²/(s_4²-s_H²) * RRstar */
  for (i = 0; i < 2*PARAM_D; i++) {
    for (j = 0; j < 2*PARAM_D; j++) {
      poly_real_from_poly_q(S->rows[i]->entries[j], RRstar[i/PARAM_D][j/PARAM_D]->rows[i%PARAM_D]->entries[j%PARAM_D]);
      poly_real_mul_scalar(S->rows[i]->entries[j], S->rows[i]->entries[j], PARAM_NEG_SHSQS4SQ_DIV_S4SQ_SHSQ);
    }
  }

  /* Add diag((s1²-sL²).I_d, s3².I_d) */
  for (i = 0; i < PARAM_D; i++)
    poly_real_add_constant(S->rows[i]->entries[i], PARAM_S1SQ_SLSQ);
  for (i = PARAM_D; i < 2*PARAM_D; i++)
    poly_real_add_constant(S->rows[i]->entries[i], PARAM_S3SQ);

  /* Schur complement elimination */
  for (i = 2*PARAM_D-1; i > 0; i--) {
    poly_real_invert(finv, S->rows[i]->entries[i]);
    for (j = 0; j < i; j++)
      poly_real_mul(S->rows[j]->entries[i], S->rows[j]->entries[i], finv);
    for (j = 0; j < i; j++) {
      for (k = 0; k < i; k++) {
        poly_real_mul(tmppoly, S->rows[j]->entries[i], S->rows[i]->entries[k]);
        poly_real_sub(S->rows[j]->entries[k], S->rows[j]->entries[k], tmppoly);
      }
    }
  }

  poly_real_clear(finv);
  poly_real_clear(tmppoly);
}

int sk_sq_spectral_norm_columns(poly_q_mat_d_d RRstar[2][2],
                                size_t columns,
                                const poly_q_mat_d_d R[2][columns],
                                double max_sq_norm) {
  size_t i,j,k;
  acb_struct num_vec[PARAM_N*2], den_vec[PARAM_N*2];
  poly_q_mat_d_d Rstar[columns][2], tmp_q;
  poly_z_mat_d_d RRstar_z, power_mat, tmp_mat;
  fmpz_t c;
  arb_t num_embedding, den_embedding;
  double sq_norm_R, test_norm;
  int is_invalid;

  for (i = 0; i < 2*PARAM_N; i++) {
    acb_init(&num_vec[i]); acb_init(&den_vec[i]);
  }
  arb_init(num_embedding); arb_init(den_embedding);
  fmpz_init(c);

  poly_z_mat_d_d_init(RRstar_z); poly_z_mat_d_d_init(power_mat);
  poly_z_mat_d_d_init(tmp_mat); poly_q_mat_d_d_init(tmp_q);

  for (i = 0; i < 2; i++)
    for (j = 0; j < columns; j++) {
      poly_q_mat_d_d_init(Rstar[j][i]);
      poly_q_mat_d_d_conjugate(Rstar[j][i], R[i][j]);
    }

  /* compute R*R^* = RRstar */
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      poly_q_mat_d_d_mul_mat_d_d(RRstar[i][j], R[i][0], Rstar[0][j]);
      for (k = 1; k < columns; k++) {
        poly_q_mat_d_d_mul_mat_d_d(tmp_q, R[i][k], Rstar[k][j]);
        poly_q_mat_d_d_add(RRstar[i][j], RRstar[i][j], tmp_q);
      }
    }
  }

  is_invalid = 0;
  for (i = 0; i < 2; i++) {
    poly_z_mat_d_d_from_poly_q_mat_d_d(RRstar_z, RRstar[i][i]);
    poly_z_mat_d_d_set(power_mat, RRstar_z);
    /* power_mat = (R_i.R_i*)^32 */
    for (k = 0; k < PARAM_IT_SPEC_NORM; k++) {
      poly_z_mat_d_d_mul_mat_d_d(tmp_mat, power_mat, power_mat);
      poly_z_mat_d_d_set(power_mat, tmp_mat);
    }
    /* power_mat = (R_i.R_i*)^33, tmp_mat = (R_i.R_i*)^32 */
    poly_z_mat_d_d_mul_mat_d_d(power_mat, tmp_mat, RRstar_z);

    for (j = 0; j < PARAM_N; j++) {
      poly_z_get_coeff(c, power_mat->rows[0]->entries[0], j);
      acb_set_fmpz(&num_vec[j], c);
      poly_z_get_coeff(c, tmp_mat->rows[0]->entries[0], j);
      acb_set_fmpz(&den_vec[j], c);
    }
    acb_dft_precomp(num_vec, num_vec, pre, SPECTRAL_NORM_PRECISION);
    acb_dft_precomp(den_vec, den_vec, pre, SPECTRAL_NORM_PRECISION);

    sq_norm_R = 0.0;
    for (j = 0; j < PARAM_N; j++) {
      acb_get_real(num_embedding, &num_vec[2*j+1]);
      acb_get_real(den_embedding, &den_vec[2*j+1]);
      arb_div(num_embedding, num_embedding, den_embedding, SPECTRAL_NORM_PRECISION);
      test_norm = arf_get_d(arb_midref(num_embedding), ARF_RND_DOWN);
      if (test_norm > sq_norm_R) sq_norm_R = test_norm;
    }
    if (sq_norm_R > max_sq_norm) {
      is_invalid = 1;
      goto spectral_cleanup;
    }
  }

spectral_cleanup:
  for (i = 0; i < 2*PARAM_N; i++) { acb_clear(&num_vec[i]); acb_clear(&den_vec[i]); }
  arb_clear(num_embedding); arb_clear(den_embedding); fmpz_clear(c);
  poly_z_mat_d_d_clear(RRstar_z); poly_z_mat_d_d_clear(power_mat);
  poly_z_mat_d_d_clear(tmp_mat); poly_q_mat_d_d_clear(tmp_q);
  for (i = 0; i < 2; i++)
    for (j = 0; j < columns; j++)
      poly_q_mat_d_d_clear(Rstar[j][i]);
  return is_invalid;
}

int sk_sq_spectral_norm(poly_q_mat_d_d RRstar[2][2],
                        const poly_q_mat_d_d R[2][PARAM_KH]) {
  return sk_sq_spectral_norm_columns(RRstar, PARAM_KH, R,
                                     PARAM_R_MAX_SQ_SPECTRAL_NORM);
}
