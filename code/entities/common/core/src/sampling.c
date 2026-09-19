/*************************************************
 * TSampler: Tag-Friendly Lattice Sampler
 * Ported from tag-friendly_sampler/implementation/code/src/sampling.c
 *
 * Adapted: poly_real_samplefz + poly_q_from_poly_real
 *       -> poly_q_samplefz (direct poly_q output)
 *************************************************/

#include "sampling.h"
#include "fips202.h"
#include "arith.h"
#include "macros.h"
#include "random.h"

void perturbation_sampler(
    poly_q_vec_d p1_p2[PARAM_KL],
    poly_q_vec_d p3,
    poly_q_vec_d p4[PARAM_KH],
    const poly_q_mat_d_d R[2][PARAM_KH],
    const poly_real_mat_2d_2d S)
{
  size_t i,j;
  poly_real tmp_sub, tmp_mul, tmp_sample;
  poly_real_vec_2d c;

  poly_real_init(tmp_sub);
  poly_real_init(tmp_mul);
  poly_real_init(tmp_sample);
  poly_real_vec_2d_init(c);

  /* sample p_2 from Gaussian with width sqrt(s_2^2 - s_L^2) */
  for (i = 1; i < PARAM_KL; i++) {
    poly_q_vec_d_gaussian_sqrt_s2sq_sLsq(p1_p2[i]);
  }

  /* sample p_4 from Gaussian with width sqrt(s_4^2 - s_H^2) */
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_vec_d_gaussian_sqrt_s4sq_sHsq(p4[i]);
  }

  /* p1_p2[0] = R1.p4, p3 = R2.p4 */
  poly_q_mat_d_d_mul_vec_d(p1_p2[0], R[0][0], p4[0]);
  poly_q_mat_d_d_mul_vec_d(p3, R[1][0], p4[0]);
  for (i = 1; i < PARAM_KH; i++) {
    poly_q_mat_d_d_muladd_vec_d(p1_p2[0], R[0][i], p4[i]);
    poly_q_mat_d_d_muladd_vec_d(p3, R[1][i], p4[i]);
  }

  /* convert to real and scale */
  for (i = 0; i < PARAM_D; i++) {
    poly_real_from_poly_q(c->entries[i],           p1_p2[0]->entries[i]);
    poly_real_from_poly_q(c->entries[i + PARAM_D], p3->entries[i]);
    poly_real_mul_scalar(c->entries[i],           c->entries[i],           PARAM_NEG_SHSQ_DIV_S4SQ_SHSQ);
    poly_real_mul_scalar(c->entries[i + PARAM_D], c->entries[i + PARAM_D], PARAM_NEG_SHSQ_DIV_S4SQ_SHSQ);
  }

  /* sample p3 with Schur complement update */
  for (i = 2*PARAM_D; i > PARAM_D; i--) {
    poly_q_samplefz(p3->entries[i-1 - PARAM_D], S->rows[i-1]->entries[i-1], c->entries[i-1]);
    poly_real_from_poly_q(tmp_sample, p3->entries[i-1 - PARAM_D]);
    poly_real_sub(tmp_sub, tmp_sample, c->entries[i-1]);
    for (j = 0; j < i-1; j++) {
      poly_real_mul(tmp_mul, tmp_sub, S->rows[j]->entries[i-1]);
      poly_real_add(c->entries[j], c->entries[j], tmp_mul);
    }
  }

  /* sample p1_p2[0] with Schur complement update */
  for (i = PARAM_D; i > 0; i--) {
    poly_q_samplefz(p1_p2[0]->entries[i-1], S->rows[i-1]->entries[i-1], c->entries[i-1]);
    poly_real_from_poly_q(tmp_sample, p1_p2[0]->entries[i-1]);
    poly_real_sub(tmp_sub, tmp_sample, c->entries[i-1]);
    for (j = 0; j < i-1; j++) {
      poly_real_mul(tmp_mul, tmp_sub, S->rows[j]->entries[i-1]);
      poly_real_add(c->entries[j], c->entries[j], tmp_mul);
    }
  }

  poly_real_vec_2d_clear(c);
  poly_real_clear(tmp_mul);
  poly_real_clear(tmp_sub);
  poly_real_clear(tmp_sample);
}

void gadget_sampler(
    poly_q_vec_d zL[PARAM_KL],
    poly_q_vec_d zH[PARAM_KH],
    const poly_q_vec_d w,
    const poly_q tag,
    const poly_q taginv)
{
  size_t i,j;
  poly_q tmp_poly, tmp_center;
  poly_q_init(tmp_poly);
  poly_q_init(tmp_center);

  for (i = 0; i < PARAM_D; i++) {
    poly_q_set(tmp_poly, w->entries[i]);
    /* sampling zL */
    for (j = 0; j < PARAM_KL; j++) {
      poly_q_mod_bL(tmp_center, tmp_poly);
      poly_q_gaussian_coset_sL(zL[j]->entries[i], tmp_center);
      poly_q_sub(tmp_center, tmp_poly, zL[j]->entries[i]);
      poly_q_div_bL(tmp_poly, tmp_center);
    }
    /* sampling zH */
    for (j = 0; j < PARAM_KH; j++) {
      poly_q_mul(tmp_center, taginv, tmp_poly);
      poly_q_mod_bH(tmp_center, tmp_center);
      poly_q_gaussian_coset_sH(zH[j]->entries[i], tmp_center);
      poly_q_mul(tmp_center, tag, zH[j]->entries[i]);
      poly_q_sub(tmp_center, tmp_poly, tmp_center);
      poly_q_div_bH(tmp_poly, tmp_center);
    }
  }
  poly_q_clear(tmp_poly);
  poly_q_clear(tmp_center);
}

void sampler(
    poly_q_vec_d v11,
    poly_q_vec_d v12,
    poly_q_vec_d v2[PARAM_KH],
    const poly_q_mat_d_d R[2][PARAM_KH],
    const poly_q_mat_d_d A,
    const poly_q_mat_d_d B[PARAM_KH],
    const poly_q_vec_d u,
    const poly_q tag,
    const poly_q taginv,
    const poly_real_mat_2d_2d S)
{
  size_t i;
  poly_q_vec_d p1_p2[PARAM_KL], p3, p4[PARAM_KH];
  poly_q_vec_d w, tmp, zL[PARAM_KL], zH[PARAM_KH];
  int64_t bexpi;

  poly_q_vec_d_init(w);
  poly_q_vec_d_init(tmp);
  for (i = 0; i < PARAM_KL; i++) {
    poly_q_vec_d_init(p1_p2[i]);
    poly_q_vec_d_init(zL[i]);
  }
  poly_q_vec_d_init(p3);
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_vec_d_init(p4[i]);
    poly_q_vec_d_init(zH[i]);
  }

  /* perturbation sampler */
  perturbation_sampler(p1_p2, p3, p4, R, S);

  /* compute G_L.p1_p2 into p1_p2[0] */
  bexpi = PARAM_BL;
  for (i = 1; i < PARAM_KL; i++) {
    poly_q_vec_d_mul_scalar(tmp, p1_p2[i], bexpi);
    poly_q_vec_d_add(p1_p2[0], p1_p2[0], tmp);
    bexpi *= PARAM_BL;
  }

  /* compute w = u - G_L.p1_p2 - A.p3 + B.p4 - qL.tag.G_H.p4 */
  poly_q_vec_d_mul_scalar(w, p4[0], PARAM_QL);
  bexpi = PARAM_QL * PARAM_BH;
  for (i = 1; i < PARAM_KH; i++) {
    poly_q_vec_d_mul_scalar(tmp, p4[i], bexpi);
    poly_q_vec_d_add(w, w, tmp);
    bexpi *= PARAM_BH;
  }
  poly_q_vec_d_mul_poly(w, w, tag);
  poly_q_vec_d_sub(w, u, w);
  poly_q_vec_d_sub(w, w, p1_p2[0]);
  poly_q_mat_d_d_mulsub_vec_d(w, A, p3);
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_mat_d_d_muladd_vec_d(w, B[i], p4[i]);
  }

  /* Gadget sampling */
  gadget_sampler(zL, zH, w, tag, taginv);

  /* v11 = G_L.p1_p2 + G_L.zL + R1.zH */
  bexpi = PARAM_BL;
  for (i = 1; i < PARAM_KL; i++) {
    poly_q_vec_d_mul_scalar(tmp, zL[i], bexpi);
    poly_q_vec_d_add(zL[0], zL[0], tmp);
    bexpi *= PARAM_BL;
  }
  poly_q_mat_d_d_mul_vec_d(v11, R[0][0], zH[0]);
  for (i = 1; i < PARAM_KH; i++) {
    poly_q_mat_d_d_muladd_vec_d(v11, R[0][i], zH[i]);
  }
  poly_q_vec_d_add(v11, v11, p1_p2[0]);
  poly_q_vec_d_add(v11, v11, zL[0]);

  /* v12 = p3 + R2.zH */
  poly_q_mat_d_d_mul_vec_d(v12, R[1][0], zH[0]);
  for (i = 1; i < PARAM_KH; i++) {
    poly_q_mat_d_d_muladd_vec_d(v12, R[1][i], zH[i]);
  }
  poly_q_vec_d_add(v12, v12, p3);

  /* v2 = p4 + zH */
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_vec_d_add(v2[i], p4[i], zH[i]);
  }

  poly_q_vec_d_clear(w);
  poly_q_vec_d_clear(tmp);
  for (i = 0; i < PARAM_KL; i++) {
    poly_q_vec_d_clear(p1_p2[i]);
    poly_q_vec_d_clear(zL[i]);
  }
  poly_q_vec_d_clear(p3);
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_vec_d_clear(p4[i]);
    poly_q_vec_d_clear(zH[i]);
  }
}
