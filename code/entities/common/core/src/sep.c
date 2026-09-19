/*************************************************
 * SEP signature scheme with TSampler (tag-friendly)
 * Ported from tag-friendly_sampler sign.c + sampling.c
 *************************************************/

#include "arith.h"
#include "randombytes.h"
#include "sep.h"
#include "sampling.h"
#include "poly_q_sampling.h"
#include "fips202.h"
#include "covariance.h"

void sep_keys_init(sep_pk_t *pk, sep_sk_t *sk) {
  for (size_t i = 0; i < PARAM_KH; i++) {
    for (size_t j = 0; j < 2; j++)
      poly_q_mat_d_d_init(sk->R[j][i]);
    poly_q_mat_d_d_init(pk->B[i]);
  }
  poly_real_mat_2d_2d_init(sk->S);
}

void sep_keys_clear(sep_pk_t *pk, sep_sk_t *sk) {
  for (size_t i = 0; i < PARAM_KH; i++) {
    for (size_t j = 0; j < 2; j++)
      poly_q_mat_d_d_clear(sk->R[j][i]);
    poly_q_mat_d_d_clear(pk->B[i]);
  }
  poly_real_mat_2d_2d_clear(sk->S);
}

void sep_sig_init(sep_sig_t *sig) {
  poly_q_init(sig->tag);
  for (size_t i = 0; i < PARAM_KH; i++)
    poly_q_vec_d_init(sig->v2[i]);
  poly_q_vec_d_init(sig->v12);
  poly_q_init(sig->v3);
}

void sep_sig_clear(sep_sig_t *sig) {
  poly_q_clear(sig->tag);
  for (size_t i = 0; i < PARAM_KH; i++)
    poly_q_vec_d_clear(sig->v2[i]);
  poly_q_vec_d_clear(sig->v12);
  poly_q_clear(sig->v3);
}

void sep_keygen(sep_pk_t *pk, sep_sk_t *sk) {
  uint8_t root_seed[SEED_BYTES], seeds[SEED_BYTES*2];
  uint8_t *public_seed = seeds, *secret_seed = &seeds[SEED_BYTES];
  poly_q_mat_d_d A, RRstar[2][2];
  uint32_t kappa;

  randombytes(root_seed, SEED_BYTES);
  sha3_512(seeds, root_seed, SEED_BYTES);

  poly_q_mat_d_d_init(A);
  poly_q_mat_d_d_init(RRstar[0][0]); poly_q_mat_d_d_init(RRstar[0][1]);
  poly_q_mat_d_d_init(RRstar[1][0]); poly_q_mat_d_d_init(RRstar[1][1]);

  poly_q_mat_d_d_uniform(A, public_seed, DOMAIN_SEPARATOR_A, 0);

  kappa = 0;
  do {
    for (size_t i = 0; i < PARAM_KH; i++)
      for (size_t j = 0; j < 2; j++)
        poly_q_mat_d_d_binomial(sk->R[j][i], secret_seed, kappa++, DOMAIN_SEPARATOR_R);
  } while(sk_sq_spectral_norm(RRstar, sk->R) > PARAM_R_MAX_SQ_SPECTRAL_NORM);

  for (size_t i = 0; i < PARAM_KH; i++) {
    poly_q_mat_d_d_mul_mat_d_d(pk->B[i], A, sk->R[1][i]);
    poly_q_mat_d_d_add(pk->B[i], pk->B[i], sk->R[0][i]);
  }
  for (size_t i = 0; i < SEED_BYTES; i++) pk->seed[i] = public_seed[i];
  compute_covariance(sk->S, RRstar);

  poly_q_mat_d_d_clear(A);
  poly_q_mat_d_d_clear(RRstar[0][0]); poly_q_mat_d_d_clear(RRstar[0][1]);
  poly_q_mat_d_d_clear(RRstar[1][0]); poly_q_mat_d_d_clear(RRstar[1][1]);
}

void _sep_sign_commitment(sep_sig_t *sig, uint8_t state[STATE_BYTES],
    const sep_sk_t *sk, const sep_pk_t *pk, const poly_q_vec_d cmt)
{
  size_t i;
  poly_q_mat_d_d A;
  poly_q_vec_d u, a3, tmp, v11;
  poly_q taginv, d_tmp;
  uint64_t n11, n12, n2, n3;

  poly_q_mat_d_d_init(A);
  poly_q_vec_d_init(u); poly_q_vec_d_init(a3);
  poly_q_vec_d_init(tmp); poly_q_vec_d_init(v11);
  poly_q_init(taginv); poly_q_init(d_tmp);

  poly_q_mat_d_d_uniform(A, pk->seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_vec_d_uniform(u, pk->seed, DOMAIN_SEPARATOR_U);
  poly_q_vec_d_uniform(a3, pk->seed, DOMAIN_SEPARATOR_A3);

  /* compute tag from state */
  poly_q_binary_fixed_weight(sig->tag, state);

  /* increment state */
  uint64_t stateinc = ((uint64_t)*(uint32_t*)state) + 1;
  *(uint32_t*)state = (uint32_t) stateinc;
  uint64_t carry = stateinc >> 32;
  for (i = 1; i < STATE_BYTES/4; i++) {
    stateinc = ((uint64_t)*(uint32_t*)&state[i*4]) + carry;
    *(uint32_t*)&state[i*4] = (uint32_t) stateinc;
    carry = stateinc >> 32;
  }

  /* invert tag modulo bH */
  poly_q_invert_mod_bH(taginv, sig->tag);

  /* add commitment to u */
  poly_q_vec_d_add(u, u, cmt);

reject_sig:
  /* sample v3 from discrete Gaussian s4 */
  poly_q_sample_gaussian_s4(sig->v3);

  /* compute w = u + cmt - a3.v3 */
  poly_q_vec_d_mul_poly(tmp, a3, sig->v3);
  poly_q_vec_d_sub(tmp, u, tmp);

  /* call TSampler */
  sampler(v11, sig->v12, sig->v2, sk->R, A, pk->B, tmp, sig->tag, taginv, sk->S);

  /* norm check */
  n11 = poly_q_vec_d_norm2(v11);
  n12 = poly_q_vec_d_norm2(sig->v12);
  n2  = poly_q_vec_d_norm2(sig->v2[0]);
  for (i = 1; i < PARAM_KH; i++) n2 += poly_q_vec_d_norm2(sig->v2[i]);
  n3  = poly_q_sq_norm2(sig->v3);

  if ((n11 > PARAM_B11SQ) || (n12 > PARAM_B12SQ) ||
      (n2  > PARAM_B2SQ)  || (n3  > PARAM_B3SQ))
    goto reject_sig;
  poly_q_mat_d_d_clear(A);
  poly_q_vec_d_clear(u); poly_q_vec_d_clear(a3);
  poly_q_vec_d_clear(tmp); poly_q_vec_d_clear(v11);
  poly_q_clear(taginv); poly_q_clear(d_tmp);
}

void sep_sign(sep_sig_t *sig, uint8_t state[STATE_BYTES],
    const sep_sk_t *sk, const sep_pk_t *pk,
    const uint8_t msg[PARAM_M*PARAM_N/8])
{
  poly_q_vec_d cmt;
  poly_q_mat_d_m D;
  poly_q_vec_m m;
  poly_q_vec_d_init(cmt);
  poly_q_mat_d_m_init(D);
  poly_q_vec_m_init(m);

  for (size_t i = 0; i < PARAM_M; i++)
    poly_q_from_bits(m->entries[i], &msg[i * PARAM_N/8]);
  poly_q_mat_d_m_uniform(D, pk->seed, DOMAIN_SEPARATOR_D);
  poly_q_mat_d_m_mul_vec_m(cmt, D, m);

  _sep_sign_commitment(sig, state, sk, pk, cmt);

  poly_q_vec_d_clear(cmt);
  poly_q_mat_d_m_clear(D);
  poly_q_vec_m_clear(m);
}

int _sep_verify_from_commitment(const sep_sig_t *sig, const poly_q_vec_d cmt, const sep_pk_t *pk)
{
  size_t i, j;
  int64_t bexpi;
  poly_q_vec_d u, a3, v11;
  poly_q_mat_d_d A;
  poly_q tag_times_bexpi;
  uint64_t n11, n12, n2, n3;
  int64_t tag_weight;

  poly_q_mat_d_d_init(A);
  poly_q_vec_d_init(u); poly_q_vec_d_init(a3);
  poly_q_vec_d_init(v11);
  poly_q_init(tag_times_bexpi);

  poly_q_mat_d_d_uniform(A, pk->seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_vec_d_uniform(u, pk->seed, DOMAIN_SEPARATOR_U);
  poly_q_vec_d_uniform(a3, pk->seed, DOMAIN_SEPARATOR_A3);

  /* v11 = u + cmt - A.v12 + (B - t.G).v2 - a3.v3 */
  poly_q_vec_d_add(v11, u, cmt);
  poly_q_vec_d_mul_poly(u, a3, sig->v3);    /* u = a3.v3 */
  poly_q_mat_d_d_muladd_vec_d(u, A, sig->v12); /* u += A.v12 */
  poly_q_vec_d_sub(v11, v11, u);

  bexpi = PARAM_QL;
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_mul_scalar(tag_times_bexpi, sig->tag, bexpi);
    poly_q_mat_d_d_set(A, pk->B[i]);
    for (j = 0; j < PARAM_D; j++)
      poly_q_sub(A->rows[j]->entries[j], A->rows[j]->entries[j], tag_times_bexpi);
    poly_q_mat_d_d_muladd_vec_d(v11, A, sig->v2[i]);
    bexpi *= PARAM_BH;
  }

  n11 = poly_q_vec_d_norm2(v11);
  n12 = poly_q_vec_d_norm2(sig->v12);
  n2  = poly_q_vec_d_norm2(sig->v2[0]);
  for (i = 1; i < PARAM_KH; i++) n2 += poly_q_vec_d_norm2(sig->v2[i]);
  n3  = poly_q_sq_norm2(sig->v3);
  tag_weight = poly_q_weight(sig->tag);

  poly_q_mat_d_d_clear(A);
  poly_q_vec_d_clear(u); poly_q_vec_d_clear(a3);
  poly_q_vec_d_clear(v11);
  poly_q_clear(tag_times_bexpi);

  return (n11 <= PARAM_B11SQ) && (n12 <= PARAM_B12SQ) &&
         (n2  <= PARAM_B2SQ)  && (n3  <= PARAM_B3SQ) &&
         (tag_weight == PARAM_W);
}

int sep_verify(const sep_sig_t *sig, const uint8_t msg[PARAM_M*PARAM_N/8], const sep_pk_t *pk)
{
  poly_q_vec_d cmt;
  poly_q_mat_d_m D;
  poly_q_vec_m m;
  poly_q_vec_d_init(cmt);
  poly_q_mat_d_m_init(D);
  poly_q_vec_m_init(m);

  for (size_t i = 0; i < PARAM_M; i++)
    poly_q_from_bits(m->entries[i], &msg[i * PARAM_N/8]);
  poly_q_mat_d_m_uniform(D, pk->seed, DOMAIN_SEPARATOR_D);
  poly_q_mat_d_m_mul_vec_m(cmt, D, m);

  int r = _sep_verify_from_commitment(sig, cmt, pk);

  poly_q_vec_d_clear(cmt);
  poly_q_mat_d_m_clear(D);
  poly_q_vec_m_clear(m);
  return r;
}
