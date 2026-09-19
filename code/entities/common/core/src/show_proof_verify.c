#include "arith.h"
#include "randombytes.h"
#include "poly_qshow_sampling.h"
#include "sep.h"
#include "show.h"
#include "macros.h"
#include "fips202.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static _Thread_local unsigned int show_verify_thread_cap = 0;

void show_verify_set_thread_cap(unsigned int cap) {
  show_verify_thread_cap = cap;
}

static int show_verify_apply_thread_cap(int threads) {
  if(show_verify_thread_cap != 0 &&
     threads > (int)show_verify_thread_cap){
    return (int)show_verify_thread_cap;
  }
  return threads;
}

static int show_verify_relation_parallel_enabled(void) {
  const char *value = getenv("REVAC_SHOW_VERIFY_RELATION_PARALLEL");

  if (value != NULL && value[0] != '\0') {
    return value[0] != '0' && show_verify_apply_thread_cap(4) > 1;
  }
#ifdef _OPENMP
  return omp_get_num_procs() >= 4 && show_verify_apply_thread_cap(4) > 1;
#else
  return 0;
#endif
}

static int show_verify_first_message_parallel_enabled(void) {
  const char *value = getenv("REVAC_SHOW_VERIFY_FIRST_MESSAGE_PARALLEL");

  if (value != NULL && value[0] != '\0') {
    return value[0] != '0' && show_verify_apply_thread_cap(2) > 1;
  }
#ifdef _OPENMP
  return omp_get_num_procs() >= 4 && show_verify_apply_thread_cap(2) > 1;
#else
  return 0;
#endif
}

static int show_verify_challenge_threads(void) {
  const char *value = getenv("REVAC_SHOW_CHALLENGE_THREADS");
  long requested = 0;

  if (value != NULL && value[0] != '\0') {
    char *end = NULL;
    requested = strtol(value, &end, 10);
    if (end != value && requested > 0 && requested <= 64) {
      return show_verify_apply_thread_cap((int)requested);
    }
  }
#ifdef _OPENMP
  {
    const int available = omp_get_num_procs();
    return show_verify_apply_thread_cap(available < 8 ? available : 8);
  }
#else
  return 1;
#endif
}

static int show_verify_profile_detail_enabled(void) {
  const char *value = getenv("REVAC_SHOW_PROFILE_DETAIL");
  return value != NULL && value[0] != '\0' && value[0] != '0';
}

static double show_verify_monotonic_seconds(void) {
  struct timespec now;
  if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return 0.0;
  return (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
}

static void show_mat_k_k_from_vec_k(poly_qshow_mat_k_k res, const poly_qshow_vec_k arg) {
  poly_qshow tmp_x[PARAM_K_SHOW];

  for (size_t i = 0; i < PARAM_K_SHOW; i++) {
    poly_qshow_init(tmp_x[i]);
  }
  for (size_t i = 1; i < PARAM_K_SHOW; i++) {
    poly_qshow_mul_x(tmp_x[i], arg->entries[i]);
  }
  for (size_t i = 0; i < PARAM_K_SHOW; i++) {
    for (size_t j = 0; j < i + 1; j++) {
      poly_qshow_set(res->rows[i]->entries[j], arg->entries[i - j]);
    }
    for (size_t j = i + 1; j < PARAM_K_SHOW; j++) {
      poly_qshow_set(res->rows[i]->entries[j], tmp_x[PARAM_K_SHOW + i - j]);
    }
  }
  for (size_t i = 0; i < PARAM_K_SHOW; i++) {
    poly_qshow_clear(tmp_x[i]);
  }
}

/*************************************************
* Name:        show_verify
*
* Description: Verification of the zero-knowledge proof of valid 
*              credentials (anonymous credentials show)   有效凭证（匿名凭证显示）的零知识证明验证
*
* Arguments:   - const osig_proof *proof: pointer to show proof structure
*              - const poly_qshow_mat_k_k *A_embed: array of polynomial matrices, subring embedding of q1.A' 多项式矩阵数组，q1.A' 的子环嵌入
*              - const poly_qshow_mat_k_k *B_embed: array of polynomial matrices, subring embedding of q1.B
*              - const poly_qshow_mat_k_k *A3_embed: array of polynomial matrices, subring embedding of q1.A3
*              - const poly_qshow_mat_k_k *Ds_embed: array of polynomial matrices, subring embedding of q1.Ds
*              - const poly_qshow_mat_k_k *D_embed: array of polynomial matrices, subring embedding of q1.D
*              - const poly_qshow_vec_k *u_embed: array of polynomial vectors, subring embedding of q1.u
*              - const uint8_t *crs_seed: pointer to byte array containing the CRS seed (allocated SEED_BYTES bytes) 指向包含 CRS 种子的字节数组的指针（分配了 SEED_BYTES 字节）
*              - const uint8_t *seed: pointer to byte array containing the seed 
*                   for public parameters (allocated SEED_BYTES bytes)指向包含种子的字节数组的指针，用于公共参数（分配 SEED_BYTES 字节）
* 
* Returns 1 if proof could be verified correctly and 0 otherwise
**************************************************/
int show_verify_handle(
    const show_proof_t       *proof, 
    const poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D*PARAM_K], 
    const poly_qshow_mat_k_k A3_embed[PARAM_D],
    const poly_qshow_mat_k_k Ds_embed[PARAM_D][2*PARAM_D], 
    const poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], 
    const poly_qshow_mat_k_k Dx_embed[PARAM_D],
    const poly_qshow_mat_k_k A_link_embed[PARAM_D][PARAM_D],
    const poly_qshow_mat_k_k D_link_embed[PARAM_D],
    const poly_qshow_vec_k   C_link_embed[PARAM_D],
    const poly_qshow_vec_k   u_embed[PARAM_D], 
    const uint8_t            crs_seed[CRS_SEED_BYTES], 
    const uint8_t            seed[SEED_BYTES]) {
  size_t i,j,k;
  int is_valid = 1;
  uint8_t buf[CHAL4_SHOW_INPUT_BYTES] = {0}, challenge_seed[SEED_BYTES];
  uint64_t sq_norm_z3;
  uint128 sq_norm_z1, sq_norm_z2;
  int64_t bexpi;
  coeff_qshow chal_2[PARAM_L_SHOW][PARAM_ARP_SHOW + 6];
  coeff_qshow tmp_coeff;
  poly_qshow tmp_poly, z1s_z1, z1s_c_one, sum_mu_gamma[7], t0, c, c_one;
  poly_qshow signature_rows[PARAM_D * PARAM_K_SHOW];
  poly_qshow link_rows[PARAM_D * PARAM_K_SHOW];
  poly_qshow_mat_d_m1 A1;
  poly_qshow_mat_d_m2 A2;
  poly_qshow_mat_256l_m2 Byg;
  poly_qshow_vec_d tmp_vec_d, a1_z1, w;
  poly_qshow_vec_k tmp_vec_k;
  poly_qshow_vec_m1 sum_gamma_r_star_i, tmp_vec_m1;
  poly_qshow_vec_k Gz1_v2[PARAM_D], chal_3_dk[PARAM_D];
  poly_qshow_vec_k chal_3_link[PARAM_D];
  poly_qshow sum_gamma_e_star_ij;
  poly_qshow_vec_m1 chal_1[PARAM_ARP_SHOW];
  poly_qshow_vec_m2 b;
  poly_qshow_vec_m2 z2_opening;
  poly_qshow_vec_256_l tmp_vec_256_l, Z;
  poly_qshow_vec_l chal_3_l;
  poly_qshow_mat_k_k ztag_mat;
  const double verify_started = show_verify_monotonic_seconds();
  double relation_started = 0.0;
  double relation_seconds = 0.0;

  // init
  // init polynomials
  poly_qshow_init(tmp_poly);
  poly_qshow_init(z1s_z1);
  poly_qshow_init(z1s_c_one);
  poly_qshow_init(c);
  poly_qshow_init(c_one);
  poly_qshow_init(t0);
  // init vectors and matrices
  poly_qshow_mat_d_m1_init(A1);
  poly_qshow_mat_d_m2_init(A2);
  poly_qshow_mat_256l_m2_init(Byg);
  poly_qshow_vec_k_init(tmp_vec_k);
  poly_qshow_mat_k_k_init(ztag_mat);
  poly_qshow_vec_m1_init(tmp_vec_m1);
  poly_qshow_vec_m1_init(sum_gamma_r_star_i);
  poly_qshow_vec_d_init(tmp_vec_d);
  poly_qshow_vec_d_init(a1_z1);
  poly_qshow_vec_d_init(w);
  poly_qshow_init(sum_gamma_e_star_ij);    
  poly_qshow_vec_m2_init(b);
  poly_qshow_vec_m2_init(z2_opening);
  poly_qshow_vec_256_l_init(tmp_vec_256_l);
  poly_qshow_vec_256_l_init(Z);
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    poly_qshow_vec_m1_init(chal_1[i]);
  }
  poly_qshow_vec_l_init(chal_3_l);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_init(chal_3_dk[i]);
    poly_qshow_vec_k_init(chal_3_link[i]);
    poly_qshow_vec_k_init(Gz1_v2[i]);
  }
  for (i = 0; i < 7; i++) {
    poly_qshow_init(sum_mu_gamma[i]);
    poly_qshow_zero(sum_mu_gamma[i]);
  }
  for (i = 0; i < PARAM_D * PARAM_K_SHOW; i++) {
    poly_qshow_init(signature_rows[i]);
    poly_qshow_init(link_rows[i]);
  }
  ////////////////////////////////////// NO INITIALIZATIONS BELOW THIS LINE

  // expanding CRS
  poly_qshow_mat_d_m1_uniform(A1, crs_seed, DOMAIN_SEPARATOR_A1_SHOW);
  poly_qshow_mat_d_m2_uniform(A2, crs_seed, DOMAIN_SEPARATOR_A2_SHOW);
  poly_qshow_mat_256l_m2_uniform(Byg, crs_seed, DOMAIN_SEPARATOR_BYG_SHOW);
  poly_qshow_vec_m2_uniform(b, crs_seed, DOMAIN_SEPARATOR_B_SHOW);

  // byte-packing the statement (only seeds)
  // the signer public key does not have to be hash (single-issuer setting)
  for (i = 0; i < CRS_SEED_BYTES; i++) {
    buf[1 + i] = crs_seed[i];
  }
  for (i = 0; i < SEED_BYTES; i++) {
    buf[1 + CRS_SEED_BYTES + i] = seed[i];
  }

  // precomputation
  for (i = 0; i < PARAM_N_SHOW; i++) {
    poly_qshow_set_coeff(c_one, i, 1); // poly with all one
  }
  poly_qshow_mul(c_one, c_one, proof->c); // holds c.one

  /****** reconstructing first message ******/
  /*
   * Compressed ABDLOP opening.  proof->tA contains tA,1 and only the first
   * m2-d entries of z2 are transmitted.  First compute
   *   v = A1*z1 + A2'*z2,1 - c*2^D*tA,1,
   * then recover w1 with the generalized hints.
   */
  poly_qshow_vec_m2_zero(z2_opening);
  for (i = 0; i < PARAM_M2_COMMIT_SHOW; i++) {
    poly_qshow_set(z2_opening->entries[i], proof->z2->entries[i]);
  }
  for (i = 0; i < PARAM_D_SHOW; i++) {
    poly_qshow_mul_scalar(tmp_vec_d->entries[i], proof->tA->entries[i],
                          (int64_t)1 << PARAM_D_ROUND_SHOW);
    poly_qshow_mul(tmp_vec_d->entries[i], tmp_vec_d->entries[i], proof->c);
  }
  #ifdef _OPENMP
  #pragma omp parallel sections num_threads(show_verify_apply_thread_cap(2)) if(show_verify_first_message_parallel_enabled())
  #endif
  {
  #ifdef _OPENMP
  #pragma omp section
  #endif
    poly_qshow_mat_d_m2_mul_vec_m2(w, A2, z2_opening);
  #ifdef _OPENMP
  #pragma omp section
  #endif
    poly_qshow_mat_d_m1_mul_vec_m1(a1_z1, A1, proof->z1);
  }
  poly_qshow_vec_d_sub(w, w, tmp_vec_d);
  poly_qshow_vec_d_add(w, w, a1_z1);
  show_use_g_hints_vec_d(tmp_vec_d, proof->hints, w);
  for (i = 0; i < PARAM_D_SHOW; i++) {
    poly_qshow_mul_scalar(z2_opening->entries[PARAM_M2_COMMIT_SHOW + i],
                          tmp_vec_d->entries[i], PARAM_GAMMA_SHOW);
    poly_qshow_sub(z2_opening->entries[PARAM_M2_COMMIT_SHOW + i],
                   w->entries[i],
                   z2_opening->entries[PARAM_M2_COMMIT_SHOW + i]);
    for (j = 0; j < PARAM_N_SHOW; j++) {
      int64_t hint = (int64_t)poly_qshow_get_coeff_centered(
          proof->hints->entries[i], j);
      const int hint_valid =
          (hint > -(PARAM_HINT_MOD_SHOW / 2)) &&
          (hint <= PARAM_HINT_MOD_SHOW / 2);
      if (!hint_valid && getenv("REVAC_SHOW_VERIFY_DEBUG") != NULL) {
        fprintf(stderr, "[show verify debug] generalized hint range failed\n");
      }
      is_valid = is_valid && hint_valid;
    }
  }
  poly_qshow_vec_d_set(w, tmp_vec_d);
  // recomputing first challenge
  buf[0] = 1;
  poly_qshow_vec_d_pack(buf + SHOW_CHALLENGE_BASE_BYTES, proof->tA);
  poly_qshow_vec_d_pack(buf + SHOW_CHALLENGE_BASE_BYTES + POLYQSHOW_VECD_PACKEDBYTES, w);
  poly_qshow_vec_256_l_pack(buf + SHOW_CHALLENGE_BASE_BYTES + 2*POLYQSHOW_VECD_PACKEDBYTES, proof->tB);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL1_SHOW_INPUT_BYTES);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(show_verify_challenge_threads()) \
    if(show_verify_challenge_threads() > 1)
#endif
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    poly_qshow_vec_m1_binomial(chal_1[i], challenge_seed, DOMAIN_SEPARATOR_CHAL1_SHOW, i, SEED_BYTES);
  }

  /****** reconstructing second message ******/
  buf[0] = 2;
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    coeff_qshow_pack(buf + CHAL1_SHOW_INPUT_BYTES + i*COEFFQSHOW_PACKEDBYTES, proof->z3[i]);
  }
  shake256(challenge_seed, SEED_BYTES, buf, CHAL2_SHOW_INPUT_BYTES);
  for (i = 0; i < PARAM_L_SHOW; i++) {
    vec_qshow_uniform(chal_2[i], challenge_seed, DOMAIN_SEPARATOR_CHAL2_SHOW, i, SEED_BYTES); // writes PARAM_ARP_SHOW + 6 uniform numbers to the first argument
  }

  /****** reconstructing third message ******/
  buf[0] = 3;
  poly_qshow_vec_l_pack(buf + CHAL2_SHOW_INPUT_BYTES, proof->h);
  poly_qshow_vec_l_uniform(chal_3_l, buf, DOMAIN_SEPARATOR_CHAL3_SHOW, CHAL3_SHOW_INPUT_BYTES);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL3_SHOW_INPUT_BYTES);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_uniform(chal_3_dk[i], challenge_seed, DOMAIN_SEPARATOR_CHAL3_SHOW, i+PARAM_L_SHOW, SEED_BYTES);
    poly_qshow_vec_k_uniform(chal_3_link[i], challenge_seed,
                             DOMAIN_SEPARATOR_CHAL3_SHOW,
                             i + PARAM_L_SHOW + PARAM_D, SEED_BYTES);
  }

  /****** reconstructing fourth message ******/
  /**********************************************
  * Recomputing commitment t0
  *   Z = c.tB - Byg.z2
  *   
  *   t0 = zFz + c.<f,z> + c^2.f + <b,z2> - c.t1
  *      = <b, z2> - c.t1 + sum_mu_gamma[0].(<z1_v1*, z1_v1> - c^2.B1^2)
  *         + sum_mu_gamma[1].(<z1_v2*, z1_v2> - c^2.B2^2)
  *         + sum_mu_gamma[2].(<z1_v3*, z1_v3> - c^2.B3^2) 
  *         + sum_mu_gamma[6].<z1_t*, z1_t>
  *         - sum_mu_gamma[3].c^2.w - sum_mu_gamma[4].<z1_t*, c.one>
  *         + sum_mu_gamma[5].<z1_m*, z1_m - c.one>
  *         + c.sum_{i < l} mu_i.(<sum_gamma_r_star[i], z1> + Z_{256/n + i}
  *         + <sum_gamma_e_star[i], Z_{:256/n}> - c.(sum_gamma_z3 + hi))
  *         + <z1_t, sum_{i < dk} mu_{l+i} G_i".z1_v2> 
  *         + c.sum_{i < dk} mu_{l+i}.([q1.I | A_embed].z1_v1 - B_embed.z1_v2 + A3_embed.z1_v3 - Ds_embed.z1_usk - D_embed.z1_m - c.u_embed)_i
  **********************************************/
  // computing <b, z2> - c.t1
  poly_qshow_vec_m2_mul_inner(t0, b, proof->z2);
  poly_qshow_mul(tmp_poly, proof->c, proof->t1);
  poly_qshow_sub(t0, t0, tmp_poly);
  
  
  
  // c holds proof->c^2 for now
  poly_qshow_mul(c, proof->c, proof->c);

  // compute sum_i mu_i gamma_{i,256 + j} for j in [6]
  for (i = 0; i < PARAM_L_SHOW; i++) {
    for (j = 0; j < 6; j++) {
      poly_qshow_mul_scalar(tmp_poly, chal_3_l->entries[i], chal_2[i][PARAM_ARP_SHOW + j]);
      poly_qshow_add(sum_mu_gamma[j], sum_mu_gamma[j], tmp_poly);
    }
  }
  poly_qshow_add(sum_mu_gamma[6], sum_mu_gamma[3], sum_mu_gamma[4]); // both sum_mu_gamma[3] and sum_mu_gamma[4] will be used later

  // computing Z = c.tB - Byg.z2
  poly_qshow_mat_256l_m2_mul_vec_m2(tmp_vec_256_l, Byg, proof->z2);
  poly_qshow_vec_256_l_mul_poly_qshow(Z, proof->tB, proof->c);
  poly_qshow_vec_256_l_sub(Z, Z, tmp_vec_256_l);

  // quadratic terms
  // <z1_v1"*,z1_v1"> - c^2 B1'^2
  poly_qshow_zero(z1s_z1);
  for (i = IDX_V1_SHOW; i < IDX_V2_SHOW; i++) {
    poly_qshow_conjugate(tmp_poly, proof->z1->entries[i]);
    poly_qshow_mul(tmp_poly, tmp_poly, proof->z1->entries[i]);
    poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
  }
  poly_qshow_mul_scalar(tmp_poly, c, PARAM_B1SQ); // c^2 B1'^2
  poly_qshow_sub(z1s_z1, z1s_z1, tmp_poly);
  poly_qshow_mul(z1s_z1, z1s_z1, sum_mu_gamma[0]); 
  poly_qshow_add(t0, t0, z1s_z1); 
  
  
  // <z1_v2"*,z1_v2"> - c^2 B2^2
  poly_qshow_zero(z1s_z1);
  for (i = IDX_V2_SHOW; i < IDX_V3_SHOW; i++) {
    poly_qshow_conjugate(tmp_poly, proof->z1->entries[i]);
    poly_qshow_mul(tmp_poly, tmp_poly, proof->z1->entries[i]);
    poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
  }
  poly_qshow_mul_scalar(tmp_poly, c, PARAM_B2SQ); // c^2 B2^2
  poly_qshow_sub(z1s_z1, z1s_z1, tmp_poly);
  poly_qshow_mul(z1s_z1, z1s_z1, sum_mu_gamma[1]); 
  poly_qshow_add(t0, t0, z1s_z1); 
  
  
  // <z1_v3"*,z1_v3"> - c^2 B3^2
  poly_qshow_zero(z1s_z1);
  for (i = IDX_V3_SHOW; i < IDX_TAG_SHOW; i++) {
    poly_qshow_conjugate(tmp_poly, proof->z1->entries[i]);
    poly_qshow_mul(tmp_poly, tmp_poly, proof->z1->entries[i]);
    poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
  }
  poly_qshow_mul_scalar(tmp_poly, c, PARAM_B3SQ); // c^2 B3^2
  poly_qshow_sub(z1s_z1, z1s_z1, tmp_poly);
  poly_qshow_mul(z1s_z1, z1s_z1, sum_mu_gamma[2]); 
  poly_qshow_add(t0, t0, z1s_z1); 
  
  
  // -c^2.w
  poly_qshow_mul_scalar(tmp_poly, c, PARAM_W); // c^2 w
  poly_qshow_mul(tmp_poly, tmp_poly, sum_mu_gamma[3]);
  poly_qshow_sub(t0, t0, tmp_poly); // substract
  
  
  // no more need for c^2, c can be used as temp variable
  // <z1_t'*,z1_t'> and <z1_t'*,c.one>
  poly_qshow_zero(z1s_z1);
  poly_qshow_zero(z1s_c_one);
  for (i = IDX_TAG_SHOW; i < IDX_USK_SHOW; i++) {
    poly_qshow_conjugate(tmp_poly, proof->z1->entries[i]);
    poly_qshow_mul(c, tmp_poly, proof->z1->entries[i]);
    poly_qshow_add(z1s_z1, z1s_z1, c);
    poly_qshow_mul(c, tmp_poly, c_one);
    poly_qshow_add(z1s_c_one, z1s_c_one, c);
  }
  poly_qshow_mul(z1s_z1, z1s_z1, sum_mu_gamma[6]);
  poly_qshow_add(t0, t0, z1s_z1);
  
    poly_qshow_mul(z1s_c_one, z1s_c_one, sum_mu_gamma[4]);
  poly_qshow_sub(t0, t0, z1s_c_one); // substract
  

  // z1s_c_one can be used as temp variable
  // <z1_m'*,z1_m' - c.one>
  poly_qshow_zero(z1s_z1);
  for (i = IDX_USK_SHOW; i < IDX_SHOW_WITNESS_END; i++) { // same for usk, all message attributes, and handle x
    poly_qshow_conjugate(tmp_poly, proof->z1->entries[i]);
    poly_qshow_sub(z1s_c_one, proof->z1->entries[i], c_one);
    poly_qshow_mul(tmp_poly, tmp_poly, z1s_c_one);
    poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
  }
  poly_qshow_mul(z1s_z1, z1s_z1, sum_mu_gamma[5]);
  poly_qshow_add(t0, t0, z1s_z1);

  // TSampler gadget quadratic part, using the same mu-weighted
  // multiplication representation as the public linear relation.
  for (i = 0; i < PARAM_K_SHOW; i++) {
    poly_qshow_set(tmp_vec_k->entries[i], proof->z1->entries[IDX_TAG_SHOW + i]);
  }
  show_mat_k_k_from_vec_k(ztag_mat, tmp_vec_k);
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_zero(Gz1_v2[i]->entries[j]);
      bexpi = (int64_t)PARAM_Q1_SHOW * PARAM_QL;
      for (k = 0; k < PARAM_KH; k++) {
        poly_qshow_mul_scalar(tmp_poly, proof->z1->entries[IDX_V2_SHOW + k*PARAM_D*PARAM_K_SHOW + i*PARAM_K_SHOW + j], bexpi);
        poly_qshow_add(Gz1_v2[i]->entries[j], Gz1_v2[i]->entries[j], tmp_poly);
        bexpi *= PARAM_BH;
      }
    }
    poly_qshow_mat_k_k_mul_vec_k(tmp_vec_k, ztag_mat, Gz1_v2[i]);
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, tmp_vec_k->entries[j], chal_3_dk[i]->entries[j]);
      poly_qshow_add(t0, t0, tmp_poly);
    }
  }
  // Linear part depending on sum_gamma_r_star and sum_gamma_e_star
  // computing Y = c.sum_{i < l} mu_i.(<sum_gamma_r_star[i], z1> + Z_{256/n + i} + <sum_gamma_e_star[i], Z_{:256/n}> - c.(sum_gamma_z3 + hi))
  // conjugating all chal_1
  for (j = 0; j < PARAM_ARP_SHOW; j++) {
    poly_qshow_vec_m1_conjugate(tmp_vec_m1, chal_1[j]);
    poly_qshow_vec_m1_set(chal_1[j], tmp_vec_m1);
  }
  
  // z1s_z1, z1s_c_one can be used as temp variables
  poly_qshow_zero(z1s_c_one); // will hold Y/c
  for (i = 0; i < PARAM_L_SHOW; i++) {
    poly_qshow_set(tmp_poly, proof->h->entries[i]); // h_i
    const int h_constant_valid = poly_qshow_get_coeff(tmp_poly, 0) == 0;
    if (!h_constant_valid && getenv("REVAC_SHOW_VERIFY_DEBUG") != NULL) {
      fprintf(stderr, "[show verify debug] h[%zu] constant check failed\n", i);
    }
    is_valid = is_valid && h_constant_valid;
    if (!is_valid) {
      goto show_verify_cleanup; // h_i does not have a zero coefficient, is_valid = 0
    }
    for (j = 0; j < PARAM_ARP_SHOW; j++) {
      poly_qshow_muladd_constant(tmp_poly, chal_2[i][j], proof->z3[j]); // adds (chal_2[i][j] * proof->z3[j]) to the constant coefficient without overflow
    }
    poly_qshow_mul(tmp_poly, tmp_poly, proof->c); 
    poly_qshow_sub(tmp_poly, Z->entries[PARAM_ARP_DIV_N_SHOW + i], tmp_poly); // Z_{256/n + i} - c.(sum_gamma_z3 + hi) 

    // sum of gamma_{ij}r_j*z1
    for (j = 0; j < PARAM_ARP_SHOW; j++) {
      if (j == 0) {
        poly_qshow_vec_m1_mul_scalar(sum_gamma_r_star_i, chal_1[j], chal_2[i][j]);
      } else {
        poly_qshow_vec_m1_mul_scalar(tmp_vec_m1, chal_1[j], chal_2[i][j]);
        poly_qshow_vec_m1_add(sum_gamma_r_star_i, sum_gamma_r_star_i, tmp_vec_m1);
      }
    }
    poly_qshow_vec_m1_mul_inner(z1s_z1, sum_gamma_r_star_i, proof->z1);
    poly_qshow_add(tmp_poly, tmp_poly, z1s_z1);

    // sum of gamma_{ij}e_j*Z = <conjugate(tau^-1([gamma_{i,0} | ... | gamma_{i,256}])), Z_{:256/n}>
    for (j = 0; j < PARAM_ARP_DIV_N_SHOW; j++) {
      poly_qshow_set_coeff(sum_gamma_e_star_ij, 0, chal_2[i][j * PARAM_N_SHOW + 0]);
      for (k = 1; k < PARAM_N_SHOW; k++) {
        poly_qshow_set_coeff(sum_gamma_e_star_ij, k, - chal_2[i][j * PARAM_N_SHOW + (PARAM_N_SHOW - k)]); // set conjugate directly
      }
      poly_qshow_mul(sum_gamma_e_star_ij, sum_gamma_e_star_ij, Z->entries[j]);
      poly_qshow_add(tmp_poly, tmp_poly, sum_gamma_e_star_ij);
    }
    poly_qshow_mul(tmp_poly, tmp_poly, chal_3_l->entries[i]);
    poly_qshow_add(z1s_c_one, z1s_c_one, tmp_poly);
  }
  poly_qshow_mul(z1s_c_one, z1s_c_one, proof->c);
  poly_qshow_add(t0, t0, z1s_c_one);
  
  
  // computing c.sum_{i < dk} mu_{l+i}.([q1.I | A_embed]z1_v1 - B_embed.z1_v2 + A3_embed.z1_v3 - Ds_embed.z1_usk - D_embed.z1_m - c.u_embed)_i
  // z1s_z1, z1s_c_one can be used as temp variables
  poly_qshow_zero(z1s_c_one); 
  relation_started = show_verify_monotonic_seconds();
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(show_verify_apply_thread_cap(4)) if(show_verify_relation_parallel_enabled()) schedule(static) private(j)
  #endif
  for (i = 0; i < PARAM_D*PARAM_K_SHOW; i++) {
    poly_qshow row_tmp;
    const size_t row_quot = i / PARAM_K_SHOW;
    const size_t row_rem = i % PARAM_K_SHOW;

    poly_qshow_init(row_tmp);
    // y1s_y1 used as temp variable to host ([A|-B|A3|-Ds|-D].[y1_v1|y1_v2|y1_v3|y1_s|y1_m])_i
    poly_qshow_mul_scalar(signature_rows[i], proof->z1->entries[IDX_V1_SHOW + i], PARAM_Q1_SHOW);
    for (j = 0; j < PARAM_D*PARAM_K_SHOW; j++) {
      poly_qshow_mul(row_tmp, A_embed[row_quot][j / PARAM_K_SHOW]->rows[row_rem]->entries[j % PARAM_K_SHOW], proof->z1->entries[IDX_V12_SHOW + j]);
      poly_qshow_add(signature_rows[i], signature_rows[i], row_tmp);
    }
    for (j = 0; j < PARAM_D*PARAM_KH*PARAM_K_SHOW; j++) {
      poly_qshow_mul(row_tmp, B_embed[row_quot][j / PARAM_K_SHOW]->rows[row_rem]->entries[j % PARAM_K_SHOW], proof->z1->entries[IDX_V2_SHOW + j]);
      poly_qshow_sub(signature_rows[i], signature_rows[i], row_tmp); // substraction
    }
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_mul(row_tmp, A3_embed[row_quot]->rows[row_rem]->entries[j], proof->z1->entries[IDX_V3_SHOW + j]);
      poly_qshow_add(signature_rows[i], signature_rows[i], row_tmp);
    }
    for (j = 0; j < 2*PARAM_D*PARAM_K_SHOW; j++) {
      poly_qshow_mul(row_tmp, Ds_embed[row_quot][j / PARAM_K_SHOW]->rows[row_rem]->entries[j % PARAM_K_SHOW], proof->z1->entries[IDX_USK_SHOW + j]);
      poly_qshow_sub(signature_rows[i], signature_rows[i], row_tmp); // substraction
    }
    for (j = 0; j < PARAM_M*PARAM_K_SHOW; j++) {
      poly_qshow_mul(row_tmp, D_embed[row_quot][j / PARAM_K_SHOW]->rows[row_rem]->entries[j % PARAM_K_SHOW], proof->z1->entries[IDX_M_SHOW + j]);
      poly_qshow_sub(signature_rows[i], signature_rows[i], row_tmp); // substraction
    }
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_mul(row_tmp, Dx_embed[row_quot]->rows[row_rem]->entries[j], proof->z1->entries[IDX_X_SHOW + j]);
      poly_qshow_sub(signature_rows[i], signature_rows[i], row_tmp);
    }
    poly_qshow_mul(row_tmp, proof->c, u_embed[row_quot]->entries[row_rem]);
    poly_qshow_sub(signature_rows[i], signature_rows[i], row_tmp);
    poly_qshow_mul(signature_rows[i], signature_rows[i], chal_3_dk[row_quot]->entries[row_rem]);
    poly_qshow_clear(row_tmp);
  }
  for (i = 0; i < PARAM_D * PARAM_K_SHOW; i++) {
    poly_qshow_add(z1s_c_one, z1s_c_one, signature_rows[i]);
  }
  poly_qshow_mul(z1s_c_one, z1s_c_one, proof->c);
  poly_qshow_add(t0, t0, z1s_c_one);
  relation_seconds = show_verify_monotonic_seconds() - relation_started;

  /* Reconstruct the separately challenged same-handle commitment relation:
   *   r0 + A_link r1 + D_link x - C_link = 0. */
  poly_qshow_zero(z1s_c_one);
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(show_verify_apply_thread_cap(4)) if(show_verify_relation_parallel_enabled()) schedule(static) private(j)
  #endif
  for (i = 0; i < PARAM_D * PARAM_K_SHOW; i++) {
    poly_qshow row_tmp;
    const size_t row_quot = i / PARAM_K_SHOW;
    const size_t row_rem = i % PARAM_K_SHOW;

    poly_qshow_init(row_tmp);
    poly_qshow_mul_scalar(link_rows[i],
                          proof->z1->entries[IDX_LINK_R0_SHOW + i],
                          PARAM_Q1_SHOW);
    for (j = 0; j < PARAM_D * PARAM_K_SHOW; j++) {
      poly_qshow_mul(
          row_tmp,
          A_link_embed[row_quot][j / PARAM_K_SHOW]
              ->rows[row_rem]->entries[j % PARAM_K_SHOW],
          proof->z1->entries[IDX_LINK_R1_SHOW + j]);
      poly_qshow_add(link_rows[i], link_rows[i], row_tmp);
    }
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_mul(row_tmp,
                     D_link_embed[row_quot]
                         ->rows[row_rem]->entries[j],
                     proof->z1->entries[IDX_X_SHOW + j]);
      poly_qshow_add(link_rows[i], link_rows[i], row_tmp);
    }
    poly_qshow_mul(row_tmp, proof->c,
                   C_link_embed[row_quot]->entries[row_rem]);
    poly_qshow_sub(link_rows[i], link_rows[i], row_tmp);
    poly_qshow_mul(link_rows[i], link_rows[i],
                   chal_3_link[row_quot]->entries[row_rem]);
    poly_qshow_clear(row_tmp);
  }
  for (i = 0; i < PARAM_D * PARAM_K_SHOW; i++) {
    poly_qshow_add(z1s_c_one, z1s_c_one, link_rows[i]);
  }
  poly_qshow_mul(z1s_c_one, z1s_c_one, proof->c);
  poly_qshow_add(t0, t0, z1s_c_one);

  // computing fourth challenge
  poly_qshow_zero(c);

  buf[0] = 4;
  /*
   * t0 is reconstructed from the responses and the statement.  Hashing a
   * prover-supplied t0 here would leave all quadratic/linear relation checks
   * outside the Fiat--Shamir transcript and make the proof unsound.
   */
  poly_qshow_pack(buf + CHAL3_SHOW_INPUT_BYTES, t0);
  poly_qshow_pack(buf + CHAL3_SHOW_INPUT_BYTES + POLYQSHOW_PACKEDBYTES, proof->t1);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL4_SHOW_INPUT_BYTES);
  
  poly_qshow_sample_challenge(c, challenge_seed, DOMAIN_SEPARATOR_CHAL4_SHOW, proof->ctr_c, SEED_BYTES);

  /* Fiat--Shamir binding: all reconstructed transcript messages, the CRS
   * seed, public seed, and session-derived statement must reproduce the
   * challenge carried by the proof.  Previously c was reconstructed but
   * never compared, so transcript and nonce modifications were accepted as
   * long as the response norms remained below their bounds. */
  {
    const int fs_valid = poly_qshow_equal(c, proof->c);
    if (!fs_valid && getenv("REVAC_SHOW_VERIFY_DEBUG") != NULL) {
      fprintf(stderr, "[show verify debug] reconstructed challenge mismatch\n");
    }
    is_valid = is_valid && fs_valid;
  }

  /****** checks ******/
  // computing z1 = y1 + c.s1, z2 = y2 + c.s2 and square l2 norms
  sq_norm_z1 = poly_qshow_vec_m1_norm2(proof->z1);
  sq_norm_z2 = poly_qshow_vec_m2_norm2(z2_opening);
  sq_norm_z3 = 0;
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    uint64_t abs_coeff;
    tmp_coeff = proof->z3[i];
    abs_coeff = tmp_coeff >= 0
                    ? (uint64_t)tmp_coeff
                    : (uint64_t)(-(int128)tmp_coeff);
    if (abs_coeff > UINT32_MAX) {
      sq_norm_z3 = UINT64_MAX;
      break;
    }
    CHK_UI_OVF_ADDITION(sq_norm_z3, abs_coeff * abs_coeff);
  }
  {
    const int z1_valid =
        sq_norm_z1 <= ((uint128)PARAM_B1SQ_SHOW_LOW64 +
                       ((uint128)PARAM_B1SQ_SHOW_HIGH64 << 64));
    const int z2_valid =
        sq_norm_z2 <= ((uint128)PARAM_B2SQ_SHOW_LOW64 +
                       ((uint128)PARAM_B2SQ_SHOW_HIGH64 << 64));
    const int z3_valid = sq_norm_z3 <= PARAM_B3SQ_SHOW;
    if ((!z1_valid || !z2_valid || !z3_valid) &&
        getenv("REVAC_SHOW_VERIFY_DEBUG") != NULL) {
      fprintf(stderr,
              "[show verify debug] response norm failed: z1=%d z2=%d z3=%d\n",
              z1_valid, z2_valid, z3_valid);
    }
    is_valid = is_valid && z1_valid && z2_valid && z3_valid;
  }
show_verify_cleanup:
  
  // clean up polynomials
  poly_qshow_clear(tmp_poly);
  poly_qshow_clear(z1s_z1);
  poly_qshow_clear(z1s_c_one);
  poly_qshow_clear(c);
  poly_qshow_clear(c_one);
  poly_qshow_clear(t0);
  // clean up vectors and matrices
  poly_qshow_mat_d_m1_clear(A1);
  poly_qshow_mat_d_m2_clear(A2);
  poly_qshow_mat_256l_m2_clear(Byg);
  poly_qshow_mat_k_k_clear(ztag_mat);
  poly_qshow_vec_k_clear(tmp_vec_k);
  poly_qshow_vec_m1_clear(tmp_vec_m1);
  poly_qshow_vec_m1_clear(sum_gamma_r_star_i);
  poly_qshow_vec_d_clear(tmp_vec_d);
  poly_qshow_vec_d_clear(a1_z1);
  poly_qshow_vec_d_clear(w);
  poly_qshow_clear(sum_gamma_e_star_ij);    
  poly_qshow_vec_m2_clear(b);
  poly_qshow_vec_m2_clear(z2_opening);
  poly_qshow_vec_256_l_clear(tmp_vec_256_l);
  poly_qshow_vec_256_l_clear(Z);
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    poly_qshow_vec_m1_clear(chal_1[i]);
  }
  poly_qshow_vec_l_clear(chal_3_l);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_clear(chal_3_dk[i]);
    poly_qshow_vec_k_clear(chal_3_link[i]);
    poly_qshow_vec_k_clear(Gz1_v2[i]);
  }
  for (i = 0; i < 7; i++) {
    poly_qshow_clear(sum_mu_gamma[i]);
  }
  for (i = 0; i < PARAM_D * PARAM_K_SHOW; i++) {
    poly_qshow_clear(link_rows[i]);
    poly_qshow_clear(signature_rows[i]);
  }

  if (show_verify_profile_detail_enabled()) {
    printf("[匿名签名 Verify 细分] 主公开关系重构 %.6f 秒；总计 %.6f 秒\n",
           relation_seconds,
           show_verify_monotonic_seconds() - verify_started);
  }
  return is_valid;
}

int show_verify(
    const show_proof_t       *proof,
    const poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
    const poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D*PARAM_K],
    const poly_qshow_mat_k_k A3_embed[PARAM_D],
    const poly_qshow_mat_k_k Ds_embed[PARAM_D][2*PARAM_D],
    const poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
    const poly_qshow_vec_k   u_embed[PARAM_D],
    const uint8_t            crs_seed[CRS_SEED_BYTES],
    const uint8_t            seed[SEED_BYTES]) {
  int ok;
  poly_qshow_mat_k_k Dx_embed[PARAM_D];
  poly_qshow_mat_k_k A_link_embed[PARAM_D][PARAM_D];
  poly_qshow_mat_k_k D_link_embed[PARAM_D];
  poly_qshow_vec_k C_link_embed[PARAM_D];

  for (size_t i = 0; i < PARAM_D; i++) {
    poly_qshow_mat_k_k_init(Dx_embed[i]);
    poly_qshow_mat_k_k_zero(Dx_embed[i]);
    poly_qshow_mat_k_k_init(D_link_embed[i]);
    poly_qshow_mat_k_k_zero(D_link_embed[i]);
    poly_qshow_vec_k_init(C_link_embed[i]);
    poly_qshow_vec_k_zero(C_link_embed[i]);
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_init(A_link_embed[i][j]);
      poly_qshow_mat_k_k_zero(A_link_embed[i][j]);
    }
  }
  ok = show_verify_handle(proof, A_embed, B_embed, A3_embed, Ds_embed, D_embed,
                          Dx_embed, A_link_embed, D_link_embed, C_link_embed,
                          u_embed, crs_seed, seed);
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_clear(A_link_embed[i][j]);
    }
    poly_qshow_vec_k_clear(C_link_embed[i]);
    poly_qshow_mat_k_k_clear(D_link_embed[i]);
    poly_qshow_mat_k_k_clear(Dx_embed[i]);
  }
  return ok;
}
