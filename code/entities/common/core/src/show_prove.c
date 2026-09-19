#include "arith.h"
#include "randombytes.h"
#include "poly_qshow_sampling.h"
#include "sep.h"
#include "show.h"
#include "macros.h"
#include "fips202.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static _Thread_local unsigned int show_user_prove_thread_cap = 0;

void show_user_prove_set_thread_cap(unsigned int cap) {
  show_user_prove_thread_cap = cap;
}

static int show_apply_thread_cap(int threads) {
  if (show_user_prove_thread_cap != 0 &&
      threads > (int)show_user_prove_thread_cap) {
    return (int)show_user_prove_thread_cap;
  }
  return threads;
}

static double show_monotonic_seconds(void) {
  struct timespec now;

  (void)clock_gettime(CLOCK_MONOTONIC, &now);
  return (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
}

static void show_print_prove_timing(const char *label, double seconds) {
  printf("[匿名签名 Prove 计时] %s：%.6f 秒\n", label, seconds);
}

static int show_profile_detail_enabled(void) {
  const char *value = getenv("REVAC_SHOW_PROFILE_DETAIL");

  return value != NULL && value[0] != '\0' && value[0] != '0';
}

static int show_round1_parallel_enabled(void) {
  const char *value = getenv("REVAC_SHOW_ROUND1_PARALLEL");

  if (value != NULL && value[0] != '\0') {
    return value[0] != '0';
  }
#ifdef _OPENMP
  /* The independent commitment products are worth splitting only on a
   * multi-core prover.  Keep an opt-out for constrained UAVs and A/B tests. */
  return omp_get_num_procs() >= 4;
#else
  return 0;
#endif
}

static int show_challenge_threads(void) {
  const char *value = getenv("REVAC_SHOW_CHALLENGE_THREADS");
  long requested = 0;

  if (value != NULL && value[0] != '\0') {
    char *end = NULL;
    requested = strtol(value, &end, 10);
    if (end != value && requested > 0 && requested <= 64) {
      return show_apply_thread_cap((int)requested);
    }
  }
#ifdef _OPENMP
  {
    const int available = omp_get_num_procs();
    return show_apply_thread_cap(available < 8 ? available : 8);
  }
#else
  return 1;
#endif
}

static int show_round1_matrix_threads(void) {
  const char *value = getenv("REVAC_SHOW_ROUND1_THREADS");

  if (value != NULL && value[0] != '\0') {
    char *end = NULL;
    long requested = strtol(value, &end, 10);
    if (end != value && requested > 0 && requested <= 64) {
      return show_apply_thread_cap((int)requested);
    }
  }
#ifdef _OPENMP
  {
    const int available = omp_get_num_procs();
    return show_apply_thread_cap(available < 12 ? available : 12);
  }
#else
  return 1;
#endif
}

static int show_round4_parallel_enabled(void) {
  const char *value = getenv("REVAC_SHOW_ROUND4_PARALLEL");

  if (value != NULL && value[0] != '\0') {
    return value[0] != '0';
  }
#ifdef _OPENMP
  return omp_get_num_procs() >= 4;
#else
  return 0;
#endif
}

static int show_round3_parallel_enabled(void) {
  const char *value = getenv("REVAC_SHOW_ROUND3_PARALLEL");

  if (value != NULL && value[0] != '\0') {
    return value[0] != '0';
  }
#ifdef _OPENMP
  return omp_get_num_procs() >= 4;
#else
  return 0;
#endif
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
* Name:        _reject_exp [static]
*
* Description: Rejection sampling given a probability threshold. Accepts
*              with probability threshold
*
* Arguments:   - double threshold: acceptance probability threshold
* 
* Returns 1 if the sample should be rejected (ie if udbl > threshold), 0 otherwise
**************************************************/
static int _reject_exp(double threshold) {
  uint64_t u = 0;
  double udbl = 0;
  randombytes((uint8_t*)&u, 6); // 48 random bits
  udbl = (double)u / (double)(1ul << 48);
  if (udbl > threshold) {
    return 1;
  }
  return 0;
}

/*************************************************
* Name:        show_user_prove_round1 [static]
*
* Description: Compute round 1 of zero-knowledge proof of valid credential
*              (anonymous credentials show)
*
* Arguments:   - show_proof_t *proof: pointer to show proof structure
*              - poly_qshow_vec_m1 *chal_1: array of polynomial vectors to host first challenge (R = R0 - R1) (initialized)
*              - uint8_t *buf: array for XOF input (allocated CHAL1_SHOW_INPUT_BYTES bytes)
*              - poly_qshow_vec_m2 s2: polynomial vector to host ABDLOP commitment randomness (initialized)
*              - poly_qshow_vec_m1 y1: polynomial vector to host mask for c.s1 (initialized)
*              - poly_qshow_vec_m2 y2: polynomial vector to host mask for c.s2 (initialized)
*              - poly_qshow_vec_256_l y3_g: polynomial vectors to host mask for R.s1 and automorphism eq. (initialized)
*              - poly_qshow_vec_d w: reusable scratch vector (initialized)
*              - const poly_qshow_vec_d tA_s1: fixed product A1.s1, computed once per proof
*              - const poly_qshow_mat_d_m1 A1: polynomial matrix for A1 (CRS)
*              - const poly_qshow_mat_d_m2 A2: polynomial matrix for A2 (CRS)
*              - const poly_qshow_mat_256l_m2 Byg: polynomial matrix for B_{y,g} (CRS)
*              - const uint8_t *randomness_seed: seed to extract randomness from (allocated SEED_BYTES bytes)
*              - const uint32_t kappa: XOF domain separator due to rejections 
**************************************************/
static void show_user_prove_round1(
    show_proof_t                 *proof,
    poly_qshow_vec_m1            chal_1[PARAM_ARP_SHOW],
    uint8_t                      buf[CHAL1_SHOW_INPUT_BYTES], 
    poly_qshow_vec_m2            s2, 
    poly_qshow_vec_m1            y1, 
    poly_qshow_vec_m2            y2, 
    poly_qshow_vec_256_l         y3_g,
    poly_qshow_vec_d             tA_low,
    poly_qshow_vec_d             w_low,
    poly_qshow_vec_d             w_high,
    poly_qshow_vec_d             w,
    const poly_qshow_vec_d       tA_s1,
    const poly_qshow_mat_d_m1    A1,
    const poly_qshow_mat_d_m2    A2,
    const poly_qshow_mat_256l_m2 Byg,
    const uint8_t                randomness_seed[SEED_BYTES],
    const uint32_t               kappa) {
  size_t i, j;
  uint32_t kpp = kappa;
  coeff_qshow tmp_coeff;
  uint8_t binomial_seed[SEED_BYTES];
  double detail_started = show_monotonic_seconds(), challenge_started;
  double sample_seconds, ta_seconds, w_seconds, tb_seconds;

  // sampling ABDLOP commitment randomness
  poly_qshow_vec_m2_binomial(s2, randomness_seed, kpp++, DOMAIN_SEPARATOR_RAND_S2_SHOW);

  // sampling Gaussian masks for c.s1, c.s2, R.s1
  poly_qshow_vec_m1_sample_gaussian_s1(y1);
  poly_qshow_vec_m2_sample_gaussian_s2(y2);
  for (i = 0; i < PARAM_ARP_DIV_N_SHOW; i++) {
    for (j = 0; j < PARAM_N_SHOW; j++) {
      tmp_coeff = SampleZ(0, PARAM_S3_SHOW);
      poly_qshow_set_coeff(y3_g->entries[i], j, tmp_coeff);
    }
  }

  // sampling uniform mask for equations with automorphisms
  for (i = PARAM_ARP_DIV_N_SHOW; i < PARAM_ARP_DIV_N_L_SHOW; i++) {
    poly_qshow_uniform_but_zero(y3_g->entries[i], randomness_seed, kpp++, DOMAIN_SEPARATOR_RAND_G_SHOW);
  }
  sample_seconds = show_monotonic_seconds() - detail_started;

  /* Compute each matrix row independently.  The previous section-level
   * parallelism left the dominant A1.y1 product serial inside one section.
   * Row parallelism keeps the exact inner-product order within every row. */
  {
    const double started = show_monotonic_seconds();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(show_round1_matrix_threads()) \
    if(show_round1_parallel_enabled())
#endif
    for (i = 0; i < PARAM_D_SHOW; i++) {
      poly_qshow row;
      poly_qshow_init(row);
      poly_qshow_vec_m2_mul_inner(row, A2->rows[i], s2);
      poly_qshow_add(proof->tA->entries[i], tA_s1->entries[i], row);
      poly_qshow_clear(row);
    }
    show_power2round_vec_d(proof->tA, tA_low, proof->tA);
    ta_seconds = show_monotonic_seconds() - started;
  }

  {
    const double started = show_monotonic_seconds();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(show_round1_matrix_threads()) \
    if(show_round1_parallel_enabled())
#endif
    for (i = 0; i < PARAM_D_SHOW; i++) {
      poly_qshow row_a1, row_a2;
      poly_qshow_init(row_a1);
      poly_qshow_init(row_a2);
      poly_qshow_vec_m1_mul_inner(row_a1, A1->rows[i], y1);
      poly_qshow_vec_m2_mul_inner(row_a2, A2->rows[i], y2);
      poly_qshow_add(w->entries[i], row_a1, row_a2);
      poly_qshow_clear(row_a2);
      poly_qshow_clear(row_a1);
    }
    show_highbits_vec_d(w_high, w);
    for (size_t row = 0; row < PARAM_D_SHOW; row++) {
      poly_qshow_mul_scalar(w_low->entries[row], w_high->entries[row],
                            PARAM_GAMMA_SHOW);
      poly_qshow_sub(w_low->entries[row], w->entries[row],
                     w_low->entries[row]);
    }
    w_seconds = show_monotonic_seconds() - started;
  }

  {
    const double started = show_monotonic_seconds();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(show_round1_matrix_threads()) \
    if(show_round1_parallel_enabled())
#endif
    for (i = 0; i < PARAM_ARP_DIV_N_L_SHOW; i++) {
      poly_qshow row;
      poly_qshow_init(row);
      poly_qshow_vec_m2_mul_inner(row, Byg->rows[i], s2);
      poly_qshow_add(proof->tB->entries[i], row, y3_g->entries[i]);
      poly_qshow_clear(row);
    }
    tb_seconds = show_monotonic_seconds() - started;
  }

  // computing first challenge  生成第一个挑战
  challenge_started = show_monotonic_seconds();
  buf[0] = 1;
  poly_qshow_vec_d_pack(buf + SHOW_CHALLENGE_BASE_BYTES, proof->tA);
  poly_qshow_vec_d_pack(buf + SHOW_CHALLENGE_BASE_BYTES + POLYQSHOW_VECD_PACKEDBYTES, w_high);
  poly_qshow_vec_256_l_pack(buf + SHOW_CHALLENGE_BASE_BYTES + 2*POLYQSHOW_VECD_PACKEDBYTES, proof->tB);
  shake256(binomial_seed, SEED_BYTES, buf, CHAL1_SHOW_INPUT_BYTES);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(show_challenge_threads()) \
    if(show_challenge_threads() > 1)
#endif
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    poly_qshow_vec_m1_binomial(chal_1[i], binomial_seed, DOMAIN_SEPARATOR_CHAL1_SHOW, i, SEED_BYTES);
  }
  if (show_profile_detail_enabled()) {
    printf("[匿名签名 Round 1 细分] 采样 %.6f 秒；tA %.6f 秒；w %.6f 秒；tB %.6f 秒；挑战 %.6f 秒\n",
           sample_seconds, ta_seconds, w_seconds, tb_seconds,
           show_monotonic_seconds() - challenge_started);
  }
}

/*************************************************
* Name:        show_user_prove_round2 [static]
*
* Description: Compute round 2 of zero-knowledge proof of valid credential
*              (anonymous credentials show) 这一轮实现了近似范围证明，证明某些值在有效范围内。
*
* Arguments:   - show_proof_t *proof: pointer to show proof structure
*              - coeff_qshow *chal_2: array of coeff_qshow to host second challenge (gamma_{i,j}) (allocated)
*              - uint8_t *buf: array for XOF input (allocated CHAL2_SHOW_INPUT_BYTES bytes)
*              - const poly_qshow_vec_m1 s1: polynomial vector, witness
*              - const poly_qshow_vec_m1 *chal_1: array of polynomial vectors, first challenge R = R0 - R1
*              - const poly_qshow_vec_256_l y3_g: polynomial vector, mask for R.s1 (ARP)
* 
* Returns 1 if round 2 passes, and 0 if it rejects
**************************************************/
static int show_user_prove_round2(
    show_proof_t               *proof,
    coeff_qshow                chal_2[PARAM_L_SHOW][PARAM_ARP_SHOW + 6],
    uint8_t                    buf[CHAL2_SHOW_INPUT_BYTES], 
    const poly_qshow_vec_m1    s1,
    const poly_qshow_vec_m1    chal_1[PARAM_ARP_SHOW],
    const poly_qshow_vec_256_l y3_g) {
  size_t i,j,k;
  int64_t tmp;
  uint64_t sq_norm_y3 = 0, sq_norm_z3 = 0;
  uint8_t challenge_seed[SEED_BYTES];
  coeff_qshow tmp_coeff;

  // computing z3 (in Z not in ring) and square norms
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    tmp = 0;
    for (j = 0; j < PARAM_M1_SHOW; j++) {
      for (k = 0; k < PARAM_N_SHOW; k++) {
          tmp += poly_qshow_get_coeff_centered(chal_1[i]->entries[j], k) * poly_qshow_get_coeff_centered(s1->entries[j], k);
      }
    }
    tmp_coeff = poly_qshow_get_coeff_centered(y3_g->entries[i / PARAM_N_SHOW], i % PARAM_N_SHOW);
    tmp += tmp_coeff;
    { uint64_t _abs = (uint64_t)(tmp_coeff >= 0 ? tmp_coeff : -tmp_coeff);
      if (_abs <= 4294967295ULL) { uint64_t _sq = _abs * _abs; CHK_UI_OVF_ADDITION(sq_norm_y3, _sq); }
      else sq_norm_y3 = UINT64_MAX; }
    { uint64_t _abs = (uint64_t)(tmp >= 0 ? tmp : -tmp);
      if (_abs <= 4294967295ULL) { uint64_t _sq = _abs * _abs; CHK_UI_OVF_ADDITION(sq_norm_z3, _sq); }
      else sq_norm_z3 = UINT64_MAX; }
    proof->z3[i] = tmp;
  }

  // rejection sampling
  // sample u3 uniform in (0,1), goto reject if u3 > exp(pi * (sq_norm_y3 - sq_norm_z3) / PARAM_S3SQ_SHOW) / PARAM_REJ3_SHOW
  if (_reject_exp(exp(M_PI *
                      ((double)sq_norm_y3 - (double)sq_norm_z3) /
                      (double)PARAM_S3SQ_SHOW) /
                  (double)PARAM_REJ3_SHOW)) {
    return 0;
  }

  // computing second challenge
  buf[0] = 2;
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    coeff_qshow_pack(buf + CHAL1_SHOW_INPUT_BYTES + i*COEFFQSHOW_PACKEDBYTES, proof->z3[i]);
  }
  shake256(challenge_seed, SEED_BYTES, buf, CHAL2_SHOW_INPUT_BYTES);
  for (i = 0; i < PARAM_L_SHOW; i++) {
    vec_qshow_uniform(chal_2[i], challenge_seed, DOMAIN_SEPARATOR_CHAL2_SHOW, i, SEED_BYTES); // writes PARAM_ARP_SHOW + 6 uniform numbers to the first argument
  }
  return 1;
}

/*************************************************
* Name:        show_user_prove_round3 [static]
*
* Description: Compute round 3 of zero-knowledge proof of valid credential
*              (anonymous credentials show)
*
* Arguments:   - show_proof_t *proof: pointer to show proof structure
*              - poly_qshow_vec_l chal_3_l: polynomial vector to host third challenge (mu_i)_{i < l} (initialized)
*              - poly_qshow_vec_k *chal_3_dk: array of polynomial vectors to host third challenge (mu_i)_{l <= i < dk+l} (initialized)
*              - uint8_t *buf: array for XOF input (allocated CHAL3_SHOW_INPUT_BYTES bytes)
*              - poly_qshow_vec_256 *sum_gamma_e_star: array polynomial vectors to host Sum_j gamma_{ij}e_j* (initialized)
*              - poly_qshow_vec_k *sum_gamma_r_star: array polynomial vectors to host Sum_j gamma_{ij}r_j* (initialized)
*              - const poly_qshow_vec_m1 s1: polynomial vector, witness
*              - const poly_qshow_vec_m1 *chal_1: array of polynomial vectors, first challenge (R = R0 - R1) 
*              - const coeff_qshow *chal_2: array of coeff_qshow, second challenge (gamma_{i,j})
*              - const poly_qshow_vec_256_l y3_g: polynomial vectors, mask for R.s1 and automorphism eq.
*              - const poly_qshow *quadratic_precomp: array of polynomials containing precomputations
**************************************************/
static void show_user_prove_round3(
    show_proof_t               *proof,         // 输出：证明结构体，存储h
    poly_qshow_vec_l           chal_3_l,                            // 输出：第三个挑战的前l个向量
    poly_qshow_vec_k           chal_3_dk[PARAM_D],                  // 输出：第三个挑战的后d×k个向量
    poly_qshow_vec_k           chal_3_link[PARAM_D],
    uint8_t                    buf[CHAL3_SHOW_INPUT_BYTES],         // 输入/输出：缓冲区
    poly_qshow_vec_256         sum_gamma_e_star[PARAM_L_SHOW],       // 输出：Σ γ_{ij}·e_j*
    poly_qshow_vec_m1          sum_gamma_r_star[PARAM_L_SHOW],       // 输出：Σ γ_{ij}·r_j*
    const poly_qshow_vec_m1    s1,                                  // 输入：见证向量
    const poly_qshow_vec_m1    chal_1[PARAM_ARP_SHOW],              // 输入：第一个挑战R
    const coeff_qshow          chal_2[PARAM_L_SHOW][PARAM_ARP_SHOW + 6],  // 输入：第二个挑战γ
    const poly_qshow_vec_256_l y3_g,                                // 输入：掩码向量
    const poly_qshow           quadratic_precomp[6])    {            // 输入：预计算的二次项 
  size_t i;
  uint8_t challenge_seed[SEED_BYTES];

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(4) if(show_round3_parallel_enabled()) schedule(static)
  #endif
  for (i = 0; i < PARAM_L_SHOW; i++) {
    poly_qshow tmp_poly;
    poly_qshow_vec_m1 tmp_vec_m1;

    poly_qshow_init(tmp_poly);
    poly_qshow_vec_m1_init(tmp_vec_m1);

    /* Since gamma is scalar,
     *   sum_j gamma_ij conjugate(r_j)
     * = conjugate(sum_j gamma_ij r_j).
     * Aggregate first and conjugate only the final m1 polynomials, instead
     * of conjugating all PARAM_ARP_SHOW*m1 challenge polynomials. */
    for (size_t col = 0; col < PARAM_ARP_SHOW; col++) {
      if (col == 0) {
        poly_qshow_vec_m1_mul_scalar(sum_gamma_r_star[i], chal_1[col],
                                     chal_2[i][col]);
      } else {
        poly_qshow_vec_m1_mul_scalar(tmp_vec_m1, chal_1[col],
                                     chal_2[i][col]);
        poly_qshow_vec_m1_add(sum_gamma_r_star[i], sum_gamma_r_star[i],
                              tmp_vec_m1);
      }
    }
    for (size_t entry = 0; entry < PARAM_M1_SHOW; entry++) {
      poly_qshow_conjugate(tmp_poly, sum_gamma_r_star[i]->entries[entry]);
      poly_qshow_set(sum_gamma_r_star[i]->entries[entry], tmp_poly);
    }

    // sum of gamma_{ij}e_j* = conjugate(tau^-1(gamma_i))
    for (size_t block = 0; block < PARAM_ARP_DIV_N_SHOW; block++) {
      poly_qshow_set_coeff(sum_gamma_e_star[i]->entries[block], 0,
                           chal_2[i][block * PARAM_N_SHOW]);
      for (size_t coeff = 1; coeff < PARAM_N_SHOW; coeff++) {
        poly_qshow_set_coeff(
            sum_gamma_e_star[i]->entries[block], coeff,
            -chal_2[i][(block + 1) * PARAM_N_SHOW - coeff]);
      }
    }

    poly_qshow_set(tmp_poly, y3_g->entries[PARAM_ARP_DIV_N_SHOW + i]); // g_i
    // sum of -gamma_{ij}z3_j
    for (size_t col = 0; col < PARAM_ARP_SHOW; col++) {
      poly_qshow_muladd_constant(tmp_poly, chal_2[i][col], -proof->z3[col]);
    }
    
    /**********************************************
    * Computing
    *   h_i = g_i - sum_j gamma_{ij}.z3_j 
    *         + sum_j gamma_{ij}.e_j*.y3 
    *         + sum_j gamma_{ij}.r_j*.s1 
    *         + gamma_{i,256}.(<v1*,v1> - B1^2)
    *         + gamma_{i,257}.(<v2*,v2> - B2^2)
    *         + gamma_{i,258}.(<v3*,v3> - B3^2)
    *         + gamma_{i,259}.(<t*,t> - w)
    *         + gamma_{i,260}.<t*,t - one>
    *         + gamma_{i,261}.<(usk|msg)*,(usk|msg) - one>
    **********************************************/
    poly_qshow_set(proof->h->entries[i], tmp_poly); // g_i - sum_j gamma_{ij}z3_j
    for (size_t block = 0; block < PARAM_ARP_DIV_N_SHOW; block++) { // <sum_e_gamma_star,y3>
      poly_qshow_mul(tmp_poly, sum_gamma_e_star[i]->entries[block],
                     y3_g->entries[block]);
      poly_qshow_add(proof->h->entries[i], proof->h->entries[i], tmp_poly);
    }
    poly_qshow_vec_m1_mul_inner(tmp_poly, sum_gamma_r_star[i], s1); // <sum_r_gamma_star,s1>
    poly_qshow_add(proof->h->entries[i], proof->h->entries[i], tmp_poly);
    for (size_t term = 0; term < 6; term++) { // precomputed quadratic terms
      poly_qshow_mul_scalar(tmp_poly, quadratic_precomp[term],
                            chal_2[i][PARAM_ARP_SHOW + term]);
      poly_qshow_add(proof->h->entries[i], proof->h->entries[i], tmp_poly);
    }

    poly_qshow_vec_m1_clear(tmp_vec_m1);
    poly_qshow_clear(tmp_poly);
  }

  // computing third challenge
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
}

/*************************************************
* Name:        show_user_prove_round4 [static]
*
* Description: Compute round 4 of zero-knowledge proof of valid credential
*              (anonymous credentials show)
*
* Arguments:   - show_proof_t *proof: pointer to show proof structure
*              - uint8_t *buf: array for XOF input (allocated CHAL4_ISS_INPUT_BYTES bytes)
*              - const poly_qshow_vec_m1 s1: polynomial vector, witness
*              - const poly_qshow_vec_m2 s2: polynomial vector, ABDLOP commitment randomness
*              - const poly_qshow_mat_256l_m2 Byg: polynomial matrix for B_{y,g} (CRS)
*              - const poly_qshow_vec_m2 b: polynomial vector for b (CRS)
*              - const poly_qshow_mat_k_k *A_embed: array of polynomial matrices, subring embedding of q1.A'
*              - const poly_qshow_mat_k_k *B_embed: array of polynomial matrices, subring embedding of q1.B
*              - const poly_qshow_mat_k_k *A3_embed: array of polynomial matrices, subring embedding of q1.A3
*              - const poly_qshow_mat_k_k *Ds_embed: array of polynomial matrices, subring embedding of q1.Ds
*              - const poly_qshow_mat_k_k *D_embed: array of polynomial matrices, subring embedding of q1.D
*              - const poly_qshow_vec_256 *sum_gamma_e_star: array polynomial vectors, Sum_j gamma_{ij}e_j*
*              - const poly_qshow_vec_k *sum_gamma_r_star: array polynomial vectors, Sum_j gamma_{ij}r_j*
*              - const poly_qshow_vec_m1 y1: polynomial vector, mask for c.s1
*              - const poly_qshow_vec_m2 y2: polynomial vector, mask for c.s2
*              - const coeff_qshow *chal_2: array of coeff_qshow, second challenge (gamma_{i,j})
*              - const poly_qshow_vec_l chal_3_l: polynomial vector, third challenge (mu_i)_{i < l} 
*              - const poly_qshow_vec_k *chal_3_dk: array of polynomial vectors, third challenge (mu_i)_{l <= i < dk+l}
*              - const poly_qshow one: polynomial with all ones
**************************************************/
static void show_user_prove_round4(
    show_proof_t                *proof,
    uint8_t                      buf[CHAL4_SHOW_INPUT_BYTES],
    const poly_qshow_vec_m1      s1,
    const poly_qshow_vec_m2      s2, 
    const poly_qshow_mat_256l_m2 Byg,
    const poly_qshow_vec_m2      b,
    const poly_qshow_mat_k_k     A_embed[PARAM_D][PARAM_D], 
    const poly_qshow_mat_k_k     B_embed[PARAM_D][PARAM_D*PARAM_K], 
    const poly_qshow_mat_k_k A3_embed[PARAM_D],
    const poly_qshow_mat_k_k     Ds_embed[PARAM_D][2*PARAM_D], 
    const poly_qshow_mat_k_k     D_embed[PARAM_D][PARAM_M], 
    const poly_qshow_mat_k_k     Dx_embed[PARAM_D],
    const poly_qshow_mat_k_k     A_link_embed[PARAM_D][PARAM_D],
    const poly_qshow_mat_k_k     D_link_embed[PARAM_D],
    const poly_qshow_vec_256     sum_gamma_e_star[PARAM_L_SHOW],
    const poly_qshow_vec_m1      sum_gamma_r_star[PARAM_L_SHOW],
    const poly_qshow_vec_m1      y1,
    const poly_qshow_vec_m2      y2,
    const coeff_qshow            chal_2[PARAM_L_SHOW][PARAM_ARP_SHOW + 6],
    const poly_qshow_vec_l       chal_3_l,
    const poly_qshow_vec_k       chal_3_dk[PARAM_D],
    const poly_qshow_vec_k       chal_3_link[PARAM_D],
    const poly_qshow             one) {
  size_t i,j,k;
  uint32_t kappa_c;
  int64_t bexpi;
  uint8_t challenge_seed[SEED_BYTES];
  poly_qshow tmp_poly, e0, e1, y1i_star, y1s_y1, y1s_s1, y1s_one, t0, sum_mu_gamma[6];
  poly_qshow signature_rows[PARAM_D * PARAM_K_SHOW];
  poly_qshow link_rows[PARAM_D * PARAM_K_SHOW];
  poly_qshow_vec_256_l tmp_vec_256_l;
  poly_qshow_vec_k tmp_vec_k, tmp_vec_k2, Gy1_v2[PARAM_D], Gs1_v2[PARAM_D];
  poly_qshow_mat_k_k ytag_mat, stag_mat;
  
  // init matrices, vectors and polynomials
  poly_qshow_init(tmp_poly);
  poly_qshow_init(e0);
  poly_qshow_init(e1);
  poly_qshow_init(y1i_star);
  poly_qshow_init(y1s_y1);
  poly_qshow_init(y1s_s1);
  poly_qshow_init(y1s_one);
  poly_qshow_init(t0);
  poly_qshow_vec_256_l_init(tmp_vec_256_l);
  poly_qshow_vec_k_init(tmp_vec_k);
  poly_qshow_vec_k_init(tmp_vec_k2);
  poly_qshow_mat_k_k_init(ytag_mat);
  poly_qshow_mat_k_k_init(stag_mat);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_init(Gy1_v2[i]);
    poly_qshow_vec_k_zero(Gy1_v2[i]);
    poly_qshow_vec_k_init(Gs1_v2[i]);
    poly_qshow_vec_k_zero(Gs1_v2[i]);
  }
  for (i = 0; i < 6; i++) {
    poly_qshow_init(sum_mu_gamma[i]);
    poly_qshow_zero(sum_mu_gamma[i]);
  }
  for (i = 0; i < PARAM_D * PARAM_K_SHOW; i++) {
    poly_qshow_init(signature_rows[i]);
    poly_qshow_init(link_rows[i]);
  }

  // compute sum_i mu_i gamma_{i,256 + j} for j in [6]
  for (i = 0; i < PARAM_L_SHOW; i++) {
    for (j = 0; j < 6; j++) {
      poly_qshow_mul_scalar(tmp_poly, chal_3_l->entries[i], chal_2[i][PARAM_ARP_SHOW + j]);
      poly_qshow_add(sum_mu_gamma[j], sum_mu_gamma[j], tmp_poly);
    }
  }
  poly_qshow_add(sum_mu_gamma[3], sum_mu_gamma[3], sum_mu_gamma[4]); // sum_mu_gamma[4] will be used later

  // quadratic terms
  // <y1_v1"*,y1_v1">, <y1_v1"*,s1_v1">
  poly_qshow_zero(y1s_y1);
  poly_qshow_zero(y1s_s1);
  for (i = IDX_V1_SHOW; i < IDX_V2_SHOW; i++) {
    poly_qshow_conjugate(y1i_star, y1->entries[i]);
    poly_qshow_mul(tmp_poly, y1i_star, y1->entries[i]);
    poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, s1->entries[i]);
    poly_qshow_add(y1s_s1, y1s_s1, tmp_poly);
  }
  poly_qshow_mul(e0, sum_mu_gamma[0], y1s_y1); 
  poly_qshow_conjugate(tmp_poly, y1s_s1);
  poly_qshow_add(y1s_s1, y1s_s1, tmp_poly); // <y1*,s1> + <y1*,s1>*
  poly_qshow_mul(e1, sum_mu_gamma[0], y1s_s1);


  // <y1_v2"*,y1_v2">, <y1_v2"*,s1_v2">
  poly_qshow_zero(y1s_y1);
  poly_qshow_zero(y1s_s1);
  for (i = IDX_V2_SHOW; i < IDX_V3_SHOW; i++) {
    poly_qshow_conjugate(y1i_star, y1->entries[i]);
    poly_qshow_mul(tmp_poly, y1i_star, y1->entries[i]);
    poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, s1->entries[i]);
    poly_qshow_add(y1s_s1, y1s_s1, tmp_poly);
  }
  poly_qshow_mul(tmp_poly, sum_mu_gamma[1], y1s_y1); 
  poly_qshow_add(e0, e0, tmp_poly); 
  poly_qshow_conjugate(tmp_poly, y1s_s1);
  poly_qshow_add(tmp_poly, tmp_poly, y1s_s1);
  poly_qshow_mul(tmp_poly, sum_mu_gamma[1], tmp_poly); 
  poly_qshow_add(e1, e1, tmp_poly);

  // <y1_v3"*,y1_v3">, <y1_v3"*,s1_v3">
  poly_qshow_zero(y1s_y1);
  poly_qshow_zero(y1s_s1);
  for (i = IDX_V3_SHOW; i < IDX_TAG_SHOW; i++) {
    poly_qshow_conjugate(y1i_star, y1->entries[i]);
    poly_qshow_mul(tmp_poly, y1i_star, y1->entries[i]);
    poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, s1->entries[i]);
    poly_qshow_add(y1s_s1, y1s_s1, tmp_poly);
  }
  poly_qshow_mul(tmp_poly, sum_mu_gamma[2], y1s_y1); 
  poly_qshow_add(e0, e0, tmp_poly); 
  poly_qshow_conjugate(tmp_poly, y1s_s1);
  poly_qshow_add(tmp_poly, tmp_poly, y1s_s1);
  poly_qshow_mul(tmp_poly, sum_mu_gamma[2], tmp_poly); 
  poly_qshow_add(e1, e1, tmp_poly);

  // <y1_t'*,y1_t'>, <y1_t'*,s1_t'>, <y1_t'*,one>
  poly_qshow_zero(y1s_y1);
  poly_qshow_zero(y1s_s1);
  poly_qshow_zero(y1s_one);
  for (i = IDX_TAG_SHOW; i < IDX_USK_SHOW; i++) {
    poly_qshow_conjugate(y1i_star, y1->entries[i]);
    poly_qshow_mul(tmp_poly, y1i_star, y1->entries[i]);
    poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, s1->entries[i]);
    poly_qshow_add(y1s_s1, y1s_s1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, one);
    poly_qshow_add(y1s_one, y1s_one, tmp_poly);
  }
  poly_qshow_mul(tmp_poly, sum_mu_gamma[3], y1s_y1); 
  poly_qshow_add(e0, e0, tmp_poly); 
  poly_qshow_conjugate(tmp_poly, y1s_s1);
  poly_qshow_add(tmp_poly, tmp_poly, y1s_s1);
  poly_qshow_mul(tmp_poly, sum_mu_gamma[3], tmp_poly); 
  poly_qshow_add(e1, e1, tmp_poly);
  poly_qshow_mul(tmp_poly, sum_mu_gamma[4], y1s_one); 
  poly_qshow_sub(e1, e1, tmp_poly); // substract linear part containing sum_mu_gamma[4].<y1*, one>

  // <y1_m'*,y1_m'>, <y1_m'*,s1_m'>, <y1_m'*,one> 
  poly_qshow_zero(y1s_y1);
  poly_qshow_zero(y1s_s1);
  poly_qshow_zero(y1s_one);
  for (i = IDX_USK_SHOW; i < IDX_SHOW_WITNESS_END; i++) {// same for usk, all message attributes, and handle x
    poly_qshow_conjugate(y1i_star, y1->entries[i]);
    poly_qshow_mul(tmp_poly, y1i_star, y1->entries[i]);
    poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, s1->entries[i]);
    poly_qshow_add(y1s_s1, y1s_s1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, one);
    poly_qshow_add(y1s_one, y1s_one, tmp_poly);
  }
  poly_qshow_mul(tmp_poly, sum_mu_gamma[5], y1s_y1); 
  poly_qshow_add(e0, e0, tmp_poly); 
  poly_qshow_conjugate(tmp_poly, y1s_s1);
  poly_qshow_add(tmp_poly, tmp_poly, y1s_s1);
  poly_qshow_sub(tmp_poly, tmp_poly, y1s_one); // contains <y1_m*, s1_m> + <y1_m*, s1_m>* - <y1_m*, one>
  poly_qshow_mul(tmp_poly, sum_mu_gamma[5], tmp_poly); 
  poly_qshow_add(e1, e1, tmp_poly);

  // Linear part depending on sum_gamma_r_star and sum_gamma_e_star
  poly_qshow_mat_256l_m2_mul_vec_m2(tmp_vec_256_l, Byg, y2);
  for (i = 0; i < PARAM_L_SHOW; i++) {
    // y1s_y1 can be used as a temp variable
    poly_qshow_vec_m1_mul_inner(y1s_y1, sum_gamma_r_star[i], y1);
    poly_qshow_sub(y1s_y1, y1s_y1, tmp_vec_256_l->entries[PARAM_ARP_DIV_N_SHOW + i]);
    for (j = 0; j < PARAM_ARP_DIV_N_SHOW; j++) {
      poly_qshow_mul(tmp_poly, sum_gamma_e_star[i]->entries[j], tmp_vec_256_l->entries[j]);
      poly_qshow_sub(y1s_y1, y1s_y1, tmp_poly);
    }
    poly_qshow_mul(y1s_y1, y1s_y1, chal_3_l->entries[i]);
    poly_qshow_add(e1, e1, y1s_y1);
  }
  
  // quadratic part depending on chal_3_quad_matrix (TSampler: KH entries, qL*bH^i scaling)
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_zero(Gy1_v2[i]->entries[j]);
      poly_qshow_zero(Gs1_v2[i]->entries[j]);
      bexpi = (int64_t)PARAM_Q1_SHOW * PARAM_QL;
      for (k = 0; k < PARAM_KH; k++) {
        poly_qshow_mul_scalar(tmp_poly, y1->entries[IDX_V2_SHOW + k*PARAM_D*PARAM_K_SHOW + i*PARAM_K_SHOW + j], bexpi);
        poly_qshow_add(Gy1_v2[i]->entries[j], Gy1_v2[i]->entries[j], tmp_poly);
        poly_qshow_mul_scalar(tmp_poly, s1->entries[IDX_V2_SHOW + k*PARAM_D*PARAM_K_SHOW + i*PARAM_K_SHOW + j], bexpi);
        poly_qshow_add(Gs1_v2[i]->entries[j], Gs1_v2[i]->entries[j], tmp_poly);
        bexpi *= PARAM_BH;
      }
    }
  }
  for (i = 0; i < PARAM_K_SHOW; i++) {
    poly_qshow_set(tmp_vec_k->entries[i], y1->entries[IDX_TAG_SHOW + i]);
    poly_qshow_set(tmp_vec_k2->entries[i], s1->entries[IDX_TAG_SHOW + i]);
  }
  show_mat_k_k_from_vec_k(ytag_mat, tmp_vec_k);
  show_mat_k_k_from_vec_k(stag_mat, tmp_vec_k2);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_mat_k_k_mul_vec_k(tmp_vec_k, ytag_mat, Gy1_v2[i]);
    poly_qshow_mat_k_k_mul_vec_k(tmp_vec_k2, ytag_mat, Gs1_v2[i]);
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, tmp_vec_k->entries[j], chal_3_dk[i]->entries[j]);
      poly_qshow_add(e0, e0, tmp_poly);
      poly_qshow_mul(tmp_poly, tmp_vec_k2->entries[j], chal_3_dk[i]->entries[j]);
      poly_qshow_add(e1, e1, tmp_poly);
    }
    poly_qshow_mat_k_k_mul_vec_k(tmp_vec_k, stag_mat, Gy1_v2[i]);
    poly_qshow_mat_k_k_mul_vec_k(tmp_vec_k2, stag_mat, Gs1_v2[i]);
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, tmp_vec_k->entries[j], chal_3_dk[i]->entries[j]);
      poly_qshow_add(e1, e1, tmp_poly);
      /*
       * The s1*s1 gadget term is part of the proved quadratic relation
       * Q(s1) + L(s1) - u = 0.  It therefore cancels with the corresponding
       * c^2 linear term when the verifier reconstructs t0; it is not a
       * separate proof element.
       */
    }
  }

  // linear part depending on embedded matrices A', B, A3, Ds, D
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(4) if(show_round4_parallel_enabled()) schedule(static) private(j)
  #endif
  for (i = 0; i < PARAM_D*PARAM_K_SHOW; i++) {
    poly_qshow row_tmp;
    const size_t row_quot = i / PARAM_K_SHOW;
    const size_t row_rem = i % PARAM_K_SHOW;

    poly_qshow_init(row_tmp);
    // y1s_y1 used as temp variable to host ([A|-B|A3|-Ds|-D].[y1_v1|y1_v2|y1_v3|y1_s|y1_m])_i
    poly_qshow_mul_scalar(signature_rows[i], y1->entries[IDX_V1_SHOW + i], PARAM_Q1_SHOW);
    for (j = 0; j < PARAM_D*PARAM_K_SHOW; j++) {
      poly_qshow_mul(row_tmp, A_embed[row_quot][j / PARAM_K_SHOW]->rows[row_rem]->entries[j % PARAM_K_SHOW], y1->entries[IDX_V12_SHOW + j]);
      poly_qshow_add(signature_rows[i], signature_rows[i], row_tmp);
    }
    for (j = 0; j < PARAM_D*PARAM_KH*PARAM_K_SHOW; j++) {
      poly_qshow_mul(row_tmp, B_embed[row_quot][j / PARAM_K_SHOW]->rows[row_rem]->entries[j % PARAM_K_SHOW], y1->entries[IDX_V2_SHOW + j]);
      poly_qshow_sub(signature_rows[i], signature_rows[i], row_tmp); // substraction
    }
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_mul(row_tmp, A3_embed[row_quot]->rows[row_rem]->entries[j], y1->entries[IDX_V3_SHOW + j]);
      poly_qshow_add(signature_rows[i], signature_rows[i], row_tmp);
    }
    for (j = 0; j < 2*PARAM_D*PARAM_K_SHOW; j++) {
      poly_qshow_mul(row_tmp, Ds_embed[row_quot][j / PARAM_K_SHOW]->rows[row_rem]->entries[j % PARAM_K_SHOW], y1->entries[IDX_USK_SHOW + j]);
      poly_qshow_sub(signature_rows[i], signature_rows[i], row_tmp); // substraction
    }
    for (j = 0; j < PARAM_M*PARAM_K_SHOW; j++) {
      poly_qshow_mul(row_tmp, D_embed[row_quot][j / PARAM_K_SHOW]->rows[row_rem]->entries[j % PARAM_K_SHOW], y1->entries[IDX_M_SHOW + j]);
      poly_qshow_sub(signature_rows[i], signature_rows[i], row_tmp); // substraction
    }
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_mul(row_tmp, Dx_embed[row_quot]->rows[row_rem]->entries[j], y1->entries[IDX_X_SHOW + j]);
      poly_qshow_sub(signature_rows[i], signature_rows[i], row_tmp);
    }
    poly_qshow_mul(signature_rows[i], signature_rows[i], chal_3_dk[row_quot]->entries[row_rem]);
    poly_qshow_clear(row_tmp);
  }
  for (i = 0; i < PARAM_D * PARAM_K_SHOW; i++) {
    poly_qshow_add(e1, e1, signature_rows[i]);
  }

  /* Independently batch the same-handle commitment equations
   *
   *   r0 + A_link r1 + D_link x - C_link = 0.
   *
   * Round 4 only needs the masked linear part L(y1).  The verifier inserts
   * the public -c*C_link term when reconstructing the transcript. */
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(4) if(show_round4_parallel_enabled()) schedule(static) private(j)
  #endif
  for (i = 0; i < PARAM_D * PARAM_K_SHOW; i++) {
    poly_qshow row_tmp;
    const size_t row_quot = i / PARAM_K_SHOW;
    const size_t row_rem = i % PARAM_K_SHOW;

    poly_qshow_init(row_tmp);
    poly_qshow_mul_scalar(link_rows[i],
                          y1->entries[IDX_LINK_R0_SHOW + i],
                          PARAM_Q1_SHOW);
    for (j = 0; j < PARAM_D * PARAM_K_SHOW; j++) {
      poly_qshow_mul(
          row_tmp,
          A_link_embed[row_quot][j / PARAM_K_SHOW]
              ->rows[row_rem]->entries[j % PARAM_K_SHOW],
          y1->entries[IDX_LINK_R1_SHOW + j]);
      poly_qshow_add(link_rows[i], link_rows[i], row_tmp);
    }
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_mul(row_tmp,
                     D_link_embed[row_quot]
                         ->rows[row_rem]->entries[j],
                     y1->entries[IDX_X_SHOW + j]);
      poly_qshow_add(link_rows[i], link_rows[i], row_tmp);
    }
    poly_qshow_mul(link_rows[i], link_rows[i],
                   chal_3_link[row_quot]->entries[row_rem]);
    poly_qshow_clear(row_tmp);
  }
  for (i = 0; i < PARAM_D * PARAM_K_SHOW; i++) {
    poly_qshow_add(e1, e1, link_rows[i]);
  }

  // committing to garbage terms
    poly_qshow_vec_m2_mul_inner(t0, b, y2);
    poly_qshow_add(t0, t0, e0);// 计算 t0 = <b, y2> + e0
    poly_qshow_vec_m2_mul_inner(tmp_poly, b, s2);
  poly_qshow_add(proof->t1, tmp_poly, e1);// 计算 t1 = <b, s2> + e1

  // computing fourth challenge
  kappa_c = 0;
  buf[0] = 4;
  poly_qshow_pack(buf + CHAL3_SHOW_INPUT_BYTES, t0);
  poly_qshow_pack(buf + CHAL3_SHOW_INPUT_BYTES + POLYQSHOW_PACKEDBYTES, proof->t1);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL4_SHOW_INPUT_BYTES);
  
  do {
    poly_qshow_sample_challenge(proof->c, challenge_seed, DOMAIN_SEPARATOR_CHAL4_SHOW, kappa_c++, SEED_BYTES);
  } while (challenge_size_show(proof->c) > PARAM_ETA_SHOW);
  proof->ctr_c = kappa_c - 1;

  
  // clean up matrices, vectors and polynomials
  poly_qshow_clear(tmp_poly);
  poly_qshow_clear(e0);
  poly_qshow_clear(e1);
  poly_qshow_clear(y1i_star);
  poly_qshow_clear(y1s_y1);
  poly_qshow_clear(y1s_s1);
  poly_qshow_clear(y1s_one);
  poly_qshow_clear(t0);
  poly_qshow_vec_256_l_clear(tmp_vec_256_l);
  poly_qshow_mat_k_k_clear(stag_mat);
  poly_qshow_mat_k_k_clear(ytag_mat);
  poly_qshow_vec_k_clear(tmp_vec_k2);
  poly_qshow_vec_k_clear(tmp_vec_k);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_clear(Gy1_v2[i]);
    poly_qshow_vec_k_clear(Gs1_v2[i]);
  }
  for (i = 0; i < 6; i++) {
    poly_qshow_clear(sum_mu_gamma[i]);
  }
  for (i = 0; i < PARAM_D * PARAM_K_SHOW; i++) {
    poly_qshow_clear(link_rows[i]);
    poly_qshow_clear(signature_rows[i]);
  }
}

/*************************************************
* Name:        show_user_prove_round5 [static]
*
* Description: Compute round 5 of zero-knowledge proof of valid credential
*              (anonymous credentials show)
*
* Arguments:   - show_proof_t *proof: pointer to showuance proof structure
*              - const poly_qshow_vec_m1 s1: polynomial vector, witness
*              - const poly_qshow_vec_m2 s2: polynomial vector, ABDLOP commitment randomness
*              - const poly_qshow_vec_m1 y1: polynomial vector, mask for c.s1
*              - const poly_qshow_vec_m2 y2: polynomial vector, mask for c.s2
* 
* Returns 1 if round 5 passes, and 0 if it rejects
**************************************************/
static int show_user_prove_round5(
    show_proof_t            *proof,
    const poly_qshow_vec_m1 s1,
    const poly_qshow_vec_m2 s2,
    const poly_qshow_vec_m1 y1,
    const poly_qshow_vec_m2 y2,
    const poly_qshow_vec_d  tA_low,
    const poly_qshow_vec_d  w_low,
    const poly_qshow_vec_d  w_high) {
  uint128 sq_norm_y1, sq_norm_z1, sq_norm_y2, sq_norm_z2;
  poly_qshow tmp;
  poly_qshow_vec_d adjusted_tail, hint_base;

  poly_qshow_init(tmp);
  poly_qshow_vec_d_init(adjusted_tail);
  poly_qshow_vec_d_init(hint_base);

  // computing z1 = y1 + c.s1, z2 = y2 + c.s2 (and square norms)
  poly_qshow_vec_m1_mul_poly_qshow(proof->z1, s1, proof->c); 
  poly_qshow_vec_m1_add(proof->z1, proof->z1, y1);
  sq_norm_y1 = poly_qshow_vec_m1_norm2(y1);
  sq_norm_z1 = poly_qshow_vec_m1_norm2(proof->z1);

  poly_qshow_vec_m2_mul_poly_qshow(proof->z2, s2, proof->c); 
  poly_qshow_vec_m2_add(proof->z2, proof->z2, y2);
  sq_norm_y2 = poly_qshow_vec_m2_norm2(y2);
  sq_norm_z2 = poly_qshow_vec_m2_norm2(proof->z2);

  // rejection sampling
  // sample u1 uniform in (0,1), goto reject if u1 > exp(pi * (sq_norm_y1 - sq_norm_z1) / PARAM_S1SQ_SHOW) / PARAM_REJ1_SHOW
  // sample u2 uniform in (0,1), goto reject if u2 > exp(pi * (sq_norm_y2 - sq_norm_z2) / PARAM_S2SQ_SHOW) / PARAM_REJ2_SHOW
  if (_reject_exp(
    exp(M_PI * ((double)sq_norm_y1 - (double)sq_norm_z1) /
        (double)PARAM_S1SQ_SHOW) / (double)PARAM_REJ1_SHOW
  ) || _reject_exp(
    exp(M_PI * ((double)sq_norm_y2 - (double)sq_norm_z2) /
        (double)PARAM_S2SQ_SHOW) / (double)PARAM_REJ2_SHOW
  )) {
    goto reject;
  }

  /*
   * Compressed ABDLOP response.  The final d components are replaced by a
   * generalized hint:
   *   z2,2 <- z2,2 - c*tA,2 - w0.
   * Only z2,1 is serialized.  The verifier recovers the norm-equivalent
   * negative of z2,2 from A1*z1 + A2'*z2,1 - c*2^D*tA,1 - gamma*w1.
   */
  for (size_t i = 0; i < PARAM_D_SHOW; i++) {
    poly_qshow_mul(tmp, proof->c, tA_low->entries[i]);
    poly_qshow_sub(proof->z2->entries[PARAM_M2_COMMIT_SHOW + i],
                   proof->z2->entries[PARAM_M2_COMMIT_SHOW + i], tmp);
    poly_qshow_sub(proof->z2->entries[PARAM_M2_COMMIT_SHOW + i],
                   proof->z2->entries[PARAM_M2_COMMIT_SHOW + i],
                   w_low->entries[i]);
    poly_qshow_set(adjusted_tail->entries[i],
                   proof->z2->entries[PARAM_M2_COMMIT_SHOW + i]);
    poly_qshow_mul_scalar(hint_base->entries[i], w_high->entries[i],
                          PARAM_GAMMA_SHOW);
    poly_qshow_sub(hint_base->entries[i], hint_base->entries[i],
                   adjusted_tail->entries[i]);
  }
  sq_norm_z2 = poly_qshow_vec_m2_norm2(proof->z2);
  if (sq_norm_z2 >
      ((uint128)PARAM_B2SQ_SHOW_LOW64 +
       ((uint128)PARAM_B2SQ_SHOW_HIGH64 << 64))) {
    goto reject;
  }
  show_make_g_hints_vec_d(proof->hints, adjusted_tail, hint_base);

  poly_qshow_vec_d_clear(hint_base);
  poly_qshow_vec_d_clear(adjusted_tail);
  poly_qshow_clear(tmp);
  return 1;

reject:
  poly_qshow_vec_d_clear(hint_base);
  poly_qshow_vec_d_clear(adjusted_tail);
  poly_qshow_clear(tmp);
  return 0;
}

/*************************************************
* Name:        show_user_prove
*
* Description: Compute zero-knowledge proof of valid credential
*              (anonymous credentials show)
*
* Arguments:   - show_proof *proof: pointer to show proof structure
*              - const poly_qshow_mat_k_k *A_embed: array of polynomial matrices, subring embedding of q1.A'
*              - const poly_qshow_mat_k_k *B_embed: array of polynomial matrices, subring embedding of q1.B
*              - const poly_qshow_mat_k_k *A3_embed: array of polynomial matrices, subring embedding of q1.A3
*              - const poly_qshow_mat_k_k *Ds_embed: array of polynomial matrices, subring embedding of q1.Ds
*              - const poly_qshow_mat_k_k *D_embed: array of polynomial matrices, subring embedding of q1.D
*              - const poly_qshow_vec_m1 s1: polynomial vector, witness
*                   s1 = [theta(v1)|a1 | theta(v2)|a2 | theta(v3)|a3 | theta(tag) | theta(usk) | theta(msg)]
*              - const uint8_t *crs_seed: pointer to byte array containing the CRS seed (allocated SEED_BYTES bytes)
*              - const uint8_t *seed: pointer to byte array containing the seed 
*                   for public parameters (allocated SEED_BYTES bytes)
**************************************************/
void show_user_prove_handle(
    show_proof_t             *proof, 
    const poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D*PARAM_K], 
    const poly_qshow_mat_k_k A3_embed[PARAM_D],
    const poly_qshow_mat_k_k Ds_embed[PARAM_D][2*PARAM_D], 
    const poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], 
    const poly_qshow_mat_k_k Dx_embed[PARAM_D],
    const poly_qshow_mat_k_k A_link_embed[PARAM_D][PARAM_D],
    const poly_qshow_mat_k_k D_link_embed[PARAM_D],
    const poly_qshow_vec_k   C_link_embed[PARAM_D],
    const poly_qshow_vec_m1  s1, 
    const uint8_t            crs_seed[CRS_SEED_BYTES], 
    const uint8_t            seed[SEED_BYTES]) {
  size_t i;
  uint8_t randomness_seed[SEED_BYTES];
  uint8_t buf[CHAL4_SHOW_INPUT_BYTES] = {0};
  uint32_t kappa;
  poly_qshow s1i_star, tmp_poly, one, quadratic_precomp[6];
  poly_qshow_mat_d_m1 A1;
  poly_qshow_mat_d_m2 A2;
  poly_qshow_mat_256l_m2 Byg;
  poly_qshow_vec_m1 y1;
  poly_qshow_vec_m1 sum_gamma_r_star[PARAM_L_SHOW];
  poly_qshow_vec_256 sum_gamma_e_star[PARAM_L_SHOW];
  poly_qshow_vec_m2 s2, y2, b;
  poly_qshow_vec_256_l y3_g;
  poly_qshow_vec_d tA_low, w_low, w_high;
  poly_qshow_vec_d tA_s1, round1_w;
  double total_started = show_monotonic_seconds();
  double step_started = show_monotonic_seconds();
  double init_seconds;
  double crs_seconds;
  double precompute_seconds;
  double round1_seconds = 0.0;
  double round2_seconds = 0.0;
  double round3_seconds = 0.0;
  double round4_seconds = 0.0;
  double round5_seconds = 0.0;
  size_t attempt_count = 0;
  size_t round2_rejections = 0;
  size_t round5_rejections = 0;

  // challenges
  poly_qshow_vec_m1 chal_1[PARAM_ARP_SHOW];
  coeff_qshow chal_2[PARAM_L_SHOW][PARAM_ARP_SHOW + 6];
  poly_qshow_vec_l chal_3_l;
  poly_qshow_vec_k chal_3_dk[PARAM_D];
  poly_qshow_vec_k chal_3_link[PARAM_D];

  // init
  // init polynomials
  poly_qshow_init(tmp_poly);
  poly_qshow_init(s1i_star);
  poly_qshow_init(one);
  for (i = 0; i < 6; i++) {
    poly_qshow_init(quadratic_precomp[i]);
    poly_qshow_zero(quadratic_precomp[i]);
  }
  // init vectors and matrices
  poly_qshow_mat_d_m1_init(A1);
  poly_qshow_vec_m1_init(y1);
  for (i = 0; i < PARAM_L_SHOW; i++) {
    poly_qshow_vec_m1_init(sum_gamma_r_star[i]);
    poly_qshow_vec_256_init(sum_gamma_e_star[i]);    
  }
  poly_qshow_mat_d_m2_init(A2);
  poly_qshow_mat_256l_m2_init(Byg);
  poly_qshow_vec_m2_init(s2);
  poly_qshow_vec_m2_init(y2);
  poly_qshow_vec_m2_init(b);
  poly_qshow_vec_256_l_init(y3_g);
  poly_qshow_vec_d_init(tA_low);
  poly_qshow_vec_d_init(w_low);
  poly_qshow_vec_d_init(w_high);
  poly_qshow_vec_d_init(tA_s1);
  poly_qshow_vec_d_init(round1_w);
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    poly_qshow_vec_m1_init(chal_1[i]);
  }
  poly_qshow_vec_l_init(chal_3_l);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_init(chal_3_dk[i]);
    poly_qshow_vec_k_init(chal_3_link[i]);
  }
  (void)C_link_embed;
  init_seconds = show_monotonic_seconds() - step_started;

  // generate random secret seed
  step_started = show_monotonic_seconds();
  randombytes(randomness_seed, SEED_BYTES);

  // expanding CRS
  poly_qshow_mat_d_m1_uniform(A1, crs_seed, DOMAIN_SEPARATOR_A1_SHOW);
  poly_qshow_mat_d_m2_uniform(A2, crs_seed, DOMAIN_SEPARATOR_A2_SHOW);
  poly_qshow_mat_256l_m2_uniform(Byg, crs_seed, DOMAIN_SEPARATOR_BYG_SHOW);
  poly_qshow_vec_m2_uniform(b, crs_seed, DOMAIN_SEPARATOR_B_SHOW);

  // byte-packing the statement
  // public parameters are derived from seed
  // the signer public key does not have to be hash (single-issuer setting)
  for (i = 0; i < CRS_SEED_BYTES; i++) {
    buf[1 + i] = crs_seed[i];
  }
  for (i = 0; i < SEED_BYTES; i++) {
    buf[1 + CRS_SEED_BYTES + i] = seed[i];
  }
  crs_seconds = show_monotonic_seconds() - step_started;

  // precomputations for quadratic equations
  step_started = show_monotonic_seconds();
  for (i = 0; i < PARAM_N_SHOW; i++) {
    poly_qshow_set_coeff(one, i, 1);
  }
  // <v1"*,v1"> - B_1'^2
  for (i = IDX_V1_SHOW; i < IDX_V2_SHOW; i++) {
    poly_qshow_conjugate(s1i_star, s1->entries[i]);
    poly_qshow_mul(tmp_poly, s1i_star, s1->entries[i]);
    poly_qshow_add(quadratic_precomp[0], quadratic_precomp[0], tmp_poly);
  }
  poly_qshow_muladd_constant(quadratic_precomp[0], -PARAM_B1SQ, 1); 

  // <v2"*,v2"> - B_2^2
  for (i = IDX_V2_SHOW; i < IDX_V3_SHOW; i++) {
    poly_qshow_conjugate(s1i_star, s1->entries[i]);
    poly_qshow_mul(tmp_poly, s1i_star, s1->entries[i]);
    poly_qshow_add(quadratic_precomp[1], quadratic_precomp[1], tmp_poly);
  }
  poly_qshow_muladd_constant(quadratic_precomp[1], -PARAM_B2SQ, 1); 

  // <v3"*,v3"> - B_3^2
  for (i = IDX_V3_SHOW; i < IDX_TAG_SHOW; i++) {
    poly_qshow_conjugate(s1i_star, s1->entries[i]);
    poly_qshow_mul(tmp_poly, s1i_star, s1->entries[i]);
    poly_qshow_add(quadratic_precomp[2], quadratic_precomp[2], tmp_poly);
  }
  poly_qshow_muladd_constant(quadratic_precomp[2], -PARAM_B3SQ, 1); 

  // <tag'*,tag'> - w and <tag'*,tag' - one>
  for (i = IDX_TAG_SHOW; i < IDX_USK_SHOW; i++) {
    poly_qshow_conjugate(s1i_star, s1->entries[i]);
    poly_qshow_mul(tmp_poly, s1i_star, s1->entries[i]);
    poly_qshow_add(quadratic_precomp[3], quadratic_precomp[3], tmp_poly); // <tag'*,tag'>
    poly_qshow_sub(tmp_poly, s1->entries[i], one);
    poly_qshow_mul(tmp_poly, s1i_star, tmp_poly);
    poly_qshow_add(quadratic_precomp[4], quadratic_precomp[4], tmp_poly); // <tag'*,tag' - one>
  }
  poly_qshow_muladd_constant(quadratic_precomp[3], -PARAM_W, 1); 

  // <m'*,m' - one> (includes usk and message)
  for (i = IDX_USK_SHOW; i < IDX_SHOW_WITNESS_END; i++) {// same for usk, all message attributes, and handle x
    poly_qshow_conjugate(s1i_star, s1->entries[i]);
    poly_qshow_sub(tmp_poly, s1->entries[i], one);
    poly_qshow_mul(tmp_poly, s1i_star, tmp_poly);
    poly_qshow_add(quadratic_precomp[5], quadratic_precomp[5], tmp_poly);
  }
  precompute_seconds = show_monotonic_seconds() - step_started;

  /* This commitment component depends only on the fixed CRS and witness.
   * Round 2/5 rejection changes the masks and s2, but never A1 or s1. */
  poly_qshow_mat_d_m1_mul_vec_m1(tA_s1, A1, s1);

  kappa = 0;
reject:
  attempt_count++;
  /****** first round ******/
  step_started = show_monotonic_seconds();
  show_user_prove_round1(proof, chal_1, buf, s2, y1, y2, y3_g,
                         tA_low, w_low, w_high, round1_w, tA_s1,
                         A1, A2, Byg,
                         randomness_seed, kappa);
  round1_seconds += show_monotonic_seconds() - step_started;
  kappa += PARAM_ARP_DIV_N_L_SHOW - PARAM_ARP_DIV_N_SHOW + 1;

  /****** second round ******/
  step_started = show_monotonic_seconds();
  if (!show_user_prove_round2(proof, chal_2, buf, s1, chal_1, y3_g)) {
    round2_seconds += show_monotonic_seconds() - step_started;
    round2_rejections++;
    goto reject;
  }
  round2_seconds += show_monotonic_seconds() - step_started;

  /****** third round ******/
  step_started = show_monotonic_seconds();
  show_user_prove_round3(proof, chal_3_l, chal_3_dk, chal_3_link,
                         buf, sum_gamma_e_star, sum_gamma_r_star, s1,
                         chal_1, chal_2, y3_g, quadratic_precomp);
  round3_seconds += show_monotonic_seconds() - step_started;
  
  /****** fourth round ******/
  step_started = show_monotonic_seconds();
  show_user_prove_round4(proof, buf, s1, s2, Byg, b, A_embed, B_embed,
                         A3_embed, Ds_embed, D_embed, Dx_embed,
                         A_link_embed, D_link_embed,
                         sum_gamma_e_star, sum_gamma_r_star, y1, y2,
                         chal_2, chal_3_l, chal_3_dk, chal_3_link, one);
  round4_seconds += show_monotonic_seconds() - step_started;

  /****** fifth round ******/
  step_started = show_monotonic_seconds();
  if (!show_user_prove_round5(proof, s1, s2, y1, y2,
                              tA_low, w_low, w_high)) {
    round5_seconds += show_monotonic_seconds() - step_started;
    round5_rejections++;
    goto reject;
  }
  round5_seconds += show_monotonic_seconds() - step_started;

  show_print_prove_timing("临时多项式/向量初始化", init_seconds);
  show_print_prove_timing("随机种子、CRS 展开与 statement 打包", crs_seconds);
  show_print_prove_timing("二次关系预计算", precompute_seconds);
  printf("[匿名签名 Prove 计时] 协议尝试 %zu 次；Round 2 拒绝 %zu 次；Round 5 拒绝 %zu 次\n",
         attempt_count, round2_rejections, round5_rejections);
  show_print_prove_timing("Round 1 累计", round1_seconds);
  show_print_prove_timing("Round 2 累计（含拒绝采样）", round2_seconds);
  show_print_prove_timing("Round 3 累计", round3_seconds);
  show_print_prove_timing("Round 4 累计", round4_seconds);
  show_print_prove_timing("Round 5 累计（含拒绝采样）", round5_seconds);

  
  // clean up
  step_started = show_monotonic_seconds();
  // clean up polynomials
  poly_qshow_clear(s1i_star);
  poly_qshow_clear(tmp_poly);
  poly_qshow_clear(one);
  for (i = 0; i < 6; i++) {
    poly_qshow_clear(quadratic_precomp[i]);
  }
  // clean up vectors and matrices
  poly_qshow_mat_d_m1_clear(A1);
  poly_qshow_vec_m1_clear(y1);
  for (i = 0; i < PARAM_L_SHOW; i++) {
    poly_qshow_vec_m1_clear(sum_gamma_r_star[i]);
    poly_qshow_vec_256_clear(sum_gamma_e_star[i]);    
  }
  poly_qshow_mat_d_m2_clear(A2);
  poly_qshow_mat_256l_m2_clear(Byg);
  poly_qshow_vec_m2_clear(s2);
  poly_qshow_vec_m2_clear(y2);
  poly_qshow_vec_m2_clear(b);
  poly_qshow_vec_256_l_clear(y3_g);
  poly_qshow_vec_d_clear(w_high);
  poly_qshow_vec_d_clear(w_low);
  poly_qshow_vec_d_clear(tA_low);
  poly_qshow_vec_d_clear(round1_w);
  poly_qshow_vec_d_clear(tA_s1);
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    poly_qshow_vec_m1_clear(chal_1[i]);
  }
  poly_qshow_vec_l_clear(chal_3_l);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_clear(chal_3_dk[i]);
    poly_qshow_vec_k_clear(chal_3_link[i]);
  }
  show_print_prove_timing("临时对象清理", show_monotonic_seconds() - step_started);
  show_print_prove_timing("show_user_prove_handle 总计",
                          show_monotonic_seconds() - total_started);
}

void show_user_prove(
    show_proof_t             *proof,
    const poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
    const poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D*PARAM_K],
    const poly_qshow_mat_k_k A3_embed[PARAM_D],
    const poly_qshow_mat_k_k Ds_embed[PARAM_D][2*PARAM_D],
    const poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
    const poly_qshow_vec_m1  s1,
    const uint8_t            crs_seed[CRS_SEED_BYTES],
    const uint8_t            seed[SEED_BYTES]) {
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
  show_user_prove_handle(proof, A_embed, B_embed, A3_embed, Ds_embed, D_embed,
                         Dx_embed, A_link_embed, D_link_embed, C_link_embed,
                         s1, crs_seed, seed);
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_clear(A_link_embed[i][j]);
    }
    poly_qshow_vec_k_clear(C_link_embed[i]);
    poly_qshow_mat_k_k_clear(D_link_embed[i]);
    poly_qshow_mat_k_k_clear(Dx_embed[i]);
  }
}
