#include <stdint.h>
#include <stddef.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "malloc.h"
#include "revac_profile.h"
#include <inttypes.h>
#include "timing.h"
#include "proofsystem.h"
#include "labrador_core.h"
#include "labradoodle.h"
#include "fips202.h"
#include "jlproj.h"
#include "poly.h"
#include "polx.h"
#include "polz.h"

static int ldd_kernel_threads(const char *variable, size_t work)
{
#ifdef _OPENMP
  const char *value = getenv(variable);
  char *end = NULL;
  long requested;

  if(value == NULL || value[0] == '\0')
    return revac_labrador_team_threads(work);
  requested = strtol(value, &end, 10);
  if(end == value || *end != '\0' || requested < 1)
    requested = 1;
  if(requested > omp_get_max_threads())
    requested = omp_get_max_threads();
  if((size_t)requested > work)
    requested = (long)work;
  return (int)requested;
#else
  (void)variable;
  (void)work;
  return 1;
#endif
}

/* Six internal CRT primes give roughly 76.697 bits of centered recovery.
 * Keep a small margin and use the exact schoolbook path only when the width
 * metadata says a vector product could approach that boundary. */
static int ldd_exact_product_required(const polxvec r,
                                      const polxvec a,
                                      const polxvec b)
{
  static size_t diagnostic_count = 0;
  size_t i;
  double max_width = 0.0;

  if(r->widths == NULL || a->widths == NULL || b->widths == NULL)
    return 1;
  for(i=0;i<b->len;i++) {
    const double product = r->widths[i*r->stride] +
                           (double)N * a->widths[0] *
                           b->widths[i*b->stride];
    if(product > max_width) max_width = product;
  }
  if(getenv("REVAC_LABRADOR_WIDTH_DIAG") != NULL &&
     diagnostic_count++ < 32)
    fprintf(stderr,
            "[Labrador width 诊断] LDD JL product "
            "log2(sqrt(width))=%.3f, exact=%d\n",
            max_width > 0.0 ? 0.5 * log2(max_width) : -INFINITY,
            max_width <= 0.0 || 0.5 * log2(max_width) >= 75.0);
  return max_width <= 0.0 || 0.5 * log2(max_width) >= 75.0;
}

static void ldd_split_mul_add(polxvec r, const polxvec a, const polxvec b)
{
  /* Decomposition digits are stored in int16_t poly coefficients. */
  const size_t bits = 13;
  const size_t limbs = (LOGQ + bits - 1) / bits;
  poly digits[limbs];
  polxvec digit, term;

  assert(a->len == 1);
  assert(r->len == b->len);
  polxvec_exact_decompose(digits, a, 1, limbs, bits);
  polxvec_init(digit, 1, 1);
  polxvec_init(term, b->len, 1);
  for(size_t limb=0;limb<limbs;limb++) {
    polxvec_frompolyvec(digit, &digits[limb], 1, 1, WIDTHMOD(bits));
    polxvec_mul(term, digit, b);
    polxvec_refresh(term);
    if(limb != 0)
      polxvec_scale(term, term, (int64_t)1 << (bits * limb));
    polxvec_add(r, r, term);
    polxvec_refresh(r);
  }
  polxvec_free(digit);
  polxvec_free(term);
}

static void ldd_collapse_jl(
  polxvec jlmat_agg[LIFTS],
  polxvec phi_jlproj[LIFTS],
  size_t nmax, 
  size_t jlbits_max,
  const uint8_t *jlmat1,
  const uint8_t *jlmat2,
  uint8_t h[16]
)
{
  size_t i;
  const int threads = ldd_kernel_threads(
      "REVAC_LABRADOR_LDD_COLLAPSE_JL_THREADS", LIFTS);
  __attribute__((aligned(64)))
  uint8_t hashbufjl[64 + QBYTES*256 + 24];
  __attribute__((aligned(64)))
  int64_t chalz_jlmat[LIFTS][256];

  for(i=0;i<LIFTS;i++){
    shake128(hashbufjl, sizeof(hashbufjl), h, 16);
    memcpy(h, hashbufjl, 16);
    jlproj_expand_challenge(chalz_jlmat[i], &hashbufjl[64]);
  }

  /* Challenge expansion is transcript ordered above.  The four resulting JL
   * collapses are independent and write disjoint outputs. */
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) if(revac_labrador_ldd_kernel_parallel_enabled("REVAC_LABRADOR_LDD_COLLAPSE_JL_PARALLEL") && threads > 1 && LIFTS >= 2 && !omp_in_parallel()) schedule(static)
  #endif
  for(i=0;i<LIFTS;i++){
    polxvec_init(jlmat_agg[i], nmax, 1);
    polxvec_init(phi_jlproj[i], jlbits_max * 256/N, 1);

    jl_aggregate_mat(jlmat_agg[i], jlmat1, jlmat2, chalz_jlmat[i]);
    jl_aggregate_proj(phi_jlproj[i], jlbits_max, chalz_jlmat[i]);
  }
}

static void ldd_aggregate_zq(
  sparsecnst *zqagg,
  polxvec phi_jlproj_full[LIFTS],
  const polxvec phi_jlproj[LIFTS],
  const size_t *jlbits,
  const lab_params pp,
  const statement ist,
  int64_t chalz_jl[LIFTS][LIFTS*pp->r],
  uint8_t h[16],
  int refresh_output
)
{
  size_t i, k, nchalz, nchalx, chalz_stride;
  int64_t *chalz;
  polxvec chalx[LIFTS];
  timing time;
  const int threads = ldd_kernel_threads(
      "REVAC_LABRADOR_LDD_ZQ_THREADS", LIFTS);
  const int profile_detail = revac_labrador_profile_detail_enabled();
  double profile_started = 0.0;

  nchalz = ist->zqcnst->sparse_nchal + ist->zqcnst->int_nchal;
  nchalx = ist->zqcnst->sigmam1_nchal;
  chalz_stride = MAX(nchalz, 1);

  if (profile_detail) profile_started = revac_labrador_profile_seconds();
  chalz = _malloc(LIFTS * chalz_stride * sizeof(int64_t));
  memset(chalx, 0, sizeof(chalx));
  for(i=0;i<LIFTS;i++){
    if(nchalx > 0){
      polxvec_init(chalx[i], nchalx, 1);
    }
    sample_chalz(&chalz[i * chalz_stride], nchalz, h);
    if(nchalx > 0){
      sample_chalx_uniform(chalx[i], h);
    }
  }
  if (profile_detail)
    fprintf(stderr,
            "[Labrador 原生细分] LDD Zq 挑战采样：nchalz=%zu nchalx=%zu %.6f 秒\n",
            nchalz, nchalx,
            revac_labrador_profile_seconds() - profile_started);

  timing_start(&time, "Set constraints to 0s");
  if (profile_detail) profile_started = revac_labrador_profile_seconds();

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) if(revac_labrador_ldd_kernel_parallel_enabled("REVAC_LABRADOR_LDD_ZQ_PARALLEL") && threads > 1 && LIFTS >= 2 && !omp_in_parallel()) schedule(static) private(i)
  #endif
  for(i=0;i<LIFTS;i++){
    sparsecnst_init(zqagg[i], 1);
    quadfunc_init(zqagg[i]->quad, 0, (pp->r_old * pp->r_old + pp->r_old)/2);
    linfunc_init(zqagg[i]->lin, 1, 1, 1);

    zqagg[i]->lin->off[0] = 0;
    polxvec_init(zqagg[i]->lin->phi[0], pp->nn, 1);
    polxvec_setzero(zqagg[i]->lin->phi[0], 0, 1, pp->nn);
  }

  timing_end(&time);
  timing_print(&time, 4);
  if (profile_detail)
    fprintf(stderr, "[Labrador 原生细分] LDD Zq 约束初始化：%.6f 秒\n",
            revac_labrador_profile_seconds() - profile_started);

  timing_start(&time, "Aggregate inputs");
  if (profile_detail) profile_started = revac_labrador_profile_seconds();

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) if(revac_labrador_ldd_kernel_parallel_enabled("REVAC_LABRADOR_LDD_ZQ_PARALLEL") && threads > 1 && LIFTS >= 2 && !omp_in_parallel()) schedule(static)
  #endif
  for(i=0;i<LIFTS;i++){
    zqcnstset_aggregate_add(zqagg[i], ist->zqcnst,
                            &chalz[i * chalz_stride], chalx[i]);

    /* Proving evaluates these constraints immediately afterwards.  Refresh
     * while the just-aggregated CRT vectors are still resident in this
     * worker's caches, instead of leaving and recreating an OpenMP team for a
     * second full pass.  Reduction does not need this extra normalization. */
    if(refresh_output){
      sparsecnst_refresh(zqagg[i]);
    }

    zqagg[i]->quad->coeffs = realloc(zqagg[i]->quad->coeffs,
                                     zqagg[i]->quad->len*sizeof(polx));
  }

  timing_end(&time);
  timing_print(&time, 4);
  if (profile_detail) {
    fprintf(stderr, "[Labrador 原生细分] LDD Zq 输入关系聚合：%.6f 秒\n",
            revac_labrador_profile_seconds() - profile_started);
  }

  timing_start(&time, "Aggregate phi for the projections");
  if (profile_detail) profile_started = revac_labrador_profile_seconds();

  sample_chalz(chalz_jl[0], LIFTS*LIFTS*pp->r, h);

  for(i=0;i<LIFTS;i++){
    polxvec_init(phi_jlproj_full[i], pp->len[LAB_PROJ], 1);
  }

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) if(revac_labrador_ldd_kernel_parallel_enabled("REVAC_LABRADOR_LDD_ZQ_PARALLEL") && threads > 1 && LIFTS >= 2 && !omp_in_parallel()) schedule(static) private(i)
  #endif
  for(k=0;k<LIFTS;k++){
    size_t j, len;
    polxvec phi;
    for(i=0;i<LIFTS;i++){
      size_t off = 0;
      for(j=0;j<pp->r;j++){
        len = (jlbits[j]-1) * 256/N;
        polxvec_init_subvec2(phi, phi_jlproj[i], 0, 1, len);
        polxvec_init_subvec(phi_jlproj_full[k], phi_jlproj_full[k], off, 1, len);
        if(i==0){
          polxvec_scale(phi_jlproj_full[k], phi, chalz_jl[k][i*pp->r + j]);
        }
        else{
          polxvec_scale_add(phi_jlproj_full[k], phi, chalz_jl[k][i*pp->r + j]);
        }
        off += len;

        // high-order bit is negated
        len = 256/N;
        polxvec_init_subvec2(phi, phi_jlproj[i], (jlbits[j]-1) * 256/N, 1, len);
        polxvec_init_subvec(phi_jlproj_full[k], phi_jlproj_full[k], off, 1, len);
        if(i==0){
          polxvec_scale(phi_jlproj_full[k], phi, -chalz_jl[k][i*pp->r + j]);
        }
        else{
          polxvec_scale_add(phi_jlproj_full[k], phi, -chalz_jl[k][i*pp->r + j]);
        }
        off += len;
      }
    }
    polxvec_init_subvec(phi_jlproj_full[k], phi_jlproj_full[k], 0, 1, 0);
    polxvec_refresh(phi_jlproj_full[k]);
  }

  timing_end(&time);
  timing_print(&time, 4);
  if (profile_detail) {
    fprintf(stderr, "[Labrador 原生细分] LDD Zq 投影 phi 聚合：%.6f 秒\n",
            revac_labrador_profile_seconds() - profile_started);
  }

  free(chalz);
  for(i=0;i<LIFTS;i++){
    if(nchalx > 0){
      polxvec_free(chalx[i]);
    }
  }
}

static void ldd_aggregate_rq(
  sparsecnst finalcnst,
  polxvec chalx_zq,
  const lab_params pp,
  const statement ist,
  const sparsecnst zqagg[LIFTS],
  uint8_t h[16]
)
{
  size_t nchalx;
  polxvec chalx, chalx_sv;
  timing time;
  const int profile_detail = revac_labrador_profile_detail_enabled();
  double profile_started = 0.0;
  double aggregate_started = 0.0;

  timing_start(&time, "Init");

  nchalx = ist->rqcnst->sparse_nchal + ist->rqcnst->com_nchal;
  polxvec_init(chalx, nchalx + LIFTS, 1);
  sample_chalx_aggregate(chalx, h);

  sparsecnst_init(finalcnst, 1);
  quadfunc_init(finalcnst->quad, 0, (pp->r_old * pp->r_old + pp->r_old)/2);
  linfunc_init(finalcnst->lin, 1, 1, 1);

  timing_end(&time);
  timing_print(&time, 4);

  timing_start(&time, "Set constraint to 0s");

  finalcnst->lin->off[0] = 0;
  polxvec_init(finalcnst->lin->phi[0], pp->nn, 1);
  polxvec_setzero(finalcnst->lin->phi[0], 0, 1, pp->nn);

  timing_end(&time);
  timing_print(&time, 4);

  timing_start(&time, "Aggregate input");
  if (profile_detail) aggregate_started = revac_labrador_profile_seconds();

  polxvec base_chalx, base_sparse, base_com, zq_chalx;
  polxvec_init_subvec2(base_chalx, chalx, 0, 1, nchalx);
  polxvec_init_subvec2(base_sparse, base_chalx, 0, 1,
                       ist->rqcnst->sparse_nchal);
  if (profile_detail) profile_started = revac_labrador_profile_seconds();
  sparsecnst_aggregate_add(finalcnst, ist->rqcnst->sparse,
                           ist->rqcnst->nsparse, base_sparse, NULL, 1);
  if (profile_detail) {
    fprintf(stderr, "[Labrador 原生细分] LDD Rq sparse 聚合：%.6f 秒\n",
            revac_labrador_profile_seconds() - profile_started);
  }
  polxvec_init_subvec2(base_com, base_chalx,
                       ist->rqcnst->sparse_nchal, 1,
                       ist->rqcnst->com_nchal);
  if (profile_detail) profile_started = revac_labrador_profile_seconds();
  comcnst_aggregate_add(finalcnst, ist->rqcnst->com,
                        ist->rqcnst->ncom, base_com);
  if (profile_detail) {
    fprintf(stderr, "[Labrador 原生细分] LDD Rq com 聚合：%.6f 秒\n",
            revac_labrador_profile_seconds() - profile_started);
  }
  polxvec_init_subvec2(zq_chalx, chalx, nchalx, 1, LIFTS);
  if (profile_detail) profile_started = revac_labrador_profile_seconds();
  sparsecnst_aggregate_add(finalcnst, zqagg, LIFTS, zq_chalx, NULL, 1);
  if (profile_detail) {
    fprintf(stderr, "[Labrador 原生细分] LDD Rq Zq 聚合：%.6f 秒\n",
            revac_labrador_profile_seconds() - profile_started);
  }
  if (profile_detail) {
    fprintf(stderr, "[Labrador 原生细分] LDD Rq 输入关系聚合：%.6f 秒\n",
            revac_labrador_profile_seconds() - aggregate_started);
  }
  
  timing_end(&time);
  timing_print(&time, 4);

  timing_start(&time, "Refresh");
  if (profile_detail) profile_started = revac_labrador_profile_seconds();

  sparsecnst_refresh(finalcnst);

  timing_end(&time);
  timing_print(&time, 4);
  if (profile_detail) {
    fprintf(stderr, "[Labrador 原生细分] LDD Rq refresh：%.6f 秒\n",
            revac_labrador_profile_seconds() - profile_started);
  }

  finalcnst->quad->coeffs = realloc(finalcnst->quad->coeffs, 
                                    finalcnst->quad->len*sizeof(polx));

  polxvec_init_subvec2(chalx_sv, chalx, nchalx, 1, LIFTS);
  polxvec_copy(chalx_zq, chalx_sv);

  polxvec_free(chalx);
}

static inline void ldd_ling_jl_add(
  polxvec hx,
  polxvec acc,
  size_t idx_phi,
  size_t idx_sx,
  size_t r,
  const int64_t chalz_jl[LIFTS][LIFTS*r],
  const polxvec chalx_zq,
  const polxvec jl_sprod_sx[LIFTS][r]
)
{
  size_t i, j;
  polxvec chalx;

  for(i=0;i<LIFTS;i++){
    polxvec_scale(acc, jl_sprod_sx[0][idx_sx], chalz_jl[i][idx_phi]);
    for(j=1;j<LIFTS;j++){
      polxvec_scale_add(acc, jl_sprod_sx[j][idx_sx], chalz_jl[i][r*j+idx_phi]);
    }
    /* The four 32-bit JL scalars can grow acc close to the six-prime CRT
     * recovery limit.  Canonicalize before multiplying by the sparse ring
     * challenge: refreshing afterwards cannot undo a wrap that already
     * happened inside polxvec_mul_add(). */
    polxvec_exact_refresh(acc);
    polxvec_init_subvec2(chalx, chalx_zq, i, 0, 1);
    polxvec_mul_add(hx, chalx, acc);

    /* Keep the accumulated h value canonical before the next lift.  The
     * vector has length one here, so the exact refresh is inexpensive. */
    polxvec_exact_refresh(hx);
  }
}

static void ldd_ling(
  poly *hy,
  size_t r, 
  size_t fu, 
  size_t bu,
  const polxvec sx[], 
  const polxvec phi[],
  const polxvec chalx_zq,
  const polxvec jl_sprod_sx[LIFTS][r],
  const int64_t chalz_jl[LIFTS][LIFTS*r]
)
{
  const int threads = ldd_kernel_threads(
      "REVAC_LABRADOR_LDD_LING_THREADS", r);
  const char *width_diag_value = getenv("REVAC_LABRADOR_WIDTH_DIAG");
  const int width_diag = width_diag_value != NULL &&
                         width_diag_value[0] != '\0' &&
                         width_diag_value[0] != '0';
  const char *ling_canonicalize_value = getenv("REVAC_LABRADOR_LING_CANONICALIZE");
  const int ling_canonicalize = ling_canonicalize_value == NULL ||
                                ling_canonicalize_value[0] == '\0' ||
                                ling_canonicalize_value[0] != '0';
  double max_width = 0.0;

  #ifdef _OPENMP
  #pragma omp parallel num_threads(threads) \
      reduction(max:max_width) \
      if(revac_labrador_ldd_kernel_parallel_enabled(\
             "REVAC_LABRADOR_LDD_LING_PARALLEL") && \
         threads > 1 && r >= 2 && !omp_in_parallel())
  #endif
  {
    polxvec hx, acc, phi_sv, sx_sv;
    size_t i, j, len, hpos;
    polxvec_init(hx, 1, 1);
    polxvec_init(acc, 1, 1);

    #ifdef _OPENMP
    #pragma omp for schedule(dynamic)
    #endif
    for(i=0;i<r;i++){
      for(j=i;j<r;j++){
        if(i == j){
          polxvec_sprod(hx, phi[i], sx[i]);
          polxvec_exact_refresh(hx);
          ldd_ling_jl_add(hx, acc, i, i, r, chalz_jl, chalx_zq,
                          jl_sprod_sx);
        }
        else{
          len = MIN(phi[i]->len, sx[j]->len);
          polxvec_init_subvec2(phi_sv, phi[i], 0, 1, len);
          polxvec_init_subvec2(sx_sv, sx[j], 0, 1, len);
          polxvec_sprod(hx, phi_sv, sx_sv);
          polxvec_exact_refresh(hx);

          len = MIN(phi[j]->len, sx[i]->len);
          polxvec_init_subvec2(phi_sv, phi[j], 0, 1, len);
          polxvec_init_subvec2(sx_sv, sx[i], 0, 1, len);
          polxvec_sprod_add(hx, phi_sv, sx_sv);
          polxvec_exact_refresh(hx);

          ldd_ling_jl_add(hx, acc, i, j, r, chalz_jl, chalx_zq,
                          jl_sprod_sx);
          ldd_ling_jl_add(hx, acc, j, i, r, chalz_jl, chalx_zq,
                          jl_sprod_sx);
        }

        if(width_diag && hx->widths != NULL && *hx->widths > max_width)
          max_width = *hx->widths;

        hpos = trimat_idx(i, j, r) * fu;
        if(ling_canonicalize)
          polxvec_exact_refresh(hx);
        polxvec_decompose(&hy[hpos], hx, 1, fu, bu);
        const char *residual_diag = getenv("REVAC_LABRADOR_RESIDUAL_DIAG");
        if(residual_diag != NULL && residual_diag[0] != '\0' &&
           residual_diag[0] != '0'){
          polxvec reconstructed;
          polxvec_init(reconstructed, 1, 1);
          polxvec_reconstruct(reconstructed, &hy[hpos], 1, fu, bu);
          polxvec_sub(reconstructed, reconstructed, hx);
          if(!polxvec_iszero(reconstructed))
            fprintf(stderr,
                    "[Labrador h 分解诊断] h[%zu,%zu] 分解重构不一致\n",
                    i, j);
          polxvec_free(reconstructed);
        }
      }
    }
    polxvec_free(hx);
    polxvec_free(acc);
  }
  if(width_diag)
    fprintf(stderr,
            "[Labrador width 诊断] LDD linear-garbage r=%zu "
            "log2(sqrt(max-width))=%.3f，limit=76.697\n",
            r, max_width > 0.0 ? 0.5 * log2(max_width) : -INFINITY);
}

static void ldd_addcheck_commit(
  comcnst c[6],
  sparsecnst czq[4*SIS1_NCOEF], 
  const lab_params pp,
  const lab_proof pi
)
{
  size_t i;

  ps_addcheck_commit_outer(c[0], pp->off[LAB_OUTCOM], pp->off[LAB_U1],
                           pp->len[LAB_U1] + pp->len[LAB_RAND1], pp->bu, pp->fu);
  ps_addcheck_commit_outer(c[1], pp->off[LAB_OUTCOM]+pp->fu, pp->off[LAB_PROJ],
                           pp->len[LAB_PROJ]+pp->len[LAB_RAND2], pp->bu, pp->fu);
  ps_addcheck_commit_outer(c[2], pp->off[LAB_OUTCOM]+2*pp->fu, pp->off[LAB_LIFT],
                           pp->len[LAB_LIFT]+pp->len[LAB_RAND3], pp->bu, pp->fu);
  ps_addcheck_commit_outer(c[3], pp->off[LAB_OUTCOM]+3*pp->fu, pp->off[LAB_U2],
                           pp->len[LAB_U2]+pp->len[LAB_RAND4], pp->bu, pp->fu);

  ps_addcheck_commit_middle(c[4], pp->kappa[1], pp->off[LAB_U1], 
                    pp->off[LAB_INCOM], pp->len[LAB_INCOM]+pp->len[LAB_QUADG]);
  ps_addcheck_commit_middle(c[5], pp->kappa[1], pp->off[LAB_U2],
                            pp->off[LAB_LING], pp->len[LAB_LING]);

  for(i=0;i<4;i++){
    ps_addcheck_commit_coeffs(&czq[i*SIS1_NCOEF], *pi->m[i], 
                              pp->off[LAB_OUTCOM] + i*pp->fu, pp->bu, pp->fu);
  }
}

static void ldd_addcheck_ling(
  sparsecnst c,
  const lab_params pp,
  const int64_t chalz_jl[LIFTS][LIFTS*pp->r],
  polx chalx_amortize[pp->r],
  const polxvec chalx_zq,
  const polxvec phi[pp->r],
  const polxvec jlmat_agg[LIFTS]
)
{
  size_t i, j, k, off;
  const char *ling_refresh_value = getenv("REVAC_LABRADOR_LING_REFRESH");
  const int ling_refresh = ling_refresh_value == NULL ||
                           ling_refresh_value[0] == '\0' ||
                           ling_refresh_value[0] != '0';
  polx cprod;
  polxvec powers, phi_sv, chalx_sv, acc1, acc2;

  sparsecnst_init(c, 1);
  linfunc_init(c->lin, 1, pp->fz + 1, pp->fz + 1);

  c->lin->off[0] = pp->off[LAB_Z];
  polxvec_init(c->lin->phi[0], pp->nmax, 1);
  polxvec_setzero(c->lin->phi[0], 0, 1, c->lin->phi[0]->len);

  for(i=0;i<pp->r;i++){
    polxvec_init_subvec(c->lin->phi[0], c->lin->phi[0], 0, 1, phi[i]->len);
    polxvec_polx_mul_add(c->lin->phi[0], chalx_amortize[i], phi[i]);
  }
  polxvec_init_subvec(c->lin->phi[0], c->lin->phi[0], 0, 1, 0);

  // add contribution of the projection matrices to phi

  polxvec_init(acc1, 1, 1);
  polxvec_init(acc2, 1, 1);

  for(k=0;k<LIFTS;k++){
    for(i=0;i<pp->r;i++){
      for(j=0;j<LIFTS;j++){
        polxvec_init_subvec2(chalx_sv, chalx_zq, j, 1, 1);
        if(j==0){
          polxvec_scale(acc1, chalx_sv, chalz_jl[j][pp->r*k + i]);
        }
        else{
          polxvec_scale_add(acc1, chalx_sv,
                            chalz_jl[j][pp->r*k + i]);
        }
      }
      if(i==0){
        polxvec_polx_mul(acc2, chalx_amortize[i], acc1);
      }
      else{
        polxvec_polx_mul_add(acc2, chalx_amortize[i], acc1);
      }
    }
    if(ling_refresh)
      polxvec_exact_refresh(acc2);
    /* Both operands are dense aggregates here.  Their polynomial product can
     * exceed the centered range of the six internal CRT primes before a
     * post-product refresh gets a chance to reduce it.  Split the 38-bit
     * scalar polynomial into three exact 13-bit limbs in that case; each limb
     * product is safely inside the CRT range and still uses the fast NTT. */
    if(ldd_exact_product_required(c->lin->phi[0], acc2, jlmat_agg[k]))
      ldd_split_mul_add(c->lin->phi[0], acc2, jlmat_agg[k]);
    else {
      polxvec_mul_add(c->lin->phi[0], acc2, jlmat_agg[k]);
      polxvec_refresh(c->lin->phi[0]);
    }
  }
  if(getenv("REVAC_LABRADOR_WIDTH_DIAG") != NULL){
    fprintf(stderr, "[Labrador width 诊断] ling phi width=%.3e\n",
            c->lin->phi[0]->widths != NULL ? *c->lin->phi[0]->widths : 0.0);
    for(k=0;k<LIFTS;k++)
      fprintf(stderr, "[Labrador width 诊断] jlmat[%zu] width=%.3e\n", k,
              jlmat_agg[k]->widths != NULL ? *jlmat_agg[k]->widths : 0.0);
  }

  if(pp->fz == 2){
    c->lin->off[1] = pp->off[LAB_Z] + pp->nmax;
    polxvec_init(c->lin->phi[1], pp->nmax, 1);
    polxvec_scale(c->lin->phi[1], c->lin->phi[0], 1<<pp->bz);
  }

  c->lin->off[pp->fz] = pp->off[LAB_LING];
  polxvec_init(c->lin->phi[pp->fz], (pp->r*pp->r + pp->r)/2 * pp->fu, 1);

  polxvec_init(powers, pp->fu, 1);
  polxvec_powers(powers, 1 << pp->bu, -1, -1);

  off = 0;
  for(i=0;i<pp->r;i++){
    for(j=i;j<pp->r;j++){
      // no need to multiply by 2 the non-diagonal since those hij are already
      // multiplied by 2 when computed
      polx_mul(cprod, chalx_amortize[i], chalx_amortize[j]);
      polxvec_init_subvec2(phi_sv, c->lin->phi[pp->fz], off, 1, pp->fu);
      polxvec_polx_mul(phi_sv, cprod, powers);
      off += pp->fu;
    }
  }
  polxvec_free(acc1);
  polxvec_free(acc2);
  polxvec_free(powers);
}

static int ldd_residual_diagnostics_enabled(void)
{
  const char *value = getenv("REVAC_LABRADOR_RESIDUAL_DIAG");
  return value != NULL && value[0] != '\0' && value[0] != '0';
}

static void ldd_diagnose_ling_components(
  const sparsecnst c,
  const lab_params pp,
  const poly *sout,
  polx chalx_amortize[pp->r],
  const polxvec sx[pp->r],
  const polxvec phi[pp->r],
  const polxvec chalx_zq,
  const polxvec jl_sprod_sx[LIFTS][pp->r],
  const int64_t chalz_jl[LIFTS][LIFTS*pp->r],
  const polxvec jlmat_agg[LIFTS]
)
{
  size_t i, residual_coefficient = N;
  polxvec base_phi, jl_phi, z_value, ling_digits;
  polxvec base_ev, jl_ev, ling_ev, total_ev;
  polz base_z, jl_z, ling_z, total_z;

  if(!ldd_residual_diagnostics_enabled()) return;

  polxvec_init(base_phi, pp->nmax, 1);
  polxvec_setzero(base_phi, 0, 1, pp->nmax);
  for(i=0;i<pp->r;i++){
    polxvec_init_subvec(base_phi, base_phi, 0, 1, phi[i]->len);
    polxvec_polx_mul_add(base_phi, chalx_amortize[i], phi[i]);
  }
  polxvec_init_subvec(base_phi, base_phi, 0, 1, 0);

  polxvec_init(jl_phi, pp->nmax, 1);
  polxvec_sub(jl_phi, c->lin->phi[0], base_phi);
  polxvec_init(z_value, pp->nmax, 1);
  polxvec_reconstruct(z_value, &sout[pp->off[LAB_Z]], pp->nmax,
                      pp->fz, pp->bz);

  polxvec_init(base_ev, 1, 1);
  polxvec_init(jl_ev, 1, 1);
  polxvec_init(ling_ev, 1, 1);
  polxvec_init(total_ev, 1, 1);
  polxvec_sprod(base_ev, base_phi, z_value);
  polxvec_sprod(jl_ev, jl_phi, z_value);

  polxvec_init(ling_digits, pp->len[LAB_LING], 1);
  polxvec_frompolyvec(ling_digits, &sout[pp->off[LAB_LING]], 1,
                      pp->len[LAB_LING], WIDTHMOD(pp->bu));
  polxvec_sprod(ling_ev, c->lin->phi[pp->fz], ling_digits);

  polxvec_add(total_ev, base_ev, jl_ev);
  polxvec_add(total_ev, total_ev, ling_ev);
  polzvec_exact_frompolxvec(&base_z, base_ev, 0, 1, 1);
  polzvec_exact_frompolxvec(&jl_z, jl_ev, 0, 1, 1);
  polzvec_exact_frompolxvec(&ling_z, ling_ev, 0, 1, 1);
  polzvec_exact_frompolxvec(&total_z, total_ev, 0, 1, 1);

  for(i=0;i<N;i++){
    zz base_coeff, jl_coeff, ling_coeff, total_coeff;
    int64_t total_value;

    polz_getcoeff(total_coeff, total_z, i);
    total_value = int64_fromzz(total_coeff);
    if(total_value == 0) continue;
    polz_getcoeff(base_coeff, base_z, i);
    polz_getcoeff(jl_coeff, jl_z, i);
    polz_getcoeff(ling_coeff, ling_z, i);
    fprintf(stderr,
            "[Labrador ling 分量诊断] coeff=%zu base=%lld jl=%lld "
            "ling=%lld total=%lld\n",
            i, (long long)int64_fromzz(base_coeff),
            (long long)int64_fromzz(jl_coeff),
            (long long)int64_fromzz(ling_coeff), (long long)total_value);
    residual_coefficient = i;
    break;
  }

  if(residual_coefficient < N){
    size_t j;
    polx cprod;
    polxvec base_h, jl_h, jl_acc, phi_sv, sx_sv;
    polxvec base_sum, jl_sum, base_delta, jl_delta;
    polxvec direct_phi, direct_acc1, direct_acc2, direct_ev, direct_delta;
    polz base_delta_z, jl_delta_z, direct_delta_z;
    zz base_delta_coeff, jl_delta_coeff, direct_delta_coeff;

    polxvec_init(base_h, 1, 1);
    polxvec_init(jl_h, 1, 1);
    polxvec_init(jl_acc, 1, 1);
    polxvec_init(base_sum, 1, 1);
    polxvec_init(jl_sum, 1, 1);
    polxvec_init(direct_phi, pp->nmax, 1);
    polxvec_setzero(direct_phi, 0, 1, pp->nmax);
    polxvec_setzero(base_sum, 0, 1, 1);
    polxvec_setzero(jl_sum, 0, 1, 1);
    for(i=0;i<pp->r;i++){
      for(j=i;j<pp->r;j++){
        if(i == j){
          polxvec_sprod(base_h, phi[i], sx[i]);
        }
        else{
          const size_t first_len = MIN(phi[i]->len, sx[j]->len);
          const size_t second_len = MIN(phi[j]->len, sx[i]->len);
          polxvec_init_subvec2(phi_sv, phi[i], 0, 1, first_len);
          polxvec_init_subvec2(sx_sv, sx[j], 0, 1, first_len);
          polxvec_sprod(base_h, phi_sv, sx_sv);
          polxvec_init_subvec2(phi_sv, phi[j], 0, 1, second_len);
          polxvec_init_subvec2(sx_sv, sx[i], 0, 1, second_len);
          polxvec_sprod_add(base_h, phi_sv, sx_sv);
        }
        polxvec_setzero(jl_h, 0, 1, 1);
        ldd_ling_jl_add(jl_h, jl_acc, i, j, pp->r, chalz_jl,
                        chalx_zq, jl_sprod_sx);
        if(i != j)
          ldd_ling_jl_add(jl_h, jl_acc, j, i, pp->r, chalz_jl,
                          chalx_zq, jl_sprod_sx);
        polx_mul(cprod, chalx_amortize[i], chalx_amortize[j]);
        polxvec_polx_mul_add(base_sum, cprod, base_h);
        polxvec_polx_mul_add(jl_sum, cprod, jl_h);
      }
    }
    /* Independent verifier-side JL recomputation: do not use jl_sprod_sx or
     * c->lin->phi.  This mirrors ldd_addcheck_ling's index order directly. */
    polxvec_init(direct_acc1, 1, 1);
    polxvec_init(direct_acc2, 1, 1);
    for(size_t k=0;k<LIFTS;k++){
      polxvec_setzero(direct_acc2, 0, 1, 1);
      for(i=0;i<pp->r;i++){
        for(size_t l=0;l<LIFTS;l++){
          polxvec chalx_sv_l;
          polxvec_init_subvec2(chalx_sv_l, chalx_zq, l, 1, 1);
          if(l == 0)
            polxvec_scale(direct_acc1, chalx_sv_l,
                          chalz_jl[l][pp->r*k + i]);
          else
            polxvec_scale_add(direct_acc1, chalx_sv_l,
                              chalz_jl[l][pp->r*k + i]);
        }
        if(i == 0)
          polxvec_polx_mul(direct_acc2, chalx_amortize[i], direct_acc1);
        else
          polxvec_polx_mul_add(direct_acc2, chalx_amortize[i], direct_acc1);
      }
      if(ldd_exact_product_required(direct_phi, direct_acc2, jlmat_agg[k]))
        ldd_split_mul_add(direct_phi, direct_acc2, jlmat_agg[k]);
      else {
        polxvec_mul_add(direct_phi, direct_acc2, jlmat_agg[k]);
        polxvec_refresh(direct_phi);
      }
    }
    polxvec_init(direct_ev, 1, 1);
    polxvec_sprod(direct_ev, direct_phi, z_value);
    polxvec_init(base_delta, 1, 1);
    polxvec_init(jl_delta, 1, 1);
    polxvec_sub(base_delta, base_ev, base_sum);
    polxvec_sub(jl_delta, jl_ev, jl_sum);
    polzvec_exact_frompolxvec(&base_delta_z, base_delta, 0, 1, 1);
    polzvec_exact_frompolxvec(&jl_delta_z, jl_delta, 0, 1, 1);
    polxvec_init(direct_delta, 1, 1);
    polxvec_sub(direct_delta, jl_ev, direct_ev);
    polzvec_exact_frompolxvec(&direct_delta_z, direct_delta, 0, 1, 1);
    polz_getcoeff(base_delta_coeff, base_delta_z, residual_coefficient);
    polz_getcoeff(jl_delta_coeff, jl_delta_z, residual_coefficient);
    polz_getcoeff(direct_delta_coeff, direct_delta_z, residual_coefficient);
    fprintf(stderr,
            "[Labrador ling 恒等式诊断] coeff=%zu base_delta=%lld "
            "jl_delta=%lld\n",
            residual_coefficient,
            (long long)int64_fromzz(base_delta_coeff),
            (long long)int64_fromzz(jl_delta_coeff));
    fprintf(stderr,
            "[Labrador ling 独立JL诊断] coeff=%zu direct_delta=%lld\n",
            residual_coefficient,
            (long long)int64_fromzz(direct_delta_coeff));
    polxvec_free(base_h);
    polxvec_free(jl_h);
    polxvec_free(jl_acc);
    polxvec_free(base_sum);
    polxvec_free(jl_sum);
    polxvec_free(base_delta);
    polxvec_free(jl_delta);
    polxvec_free(direct_phi);
    polxvec_free(direct_acc1);
    polxvec_free(direct_acc2);
    polxvec_free(direct_ev);
    polxvec_free(direct_delta);
  }

  polxvec_free(base_phi);
  polxvec_free(jl_phi);
  polxvec_free(z_value);
  polxvec_free(ling_digits);
  polxvec_free(base_ev);
  polxvec_free(jl_ev);
  polxvec_free(ling_ev);
  polxvec_free(total_ev);
}

static void ldd_addcheck_system(
  sparsecnst c,
  const lab_params pp,
  const sparsecnst finalcnst,
  const polxvec chalx,
  const polxvec phi_jl[LIFTS]
)
{
  size_t i, j, k, l, idx, off, ncoeffs;
  size_t rmap[pp->r_old];
  polxvec powers, chalx_sv, phi_sv, gphi;

  sparsecnst_init(c, 1);

  ncoeffs = 0;
  for(i=0;i<finalcnst->quad->len;i++){
    ncoeffs += pp->rr[finalcnst->quad->rows[i]];
  }
  off = 0;
  for(i=0;i<pp->r_old;i++){
    rmap[i] = off;
    off += pp->rr[i];
  }

  linfunc_init(c->lin, 1, ncoeffs + pp->r + 2, ncoeffs + pp->r + 2);

  polxvec_copy(c->b, finalcnst->b);

  idx = 0;

  // gij
  polxvec_init(gphi, pp->fg, 1);
  polxvec_init(powers, pp->fg, 1);
  polxvec_powers(powers, 1 << pp->bg, 1, 1);
  for(k=0;k<finalcnst->quad->len;k++){
    i = finalcnst->quad->rows[k];
    j = finalcnst->quad->cols[k];
    polxvec_polx_mul(gphi, finalcnst->quad->coeffs[k], powers);
    polxvec_refresh(gphi);
    for(l=0;l<pp->rr[i];l++){
      // a[i][j] * g[rmap[i]+l][rmap[j]+l]
      c->lin->off[idx] = pp->off[LAB_QUADG] 
                         + pp->fg * trimat_idx(rmap[i]+l, rmap[j]+l, pp->r);
      polxvec_init(c->lin->phi[idx], pp->fg, 1);
      polxvec_copy(c->lin->phi[idx], gphi);
      idx++;
    }
  }
  polxvec_free(gphi);
  polxvec_free(powers);

  // hii
  polxvec_init(powers, pp->fu, 1);
  polxvec_powers(powers, 1 << pp->bu, 1, 1);
  off = 0;
  for(i=0;i<pp->r;i++){
    c->lin->off[idx] = pp->off[LAB_LING] + off;
    polxvec_init(c->lin->phi[idx], pp->fu, 1);
    polxvec_copy(c->lin->phi[idx], powers);
    off += pp->fu * (pp->r - i);
    idx++;
  }
  polxvec_free(powers);

  // liftings
  c->lin->off[idx] = pp->off[LAB_LIFT];
  polxvec_init(c->lin->phi[idx], pp->len[LAB_LIFT], 1);
  
  polxvec_init(powers, LOGQ, 1);
  polxvec_powers(powers, 2, -1, 1);
  for(i=0;i<LIFTS;i++){
    polxvec_init_subvec2(chalx_sv, chalx, i, 1, 1);
    polxvec_init_subvec2(phi_sv, c->lin->phi[idx], i * LOGQ, 1, LOGQ);
    polxvec_mul(phi_sv, chalx_sv, powers);
  }
  polxvec_free(powers);
  idx++;

  // projections
  c->lin->off[idx] = pp->off[LAB_PROJ];
  polxvec_init(c->lin->phi[idx], pp->len[LAB_PROJ], 1);

  for(i=0;i<LIFTS;i++){
    polxvec_init_subvec2(chalx_sv, chalx, i, 1, 1);
    if(i==0){
      polxvec_mul(c->lin->phi[idx], chalx_sv, phi_jl[i]);
    }
    else{
      polxvec_mul_add(c->lin->phi[idx], chalx_sv, phi_jl[i]);
    }
  }
}

static void ldd_addchecks(
  statement ost,
  const lab_params pp,
  const lab_proof pi,
  const int64_t chalz_jl[LIFTS][LIFTS*pp->r],
  polx *chalx_amortize,
  const polxvec chalx_zq,
  const polxvec phi[pp->r],
  const polxvec jlmat_agg[LIFTS],
  const polxvec phi_jl[LIFTS],
  const sparsecnst finalcnst
)
{
  const int threads = ldd_kernel_threads(
      "REVAC_LABRADOR_LDD_ADDCHECK_THREADS", 6);

  rqcnstset_init(ost->rqcnst, 3, 7);
  zqcnstset_init(ost->zqcnst, LIFTS+4*SIS1_NCOEF, LIFTS+4*SIS1_NCOEF + 1, 0, 1, 0);

  ost->rqcnst->sparse_nchal = 3;
  ost->rqcnst->com_nchal = pp->kappa[0] + 2*pp->kappa[1] + 4*pp->kappa[2];
  ost->zqcnst->sparse_nchal = LIFTS + 4*SIS1_NCOEF;

  #ifdef _OPENMP
  #pragma omp parallel sections num_threads(threads) \
      if(threads > 1 && !omp_in_parallel())
  #endif
  {
    #ifdef _OPENMP
    #pragma omp section
    #endif
    {
      const double started = revac_labrador_profile_seconds();
      ps_addcheck_lift_zero_coeff(ost->zqcnst->sparse, pp->off[LAB_LIFT]);
      if(revac_labrador_profile_detail_enabled())
        fprintf(stderr, "[Labrador Add checks] lift zero: %.6f s\n",
                revac_labrador_profile_seconds() - started);
    }

    #ifdef _OPENMP
    #pragma omp section
    #endif
    {
      const double started = revac_labrador_profile_seconds();
      ldd_addcheck_commit(ost->rqcnst->com,
                          &ost->zqcnst->sparse[LIFTS], pp, pi);
      if(revac_labrador_profile_detail_enabled())
        fprintf(stderr, "[Labrador Add checks] commit: %.6f s\n",
                revac_labrador_profile_seconds() - started);
    }

    #ifdef _OPENMP
    #pragma omp section
    #endif
    {
      const double started = revac_labrador_profile_seconds();
      lab_addcheck_amortization(ost->rqcnst->com[6], pp->kappa[0], pp->r,
                                pp->nmax, pp->off[LAB_Z], pp->off[LAB_INCOM],
                                pp->fz, pp->bz, pp->fu, pp->bu,
                                chalx_amortize);
      if(revac_labrador_profile_detail_enabled())
        fprintf(stderr, "[Labrador Add checks] amortization: %.6f s\n",
                revac_labrador_profile_seconds() - started);
    }

    #ifdef _OPENMP
    #pragma omp section
    #endif
    {
      const double started = revac_labrador_profile_seconds();
      lab_addcheck_quadg(ost->rqcnst->sparse[0], pp->r,
                         pp->off[LAB_QUADG], pp->fz, pp->bz, pp->fg, pp->bg,
                         chalx_amortize);
      if(revac_labrador_profile_detail_enabled())
        fprintf(stderr, "[Labrador Add checks] quadg: %.6f s\n",
                revac_labrador_profile_seconds() - started);
    }

    #ifdef _OPENMP
    #pragma omp section
    #endif
    {
      const double started = revac_labrador_profile_seconds();
      ldd_addcheck_ling(ost->rqcnst->sparse[1], pp, chalz_jl,
                        chalx_amortize, chalx_zq, phi, jlmat_agg);
      if(revac_labrador_profile_detail_enabled())
        fprintf(stderr, "[Labrador Add checks] ling: %.6f s\n",
                revac_labrador_profile_seconds() - started);
    }

    #ifdef _OPENMP
    #pragma omp section
    #endif
    {
      const double started = revac_labrador_profile_seconds();
      ldd_addcheck_system(ost->rqcnst->sparse[2], pp, finalcnst, chalx_zq,
                          phi_jl);
      if(revac_labrador_profile_detail_enabled())
        fprintf(stderr, "[Labrador Add checks] system: %.6f s\n",
                revac_labrador_profile_seconds() - started);
    }
  }
}

void ldd_prove(
  lab_proof pi, 
  statement ost, 
  witness owt, 
  const statement ist, 
  const witness iwt, 
  const lab_params pp
)
{
  size_t i, j, tpos, projpos, off, projoff[pp->r];
  size_t jlbits[pp->r], jlbits_max;
  uint8_t *jlmat1, *jlmat2;
  int64_t bufbits[4*pp->randlen*N], chalz_jl[LIFTS][LIFTS*pp->r];
  poly *sy[pp->r], *sout;
  polx *chalx_amortize;
  polxvec sxl, sxl_jl, sxq[pp->r], sxq_old[ist->r], phi[pp->r];
  polxvec tx, u1x, tgd, U1, u1d, Up, projd, liftd, Ub, u2x, hd;
  polxvec chalx_zq, U2, u2d, zx;
  polxvec jl_sprod_sx[LIFTS][pp->r], jlmat_agg[LIFTS], phi_jlproj[LIFTS];
  polxvec phi_jlproj_full[LIFTS];
  sparsecnst zqagg[LIFTS], finalcnst;
  timing time;

  timing_start(&time, "Process input witness");

  // View input witness with new split

  sy[0] = iwt->s[0];
  for(i=1;i<pp->r;i++){
    sy[i] = &sy[i-1][pp->n[i-1]];
  }

  // Input witness to polx

  polxvec_init(sxl, pp->nn, 1);
  polxvec_frompolyvec(sxl, sy[0], 1, pp->nn, 0);
  
  off = 0;
  for(i=0;i<pp->r;i++){
    polxvec_setwidths1(sxl, off, 1, pp->n[i], pp->normsq[i]/(pp->n[i]*N));
    polxvec_init_subvec2(sxq[i], sxl, off, 1, pp->n[i]);
    off += pp->n[i];
  }
  off = 0;
  for(i=0;i<ist->r;i++){
    polxvec_init_subvec2(sxq_old[i], sxl, off, 1, ist->n[i]);
    off += ist->n[i];
  }

  timing_end(&time);
  timing_print(&time, 3);

  // Init

  timing_start(&time, "Init output");
  
  lab_witness_init(owt, pp);
  lab_proof_init(pi, pp);
  lab_statement_init(ost, ist, pp);
  sout = owt->s[0];

  jlbits_max = 0;
  for(i=0;i<pp->r;i++){
    jlbits[i] = ceil(log2(JL_INF_MULT*sqrt(pp->normsq[i])));
    jlbits_max = MAX(jlbits_max, jlbits[i]);
  }

  timing_end(&time);
  timing_print(&time, 3);

  timing_start(&time, "Init comkey");

  lab_comkey_init(pp);

  timing_end(&time);
  timing_print(&time, 3);

  // Sample outer commitment randomness

  if(pp->randlen > 0){
    randombits64(bufbits, 4 * pp->randlen * N);
    polyvec_fromint64vec(&sout[pp->off[LAB_RAND1]], bufbits, 1, pp->randlen, 1, 
                        NULL);
    polyvec_fromint64vec(&sout[pp->off[LAB_RAND2]], &bufbits[pp->randlen*N], 
                        1, pp->randlen, 1, NULL);
    polyvec_fromint64vec(&sout[pp->off[LAB_RAND3]], &bufbits[2*pp->randlen*N], 
                        1, pp->randlen, 1, NULL);
    polyvec_fromint64vec(&sout[pp->off[LAB_RAND4]], &bufbits[3*pp->randlen*N], 
                        1, pp->randlen, 1, NULL);
  }

  // Inner commitments

  timing_start(&time, "Inner commitments");

  polxvec_init(tx, pp->kappa[0], 1);
  tpos = pp->off[LAB_INCOM];
  for(i=0;i<pp->r;i++){
    commit(tx, sxq[i]);
    polxvec_decompose(&sout[tpos], tx, tx->len, pp->fu, pp->bu);
    tpos += pp->kappa[0] * pp->fu;
  }

  timing_end(&time);
  timing_print(&time, 3);

  // Quadratic garbage

  timing_start(&time, "Quadratic garbage");

  lab_quadg(&sout[pp->off[LAB_QUADG]], NULL, sxq, pp->r, pp->fg, pp->bg, 0);

  timing_end(&time);
  timing_print(&time, 3);

  // Middle commitment u1

  polxvec_init(u1x, pp->kappa[1], 1);
  polxvec_init(tgd, pp->len[LAB_INCOM] + pp->len[LAB_QUADG], 1);
  polxvec_frompolyvec(tgd, &sout[pp->off[LAB_INCOM]], 1, pp->len[LAB_INCOM],
                      WIDTHMOD(pp->bu));
  polxvec_init_subvec(tgd, tgd, pp->len[LAB_INCOM], 1, pp->len[LAB_QUADG]);
  polxvec_frompolyvec(tgd, &sout[pp->off[LAB_QUADG]], 1, tgd->len, 
                      WIDTHMOD(pp->bg));
  polxvec_init_subvec(tgd, tgd, 0, 1, 0);
  commit(u1x, tgd);
  polxvec_bindec(&sout[pp->off[LAB_U1]], u1x, u1x->len, LOGQ);

  // Outer commitment U1

  polxvec_init(U1, pp->kappa[2], 1);
  polxvec_init(u1d, pp->len[LAB_U1] + pp->len[LAB_RAND1], 1);
  polxvec_frompolyvec(u1d, &sout[pp->off[LAB_U1]], 1, pp->len[LAB_U1], 0.5);
  if(pp->len[LAB_RAND1] > 0){
    polxvec_init_subvec(u1d, u1d, pp->len[LAB_U1], 1, pp->len[LAB_RAND1]);
    polxvec_frompolyvec(u1d, &sout[pp->off[LAB_RAND1]], 1, u1d->len, 0.5);
    polxvec_init_subvec(u1d, u1d, 0, 1, 0);
  }
  commit(U1, u1d);
  polxvec_decompose(&sout[pp->off[LAB_OUTCOM]], U1, 1, pp->fu, pp->bu);
  polzvec_frompolxvec(pi->m[0], U1, 0, 1, pp->kappa[2]);
  outcom_clear(*pi->m[0]);

  update_hash_polz(ost->h, pi->m[0], pp->kappa[2]);

  // JL projections

  timing_start(&time, "Compute JL projections");
  
  jl_sample_mat(&jlmat1, &jlmat2, ost->h, pp->nmax);
  const int projection_threads = ldd_kernel_threads(
      "REVAC_LABRADOR_LDD_PROJECTION_THREADS", pp->r);

  projpos = pp->off[LAB_PROJ];
  for(i=0;i<pp->r;i++){
    projoff[i] = projpos;
    projpos += 256/N * jlbits[i];
  }

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(projection_threads) \
      if(revac_labrador_ldd_kernel_parallel_enabled(\
             "REVAC_LABRADOR_LDD_PROJECTION_PARALLEL") && \
         projection_threads > 1 && pp->r >= 2 && !omp_in_parallel()) \
      schedule(static)
  #endif
  for(i=0;i<pp->r;i++){
    int32_t proj32_local[256];
    int64_t proj64_local[256];
    polxvec projx_local;
    jl_project(proj32_local, sy[i], pp->n[i], jlmat1, jlmat2);
    for(size_t coeff=0;coeff<256;coeff++){
      proj64_local[coeff] = proj32_local[coeff];
    }
    polxvec_init(projx_local, 256/N, 0);
    polxvec_fromint64vec(projx_local, proj64_local, 256/N, 1, 0);
    polxvec_bindec(&sout[projoff[i]], projx_local, projx_local->len,
                   jlbits[i]);
    polxvec_free(projx_local);
  }
  polxvec_init(sxl_jl, pp->len[LAB_PROJ], 1);
  polxvec_frompolyvec(sxl_jl, &sout[pp->off[LAB_PROJ]], 1, sxl_jl->len, 0.5);

  timing_end(&time);
  timing_print(&time, 3);

  // Outer commitment Up

  polxvec_init(Up, pp->kappa[2], 1);
  polxvec_init(projd, pp->len[LAB_PROJ] + pp->len[LAB_RAND2], 1);
  polxvec_frompolyvec(projd, &sout[pp->off[LAB_PROJ]], 1,pp->len[LAB_PROJ],0.5);
  if(pp->len[LAB_RAND2] > 0){
    polxvec_init_subvec(projd, projd, pp->len[LAB_PROJ], 1, pp->len[LAB_RAND2]);
    polxvec_frompolyvec(projd, &sout[pp->off[LAB_RAND2]], 1, projd->len, 0.5);
    polxvec_init_subvec(projd, projd, 0, 1, 0);
  }
  commit(Up, projd);
  polxvec_decompose(&sout[pp->off[LAB_OUTCOM]+pp->fu], Up, 1, pp->fu, pp->bu);
  polzvec_frompolxvec(pi->m[1], Up, 0, 1, pp->kappa[2]);
  outcom_clear(*pi->m[1]);

  update_hash_polz(ost->h, pi->m[1], pp->kappa[2]);
  
  // Collapse projection matrices
  
  timing_start(&time, "Collapse JL");

  ldd_collapse_jl(jlmat_agg, phi_jlproj, pp->nmax, jlbits_max, jlmat1, jlmat2,
                  ost->h);

  timing_end(&time);
  timing_print(&time, 3);
  
  // Aggregate Zq constraints

  timing_start(&time, "Zq aggregation");

  ldd_aggregate_zq(zqagg, phi_jlproj_full, phi_jlproj, jlbits, pp, ist,
                   chalz_jl, ost->h, 1);

  timing_end(&time);
  timing_print(&time, 3);

  // Collapsed JL matrices times witness

  timing_start(&time, "Collapsed JL matrices times witness");
  const char *jl_sprod_refresh_value = getenv("REVAC_LABRADOR_JL_SPROD_REFRESH");
  const int jl_sprod_refresh = jl_sprod_refresh_value == NULL ||
                               jl_sprod_refresh_value[0] == '\0' ||
                               jl_sprod_refresh_value[0] != '0';
  const int jl_sprod_threads = ldd_kernel_threads(
      "REVAC_LABRADOR_LDD_JL_SPROD_THREADS", LIFTS * pp->r);

  #ifdef _OPENMP
  #pragma omp parallel for collapse(2) num_threads(jl_sprod_threads) if(revac_labrador_ldd_kernel_parallel_enabled("REVAC_LABRADOR_LDD_JL_SPROD_PARALLEL") && jl_sprod_threads > 1 && LIFTS * pp->r >= 2 && !omp_in_parallel()) schedule(static)
  #endif
  for(i=0;i<LIFTS;i++){
    for(j=0;j<pp->r;j++){
      polxvec jl_view;
      polxvec_init_subvec2(jl_view, jlmat_agg[i], 0, 1, pp->n[j]);
      polxvec_init(jl_sprod_sx[i][j], 1, 1);
      polxvec_sprod(jl_sprod_sx[i][j], jl_view, sxq[j]);
      if(jl_sprod_refresh)
        polxvec_exact_refresh(jl_sprod_sx[i][j]);
    }
  }

  timing_end(&time);
  timing_print(&time, 3);

  // Lift constraints

  timing_start(&time, "Lift constraints");
  const int lift_threads = ldd_kernel_threads(
      "REVAC_LABRADOR_LDD_LIFT_THREADS", LIFTS);

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(lift_threads) if(revac_labrador_ldd_kernel_parallel_enabled("REVAC_LABRADOR_LDD_LIFT_PARALLEL") && lift_threads > 1 && LIFTS >= 2 && !omp_in_parallel()) schedule(static)
  #endif
  for(i=0;i<LIFTS;i++){
    polxvec liftx_local;
    size_t jj, kk;
    polxvec_init(liftx_local, 1, 1);
    sparsecnst_eval(liftx_local, zqagg[i], sxq_old, sxl);
    polxvec_sprod_add(liftx_local, phi_jlproj_full[i], sxl_jl);

    for(jj=0;jj<LIFTS;jj++){
      for(kk=0;kk<pp->r;kk++){
        polxvec_scale_add(liftx_local, jl_sprod_sx[jj][kk],
                          chalz_jl[i][jj*pp->r + kk]);
      }
    }
    polxvec_bindec(&sout[pp->off[LAB_LIFT] + i * LOGQ], liftx_local,
                    liftx_local->len, LOGQ);
    polxvec_free(liftx_local);
  }

  timing_end(&time);
  timing_print(&time, 3);

  // Outer commitment Ub

  polxvec_init(Ub, pp->kappa[2], 1);
  polxvec_init(liftd, pp->len[LAB_LIFT] + pp->len[LAB_RAND3], 1);
  polxvec_frompolyvec(liftd, &sout[pp->off[LAB_LIFT]], 1,pp->len[LAB_LIFT],0.5);
  if(pp->len[LAB_RAND3] > 0){
    polxvec_init_subvec(liftd, liftd, pp->len[LAB_LIFT], 1, pp->len[LAB_RAND3]);
    polxvec_frompolyvec(liftd, &sout[pp->off[LAB_RAND3]], 1, liftd->len, 0.5);
    polxvec_init_subvec(liftd, liftd, 0, 1, 0);
  }
  commit(Ub, liftd);
  polxvec_decompose(&sout[pp->off[LAB_OUTCOM]+2*pp->fu], Ub, 1, pp->fu, pp->bu);
  polzvec_frompolxvec(pi->m[2], Ub, 0, 1, pp->kappa[2]);
  outcom_clear(*pi->m[2]);

  update_hash_polz(ost->h, pi->m[2], pp->kappa[2]);

  // Aggregate Rq constraints

  timing_start(&time, "Rq aggregation");

  polxvec_init(chalx_zq, LIFTS, 1);
  ldd_aggregate_rq(finalcnst, chalx_zq, pp, ist, zqagg, ost->h);

  timing_end(&time);
  timing_print(&time, 3);

  off = 0;
  for(i=0;i<pp->r;i++){
    polxvec_init_subvec2(phi[i], finalcnst->lin->phi[0], off, 1, pp->n[i]);
    off += pp->n[i];
  }

  // Linear garbage hij

  timing_start(&time, "Linear garbage");

  ldd_ling(&sout[pp->off[LAB_LING]], pp->r, pp->fu, pp->bu, sxq, phi, chalx_zq, 
           jl_sprod_sx, chalz_jl);

  timing_end(&time);
  timing_print(&time, 3);

  // Middle commitment u2

  polxvec_init(u2x, pp->kappa[1], 1);
  polxvec_init(hd, pp->len[LAB_LING], 1);
  polxvec_frompolyvec(hd, &sout[pp->off[LAB_LING]], 1, pp->len[LAB_LING],
                      WIDTHMOD(pp->bu));
  commit(u2x, hd);
  polxvec_bindec(&sout[pp->off[LAB_U2]], u2x, u2x->len, LOGQ);

  // Outer commitment U2

  polxvec_init(U2, pp->kappa[2], 1);
  polxvec_init(u2d, pp->len[LAB_U2] + pp->len[LAB_RAND4], 1);
  polxvec_frompolyvec(u2d, &sout[pp->off[LAB_U2]], 1, pp->len[LAB_U2], 0.5);
  if(pp->len[LAB_RAND4] > 0){
    polxvec_init_subvec(u2d, u2d, pp->len[LAB_U2], 1, pp->len[LAB_RAND4]);
    polxvec_frompolyvec(u2d, &sout[pp->off[LAB_RAND4]], 1, u2d->len, 0.5);
    polxvec_init_subvec(u2d, u2d, 0, 1, 0);
  }
  commit(U2, u2d);
  polxvec_decompose(&sout[pp->off[LAB_OUTCOM]+3*pp->fu], U2, 1, pp->fu, pp->bu);
  polzvec_frompolxvec(pi->m[3], U2, 0, 1, pp->kappa[2]);
  outcom_clear(*pi->m[3]);

  update_hash_polz(ost->h, pi->m[3], pp->kappa[2]);

  // Amortization

  timing_start(&time, "Amortization");

  chalx_amortize = _aligned_alloc(64, pp->r * sizeof(polx));
  sample_chalx_amortize(chalx_amortize, pp->r, ost->h);
  
  polxvec_init(zx, pp->nmax, 1);
  polxvec_setzero(zx, 0, 1, pp->nmax);
  for(i=0;i<pp->r;i++){
    polxvec_init_subvec(zx, zx, 0, 1, pp->n[i]);
    polxvec_polx_mul_add(zx, chalx_amortize[i], sxq[i]);
  }
  polxvec_init_subvec(zx, zx, 0, 1, 0);
  polxvec_decompose(&sout[pp->off[LAB_Z]], zx, zx->len, pp->fz, pp->bz);
  if(ldd_residual_diagnostics_enabled()){
    polxvec z_reconstructed;
    polxvec_init(z_reconstructed, zx->len, 1);
    polxvec_reconstruct(z_reconstructed, &sout[pp->off[LAB_Z]], zx->len,
                        pp->fz, pp->bz);
    polxvec_sub(z_reconstructed, z_reconstructed, zx);
    if(!polxvec_iszero(z_reconstructed))
      fprintf(stderr,
              "[Labrador z 分解诊断] amortized z 分解重构不一致\n");
    polxvec_free(z_reconstructed);
  }

  timing_end(&time);
  timing_print(&time, 3);

  // Generate constraints for the verifier's checks

  timing_start(&time, "Add checks");

  ldd_addchecks(ost, pp, pi, chalz_jl, chalx_amortize, chalx_zq, phi, jlmat_agg, 
                phi_jlproj_full, finalcnst);
    ldd_diagnose_ling_components(ost->rqcnst->sparse[1], pp, sout,
                               chalx_amortize, sxq, phi, chalx_zq,
                               jl_sprod_sx, chalz_jl, jlmat_agg);

  timing_end(&time);
  timing_print(&time, 3);

  timing_start(&time, "Free");

  for(i=0;i<LIFTS;i++){
    sparsecnst_free(zqagg[i]);
    polxvec_free(jlmat_agg[i]);
    polxvec_free(phi_jlproj[i]);
    polxvec_free(phi_jlproj_full[i]);
    for(j=0;j<pp->r;j++){
      polxvec_free(jl_sprod_sx[i][j]);
    }
  }
  sparsecnst_free(finalcnst);

  polxvec_free(sxl);
  polxvec_free(sxl_jl);
  polxvec_free(u1x);
  polxvec_free(tgd);
  polxvec_free(U1);
  polxvec_free(u1d);
  polxvec_free(Up);
  polxvec_free(projd);
  polxvec_free(liftd);
  polxvec_free(Ub);
  polxvec_free(u2x);
  polxvec_free(hd);
  polxvec_free(chalx_zq);
  polxvec_free(U2);
  polxvec_free(u2d);
  polxvec_free(zx);
  free(jlmat1);
  free(chalx_amortize);

  timing_end(&time);
  timing_print(&time, 3);
}

void ldd_reduce(
  statement ost,
  const statement ist, 
  const lab_proof pi,
  const lab_params pp
)
{
  size_t i, jlbits[pp->r], off, jlbits_max;
  uint8_t *jlmat1, *jlmat2;
  int64_t chalz_jl[LIFTS][LIFTS*pp->r];
  polx *chalx_amortize;
  polxvec phi[pp->r], chalx_zq, jlmat_agg[LIFTS], phi_jlproj[LIFTS];
  polxvec phi_jlproj_full[LIFTS];
  sparsecnst zqagg[LIFTS], finalcnst;
  timing time;
  const int profile_detail = revac_labrador_profile_detail_enabled();
  double profile_started = 0.0;

  if (profile_detail) profile_started = revac_labrador_profile_seconds();
  lab_statement_init(ost, ist, pp);
  lab_comkey_init(pp);

  jlbits_max = 0;
  for(i=0;i<pp->r;i++){
    jlbits[i] = ceil(log2(JL_INF_MULT*sqrt(pp->normsq[i])));
    jlbits_max = MAX(jlbits_max, jlbits[i]);
  }
  if (profile_detail) {
    fprintf(stderr, "[Labrador 原生细分] LDD 输出 statement 初始化：%.6f 秒\n",
            revac_labrador_profile_seconds() - profile_started);
  }

  // Hash first message

  if (profile_detail) profile_started = revac_labrador_profile_seconds();
  update_hash_polz(ost->h, pi->m[0], pp->kappa[2]);

  // Sample JL matrices

  jl_sample_mat(&jlmat1, &jlmat2, ost->h, pp->nmax);

  // Hash second message

  update_hash_polz(ost->h, pi->m[1], pp->kappa[2]);
  if (profile_detail) {
    fprintf(stderr, "[Labrador 原生细分] LDD transcript + JL 采样：%.6f 秒\n",
            revac_labrador_profile_seconds() - profile_started);
  }

  // Collapse projection matrices
  
  timing_start(&time, "Collapse JL");
  if (profile_detail) profile_started = revac_labrador_profile_seconds();

  ldd_collapse_jl(jlmat_agg, phi_jlproj, pp->nmax, jlbits_max, jlmat1, jlmat2,
                  ost->h);
  if (profile_detail) {
    fprintf(stderr, "[Labrador 原生细分] LDD Collapse JL：%.6f 秒\n",
            revac_labrador_profile_seconds() - profile_started);
  }

  timing_end(&time);
  timing_print(&time, 3);

  // Aggregate Zq constraints

  timing_start(&time, "Zq aggregation");
  if (profile_detail) profile_started = revac_labrador_profile_seconds();

  ldd_aggregate_zq(zqagg, phi_jlproj_full, phi_jlproj, jlbits, pp, ist,
                   chalz_jl, ost->h, 0);
  if (profile_detail) {
    fprintf(stderr, "[Labrador 原生细分] LDD Zq aggregation：%.6f 秒\n",
            revac_labrador_profile_seconds() - profile_started);
  }

  timing_end(&time);
  timing_print(&time, 3);

  // Hash third message

  if (profile_detail) profile_started = revac_labrador_profile_seconds();
  update_hash_polz(ost->h, pi->m[2], pp->kappa[2]);
  if (profile_detail) {
    fprintf(stderr, "[Labrador 原生细分] LDD 第三消息哈希：%.6f 秒\n",
            revac_labrador_profile_seconds() - profile_started);
  }

  // Aggregate Rq constraints

  timing_start(&time, "Rq aggregation");
  if (profile_detail) profile_started = revac_labrador_profile_seconds();

  polxvec_init(chalx_zq, LIFTS, 1);
  ldd_aggregate_rq(finalcnst, chalx_zq, pp, ist, zqagg, ost->h);
  if (profile_detail) {
    fprintf(stderr, "[Labrador 原生细分] LDD Rq aggregation：%.6f 秒\n",
            revac_labrador_profile_seconds() - profile_started);
  }

  timing_end(&time);
  timing_print(&time, 3);

  off = 0;
  for(i=0;i<pp->r;i++){
    polxvec_init_subvec2(phi[i], finalcnst->lin->phi[0], off, 1, pp->n[i]);
    off += pp->n[i];
  }

  // Hash fourth message

  if (profile_detail) profile_started = revac_labrador_profile_seconds();
  update_hash_polz(ost->h, pi->m[3], pp->kappa[2]);

  // Sample challenges for amortization

  chalx_amortize = _aligned_alloc(64, pp->r * sizeof(polx));
  sample_chalx_amortize(chalx_amortize, pp->r, ost->h);
  if (profile_detail) {
    fprintf(stderr, "[Labrador 原生细分] LDD 第四消息哈希与挑战：%.6f 秒\n",
            revac_labrador_profile_seconds() - profile_started);
  }

  // Generate constraints for the verifier's checks

  timing_start(&time, "Add checks");
  if (profile_detail) profile_started = revac_labrador_profile_seconds();

  ldd_addchecks(ost, pp, pi, chalz_jl, chalx_amortize, chalx_zq, phi, jlmat_agg, 
                phi_jlproj_full, finalcnst);
  if (profile_detail) {
    fprintf(stderr, "[Labrador 原生细分] LDD Add checks：%.6f 秒\n",
            revac_labrador_profile_seconds() - profile_started);
  }

  timing_end(&time);
  timing_print(&time, 3);

  timing_start(&time, "Free");
  if (profile_detail) profile_started = revac_labrador_profile_seconds();

  for(i=0;i<LIFTS;i++){
    sparsecnst_free(zqagg[i]);
    polxvec_free(jlmat_agg[i]);
    polxvec_free(phi_jlproj[i]);
    polxvec_free(phi_jlproj_full[i]);
  }
  sparsecnst_free(finalcnst);

  polxvec_free(chalx_zq);
  free(jlmat1);
  free(chalx_amortize);
  if (profile_detail) {
    fprintf(stderr, "[Labrador 原生细分] LDD Free：%.6f 秒\n",
            revac_labrador_profile_seconds() - profile_started);
  }

  timing_end(&time);
  timing_print(&time, 3);
}
