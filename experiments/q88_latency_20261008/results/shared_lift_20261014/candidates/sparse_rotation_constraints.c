#include <stdio.h>
#include "malloc.h"
#include "polx.h"
#include "polz.h"
#include "constraints.h"
#include "comkey.h"
#include "revac_profile.h"

static size_t q88_sparse_rotation_blocks;
__attribute__((destructor)) static void q88_sparse_rotation_report(void) {
  const char *mode=getenv("REVAC_Q88_ROTATION_SPARSE");
  if(mode && mode[0]!='0') fprintf(stderr,"[q88 sparse rotation] optimized_blocks=%zu\n",q88_sparse_rotation_blocks);
}

static int residual_diagnostics_enabled(void){
  const char *value = getenv("REVAC_LABRADOR_RESIDUAL_DIAG");
  return value != NULL && value[0] != '\0' && value[0] != '0';
}

static void sparsecnst_print_residual(const polxvec ev){
  const int64_t q = ((int64_t)1 << LOGQ) - QOFF;

  for(size_t vector_index=0;vector_index<ev->len;vector_index++){
    polxvec view;
    polz approximate, exact;
    zz approximate_coefficient, exact_coefficient;

    polxvec_init_subvec2(view, ev, vector_index, 1, 1);
    polzvec_frompolxvec(&approximate, view, 0, 1, 1);
    polz_center(approximate);
    polzvec_exact_frompolxvec(&exact, view, 0, 1, 1);
    for(size_t coefficient=0;coefficient<N;coefficient++){
      int64_t approximate_value, exact_value, delta;

      polz_getcoeff(approximate_coefficient, approximate, coefficient);
      polz_getcoeff(exact_coefficient, exact, coefficient);
      approximate_value = int64_fromzz(approximate_coefficient);
      exact_value = int64_fromzz(exact_coefficient);
      if(approximate_value == 0 && exact_value == 0)
        continue;
      delta = approximate_value - exact_value;
      fprintf(stderr,
              "[Labrador residual 诊断] vector=%zu coeff=%zu "
              "approx=%lld exact=%lld delta=%lld delta/q=%lld "
              "delta%%q=%lld\n",
              vector_index, coefficient,
              (long long)approximate_value, (long long)exact_value,
              (long long)delta, (long long)(delta / q),
              (long long)(delta % q));
      return;
    }
  }
  fprintf(stderr,
          "[Labrador residual 诊断] approximate/exact residual 均为零\n");
}

/*
  Initializes quad to have length len, while allocating enough memory to support
  up to maxlen elements.
*/
void quadfunc_init(quadfunc quad, size_t len, size_t maxlen){
  quad->len = len;
  if(maxlen == 0){
    quad->rows = NULL;
    quad->cols = NULL;
    quad->coeffs = NULL;
  }
  else{
    quad->rows = _malloc(2*maxlen*sizeof(size_t));
    quad->cols = &quad->rows[maxlen];
    quad->coeffs = _aligned_alloc(64, maxlen*sizeof(polx));
  }
}

static void quadfunc_copy(quadfunc out, const quadfunc in){
  size_t i;

  quadfunc_init(out, in->len, in->len);
  for(i=0;i<in->len;i++){
    out->rows[i] = in->rows[i];
    out->cols[i] = in->cols[i];
    polx_copy(out->coeffs[i], in->coeffs[i]);
  }
}

static void quadfunc_copy2(quadfunc out, const quadfunc in, size_t maxlen){
  size_t i;

  quadfunc_init(out, in->len, maxlen);
  for(i=0;i<in->len;i++){
    out->rows[i] = in->rows[i];
    out->cols[i] = in->cols[i];
    polx_copy(out->coeffs[i], in->coeffs[i]);
  }
}

void quadfunc_free(quadfunc quad){
  free(quad->rows);
  free(quad->coeffs);
  quad->len = 0;
}

/* ev += quad(sx) */
void quadfunc_eval_add(
  polxvec ev, 
  const quadfunc quad, 
  const polxvec sx[]
)
{
  size_t i, j, k, len, off, chunk;
  polxvec sxj, sxk, tmp;
  polxvec_init(tmp, ev->len, 1);

  for(i=0;i<quad->len;i++) {
    j = quad->rows[i];
    k = quad->cols[i];
    len = MIN(sx[j]->len, sx[k]->len);

    off = 0;
    while(len > 0){
      chunk = MIN(1000 * ev->len, len);
      polxvec_init_subvec2(sxj, sx[j], off, 1, chunk);
      polxvec_init_subvec2(sxk, sx[k], off, 1, chunk);
      polxvec_sprod_extension(tmp, sxj, sxk);
      polxvec_refresh(tmp);
      polxvec_polx_mul_add(ev, quad->coeffs[i], tmp);
      polxvec_refresh(ev); //XXX required ?
      off += chunk;
      len -= chunk;
    }
  }
  polxvec_free(tmp);
}

/*
  r += c*a, where c->len == 1. The quadfunc r is assumed to have enough 
  allocated memory for the update.
*/
static void quadfunc_polxvec_mul_add(
  quadfunc r, 
  const polxvec c, 
  const quadfunc a
)
{
  size_t i,j;
  polxvec rj, ai;

  for(i=0;i<a->len;i++) {
    for(j=0;j<r->len;j++) {
      if(r->rows[j] == a->rows[i] && r->cols[j] == a->cols[i]) {
        polxvec_init_frompolx(rj, r->coeffs[j]);
        polxvec_init_frompolx(ai, a->coeffs[i]);
        polxvec_mul_add(rj, c, ai);
        break;
      }
    }
    if(j == r->len) {
      r->rows[j] = a->rows[i];
      r->cols[j] = a->cols[i];
      r->coeffs[j]->width = 0;
      polxvec_init_frompolx(rj, r->coeffs[j]);
      polxvec_init_frompolx(ai, a->coeffs[i]);
      polxvec_mul(rj, c, ai);
      r->len += 1;
    }
  }
}

/*
  r += s*a. The quadfunc r is assumed to have enough allocated memory for the 
  update.
*/
static void quadfunc_scale_add(quadfunc r, const quadfunc a, int64_t s) {
  size_t i,j;

  for(i=0;i<a->len;i++) {
    for(j=0;j<r->len;j++) {
      if(r->rows[j] == a->rows[i] && r->cols[j] == a->cols[i]) {
        polx_scale_add(r->coeffs[j],a->coeffs[i],s);
        break;
      }
    }
    if(j == r->len) {
      r->rows[j] = a->rows[i];
      r->cols[j] = a->cols[i];
      polx_scale(r->coeffs[j],a->coeffs[i],s);
      r->len += 1;
    }
  }
}

/*
  Initializes lin with the specified rank and nparts, while allocating enough
  memory to support up to maxnparts.
*/
void linfunc_init(linfunc lin, size_t rank, size_t nparts, size_t maxnparts){
  lin->rank = rank;
  lin->nparts = nparts;
  if(maxnparts > 0){
    lin->off = _malloc(maxnparts*(sizeof(size_t) + sizeof(polxvec)));
    lin->phi = (polxvec *) &lin->off[maxnparts];
  }
  else{
    lin->off = NULL;
    lin->phi = NULL;
  }
}

static void linfunc_copy(linfunc out, const linfunc in){
  size_t i;

  linfunc_init(out, in->rank, in->nparts, in->nparts);
  for(i=0;i<in->nparts;i++){
    out->off[i] = in->off[i];
    polxvec_init(out->phi[i], in->phi[i]->len, 1);
    polxvec_copy(out->phi[i], in->phi[i]);
  }
}

void linfunc_free(linfunc lin){
  size_t i;
  for(i=0;i<lin->nparts;i++){
    polxvec_free(lin->phi[i]);
  }
  free(lin->off);
  lin->nparts = 0;
}

/* ev += lin(sx) */
void linfunc_eval_add(polxvec ev, const linfunc lin, const polxvec sx){
  size_t i;
  polxvec sx_sv;

  for(i=0;i<lin->nparts;i++){
    polxvec_init_subvec2(sx_sv, sx, lin->off[i], 1, lin->phi[i]->len);
    polxvec_sprod_extension_add(ev, lin->phi[i], sx_sv);
  }
}

/*
  Initializes a zero sparse constraint, the vector b is allocated and set to 0.
*/
void sparsecnst_init(sparsecnst cnst, size_t rank){
  cnst->borrowed = 0;
  cnst->transcript_lin = NULL;
  cnst->transcript_lin_len = 0;
  cnst->transcript_quad = NULL;
  cnst->transcript_quad_len = 0;
  quadfunc_init(cnst->quad, 0, 0);
  linfunc_init(cnst->lin, rank, 0, 0);
  polxvec_init(cnst->b, rank, 1);
  polxvec_setzero(cnst->b, 0, 1, rank);
}

/* Create an immutable, non-owning view.  This is used only while the source
 * statement is alive and avoids deep-copying a large static constraint
 * prefix into a short-lived Dachshund statement. */
void sparsecnst_borrow(sparsecnst out, const sparsecnst in){
  memcpy(out, in, sizeof(*out));
  out->borrowed = 1;
}

void sparsecnst_copy(sparsecnst out, const sparsecnst in){
  sparsecnst_init(out, in->b->len);
  quadfunc_copy(out->quad, in->quad);
  linfunc_copy(out->lin, in->lin);
  polxvec_copy(out->b, in->b);
}

void sparsecnst_copy2(sparsecnst out, const sparsecnst in, size_t maxlen){
  sparsecnst_init(out, in->b->len);
  quadfunc_copy2(out->quad, in->quad, maxlen);
  linfunc_copy(out->lin, in->lin);
  polxvec_copy(out->b, in->b);
}

void sparsecnst_free(sparsecnst cnst){
  if(cnst->borrowed){
    cnst->borrowed = 0;
    return;
  }
  quadfunc_free(cnst->quad);
  linfunc_free(cnst->lin);
  polxvec_free(cnst->b);
  free(cnst->transcript_lin);
  free(cnst->transcript_quad);
  cnst->transcript_lin = NULL;
  cnst->transcript_lin_len = 0;
  cnst->transcript_quad = NULL;
  cnst->transcript_quad_len = 0;
}

void sparsecnst_refresh(sparsecnst cnst){
  size_t i;

  for(i=0;i<cnst->quad->len;i++){
    polx_refresh(cnst->quad->coeffs[i]);
  }
  for(i=0;i<cnst->lin->nparts;i++){
    polxvec_refresh(cnst->lin->phi[i]);
  }
  polxvec_refresh(cnst->b);
}

/* 
  ev = cnst(sx)
  The parameters sxq and sxl are two views of the same witness vector sx:
    - sxq: array of the r polxvec that form the witness sx.
    - sxl: concatenation of all the polxvec that form the witness sx.
*/
void sparsecnst_eval(
  polxvec ev,
  const sparsecnst cnst,
  const polxvec sxq[],
  const polxvec sxl
)
{
  polxvec_setzero(ev, 0, 1, ev->len);
  
  quadfunc_eval_add(ev, cnst->quad, sxq);
  linfunc_eval_add(ev, cnst->lin, sxl);
  polxvec_sub(ev, ev, cnst->b);
}

/* 
  When full=1, sparsecnst_check returns 1 if cnst(sx) = 0, and 0 
  otherwise. When full=0, sparse_check returns 1 if ctcoef(cnst(sx)) = 0,
  and 0 otherwise.
  The parameters sxq and sxl are two views of the same witness vector sx:
    - sxq: array of the r polxvec that form the witness sx.
    - sxl: concatenation of all the polxvec that form the witness sx.
*/
int sparsecnst_check(
  const sparsecnst cnst,
  const polxvec sxq[],
  const polxvec sxl,
  int full
)
{
  int check;
  polxvec ev;

  polxvec_init(ev, cnst->lin->rank, 1);
  sparsecnst_eval(ev, cnst, sxq, sxl);
  if(full){
    check = polxvec_iszero(ev);
  }
  else{
    check = polxvec_iszero_constcoeff(ev, 0);
  }
  if(!check && residual_diagnostics_enabled())
    sparsecnst_print_residual(ev);
  polxvec_free(ev);
  return check;
}

/*
  The polxvec b is seen as a vector of polynomials in a higher degree ring,
  where the degree is the smallest power of 2 greater or equal to the rank,
  and the rank is a->len. Each such polynomial is expanded into a rotation
  matrix, where the number of rows is equal to the rank, and the number of
  columns is equal to the degree, except for the last polynomial, where the
  number of columns will be truncated to match the length of r. Each such
  rotation matrix is aggregated using the polynomial challenges in a and added
  to the corresponding subvector in r.
*/
static void rotation_aggregate_internal(
  polxvec r, 
  const polxvec a, 
  const polxvec b,
  int batch, int sparse
)
{
  size_t rank, deg, len, nblocks;
  int threads;
  const char *lazy_value = getenv("REVAC_LABRADOR_ROTATION_LAZY_SCRATCH");
  const int lazy_scratch = lazy_value == NULL || lazy_value[0] != '0';

  rank = a->len;
  deg = next2power(rank);
  len = r->len;
  nblocks = (len + deg - 1) / deg;
  threads = revac_labrador_rotation_threads(nblocks);

  /* Each block reads the same challenge vector but disjoint slices of b and
   * r.  All mutable scratch objects live inside the parallel region. */
  #ifdef _OPENMP
  #pragma omp parallel \
      num_threads(threads) \
      if(revac_labrador_rotation_parallel_enabled() && \
         threads > 1 && nblocks >= 4 && !omp_in_parallel())
  #endif
  {
    polxvec r_sv, a_sv, b_sv, b0_sv, tutl, tu, tl, tu_sv, tl_sv, tmp, tmp_sv;
    polx x;
    size_t block, row, col;
    int scratch_initialized = 0;

    polx_monomial(x, 1, 1);
    if(!lazy_scratch){
      polxvec_init(tutl, 2*deg, 1);
      polxvec_init(tmp, deg, 1);
      scratch_initialized = 1;
    }

    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for(block=0;block<nblocks;block++){
      const size_t off = block * deg;
      const size_t ncol = MIN(deg, len - off);
      int scalar_block = b->widths != NULL;

      polxvec_init_subvec2(b_sv, b, off, 1, deg);
      polxvec_init_subvec2(r_sv, r, off, 1, ncol);

      /* A surprisingly common case in generated statements is an
       * extension-ring coefficient which is a scalar: only coordinate zero
       * of this degree-sized block is non-zero.  Expanding its full rotation
       * matrix performs deg^2 pointwise products even though the result is
       * simply r[j] += b[0] * a[j].  Width zero is an exact invariant for a
       * zero polx in this representation, so this test does not inspect or
       * approximate coefficients and does not alter the transcript.  The
       * q96 rank-32 aggregation equalities use precisely this layout. */
      for(col=1;scalar_block && col<deg;col++){
        scalar_block = b->widths[(off + col) * b->stride] == 0;
      }
      if(scalar_block){
        polxvec_init_subvec2(b0_sv, b, off, 1, 1);
        polxvec_init_subvec2(a_sv, a, 0, 1, ncol);
        polxvec_mul_add(r_sv, b0_sv, a_sv);
        continue;
      }

      if(!scratch_initialized){
        polxvec_init(tutl, 2*deg, 1);
        polxvec_init(tmp, deg, 1);
        scratch_initialized = 1;
      }

      polxvec_setzero(tutl, 0, 1, 2*ncol);
      polxvec_init_subvec2(tu, tutl, 0, 1, ncol);
      polxvec_init_subvec2(tl, tutl, ncol, 1, ncol);

      size_t positions[4], nonzero=0;
      if(sparse && b->widths && deg>=16) {
        for(size_t position=0;position<deg;position++) {
          if(b->widths[(off+position)*b->stride]==0) continue;
          if(nonzero==4) { nonzero=5;break; }
          positions[nonzero++]=position;
        }
      } else nonzero=5;
      if(nonzero<=4) {
        /* Preserve the original per-output row order, and its exact width
         * arithmetic. Only products with publicly known zero coefficients
         * are omitted; no private witness values are inspected. */
        for(row=0;row<rank;row++) {
          polxvec_init_subvec2(a_sv,a,row,1,1);
          for(size_t entry=0;entry<nonzero;entry++) {
            const size_t position=positions[entry];
            const size_t output=position<=row ? row-position : deg+row-position;
            if(output>=ncol) continue;
            polxvec_init_subvec2(b0_sv,b_sv,position,1,1);
            polxvec_init_subvec2(tmp_sv,tmp,0,1,1);
            polxvec_mul(tmp_sv,a_sv,b0_sv);
            if(position<=row) polxvec_init_subvec2(tl_sv,tl,output,1,1);
            else polxvec_init_subvec2(tl_sv,tu,output,1,1);
            polxvec_add(tl_sv,tl_sv,tmp_sv);
          }
        }
        polxvec_polx_mul_add(r_sv,x,tu);
        polxvec_add(r_sv,r_sv,tl);
        __atomic_fetch_add(&q88_sparse_rotation_blocks,1,__ATOMIC_RELAXED);
        continue;
      }

      for(row=0;row<rank;row++){
        polxvec_init_subvec2(tmp_sv, tmp, 0, 1, deg);
        polxvec_init_subvec2(a_sv, a, row, 1, 1);
        polxvec_mul(tmp_sv, a_sv, b_sv);

        if(batch){
          /* Same per-coordinate additions and row order as the scalar
           * loop, using two reversed views instead of ncol tiny calls. */
          const size_t lower = MIN(ncol, row + 1);
          polxvec_init_subvec2(tl_sv, tl, 0, 1, lower);
          polxvec_init_subvec2(tmp_sv, tmp, row, -1, lower);
          polxvec_add(tl_sv, tl_sv, tmp_sv);
          if(lower < ncol){
            polxvec_init_subvec2(tu_sv, tu, lower, 1, ncol-lower);
            polxvec_init_subvec2(tmp_sv, tmp, deg-1, -1, ncol-lower);
            polxvec_add(tu_sv, tu_sv, tmp_sv);
          }
        }
        else for(col=0;col<ncol;col++){
          if(col<=row){
            // tl[col] += tmp[row-col]
            polxvec_init_subvec2(tl_sv, tl, col, 1, 1);
            polxvec_init_subvec2(tmp_sv, tmp, row-col, 1, 1);
            polxvec_add(tl_sv, tl_sv, tmp_sv);
          }
          else{
            // tu[col] += tmp[deg+row-col]
            polxvec_init_subvec2(tu_sv, tu, col, 1, 1);
            polxvec_init_subvec2(tmp_sv, tmp, deg+row-col, 1, 1);
            polxvec_add(tu_sv, tu_sv, tmp_sv);
          }
        }
      }
      polxvec_polx_mul_add(r_sv, x, tu);
      polxvec_add(r_sv, r_sv, tl);
    }

    if(scratch_initialized){
      polxvec_free(tutl);
      polxvec_free(tmp);
    }
  }
}


void polxvec_rotation_aggregate_add(polxvec r, const polxvec a,
                                    const polxvec b)
{
  const char *mode = getenv("REVAC_LABRADOR_ROTATION_BATCH");
  const char *audit = getenv("REVAC_LABRADOR_ROTATION_CHECK");
  const int batch = mode != NULL && mode[0] == '1';
  const char *sparse_mode=getenv("REVAC_Q88_ROTATION_SPARSE");
  const int sparse=sparse_mode && (sparse_mode[0]=='1' || sparse_mode[0]=='2');
  if((audit != NULL && audit[0] == '1') || (sparse_mode && sparse_mode[0]=='2')){
    polxvec reference;
    polxvec_init(reference, r->len, r->widths != NULL);
    polxvec_copy(reference, r);
    rotation_aggregate_internal(reference, a, b, 0, 0);
    rotation_aggregate_internal(r, a, b, 1, sparse);
    for(size_t j=0;j<r->len;j++){
      for(size_t k=0;k<K;k++){
        for(size_t i=0;i<N;i++){
          if(r->proj[k][j*r->stride]->c[i] !=
             reference->proj[k][j]->c[i]){
            fprintf(stderr, "rotation batch coefficient mismatch\n");
            abort();
          }
        }
      }
      if(r->widths && r->widths[j*r->stride] != reference->widths[j]){
        fprintf(stderr, "rotation batch width mismatch\n");
        abort();
      }
    }
    polxvec_free(reference);
    return;
  }
  rotation_aggregate_internal(r, a, b, batch, sparse);
}

/*
  out += c*in, where c->len == in->lin->rank.
*/
static void sparsecnst_polxvec_mul_add(
  sparsecnst out,
  const polxvec c,
  const sparsecnst in,
  double profile_parts[3]
)
{
  size_t i;
  polxvec outphi;
  const int profile = profile_parts != NULL;
  double started = 0.0;

  if(profile) started = revac_labrador_profile_seconds();
  if(in->lin->rank == 1){
    quadfunc_polxvec_mul_add(out->quad, c, in->quad);
  }
  if(profile) {
    profile_parts[0] += revac_labrador_profile_seconds() - started;
    started = revac_labrador_profile_seconds();
  }

  polxvec_sprod_add(out->b, c, in->b);
  if(profile) {
    profile_parts[1] += revac_labrador_profile_seconds() - started;
    started = revac_labrador_profile_seconds();
  }

  for(i=0;i<in->lin->nparts;i++){
    polxvec_init_subvec2(outphi, out->lin->phi[0], in->lin->off[i], 1,
                        in->lin->phi[i]->len);
    if(in->lin->rank == 1){
      polxvec_mul_add(outphi, c, in->lin->phi[i]);
    }
    else{
      polxvec_rotation_aggregate_add(outphi, c, in->lin->phi[i]);
    }
  }
  if(profile) profile_parts[2] += revac_labrador_profile_seconds() - started;
}

/*
  out += s*in
*/
static void sparsecnst_scale_add(
  sparsecnst out,
  const sparsecnst in,
  int64_t s
)
{
  size_t i;
  polxvec outphi;
  
  quadfunc_scale_add(out->quad, in->quad, s);

  polxvec_scale_add(out->b, in->b, s);

  for(i=0;i<in->lin->nparts;i++){
    polxvec_init_subvec2(outphi, out->lin->phi[0], in->lin->off[i], 1,
                        in->lin->phi[i]->len);
    polxvec_scale_add(outphi, in->lin->phi[i], s);
  }
}

/* 
  Add aggregation of all the input constraints into a single constraint.
    - out: single output constraint, over the base ring, with rank 1, and one
      part such that off[0] is 0 and phi[0]->len is the length of the witness.
    - in: array of input constraints.
    - ncnst: number of input constraints.
    - chal: polynomial challenges, used when full=1. In that case, there are
      as many challenges as the sum of the ranks of the input constraints.
    - chalz: integer scalar challenges, used when full=0. In that case, there 
      are as many challenges as the number of input constraints.
    - full: when set to 1, the constraints hold over the full ring, and are
      aggregated with the polynomial challenges in chal. Otherwise, the
      constraints only need to hold when looking at the constant coefficient, 
      and the aggregation is done with the integer scalars in chalz.
*/
void sparsecnst_aggregate_add(
  sparsecnst out, 
  const sparsecnst *in, 
  size_t ncnst, 
  const polxvec chal, 
  const int64_t *chalz, 
  int full
)
{
  size_t i, off;
  polxvec chal_sv;
  const int profile = revac_labrador_profile_detail_enabled() && full;
  double profile_parts[3] = {0.0, 0.0, 0.0};

  off = 0;
  for(i=0;i<ncnst;i++){
    if(full){
      polxvec_init_subvec2(chal_sv, chal, off, 1, in[i]->lin->rank);
      sparsecnst_polxvec_mul_add(out, chal_sv, in[i],
                                 profile ? profile_parts : NULL);
    }
    else{
      sparsecnst_scale_add(out, in[i], chalz[off]);
    }
    off += in[i]->lin->rank;

  }
  if(profile && ncnst > 0){
    fprintf(stderr,
            "[Labrador sparse 细分] n=%zu quad=%.6f b=%.6f phi=%.6f 秒\n",
            ncnst, profile_parts[0], profile_parts[1], profile_parts[2]);
  }
}

/*
  Initializes cnst with the specified rank, ncom and nphi, while allocating
  enough memory to support up to maxnphi parts.
*/
void comcnst_init(
  comcnst cnst, 
  size_t rank, 
  size_t ncom, 
  size_t nphi, 
  size_t maxnphi
)
{
  size_t i;

  cnst->borrowed = 0;
  cnst->rank = rank;
  polxvec_init(cnst->b, rank, 1);

  cnst->ncom = ncom;
  cnst->comk_off = _malloc(ncom*(3*sizeof(size_t) + sizeof(int64_t)));
  cnst->comw_off = &cnst->comk_off[ncom];
  cnst->comw_len = &cnst->comw_off[ncom];
  cnst->scalar = (int64_t *) &cnst->comw_len[ncom];

  cnst->nphi = nphi;
  cnst->phiw_off = _malloc(maxnphi*(sizeof(size_t) + sizeof(polxvec)));
  cnst->phi = (polxvec *) &cnst->phiw_off[maxnphi];

  for(i=0;i<ncom;i++){
    cnst->scalar[i] = 1;
  }
}

void comcnst_borrow(comcnst out, const comcnst in){
  *out = *in;
  out->borrowed = 1;
}

void comcnst_copy(comcnst out, const comcnst in){
  size_t i;

  comcnst_init(out, in->rank, in->ncom, in->nphi, in->nphi);
  for(i=0;i<in->ncom;i++){
    out->comk_off[i] = in->comk_off[i];
    out->comw_off[i] = in->comw_off[i];
    out->comw_len[i] = in->comw_len[i];
    out->scalar[i] = in->scalar[i];
  }
  for(i=0;i<in->nphi;i++){
    out->phiw_off[i] = in->phiw_off[i];
    polxvec_init(out->phi[i], in->phi[i]->len, 1);
    polxvec_copy(out->phi[i], in->phi[i]);
  }
}

void comcnst_free(comcnst cnst){
  size_t i;
  if(cnst->borrowed){
    memset(cnst, 0, sizeof(*cnst));
    return;
  }
  polxvec_free(cnst->b);
  free(cnst->comk_off);
  
  for(i=0;i<cnst->nphi;i++){
    polxvec_free(cnst->phi[i]);
  }
  free(cnst->phiw_off);

  cnst->ncom = 0;
  cnst->nphi = 0;
}

/*
  ev = cnst(sx)
*/
void comcnst_eval(polxvec ev, const comcnst cnst, const polxvec sx){
  size_t i, j, rank;
  polxvec tmp, comkey_sv, sx_sv, phi_sv;

  rank = cnst->rank;
  polxvec_setzero(ev, 0, 1, rank);
  polxvec_init(tmp, rank, 1);

  for(i=0;i<cnst->ncom;i++){
    polxvec_init_subvec2(comkey_sv, comkey, cnst->comk_off[i], 1,
                        comkey->len - cnst->comk_off[i]);
    polxvec_init_subvec2(sx_sv, sx, cnst->comw_off[i], 1, cnst->comw_len[i]);
    if(cnst->scalar[i] == 1){
      polxvec_sprod_extension_add(ev, comkey_sv, sx_sv);
    }
    else{
      polxvec_sprod_extension(tmp, comkey_sv, sx_sv);
      polxvec_scale_add(ev, tmp, cnst->scalar[i]);
    }
  }

  for(i=0;i<cnst->nphi;i++){
    for(j=0;j<cnst->phi[i]->len;j++){
      polxvec_init_subvec2(sx_sv, sx, cnst->phiw_off[i] + j*rank, 1, rank);
      polxvec_init_subvec2(phi_sv, cnst->phi[i], j, 1, 1);
      polxvec_mul_add(ev, phi_sv, sx_sv);
    }
  }

  polxvec_sub(ev, ev, cnst->b);

  polxvec_free(tmp);
}

/*
  Returns 1 if cnst(sx) = 0; otherwise returns 0.
*/
int comcnst_check(const comcnst cnst, const polxvec sx){
  int check;
  polxvec ev;

  polxvec_init(ev, cnst->rank, 1);
  comcnst_eval(ev, cnst, sx);
  check = polxvec_iszero(ev);
  polxvec_free(ev);

  return check;
}

/*
  out += c*in
*/
static void comcnst_polxvec_mul_add(
  sparsecnst out,
  polxvec c,
  const comcnst in
)
{
  size_t i, j;
  polxvec inphi, outphi, tmp;

  polxvec_init(tmp, c->len, 1);

  polxvec_sprod_add(out->b, c, in->b);
  polxvec_refresh(out->b);

  for(i=0;i<in->ncom;i++){
    polxvec_init_subvec2(inphi, comkey, in->comk_off[i], 1,
                        comkey->len - in->comk_off[i]);
    polxvec_init_subvec2(outphi, out->lin->phi[0], in->comw_off[i], 1,
                        in->comw_len[i]);
    if(in->scalar[i] == 1){
      polxvec_rotation_aggregate_add(outphi, c, inphi);
    }
    else{
      polxvec_scale(tmp, c, in->scalar[i]);
      polxvec_refresh(tmp); // FIX: needed but assertion not triggered
      polxvec_rotation_aggregate_add(outphi, tmp, inphi);
    }
    polxvec_refresh(outphi);
  }

  for(i=0;i<in->nphi;i++){
    for(j=0;j<in->phi[i]->len;j++){
      polxvec_init_subvec2(inphi, in->phi[i], j, 1, 1);
      polxvec_init_subvec2(outphi, out->lin->phi[0], in->phiw_off[i]+j*in->rank,
                          1, in->rank);
      polxvec_mul_add(outphi, inphi, c);
      polxvec_refresh(outphi);
    }
  }

  polxvec_free(tmp);
}

/* 
  Add aggregation of all the input commitment constraints into a single sparse
  constraint.
    - out: output sparse constraint, over the base ring, with rank 1, and one
      part such that off[0] is 0 and phi[0]->len is the length of the witness.
    - in: array of input commitment constraints.
    - ncnst: number of input constraints.
    - chal: vecor of polynomial challenges used to aggregate. There are as many 
      challenges as the sum of the ranks of the input constraints.
*/
void comcnst_aggregate_add(
  sparsecnst out, 
  const comcnst *in, 
  size_t ncnst, 
  const polxvec chal
)
{
  size_t i, off;
  polxvec chal_sv;

  off = 0;
  for(i=0;i<ncnst;i++){
    polxvec_init_subvec2(chal_sv, chal, off, 1, in[i]->rank);
    comcnst_polxvec_mul_add(out, chal_sv, in[i]);
    off += in[i]->rank;
  }
}

void sigmam1cnst_init(
  sigmam1cnst cnst, 
  size_t off1, 
  size_t off2, 
  size_t len,
  int mul
)
{
  cnst->borrowed = 0;
  cnst->off1 = off1;
  cnst->off2 = off2;
  cnst->len = len;
  cnst->mul = mul;
  if(mul){
    polxvec_init(cnst->c, 1, 1);
  }
}

void sigmam1cnst_borrow(sigmam1cnst out, const sigmam1cnst in){
  *out = *in;
  out->borrowed = 1;
}

void sigmam1cnst_copy(sigmam1cnst out, const sigmam1cnst in){
  sigmam1cnst_init(out, in->off1, in->off2, in->len, in->mul);
  if(in->mul){
    polxvec_copy(out->c, in->c);
  }
}

void sigmam1cnst_free(sigmam1cnst cnst){
  if(cnst->borrowed){
    memset(cnst, 0, sizeof(*cnst));
    return;
  }
  if(cnst->mul){
    polxvec_free(cnst->c);
    cnst->mul = 0;
  }
}

/*
  Returns 1 if cnst(sx)=0; otherwise returns 0.
*/
int sigmam1cnst_check(
  const sigmam1cnst cnst, 
  const polxvec sx
)
{
  int check;
  polxvec sx1, sx2, eq;

  polxvec_init_subvec2(sx1, sx, cnst->off1, 1, cnst->len);
  polxvec_init_subvec2(sx2, sx, cnst->off2, 1, cnst->len);
  polxvec_init(eq, cnst->len, 1);
  polxvec_sigmam1(eq, sx1);
  if(cnst->mul){
    polxvec_mul(eq, cnst->c, eq);
  }
  polxvec_sub(eq, eq, sx2);
  check = polxvec_iszero(eq);
  polxvec_free(eq);
  return check;
}

/* 
  Add aggregation of all the input sigmam1 constraints into a single sparse
  constraint.
    - out: output sparse constraint, over the base ring, with rank 1, and one
      part such that off[0] is 0 and phi[0]->len is the length of the witness.
    - in: array of input sigmam1 constraints.
    - ncnst: number of input constraints.
    - chalx: uniform polynomial challenges used to aggregate. The number of
      challenges is (sum_i in[i]->len)
*/
void sigmam1cnst_aggregate_add(
  sparsecnst out, 
  const sigmam1cnst *in, 
  size_t ncnst, 
  const polxvec chalx
)
{
  size_t i, off, maxlen;
  polxvec phi1, phi2, chalx_phi1, chalx_sv;

  if(ncnst == 0){
    return;
  }

  maxlen = 0;
  for(i=0;i<ncnst;i++){
    maxlen = MAX(maxlen, in[i]->len);
  }

  polxvec_init(chalx_phi1, maxlen, 1);

  off = 0;
  for(i=0;i<ncnst;i++){
    polxvec_init_subvec(chalx_phi1, chalx_phi1, 0, 1, in[i]->len);
    polxvec_init_subvec2(chalx_sv, chalx, off, 1, in[i]->len);
    polxvec_init_subvec2(phi1, out->lin->phi[0], in[i]->off1, 1, in[i]->len);
    polxvec_init_subvec2(phi2, out->lin->phi[0], in[i]->off2, 1, in[i]->len);

    polxvec_add(phi2, phi2, chalx_sv);

    if(in[i]->mul){
      polxvec_mul(chalx_phi1, in[i]->c, chalx_sv);
      polxvec_sigmam1(chalx_phi1, chalx_phi1);
    }
    else{
      polxvec_sigmam1(chalx_phi1, chalx_sv);
    }
    polxvec_sub(phi1, phi1, chalx_phi1);

    off += in[i]->len;
  }

  polxvec_free(chalx_phi1);
}

void intcnst_copy(intcnst out, const intcnst in){
  out->off = in->off;
  out->rank = in->rank;
}

int intcnst_check(const intcnst cnst, const polxvec sx){
  size_t i, j;
  polz sxz;
  polxvec pol;
  int check;

  check = 1;

  polzvec_frompolxvec(&sxz, sx, cnst->off, 1, 1);
  for(i=0;i<L;i++){
    for(j=1;j<N;j++){
      if(sxz->limbs[i]->c[j] != 0){
        check = 0;
      }
    }
  }

  if(cnst->rank > 1){
    polxvec_init_subvec2(pol, sx, cnst->off+1, 1, cnst->rank-1);
    check &= polxvec_iszero(pol);
  }

  return check;
}

void intcnst_aggregate_add(
  sparsecnst out, 
  const intcnst *in, 
  size_t ncnst,
  const int64_t *chalz
)
{
  size_t i, j, k, chalz_off;
  int64_t phi64[N];
  polxvec phi_in, phi_out;
  double width = 0x1p64; // 2**64

  polxvec_init(phi_in, 1, 1);

  chalz_off = 0;
  for(i=0;i<ncnst;i++){
    for(j=0;j<in[i]->rank;j++){
      phi64[0] = (j == 0) ? 0 : chalz[chalz_off++];
      for(k=1;k<N;k++){
        phi64[k] = chalz[chalz_off++];
      }
      polxvec_fromint64vec2(phi_in, phi64, 1, 1, width);
      polxvec_init_subvec2(phi_out, out->lin->phi[0], in[i]->off+j, 1, 1);
      polxvec_add(phi_out, phi_out, phi_in);
    }
  }

  polxvec_free(phi_in);
}

void rqcnstset_init(rqcnstset rqc, size_t nsparse, size_t ncommit){
  rqc->nsparse = nsparse;
  rqc->ncom = ncommit;
  if(nsparse + ncommit == 0){
    rqc->sparse = NULL;
  }
  else{
    rqc->sparse = _malloc(nsparse*sizeof(sparsecnst) + ncommit*sizeof(comcnst));
    rqc->com = (comcnst *) &rqc->sparse[nsparse];
  }
  rqc->sparse_nchal = 0;
  rqc->com_nchal = 0;
}

void rqcnstset_free(rqcnstset rqc){
  size_t i;

  for(i=0;i<rqc->nsparse;i++){
    sparsecnst_free(rqc->sparse[i]);
  }
  for(i=0;i<rqc->ncom;i++){
    comcnst_free(rqc->com[i]);
  }
  free(rqc->sparse);

  rqc->nsparse = 0;
  rqc->sparse_nchal = 0;
  rqc->ncom = 0;
  rqc->com_nchal = 0;
}

/*
  Returns 1 if all the constraints in rqc are satisfied; otherwise returns 0.
  The parameters sxq and sxl are two views of the same witness vector sx:
    - sxq: array of the r polxvec that form the witness sx.
    - sxl: concatenation of all the polxvec that form the witness sx.
*/
int rqcnstset_check(
  const rqcnstset rqc, 
  const polxvec sxq[], 
  const polxvec sxl
)
{
  int check;
  size_t i;

  check = 1;
  for(i=0;i<rqc->nsparse;i++){
    const int valid = sparsecnst_check(rqc->sparse[i], sxq, sxl, 1);
    if(!valid && revac_labrador_profile_detail_enabled())
      fprintf(stderr,
              "[Labrador Rq 诊断] sparse[%zu]/%zu 失败：rank=%zu, quad=%zu, lin=%zu\n",
              i, rqc->nsparse, rqc->sparse[i]->b->len,
              rqc->sparse[i]->quad->len,
              rqc->sparse[i]->lin->nparts);
    check &= valid;
  }
  for(i=0;i<rqc->ncom;i++){
    const int valid = comcnst_check(rqc->com[i], sxl);
    if(!valid && revac_labrador_profile_detail_enabled())
      fprintf(stderr, "[Labrador Rq 诊断] commitment[%zu]/%zu 失败\n",
              i, rqc->ncom);
    check &= valid;
  }
  return check;
}

/* 
  Add aggregation of all the constraints in rqc into a single sparse constraint.
    - out: output sparse constraint, over the base ring, with rank 1, and one
      part such that off[0] is 0 and phi[0]->len is the length of the witness.
    - rqc: input set of Rq constraints.
    - chal: vector  of polynomial challenges used to aggregate. The number of
      challenges is rqc->sparse_nchal + rqc->com_nchal.
*/
void rqcnstset_aggregate_add(
  sparsecnst out, 
  const rqcnstset rqc, 
  const polxvec chal
)
{
  polxvec chal_sv;

  polxvec_init_subvec2(chal_sv, chal, 0, 1, rqc->sparse_nchal);
  sparsecnst_aggregate_add(out, rqc->sparse, rqc->nsparse, chal_sv, NULL, 1);

  polxvec_init_subvec2(chal_sv, chal, rqc->sparse_nchal, 1, rqc->com_nchal);
  comcnst_aggregate_add(out, rqc->com, rqc->ncom, chal_sv);
}

void zqcnstset_init(
  zqcnstset zqc, 
  size_t nsparse, 
  size_t maxnsparse, 
  size_t nsigmam1, 
  size_t maxnsigmam1,
  size_t nint
)
{
  zqc->nsparse = nsparse;
  zqc->nsigmam1 = nsigmam1;
  zqc->nint = nint;
  
  if(maxnsparse + maxnsigmam1 + nint == 0){
    zqc->sparse = NULL;
  }
  else{
    zqc->sparse = _malloc(maxnsparse * sizeof(sparsecnst)
                        + maxnsigmam1 * sizeof(sigmam1cnst)
                        + nint * sizeof(intcnst));
    zqc->sigmam1 = (sigmam1cnst *) &zqc->sparse[maxnsparse];
    zqc->intc = (intcnst *) &zqc->sigmam1[maxnsigmam1];
  }
  zqc->sparse_nchal = 0;
  zqc->sigmam1_nchal = 0;
  zqc->int_nchal = 0;
#ifndef NDEBUG
  zqc->maxnsparse = maxnsparse;
  zqc->maxnsigmam1 = maxnsigmam1;
#endif
}

void zqcnstset_free(zqcnstset zqc){
  size_t i;

  for(i=0;i<zqc->nsparse;i++){
    sparsecnst_free(zqc->sparse[i]);
  }
  for(i=0;i<zqc->nsigmam1;i++){
    sigmam1cnst_free(zqc->sigmam1[i]);
  }
  free(zqc->sparse);

  zqc->nsparse = 0;
  zqc->sparse_nchal = 0;
  zqc->nsigmam1 = 0;
  zqc->sigmam1_nchal = 0;
#ifndef NDEBUG
  zqc->maxnsparse = 0;
  zqc->maxnsigmam1 = 0;
#endif
}

/*
  Returns 1 if all the constraints in zqc are satisfied; otherwise returns 0.
  The parameters sxq and sxl are two views of the same witness vector sx:
    - sxq: array of the r polxvec that form the witness sx.
    - sxl: concatenation of all the polxvec that form the witness sx.
*/
int zqcnstset_check(
  const zqcnstset zqc, 
  const polxvec sxq[], 
  const polxvec sxl
)
{
  int check;
  size_t i;

  check = 1;
  for(i=0;i<zqc->nsparse;i++){
    check &= sparsecnst_check(zqc->sparse[i], sxq, sxl, 0);
  }
  for(i=0;i<zqc->nsigmam1;i++){
    check &= sigmam1cnst_check(zqc->sigmam1[i], sxl);
  }
  for(i=0;i<zqc->nint;i++){
    check &= intcnst_check(zqc->intc[i], sxl);
  }
  return check;
}

/* 
  Add aggregation of all the constraints in zqc into a single sparse constraint.
    - out: output sparse constraint, over the base ring, with rank 1, and one
      part such that off[0] is 0 and phi[0]->len is the length of the witness.
    - zqc: input set of Zq constraints.
    - chalz: integer scalar challenges used to aggregate zqc->sparse. The number 
      of challenges is zqc->sparse_nchal.
    - chalx: polynomial challenges used to aggregate zqc->sigmam1. The length of 
      chalx is zqc->sigmam1_nchal.
*/
void zqcnstset_aggregate_add(
  sparsecnst out, 
  const zqcnstset zqc, 
  const int64_t *chalz,
  const polxvec chalx
)
{
  sparsecnst_aggregate_add(out, zqc->sparse, zqc->nsparse, NULL, chalz, 0);

  sigmam1cnst_aggregate_add(out, zqc->sigmam1, zqc->nsigmam1, chalx);

  intcnst_aggregate_add(out, zqc->intc, zqc->nint, &chalz[zqc->sparse_nchal]);
}
