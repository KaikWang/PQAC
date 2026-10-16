#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <dlfcn.h>
#include "data.h"
#include "malloc.h"
#include "polx.h"
#include "poly.h"
#include "polz.h"
#include "revac_profile.h"
#include "scratch_pool.h"

static int64_t cmodq(int64_t a) {
  int64_t t;
  const int64_t mask = ((int64_t)1 << LOGQ) - 1;
  const int64_t q = ((int64_t)1 << LOGQ) - QOFF;

  t = a >> LOGQ;
  a &= mask;
  a += t*QOFF;
  t = q/2 - a;
  a -= (t >> 63)&q;
  return a;
}

static int revac_polxvec_parallel_enabled(void) {
  const char *value = getenv("REVAC_LABRADOR_POLXVEC_PARALLEL");
  return value != NULL && value[0] != '\0' && value[0] != '0';
}

/* Optional CUDA sprod backend.  The shared object is deliberately loaded at
 * runtime: CPU-only deployments and future hosts without CUDA keep the exact
 * AVX512 path, while a matching GPU can opt in without changing the proof
 * relation or the wire format. */
typedef int (*revac_cuda_sprod_fn)(const int16_t *const *a,
                                   const int16_t *const *b,
                                   int16_t *const *r, size_t len, size_t n,
                                   size_t k, const int16_t *p,
                                   const int16_t *inv_r);

static size_t revac_sprod_calls, revac_sprod_max_len;
static size_t revac_sprod_eligible, revac_sprod_gpu_calls;

__attribute__((destructor)) static void revac_sprod_report(void) {
  if(getenv("REVAC_LABRADOR_CUDA_AUDIT") != NULL)
    fprintf(stderr, "[CUDA audit] calls=%zu max_len=%zu eligible=%zu gpu_success=%zu\n",
            revac_sprod_calls, revac_sprod_max_len,
            revac_sprod_eligible, revac_sprod_gpu_calls);
}

static revac_cuda_sprod_fn revac_cuda_sprod_loader(void) {
  static void *handle = NULL;
  static revac_cuda_sprod_fn fn = NULL;
  static int state = 0;
  const char *path;

  if(state != 0) return fn;
  state = -1;
  path = getenv("REVAC_LABRADOR_CUDA_LIBRARY");
  if(path == NULL || path[0] == '\0') path = "librevac_labrador_cuda.so";
  handle = dlopen(path, RTLD_NOW | RTLD_LOCAL);
  if(handle == NULL) return NULL;
  fn = (revac_cuda_sprod_fn)dlsym(handle, "revac_cuda_sprod");
  if(fn == NULL) {
    dlclose(handle);
    handle = NULL;
    return NULL;
  }
  state = 1;
  return fn;
}

static int revac_cuda_sprod(polxvec r, const polxvec a, const polxvec b) {
  const char *enabled = getenv("REVAC_LABRADOR_CUDA_SPROD");
  const char *min_value = getenv("REVAC_LABRADOR_CUDA_MIN_LEN");
  size_t min_len = 8192;
  char *end = NULL;
  revac_cuda_sprod_fn fn;
  const int16_t *ap[K], *bp[K];
  int16_t *rp[K], primes_local[K], inv_r[K];
  size_t i;

  if(getenv("REVAC_LABRADOR_CUDA_AUDIT") != NULL) {
    size_t previous = __atomic_load_n(&revac_sprod_max_len, __ATOMIC_RELAXED);
    __atomic_fetch_add(&revac_sprod_calls, 1, __ATOMIC_RELAXED);
    while(previous < a->len && !__atomic_compare_exchange_n(
        &revac_sprod_max_len, &previous, a->len, 1,
        __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {}
  }

  if(enabled == NULL || enabled[0] == '\0' || enabled[0] == '0' ||
     a->len < 2 || a->len != b->len || r->len != 1 ||
     a->stride != 1 || b->stride != 1 || r->stride != 1)
    return 0;
  if(min_value != NULL && min_value[0] != '\0') {
    unsigned long parsed = strtoul(min_value, &end, 10);
    if(end != min_value && *end == '\0') min_len = (size_t)parsed;
  }
  if(a->len < min_len) return 0;
  if(getenv("REVAC_LABRADOR_CUDA_AUDIT") != NULL)
    __atomic_fetch_add(&revac_sprod_eligible, 1, __ATOMIC_RELAXED);
  fn = revac_cuda_sprod_loader();
  if(fn == NULL) return 0;
  for(i = 0; i < K; i++) {
    int64_t x = 1, y = ((int64_t)1 << 16) % primes[i]->p;
    int64_t q = primes[i]->p;
    /* Extended Euclid for R^-1 mod p; this is public parameter data. */
    int64_t old_r = y, new_r = q, old_s = 1, new_s = 0;
    while(new_r != 0) {
      int64_t quotient = old_r / new_r;
      int64_t tmp = old_r - quotient * new_r;
      old_r = new_r; new_r = tmp;
      tmp = old_s - quotient * new_s;
      old_s = new_s; new_s = tmp;
    }
    (void)x;
    if(old_s < 0) old_s += q;
    primes_local[i] = primes[i]->p;
    inv_r[i] = (int16_t)old_s;
    ap[i] = (const int16_t *)a->proj[i];
    bp[i] = (const int16_t *)b->proj[i];
    rp[i] = (int16_t *)r->proj[i];
  }
  int success = fn(ap, bp, rp, a->len, N, K, primes_local, inv_r) == 1;
  if(success && getenv("REVAC_LABRADOR_CUDA_AUDIT") != NULL)
    __atomic_fetch_add(&revac_sprod_gpu_calls, 1, __ATOMIC_RELAXED);
  return success;
}

static int revac_polxvec_scalar_threads(size_t work) {
#ifdef _OPENMP
  const char *value = getenv("REVAC_LABRADOR_POLXVEC_SCALAR_THREADS");
  char *end = NULL;
  long requested;

  if(value == NULL || value[0] == '\0') return 1;
  requested = strtol(value, &end, 10);
  if(end == value || *end != '\0' || requested < 1) requested = 1;
  if(requested > omp_get_max_threads()) requested = omp_get_max_threads();
  if((size_t)requested > work) requested = (long)work;
  return (int)requested;
#else
  (void)work;
  return 1;
#endif
}

static int revac_polxvec_kernel_threads(const char *variable, size_t work) {
#ifdef _OPENMP
  const char *value = getenv(variable);
  char *end = NULL;
  long requested;

  if(value == NULL || value[0] == '\0') return 1;
  requested = strtol(value, &end, 10);
  if(end == value || *end != '\0' || requested < 1) requested = 1;
  if(requested > omp_get_max_threads()) requested = omp_get_max_threads();
  if((size_t)requested > work) requested = (long)work;
  return (int)requested;
#else
  (void)variable;
  (void)work;
  return 1;
#endif
}

static int revac_polxvec_linear_threads(size_t work) {
  return revac_polxvec_kernel_threads(
      "REVAC_LABRADOR_LINEAR_THREADS", work);
}

static int revac_exact_refresh_threads(size_t work) {
#ifdef _OPENMP
  const char *value = getenv("REVAC_LABRADOR_EXACT_REFRESH_THREADS");
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
  (void)work;
  return 1;
#endif
}

void polxvec_init(polxvec r, size_t len, int widths) {
  size_t i;
  void *buf;

  r->alloc = r->len = len;
  r->off = 0;
  r->stride = 1;
  buf = q88_pool_take(len*K*sizeof(poly));
  for(i=0;i<K;i++) {
    r->proj[i] = buf;
    buf = &r->proj[i][len];
  }
  if(widths)
    r->widths = _malloc(len*sizeof(double));
  else
    r->widths = NULL;
}

void polxvec_init_subvec(polxvec r, const polxvec a, size_t off, ssize_t stride, size_t len) {
  size_t i;

  assert(off < a->alloc);

  for(i=0;i<K;i++)
    r->proj[i] = &a->proj[i][off-a->off];
  r->widths = (a->widths) ? &a->widths[off-a->off] : NULL;

  r->alloc = a->alloc;
  r->off = off;
  r->stride = stride;
  if(len)
    r->len = len;
  else if(stride >= 0)
    r->len = (r->alloc-off)/stride;
  else
    r->len = off/-stride;

  assert(off + stride*(r->len-1) < r->alloc);
}

void polxvec_init_subvec2(polxvec r, const polxvec a, size_t off, ssize_t stride, size_t len) {
  size_t i;

  assert(off < a->len);

  for(i=0;i<K;i++)
    r->proj[i] = &a->proj[i][a->stride*off];
  r->widths = (a->widths) ? &a->widths[a->stride*off] : NULL;

  r->alloc = a->alloc;
  r->off = off = a->off + a->stride*off;
  r->stride = stride = stride*a->stride;
  if(len)
    r->len = len;
  else if(stride >= 0)
    r->len = (r->alloc-off)/stride;
  else
    r->len = off/-stride;
  assert(off + stride*(r->len-1) < r->alloc);
}

void polxvec_init_frompolx(polxvec r, polx a) {
  size_t i;

  r->alloc = 0;
  r->off = 0;
  r->stride = 0;
  r->len = 1;
  r->widths = (a->width == MAXWIDTH) ? NULL : &a->width;
  for(i=0;i<K;i++)
    r->proj[i] = &a->proj[i];
}

void polxvec_free(polxvec r) {
  if(!r->alloc) return;
  q88_pool_give(&r->proj[0][-r->off],r->alloc*K*sizeof(poly));
  free(&r->widths[-r->off]);
  r->alloc = 0;
}

void polxvec_setwidths1(polxvec r, size_t off, ssize_t stride, size_t len, double width) {
  size_t i;

  if(!r->widths)
    return;

  assert(width < MAXWIDTH);
  assert(off < r->len);
  if(!len && stride >= 0)
    len = (r->len-off)/stride;
  else if(!len)
    len = off/-stride;
  assert(off + stride*(len-1) < r->len);
  off *= r->stride;
  stride *= r->stride;
  for(i=0;i<len;i++)
    r->widths[off+i*stride] = width;
}

static void polxvec_copywidths(polxvec r, const polxvec a) {
  size_t i;

  assert(r->len >= a->len);
  assert(!r->widths || a->widths);
  if(!r->widths)
    return;
  for(i=0;i<a->len;i++)
    r->widths[i*r->stride] = a->widths[i*a->stride];
}

void polx_print(const polx a) {
  polz t;

  polz_frompolx(t,a);
  polz_center(t);
  polz_print(t);
}

void polxvec_print(const polxvec a, size_t off, ssize_t stride, size_t len) {
  size_t i;
  polz t;

  assert(off < a->len);
  if(!len && stride >= 0)
    len = (a->len-off)/stride;
  else if(!len)
    len = off/-stride;
  assert(off + stride*(len-1) < a->len);
  for(i=0;i<len;i++) {
    polzvec_frompolxvec(&t,a,off+i*stride,0,1);
    polz_center(t);
    polz_print(t);
  }
}

void polx_getcoeff(zz r, const polx a, int k) {
  polz t;

  polz_frompolx(t,a);
  polz_getcoeff(r,t,k);
}

void polxvec_getcoeff(zz r, const polxvec a, size_t j, int k) {
  polz t;

  assert(a->len > j);
  polzvec_frompolxvec(&t,a,j,1,1);
  polz_getcoeff(r,t,k);
}

void polx_setzero(polx r) {
  r->width = 0;
  polyvec_setzero(&r->proj[0],1,K);
}

void polxvec_setzero(polxvec r, size_t off, ssize_t stride, size_t len) {
  assert(off < r->len);
  if(!len && stride >= 0)
    len = (r->len-off)/stride;
  else if(!len)
    len = off/-stride;
  assert(off + stride*(len-1) < r->len);
  polxvec_setwidths1(r,off,stride,len,0);
  off *= r->stride;
  stride *= r->stride;

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(revac_labrador_team_threads(K)) \
      if(revac_labrador_zero_parallel_enabled() && len >= 256 && \
         revac_labrador_parallel_work(K, 2)) schedule(static)
  #endif
  for(size_t i=0;i<K;i++)
    polyvec_setzero(&r->proj[i][off],stride,len);
}

int polx_iszero(const polx a) {
  polz t;

  polz_frompolx(t,a);
  polz_center(t);
  return polz_iszero(t);
}

int polxvec_iszero(const polxvec a) {
  size_t j,k;
  int64_t r = 0;
  polz t[32];

  for(j=0;(k=MIN(32,a->len-j));j+=k) {
    polzvec_frompolxvec(t,a,j,1,k);
    polzvec_center(t,k);
    r += polzvec_iszero(t,k);
  }

  r -= (a->len+31)/32;
  r >>= 63;
  r += 1;
  return r;
}

int polx_iszero_constcoeff(const polx a) {
  polz t;

  polz_frompolx(t,a);
  polz_center(t);
  return polz_iszero_constcoeff(t);
}

int polxvec_iszero_constcoeff(const polxvec a, size_t j) {
  polz t;

  assert(a->len > j);
  polzvec_frompolxvec(&t,a,j,0,1);
  polz_center(t);
  return polz_iszero_constcoeff(t);
}

void polx_copy(polx r, const polx a){
  size_t i;

  r->width = a->width;
  for(i=0;i<K;i++)
    polyvec_copy(&r->proj[i],&a->proj[i],0,0,1);
}

void polxvec_copy(polxvec r, const polxvec a) {
  size_t i;

  assert(r->len >= a->len);
  polxvec_copywidths(r,a);
  for(i=0;i<K;i++)
    polyvec_copy(r->proj[i],a->proj[i],r->stride,a->stride,a->len);
}

void polx_monomial(polx r, int k, int64_t v) {
  size_t i;

  v = cmodq(v);
  r->width = v/sqrt(N);
  for(i=0;i<K;i++)
    poly_monomial_ntt(r->proj[i],modp(v,primes[i]),k,primes[i]);
}

void polxvec_monomial(polxvec r, size_t j, int k, int64_t v) {
  size_t i;

  assert(r->len > j);
  v = cmodq(v);
  r->widths[j] = (double)v*v/N;
  for(i=0;i<K;i++)
    poly_monomial_ntt(r->proj[i][j*r->stride],modp(v,primes[i]),k,primes[i]);
}

void polx_fromint64vec(polx r, const int64_t *a, double width) {
  size_t i;

  assert(width <= MAXWIDTH);
  r->width = width;
  for(i=0;i<K;i++) {
    polyvec_fromint64vec(&r->proj[i],a,0,1,1,primes[i]);
    poly_scale(r->proj[i],r->proj[i],primes[i]->s,primes[i]);
    poly_ntt(r->proj[i],r->proj[i],primes[i]);
  }
}

void polxvec_fromint64vec(polxvec r, const int64_t *a, size_t len, size_t deg, double width) {
  size_t i,k;
  const size_t nblocks = (len + 31) / 32;
  const int threads = revac_labrador_team_threads(nblocks);

  assert(r->len >= len*deg);
  polxvec_setwidths1(r,0,1,len*deg,width);
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) \
      if(revac_labrador_conversion_parallel_enabled(\
             "REVAC_LABRADOR_FROMINT64_PARALLEL") && \
         revac_labrador_parallel_work(nblocks, 2)) \
      schedule(static) private(i,k)
  #endif
  for(size_t block=0;block<nblocks;block++) {
    const size_t j = block * 32;
    k = MIN(32,len-j);
    for(i=0;i<K;i++) {
      polyvec_fromint64vec(&r->proj[i][r->stride*j*deg],&a[j*deg*N],r->stride,k,deg,primes[i]);
      polyvec_scale(&r->proj[i][r->stride*j*deg],&r->proj[i][r->stride*j*deg],
                    r->stride,r->stride,k*deg,primes[i]->s,primes[i]);
      polyvec_ntt(&r->proj[i][r->stride*j*deg],&r->proj[i][r->stride*j*deg],r->stride,r->stride,k*deg,primes[i]);
    }
  }
}

void polxvec_fromint64vec2(polxvec r, const int64_t *a, size_t len, size_t deg, double width) {
  polz t[len*deg];

  polzvec_fromint64vec(t,len,deg,a);
  polzvec_topolxvec(r,t,0,1,len*deg);
  polxvec_setwidths1(r,0,1,len,width);
}

void polx_frompoly(polx r, const poly a, double width) {
  size_t i;

  assert(width <= MAXWIDTH);
  r->width = width;
  for(i=0;i<K;i++)
    poly_scale(r->proj[i],a,primes[i]->s,primes[i]);
  polx_ntt(r,r);
}

void polxvec_frompolyvec(polxvec r, const poly *a, ssize_t stride, size_t len, double width) {
  size_t i,k;
  const size_t nblocks = (len + 31) / 32;
  int threads = revac_polxvec_kernel_threads(
      "REVAC_LABRADOR_FROMPOLY_THREADS", nblocks);
  const int explicit_threads = threads > 1;

  if(!explicit_threads) threads = revac_labrador_team_threads(nblocks);

  assert(r->len >= len);
  polxvec_setwidths1(r,0,1,len,width);
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) \
      if((explicit_threads || revac_labrador_conversion_parallel_enabled(\
             "REVAC_LABRADOR_FROMPOLY_PARALLEL")) && \
         threads > 1 && nblocks >= 2 && !omp_in_parallel()) \
      schedule(static) private(i,k)
  #endif
  for(size_t block=0;block<nblocks;block++) {
    const size_t j = block * 32;
    k = MIN(32,len-j);
    for(i=0;i<K;i++) {
      polyvec_scale(&r->proj[i][j*r->stride],&a[j*stride],r->stride,stride,k,primes[i]->s,primes[i]);
      polyvec_ntt(&r->proj[i][j*r->stride],&r->proj[i][j*r->stride],r->stride,r->stride,k,primes[i]);
    }
  }
}

void polx_frompolxvec(polx r, const polxvec a, size_t j) {
  size_t i;

  assert(j < a->len);
  r->width = a->widths ? a->widths[j] : MAXWIDTH;
  for(i=0;i<K;i++)
    polyvec_copy(&r->proj[i],&a->proj[i][j*a->stride],0,0,1);
}

void polx_almostuniform(polx r, const uint8_t seed[16], uint64_t nonce) {
  size_t i;
  polz t;

  r->width = ldexp(1,2*LOGQ)/12.0;
  polzvec_almostuniform(&t,1,seed,nonce);
  polz_center(t);  // FIXME: Performance?
  for(i=0;i<K;i++) {
    polz_topoly_montgomery(r->proj[i],t,primes[i]);
    poly_ntt(r->proj[i],r->proj[i],primes[i]);
  }
}

void polxvec_almostuniform(polxvec r, const uint8_t seed[16], uint64_t nonce) {
  size_t i,j,k;
  polz t[32];

  polxvec_setwidths1(r,0,1,r->len,ldexp(1,2*LOGQ)/12.0);
  for(j=0;(k=MIN(32,r->len-j));j+=k) {
    polzvec_almostuniform(t,k,seed,nonce);
    polzvec_center(t,k);  // FIXME: Performance?
    nonce += (uint64_t)1 << 32;
    for(i=0;i<K;i++) {
      polzvec_topolyvec_montgomery(&r->proj[i][j*r->stride],t,r->stride,k,primes[i]);
      polyvec_ntt(&r->proj[i][j*r->stride],&r->proj[i][j*r->stride],r->stride,r->stride,k,primes[i]);
    }
  }
}

void polx_ternary(polx r, const uint8_t seed[16], uint64_t nonce) {
  polyvec_ternary(&r->proj[K-1],0,1,seed,nonce);
  polx_frompoly(r,r->proj[K-1],10/16.0);
}

void polxvec_ternary(polxvec r, const uint8_t seed[16], uint64_t nonce) {
  polyvec_ternary(r->proj[K-1],r->stride,r->len,seed,nonce);
  polxvec_frompolyvec(r,r->proj[K-1],r->stride,r->len,10/16.0);
}

void polx_quarternary(polx r, const uint8_t seed[16], uint64_t nonce) {
  polyvec_quarternary(&r->proj[K-1],0,1,seed,nonce);
  polx_frompoly(r,r->proj[K-1],1.5);
}

void polxvec_quarternary(polxvec r, const uint8_t seed[16], uint64_t nonce) {
  polyvec_quarternary(r->proj[K-1],r->stride,r->len,seed,nonce);
  polxvec_frompolyvec(r,r->proj[K-1],r->stride,r->len,1.5);
}

void polx_challenge(polx r, const uint8_t seed[16], uint64_t nonce) {
  polyvec_challenge(&r->proj[K-1],0,1,seed,nonce);
  polx_frompoly(r,r->proj[K-1],(double)(TAU1+4*TAU2)/N);
}

void polxvec_challenge(polxvec r, const uint8_t seed[16], uint64_t nonce) {
  polyvec_challenge(r->proj[K-1],r->stride,r->len,seed,nonce);
  polxvec_frompolyvec(r,r->proj[K-1],r->stride,r->len,(double)(TAU1+4*TAU2)/N);
}

void polx_refresh(polx r) {
  polz t;

  polz_frompolx(t,r);
  polz_center(t);
  polz_topolx(r,t);
}

static int q88_zero_refresh_mode(void) {
  static _Thread_local int mode=-1;
  if(mode<0) {
    const char *env=getenv("REVAC_LABRADOR_ZERO_REFRESH");
    mode=env ? atoi(env) : 0;
  }
  return mode;
}
static int q88_polxvec_range_zero(const polxvec r,size_t off,size_t len) {
  if(r->widths==NULL) return 0;
  for(size_t j=off;j<off+len;j++) if(r->widths[j*r->stride]!=0.0) return 0;
  for(size_t k=0;k<K;k++) for(size_t j=off;j<off+len;j++)
    for(size_t c=0;c<N;c++) if(r->proj[k][j*r->stride]->c[c]!=0) return 0;
  return 1;
}

void polxvec_refresh(polxvec r) {
  const size_t nblocks = (r->len + 31) / 32;
  int threads = revac_polxvec_kernel_threads(
      "REVAC_LABRADOR_REFRESH_THREADS", nblocks);
  const int explicit_threads = threads > 1;

  if(!explicit_threads) threads = revac_labrador_team_threads(nblocks);

  /* One team processes the entire long vector.  The former loop invoked
   * polzvec_center() once per 32-polynomial block, which repeatedly created
   * small OpenMP teams and dominated fixed LaBRADOR folding. */
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) \
      if((explicit_threads || revac_labrador_conversion_parallel_enabled(\
             "REVAC_LABRADOR_REFRESH_PARALLEL")) && \
         threads > 1 && nblocks >= 2 && !omp_in_parallel()) schedule(static)
  #endif
  for(size_t block=0;block<nblocks;block++) {
    const size_t j = block * 32;
    const size_t k = MIN(32,r->len-j);
    if(q88_zero_refresh_mode() && q88_polxvec_range_zero(r,j,k)) continue;
    polz t[32];
    polzvec_frompolxvec(t,r,j,1,k);
    polzvec_center(t,k);
    polzvec_topolxvec(r,t,j,1,k);
  }
}

void polx_reduce(polx r) {
  size_t i;

  for(i=0;i<K;i++)
    poly_reduce(r->proj[i],primes[i]);
}

void polxvec_reduce(polxvec r) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_reduce(r->proj[i],r->stride,r->len,primes[i]);
}

void polx_neg(polx r, const polx a) {
  size_t i;

  r->width = a->width;
  for(i=0;i<K;i++)
    poly_neg(r->proj[i],a->proj[i]);
}

void polxvec_neg(polxvec r, const polxvec a) {
  size_t i;

  assert(r->len == a->len);
  polxvec_copywidths(r,a);
  for(i=0;i<K;i++)
    polyvec_neg(r->proj[i],a->proj[i],r->stride,a->stride,a->len);
}

void polx_add(polx r, const polx a, const polx b) {
  size_t i;

  r->width = a->width + b->width;
  assert(r->width <= MAXWIDTH);
  for(i=0;i<K;i++)
    poly_add(r->proj[i],a->proj[i],b->proj[i]);
  polx_reduce(r);
}

void polxvec_add(polxvec r, const polxvec a, const polxvec b) {
  size_t i;
  const size_t nblocks = (r->len + 31) / 32;
  const int threads = revac_polxvec_linear_threads(nblocks);

  assert(r->len == a->len && a->len == b->len);
  assert(!r->widths || (a->widths && b->widths));
  if(r->widths) {
    for(i=0;i<r->len;i++) {
      r->widths[i*r->stride] = a->widths[i*a->stride] + b->widths[i*b->stride];
      assert(r->widths[i*r->stride] <= MAXWIDTH);
    }
  }
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) \
      if(threads > 1 && nblocks >= 4 && !omp_in_parallel()) schedule(static)
  #endif
  for(size_t block=0;block<nblocks;block++) {
    const size_t j = block * 32;
    const size_t k = MIN(32,r->len-j);
    for(size_t prime_idx=0;prime_idx<K;prime_idx++) {
      polyvec_add(&r->proj[prime_idx][j*r->stride],
                  &a->proj[prime_idx][j*a->stride],
                  &b->proj[prime_idx][j*b->stride],
                  r->stride,a->stride,b->stride,k);
      polyvec_reduce(&r->proj[prime_idx][j*r->stride],r->stride,k,
                     primes[prime_idx]);
    }
  }
}

void polx_sub(polx r, const polx a, const polx b) {
  size_t i;

  r->width = a->width + b->width;
  assert(r->width <= MAXWIDTH);
  for(i=0;i<K;i++)
    poly_sub(r->proj[i],a->proj[i],b->proj[i]);
  polx_reduce(r);
}

void polxvec_sub(polxvec r, const polxvec a, const polxvec b) {
  size_t i;
  const size_t nblocks = (r->len + 31) / 32;
  const int threads = revac_polxvec_linear_threads(nblocks);

  assert(r->len == a->len && a->len == b->len);
  assert(!r->widths || (a->widths && b->widths));
  if(r->widths) {
    for(i=0;i<r->len;i++) {
      r->widths[i*r->stride] = a->widths[i*a->stride] + b->widths[i*b->stride];
      assert(r->widths[i*r->stride] <= MAXWIDTH);
    }
  }
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) \
      if(threads > 1 && nblocks >= 4 && !omp_in_parallel()) schedule(static)
  #endif
  for(size_t block=0;block<nblocks;block++) {
    const size_t j = block * 32;
    const size_t k = MIN(32,r->len-j);
    for(size_t prime_idx=0;prime_idx<K;prime_idx++) {
      polyvec_sub(&r->proj[prime_idx][j*r->stride],
                  &a->proj[prime_idx][j*a->stride],
                  &b->proj[prime_idx][j*b->stride],
                  r->stride,a->stride,b->stride,k);
      polyvec_reduce(&r->proj[prime_idx][j*r->stride],r->stride,k,
                     primes[prime_idx]);
    }
  }
}

void polx_ntt(polx r, const polx a) {
  size_t i;

  for(i=0;i<K;i++)
    poly_ntt(r->proj[i],a->proj[i],primes[i]);
}

void polxvec_ntt(polxvec r, const polxvec a) {
  size_t i;

  assert(r->len == a->len);
  for(i=0;i<K;i++)
    polyvec_ntt(r->proj[i],a->proj[i],r->stride,a->stride,r->len,primes[i]);
}

void polx_invntt(polx r, const polx a) {
  size_t i;

  for(i=0;i<K;i++)
    poly_invntt(r->proj[i],a->proj[i],primes[i]);
}

void polxvec_invntt(polxvec r, const polxvec a) {
  size_t i;

  assert(r->len == a->len);
  for(i=0;i<K;i++)
    polyvec_invntt(r->proj[i],a->proj[i],r->stride,a->stride,r->len,primes[i]);
}

void polxvec_ntt_interleaved_half(polxvec r, const polxvec a) {
  size_t i;

  assert(r->len == 2*a->len);
  if(r->widths) {
    for(i=0;i<r->len;i++) {
      r->widths[i*r->stride] = 0;  // FIXME
      assert(r->widths[i*r->stride] <= MAXWIDTH);
    }
  }
  for(i=0;i<K;i++)
    polyvec_ntt_interleaved_half(r->proj[i],a->proj[i],r->stride,a->stride,r->len,primes[i]);
}

void polxvec_invntt_interleaved(polxvec r, const polxvec a) {
  size_t i;

  assert(r->len == a->len);
  if(r->widths) {
    for(i=0;i<r->len;i++) {
      r->widths[i*r->stride] = 0;  // FIXME
      assert(r->widths[i*r->stride] <= MAXWIDTH);
    }
  }
  for(i=0;i<K;i++)
    polyvec_invntt_interleaved(r->proj[i],a->proj[i],r->stride,a->stride,r->len,primes[i]);
}

void polx_mul(polx r, const polx a, const polx b) {
  size_t i;

  r->width = N*a->width*b->width;
  assert(r->width <= MAXWIDTH);
  for(i=0;i<K;i++)
    poly_pointwise(r->proj[i],a->proj[i],b->proj[i],primes[i]);
}

void polx_poly_mul(polx r, const polx a, const poly b, double width) {
  polx t;

  polx_frompoly(t,b,width);
  polx_mul(r,a,t);
}

void polxvec_mul(polxvec r, const polxvec a, const polxvec b) {
  size_t i;

  if(a->len == 1) {
    const size_t nchunks = (b->len + 31) / 32;
    const int threads = revac_polxvec_scalar_threads(nchunks);
    assert(r->len == b->len);
    assert(!r->widths || (a->widths && b->widths));
    if(r->widths) {
      for(i=0;i<r->len;i++) {
        r->widths[i*r->stride] = N*a->widths[0]*b->widths[i*b->stride];
        assert(r->widths[i*r->stride] <= MAXWIDTH);
      }
    }
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(threads) \
        if(threads > 1 && nchunks >= 2 && !omp_in_parallel()) schedule(static)
    #endif
    for(size_t chunk_idx=0;chunk_idx<nchunks;chunk_idx++) {
      const size_t chunk_off = 32 * chunk_idx;
      const size_t chunk_len = MIN(32, b->len - chunk_off);
      for(size_t prime_idx=0;prime_idx<K;prime_idx++)
        polyvec_poly_pointwise(
            &r->proj[prime_idx][chunk_off*r->stride], *a->proj[prime_idx],
            &b->proj[prime_idx][chunk_off*b->stride], r->stride, b->stride,
            chunk_len, primes[prime_idx]);
    }
  }
  else {
    assert(r->len == a->len && a->len == b->len);
    assert(!r->widths || (a->widths && b->widths));
    if(r->widths) {
      for(i=0;i<r->len;i++) {
        r->widths[i*r->stride] = N*a->widths[i*a->stride]*b->widths[i*b->stride];
        assert(r->widths[i*r->stride] <= MAXWIDTH);
      }
    }
    for(i=0;i<K;i++)
      polyvec_pointwise(r->proj[i],a->proj[i],b->proj[i],r->stride,a->stride,b->stride,r->len,primes[i]);
  }
}

void polxvec_polx_mul(polxvec r, polx a, const polxvec b) {
  polxvec aa;

  polxvec_init_frompolx(aa,a);
  polxvec_mul(r,aa,b);
}

void polx_mul_add(polx r, const polx a, const polx b) {
  size_t i;

  r->width += N*a->width*b->width;
  assert(r->width <= MAXWIDTH);
  for(i=0;i<K;i++)
    poly_pointwise_add(r->proj[i],a->proj[i],b->proj[i],primes[i]);
  polx_reduce(r);
}

void polxvec_mul_add(polxvec r, const polxvec a, const polxvec b) {
  size_t i,j,k;

  if(a->len == 1) {
    const size_t nchunks = (b->len + 31) / 32;
    int threads = revac_polxvec_scalar_threads(nchunks);
    assert(r->len == b->len);
    assert(!r->widths || (a->widths && b->widths));
    if(r->widths) {
      for(i=0;i<b->len;i++) {
        r->widths[i*r->stride] += N*a->widths[0]*b->widths[i*b->stride];
        assert(r->widths[i*r->stride] <= MAXWIDTH);
      }
    }
    if(threads == 1 && revac_polxvec_parallel_enabled())
      threads = revac_labrador_team_threads(nchunks);
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(threads) \
        if(threads > 1 && nchunks >= 2 && !omp_in_parallel()) schedule(static)
    #endif
    for(size_t chunk_idx=0;chunk_idx<nchunks;chunk_idx++) {
      const size_t chunk_off = 32 * chunk_idx;
      const size_t chunk_len = MIN(32, b->len - chunk_off);
      for(size_t prime_idx=0;prime_idx<K;prime_idx++) {
        polyvec_poly_pointwise_add(
            &r->proj[prime_idx][chunk_off*r->stride], *a->proj[prime_idx],
            &b->proj[prime_idx][chunk_off*b->stride], r->stride, b->stride,
            chunk_len, primes[prime_idx]);
        polyvec_reduce(&r->proj[prime_idx][chunk_off*r->stride], r->stride,
                       chunk_len, primes[prime_idx]);
      }
    }
  }
  else {
    assert(r->len == a->len && a->len == b->len);
    assert(!r->widths || (a->widths && b->widths));
    if(r->widths) {
      for(i=0;i<b->len;i++) {
        r->widths[i*r->stride] += N*a->widths[i*a->stride]*b->widths[i*b->stride];
        assert(r->widths[i*r->stride] <= MAXWIDTH);
      }
    }
    for(j=0;(k=MIN(32,a->len-j));j+=k) {
      for(i=0;i<K;i++) {
        polyvec_pointwise_add(&r->proj[i][j*r->stride],&a->proj[i][j*a->stride],&b->proj[i][j*b->stride],
                              r->stride,a->stride,b->stride,k,primes[i]);
        polyvec_reduce(&r->proj[i][j*r->stride],r->stride,k,primes[i]);
      }
    }
  }
}

void polxvec_polx_mul_add(polxvec r, polx a, const polxvec b) {
  polxvec aa;

  polxvec_init_frompolx(aa,a);
  polxvec_mul_add(r,aa,b);
}

void polxvec_sprod(polxvec r, const polxvec a, const polxvec b) {
  size_t i;

  assert(r->len == 1);
  assert(a->len == b->len);
  assert(r->widths || (a->widths && b->widths));
  if(a->len == 0) {
    polxvec_setzero(r,0,0,1);
    return;
  }
  if(r->widths) {
    *r->widths = 0;
    for(i=0;i<a->len;i++)
      *r->widths += N*a->widths[i*a->stride]*b->widths[i*b->stride];
    assert(*r->widths <= MAXWIDTH);
  }
  if(revac_cuda_sprod(r, a, b)) return;
  /* Each prime projection owns a separate output. Preserve the original
     reduction order within each projection; avoid nested OpenMP teams. */
  const int sprod_threads = revac_polxvec_kernel_threads(
      "REVAC_Q88_SPROD_THREADS", K);
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(sprod_threads) \
      if(sprod_threads > 1 && a->len >= 256 && !omp_in_parallel()) schedule(static)
  #endif
  for(i=0;i<K;i++)
    polyvec_sprod_pointwise(*r->proj[i],a->proj[i],b->proj[i],a->stride,b->stride,a->len,primes[i]);
}

void polxvec_sprod_add(polxvec r, const polxvec a, const polxvec b) {
  size_t i;

  assert(r->len == 1);
  assert(a->len == b->len);
  assert(r->widths || (a->widths && b->widths));
  if(r->widths) {
    for(i=0;i<a->len;i++)
      *r->widths += N*a->widths[i*a->stride]*b->widths[i*b->stride];
    assert(*r->widths <= MAXWIDTH);
  }
  /* Each prime projection owns a separate output. Preserve the original
     reduction order within each projection; avoid nested OpenMP teams. */
  const int sprod_threads = revac_polxvec_kernel_threads(
      "REVAC_Q88_SPROD_THREADS", K);
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(sprod_threads) \
      if(sprod_threads > 1 && a->len >= 256 && !omp_in_parallel()) schedule(static)
  #endif
  for(i=0;i<K;i++)
    polyvec_sprod_pointwise_add(*r->proj[i],a->proj[i],b->proj[i],a->stride,b->stride,a->len,primes[i]);
  polxvec_reduce(r);
}

static int64_t polx_exact_center(__int128 value)
{
  const int64_t q = ((int64_t)1 << LOGQ) - QOFF;

  value %= q;
  if(value > q/2)
    value -= q;
  else if(value < -q/2)
    value += q;
  return (int64_t)value;
}

static int64_t polx_exact_positive_mod(__int128 value, int64_t mod)
{
  int64_t reduced = (int64_t)(value % mod);
  return reduced < 0 ? reduced + mod : reduced;
}

/* Recover coefficients exactly from the six internal CRT projections. */
static void polz_exact_frompolx_crt(polz r, const polx a)
{
  static int initialized = 0;
  static int64_t inv_scale[K];
  static int64_t crt_inv[K];
  static int64_t scaled_crt_inv[K];
  static __int128 crt_factor[K];
  static __int128 crt_product;
  polx coeffs;
  size_t i, k;
  const char *check_value = getenv("REVAC_LABRADOR_CRT_SCALE_CHECK");
  const int check_scale = check_value != NULL && check_value[0] == '1';
  const char *reduce_value = getenv("REVAC_LABRADOR_CRT_BOUNDED_REDUCE");
  const int bounded_reduce = reduce_value != NULL && reduce_value[0] != '0';
  const int check_reduce = reduce_value != NULL && reduce_value[0] == '2';

  if(!__atomic_load_n(&initialized, __ATOMIC_ACQUIRE)){
    #ifdef _OPENMP
    #pragma omp critical(polx_exact_crt_init)
    #endif
    {
      if(!initialized){
        polx one, one_coeffs;

        polx_monomial(one, 0, 1);
        polx_invntt(one_coeffs, one);
        crt_product = 1;
        for(k=0;k<K;k++)
          crt_product *= primes[k]->p;
        for(k=0;k<K;k++){
          int64_t scale = polx_exact_positive_mod(
              one_coeffs->proj[k]->c[0], primes[k]->p);
          inv_scale[k] = polx_exact_positive_mod(
              pow_simple((int16_t)scale, primes[k]->p-2, primes[k]),
              primes[k]->p);
          crt_factor[k] = crt_product / primes[k]->p;
          crt_inv[k] = polx_exact_positive_mod(
              pow_simple((int16_t)polx_exact_positive_mod(
                             crt_factor[k], primes[k]->p),
                         primes[k]->p-2, primes[k]),
              primes[k]->p);
          scaled_crt_inv[k] = (inv_scale[k] * crt_inv[k]) % primes[k]->p;
        }
        __atomic_store_n(&initialized, 1, __ATOMIC_RELEASE);
      }
    }
  }

  polx_invntt(coeffs, a);
  for(i=0;i<N;i++){
    __int128 value = 0;

    for(k=0;k<K;k++){
      const int64_t prime = primes[k]->p;
      /* Combine the two fixed CRT scale factors modulo p.  An int16
       * coefficient times a residue below the int16 prime fits in int64;
       * preserve the canonical nonnegative representative exactly. */
      int64_t alpha = ((int64_t)coeffs->proj[k]->c[i] *
                       scaled_crt_inv[k]) % prime;
      if(alpha < 0) alpha += prime;
      if(check_scale){
        const int64_t residue = polx_exact_positive_mod(
            (__int128)coeffs->proj[k]->c[i] * inv_scale[k], prime);
        const int64_t reference = polx_exact_positive_mod(
            (__int128)residue * crt_inv[k], prime);
        if(alpha != reference){
          fprintf(stderr, "CRT scale equivalence failure\n");
          abort();
        }
      }
      value += crt_factor[k] * alpha;
    }
    if(bounded_reduce){
      /* Each nonnegative CRT term is less than P, so 0 <= value < K*P.
       * Binary conditional subtraction replaces a variable 128-bit divide.
       * Centering uses the same floor((value+P/2)/P) as the reference. */
      const __int128 original = value;
      __int128 shifted = value + crt_product/2;
      for(size_t step=next2power(K);step;step>>=1){
        const __int128 multiple = (__int128)step * crt_product;
        if(shifted >= multiple) shifted -= multiple;
      }
      value = shifted - crt_product/2;
      if(check_reduce && value != original -
         ((original + crt_product/2) / crt_product) * crt_product){
        fprintf(stderr, "CRT bounded reduction mismatch\n");
        abort();
      }
    }
    else
      value -= ((value + crt_product/2) / crt_product) * crt_product;
    polz_setcoeff_fromint64(r, polx_exact_center(value), i);
  }
}

static void polz_exact_frompolx(polz r, const polx a)
{
  static int use_exact_crt = -1;

  if(use_exact_crt < 0){
    const char *value = getenv("REVAC_LABRADOR_EXACT_CRT");
    use_exact_crt = value != NULL && value[0] != '\0' && value[0] != '0';
  }
  if(use_exact_crt){
    polz_exact_frompolx_crt(r, a);
    return;
  }
  polz_frompolx(r, a);
  polz_center(r);
}

void polxvec_exact_refresh(polxvec r)
{
  const int threads = revac_exact_refresh_threads(r->len);

  /* Products and sums are exact in every internal CRT projection.  Recover
   * each final coefficient once, reduce it modulo q, then re-encode it.  This
   * preserves the exact centered-q result while avoiding an inverse NTT/CRT
   * reconstruction after every individual multiply-add. */
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) \
      if(threads > 1 && r->len >= 32 && !omp_in_parallel()) schedule(static)
  #endif
  for(size_t i=0;i<r->len;i++){
    if(q88_zero_refresh_mode() && q88_polxvec_range_zero(r,i,1)) continue;
    polx value;
    polz exact;
    polxvec view, result;

    polx_frompolxvec(value, r, i);
    polz_exact_frompolx_crt(exact, value);
    polz_topolx(value, exact);
    polxvec_init_subvec2(view, r, i, 1, 1);
    polxvec_init_frompolx(result, value);
    polxvec_copy(view, result);
  }
}

static void polz_exact_load(int64_t out[N], const polz in)
{
  size_t i;
  zz coeff;

  for(i=0;i<N;i++){
    polz_getcoeff(coeff, in, i);
    out[i] = int64_fromzz(coeff);
  }
}

static void polx_exact_accumulate(__int128 acc[N], const polx a, const polx b)
{
  size_t i, j, anz = 0, bnz = 0;
  int64_t ac[N], bc[N];
  polz az, bz;

  polz_exact_frompolx(az, a);
  polz_exact_frompolx(bz, b);
  polz_exact_load(ac, az);
  polz_exact_load(bc, bz);
  for(i=0;i<N;i++){
    anz += ac[i] != 0;
    bnz += bc[i] != 0;
  }
  if(bnz < anz){
    for(i=0;i<N;i++){
      const int64_t tmp = ac[i];
      ac[i] = bc[i];
      bc[i] = tmp;
    }
  }
  for(i=0;i<N;i++){
    if(ac[i] == 0)
      continue;
    for(j=0;j<N-i;j++)
      acc[i+j] += (__int128)ac[i] * bc[j];
    for(;j<N;j++)
      acc[i+j-N] -= (__int128)ac[i] * bc[j];
  }
}

static void polx_exact_accumulate_crt(__int128 acc[N], const polx a,
                                      const polx b)
{
  size_t i, j;
  int64_t ac[N], bc[N];
  polz az, bz;

  polz_exact_frompolx_crt(az, a);
  polz_exact_frompolx_crt(bz, b);
  polz_exact_load(ac, az);
  polz_exact_load(bc, bz);
  for(i=0;i<N;i++){
    if(ac[i] == 0)
      continue;
    for(j=0;j<N-i;j++)
      acc[i+j] += (__int128)ac[i] * bc[j];
    for(;j<N;j++)
      acc[i+j-N] -= (__int128)ac[i] * bc[j];
  }
}

static void polx_exact_accumulate_sparse(__int128 acc[N],
                                         const uint16_t pos[N],
                                         const int64_t coeffs[N], size_t len,
                                         const polz b)
{
  size_t k, i, j;
  int64_t bc[N];

  polz_exact_load(bc, b);
  for(k=0;k<len;k++){
    i = pos[k];
    for(j=0;j<N-i;j++)
      acc[i+j] += (__int128)coeffs[k] * bc[j];
    for(;j<N;j++)
      acc[i+j-N] -= (__int128)coeffs[k] * bc[j];
  }
}

static size_t polx_exact_sparse_coeffs(uint16_t pos[N], int64_t coeffs[N],
                                       const int64_t dense[N])
{
  size_t i, len = 0;

  for(i=0;i<N;i++){
    if(dense[i] == 0)
      continue;
    pos[len] = (uint16_t)i;
    coeffs[len] = dense[i];
    len++;
  }
  return len;
}

static void polx_exact_store(polx r, __int128 acc[N])
{
  size_t i;
  polz rz;

  for(i=0;i<N;i++)
    polz_setcoeff_fromint64(rz, polx_exact_center(acc[i]), i);
  polz_topolx(r, rz);
}

static void polx_exact_add_current(__int128 acc[N], const polx r)
{
  size_t i;
  polz rz;

  polz_exact_frompolx(rz, r);
  for(i=0;i<N;i++){
    zz coeff;
    polz_getcoeff(coeff, rz, i);
    acc[i] += int64_fromzz(coeff);
  }
}

void polxvec_exact_mul(polxvec r, const polxvec a, const polxvec b)
{
  size_t i;

  assert(a->len == 1 || a->len == b->len);
  assert(r->len == b->len);
  if(a->len == 1){
    int64_t ac[N], acoeffs[N];
    uint16_t apos[N];
    size_t alen;
    polx multiplier;
    polz multiplierz;

    /* Quadratic constraints commonly multiply a vector by one public
     * polynomial.  Recover that polynomial once and retain only its nonzero
     * coefficients.  In the q88 accumulator relation it is a degree-zero
     * scalar, so this replaces one CRT recovery and an N^2 convolution per
     * output polynomial with one shared recovery and an O(N) product. */
    polx_frompolxvec(multiplier, a, 0);
    polz_exact_frompolx(multiplierz, multiplier);
    polz_exact_load(ac, multiplierz);
    alen = polx_exact_sparse_coeffs(apos, acoeffs, ac);
    for(i=0;i<r->len;i++){
      __int128 acc[N] = {0};
      polx rr, rb;
      polz bz;
      polxvec r_sv, result;

      polx_frompolxvec(rb, b, i);
      polz_exact_frompolx(bz, rb);
      polx_exact_accumulate_sparse(acc, apos, acoeffs, alen, bz);
      polx_exact_store(rr, acc);
      polxvec_init_subvec2(r_sv, r, i, 1, 1);
      polxvec_init_frompolx(result, rr);
      polxvec_copy(r_sv, result);
    }
    return;
  }
  for(i=0;i<r->len;i++){
    __int128 acc[N] = {0};
    polx ra, rb;
    polxvec r_sv, result;

    polx_frompolxvec(ra, a, a->len == 1 ? 0 : i);
    polx_frompolxvec(rb, b, i);
    polx_exact_accumulate(acc, ra, rb);
    polx_exact_store(ra, acc);
    polxvec_init_subvec2(r_sv, r, i, 1, 1);
    polxvec_init_frompolx(result, ra);
    polxvec_copy(r_sv, result);
  }
}

/* The rank-1 quadratic-left path must use the exact CRT conversion for both
 * operands; keeping this separate avoids imposing the cost on every vector. */
void polxvec_exact_crt_mul(polxvec r, const polxvec a, const polxvec b)
{
  size_t i;

  assert(a->len == 1 || a->len == b->len);
  assert(r->len == b->len);
  for(i=0;i<r->len;i++){
    __int128 acc[N] = {0};
    polx ra, rb;
    polxvec r_sv, result;

    polx_frompolxvec(ra, a, a->len == 1 ? 0 : i);
    polx_frompolxvec(rb, b, i);
    polx_exact_accumulate_crt(acc, ra, rb);
    polx_exact_store(ra, acc);
    polxvec_init_subvec2(r_sv, r, i, 1, 1);
    polxvec_init_frompolx(result, ra);
    polxvec_copy(r_sv, result);
  }
}

void polxvec_exact_mul_add(polxvec r, const polxvec a, const polxvec b)
{
  size_t i;

  assert(a->len == 1 || a->len == b->len);
  assert(r->len == b->len);
  if(a->len == 1){
    int64_t ac[N], acoeffs[N];
    uint16_t apos[N];
    size_t alen;
    polx challenge;
    polz challengez;

    polx_frompolxvec(challenge, a, 0);
    polz_exact_frompolx(challengez, challenge);
    polz_exact_load(ac, challengez);
    alen = polx_exact_sparse_coeffs(apos, acoeffs, ac);
    for(i=0;i<r->len;i++){
      __int128 acc[N] = {0};
      polx rr, rb;
      polz bz;
      polxvec r_sv, result;

      polx_frompolxvec(rr, r, i);
      polx_frompolxvec(rb, b, i);
      polx_exact_add_current(acc, rr);
      polz_exact_frompolx(bz, rb);
      polx_exact_accumulate_sparse(acc, apos, acoeffs, alen, bz);
      polx_exact_store(rr, acc);
      polxvec_init_subvec2(r_sv, r, i, 1, 1);
      polxvec_init_frompolx(result, rr);
      polxvec_copy(r_sv, result);
    }
    return;
  }
  for(i=0;i<r->len;i++){
    __int128 acc[N] = {0};
    polx rr, ra, rb;
    polxvec r_sv, result;

    polx_frompolxvec(rr, r, i);
    polx_frompolxvec(ra, a, i);
    polx_frompolxvec(rb, b, i);
    polx_exact_add_current(acc, rr);
    polx_exact_accumulate(acc, ra, rb);
    polx_exact_store(rr, acc);
    polxvec_init_subvec2(r_sv, r, i, 1, 1);
    polxvec_init_frompolx(result, rr);
    polxvec_copy(r_sv, result);
  }
}

void polxvec_exact_sprod_add(polxvec r, const polxvec a, const polxvec b)
{
  size_t i;
  __int128 acc[N] = {0};
  polx rr, ra, rb;
  polxvec result;

  assert(r->len == 1);
  assert(a->len == b->len);
  polx_frompolxvec(rr, r, 0);
  polx_exact_add_current(acc, rr);
  for(i=0;i<a->len;i++){
    polx_frompolxvec(ra, a, i);
    polx_frompolxvec(rb, b, i);
    polx_exact_accumulate(acc, ra, rb);
  }
  polx_exact_store(rr, acc);
  polxvec_init_frompolx(result, rr);
  polxvec_copy(r, result);
}

void polx_scale(polx r, const polx a, int64_t s) {
  size_t i;

  s = cmodq(s);
  r->width = (double)s*s*a->width;
  assert(r->width <= MAXWIDTH);
  s <<= 16;  // Montgomery factor
  for(i=0;i<K;i++)
    poly_scale(r->proj[i],a->proj[i],modp(s,primes[i]),primes[i]);
}

void polx_scale_frompoly(polx r, const poly a, double width, int64_t s) {
  size_t i;

  s = cmodq(s);
  r->width = (double)s*s*width;
  assert(r->width <= MAXWIDTH);
  for(i=0;i<K;i++)
    poly_scale(r->proj[i],a,modp(s*primes[i]->s,primes[i]),primes[i]);
  polx_ntt(r,r);
}

void polxvec_scale(polxvec r, const polxvec a, int64_t s) {
  size_t i;
  const size_t nblocks = (r->len + 31) / 32;
  const int threads = revac_polxvec_linear_threads(nblocks);
  int16_t sp[K];

  s = cmodq(s);
  assert(r->len == a->len);
  assert(!r->widths || a->widths);
  if(r->widths) {
    for(i=0;i<r->len;i++) {
      r->widths[i*r->stride] = (double)s*s*a->widths[i*a->stride];
      assert(r->widths[i*r->stride] <= MAXWIDTH);
    }
  }
  s <<= 16;
  for(i=0;i<K;i++)
    sp[i] = modp(s,primes[i]);
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) \
      if(threads > 1 && nblocks >= 4 && !omp_in_parallel()) schedule(static)
  #endif
  for(size_t block=0;block<nblocks;block++) {
    const size_t j = block * 32;
    const size_t k = MIN(32,r->len-j);
    for(size_t prime_idx=0;prime_idx<K;prime_idx++)
      polyvec_scale(&r->proj[prime_idx][j*r->stride],
                    &a->proj[prime_idx][j*a->stride],
                    r->stride,a->stride,k,sp[prime_idx],primes[prime_idx]);
  }
}

void polxvec_scale_frompolyvec(polxvec r, const poly *a, ssize_t stride, size_t len, double width, int64_t s) {
  size_t i,j,k;
  int16_t sp[K];

  assert(r->len == len);
  polxvec_setwidths1(r,0,1,len,(double)s*s*width);
  s = cmodq(s);
  for(i=0;i<K;i++)
    sp[i] = modp(s*primes[i]->s,primes[i]);
  for(j=0;(k=MIN(32,len-j));j+=k) {
    for(i=0;i<K;i++) {
      polyvec_scale(&r->proj[i][j*r->stride],&a[j*stride],r->stride,stride,k,sp[i],primes[i]);
      polyvec_ntt(&r->proj[i][j*r->stride],&r->proj[i][j*r->stride],r->stride,r->stride,k,primes[i]);
    }
  }
}

void polx_scale_add(polx r, const polx a, int64_t s) {
  size_t i;

  s = cmodq(s);
  r->width += (double)s*s*a->width;
  assert(r->width <= MAXWIDTH);
  s <<= 16;
  for(i=0;i<K;i++)
    poly_scale_add(r->proj[i],a->proj[i],modp(s,primes[i]),primes[i]);
  polx_reduce(r);
}

void polxvec_scale_add(polxvec r, const polxvec a, int64_t s) {
  size_t i;
  const size_t nblocks = (r->len + 31) / 32;
  const int threads = revac_polxvec_linear_threads(nblocks);
  int16_t sp[K];

  assert(r->len == a->len);
  assert(!r->widths || a->widths);
  s = cmodq(s);
  if(r->widths) {
    for(i=0;i<r->len;i++) {
      r->widths[i*r->stride] += (double)s*s*a->widths[i*a->stride];
      assert(r->widths[i*r->stride] <= MAXWIDTH);
    }
  }
  s <<= 16;
  for(i=0;i<K;i++)
    sp[i] = modp(s,primes[i]);
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) \
      if(threads > 1 && nblocks >= 4 && !omp_in_parallel()) schedule(static)
  #endif
  for(size_t block=0;block<nblocks;block++) {
    const size_t j = block * 32;
    const size_t k = MIN(32,r->len-j);
    for(size_t prime_idx=0;prime_idx<K;prime_idx++) {
      polyvec_scale_add(&r->proj[prime_idx][j*r->stride],
                        &a->proj[prime_idx][j*a->stride],
                        r->stride,a->stride,k,sp[prime_idx],primes[prime_idx]);
      polyvec_reduce(&r->proj[prime_idx][j*r->stride],r->stride,k,
                     primes[prime_idx]);
    }
  }
}

static size_t polxvec_sprod_extension_internal(polxvec c, const polxvec a, const polxvec b, int add) {
  size_t i,j,k;
  size_t deg;
  size_t off_a = 0;
  double width_a,width_b,width_c;
  poly tmp[c->len];

  deg = next2power(c->len);
  k  = extlen(b->len,deg);
  assert(a->len >= off_a+k);
  assert(!c->widths || (a->widths && b->widths));
  if(c->widths) {
    width_c = 0;
    for(i=0;i<b->len;i+=deg) {
      width_a = width_b = 0;
      for(j=0;j<deg;j++) {
         width_a += a->widths[(off_a+i+j)*a->stride];
         width_b += (i+j < b->len) ? b->widths[(i+j)*b->stride] : 0;
      }
      width_c += width_a*width_b;
    }
    width_c *= N/deg;
    if(!add) {
      assert(width_c <= MAXWIDTH);
      for(i=0;i<c->len;i++)
        c->widths[i*c->stride] = width_c;
    }
    else {
      for(i=0;i<c->len;i++) {
        c->widths[i*c->stride] += width_c;
        assert(c->widths[i*c->stride] <= MAXWIDTH);
      }
    }
  }
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(revac_labrador_team_threads(K)) \
      if(!add && revac_labrador_commit_parallel_enabled() && b->len >= 512 && \
         revac_labrador_parallel_work(K, 2)) schedule(static)
  #endif
  for(i=0;i<K;i++) {
    if(!add)
      polyvec_sprod_extension(c->proj[i],a->proj[i],b->proj[i],
                              c->stride,a->stride,b->stride,deg,c->len,b->len,primes[i]);
    else {
      polyvec_sprod_extension(tmp,&a->proj[i][off_a*a->stride],b->proj[i],
                              1,a->stride,b->stride,deg,c->len,b->len,primes[i]);
      polyvec_add(c->proj[i],c->proj[i],tmp,c->stride,c->stride,1,c->len);  // TODO: better rely on interleave_reduce_add in sprod_extension
      polyvec_reduce(c->proj[i],c->stride,c->len,primes[i]);
    }
  }
  return k;
}

size_t polxvec_sprod_extension(polxvec c, const polxvec a, const polxvec b){
  return polxvec_sprod_extension_internal(c, a, b, 0);
}

size_t polxvec_sprod_extension_add(polxvec c, const polxvec a, const polxvec b){
  return polxvec_sprod_extension_internal(c, a, b, 1);
}

size_t polxvec_collaps_add_extension(polxvec c, const polxvec b, const polxvec a, size_t off_a) {
  size_t i,j,k;
  size_t deg;
  double width_a,width_b,width_c;

  deg = next2power(b->len);
  k  = extlen(c->len,deg);
  assert(a->len >= off_a+k);
  assert(!c->widths || (a->widths && b->widths));
  if(c->widths) {
    width_b = 0;
    for(i=0;i<b->len;i++)
      width_b += b->widths[i*b->stride];
    for(i=0;i<c->len;i+=deg) {
      width_a = 0;
      for(j=0;j<deg;j++)
         width_a += a->widths[(off_a+i+j)*a->stride];
      width_c = width_a*width_b*N/deg;
      for(j=0;j<MIN(deg,c->len-i);j++) {
        c->widths[(i+j)*c->stride] += width_c;
        assert(c->widths[(i+j)*c->stride] <= MAXWIDTH);
      }
    }
  }
  for(i=0;i<K;i++)
    polyvec_collaps_add_extension(c->proj[i],b->proj[i],&a->proj[i][off_a*a->stride],
                                  c->stride,b->stride,a->stride,deg,c->len,b->len,primes[i]);
  return k;
}

size_t polxvec_pairwise_mul(polxvec h, const polxvec a) {
  size_t i;
  const size_t r = a->len;
  const size_t m = (r*r+r)/2;

  assert(h->len == m);
  if(h->widths) {
    for(i=0;i<m;i++) {
      h->widths[i*h->stride] = 0;  //FIXME
      assert(h->widths[i*h->stride] <= MAXWIDTH);
    }
  }
  for(i=0;i<K;i++)
    polyvec_pairwise_pointwise(h->proj[i],a->proj[i],a->proj[i],h->stride,a->stride,a->stride,r,r,primes[i]);
  return m;
}

size_t polxvec_pairwise_sprod(polxvec g, const polxvec a, const polxvec b, size_t r, size_t n) {
  size_t i;
  const size_t m = (r*r+r)/2;

  assert(g->len == m);
  assert(a->len == r*n && b->len == r*n);
  if(g->widths) {
    for(i=0;i<m;i++) {
      g->widths[i*g->stride] = 0;  //FIXME
      assert(g->widths[i*g->stride] <= MAXWIDTH);
    }
  }
  for(i=0;i<K;i++)
    polyvec_pairwise_sprod(g->proj[i],a->proj[i],b->proj[i],r,n,primes[i]);
  return m;
}

void polx_bindec(poly *r, const polx a, ssize_t stride, size_t t) {
  polz b;

  polz_frompolx(b,a);
  polz_center(b);
  polz_bindec(r,b,stride,t);
}

void polxvec_bindec(poly *r, const polxvec a, ssize_t stride, size_t t) {
  size_t j,k;
  polz b[32];

  assert(stride >= (ssize_t)a->len);
  for(j=0;(k=MIN(32,a->len-j));j+=k) {
    polzvec_frompolxvec(b,a,j,1,k);
    polzvec_center(b,k);
    polzvec_bindec(&r[j],b,k,stride,t);
  }
}

void polx_decompose(poly *r, const polx a, ssize_t stride, size_t t, size_t d) {
  polz b;

  polz_frompolx(b,a);
  polz_center(b);
  polz_decompose(r,b,stride,t,d);
}

void polxvec_decompose(poly *r, const polxvec a, ssize_t stride, size_t t, size_t d) {
  size_t j,k;
  polz b[32];

  assert(stride >= (ssize_t)a->len);
  for(j=0;(k=MIN(32,a->len-j));j+=k) {
    polzvec_frompolxvec(b,a,j,1,k);
    polzvec_center(b,k);
    polzvec_decompose(&r[j],b,k,stride,t,d);
  }
}

void polxvec_exact_decompose(poly *r, const polxvec a, ssize_t stride,
                             size_t t, size_t d) {
  polz b;

  assert(a->len == 1);
  assert(stride >= 1);
  {
    polx value;
    polx_frompolxvec(value, a, 0);
    polz_exact_frompolx_crt(b, value);
  }
  polz_decompose(r, b, stride, t, d);
}

void polx_reconstruct(polx r, const poly *a, ssize_t stride, size_t t, size_t d) {
  polz b;

  polz_reconstruct(b,a,stride,t,d);
  polz_topolx(r,b);
}

void polxvec_reconstruct(polxvec r, const poly *a, ssize_t stride, size_t t, size_t d) {
  size_t j,k;
  polz b[32];

  //polxvec_setwidths1(r,0,r->len,ldexp(1,2*d*t)/12.0);
  for(j=0;(k=MIN(32,r->len-j));j+=k) {
    polzvec_reconstruct(b,&a[j],k,stride,t,d);
    polzvec_topolxvec(r,b,j,1,k);
  }
}

void polx_sigmam1(polx r, const polx a) {
  size_t i;

  r->width = a->width;
  for(i=0;i<K;i++)
    poly_sigmam1_ntt(r->proj[i],a->proj[i]);
}

void polxvec_sigmam1(polxvec r, const polxvec a) {
  size_t i;

  assert(r->len == a->len);
  polxvec_copywidths(r,a);
  for(i=0;i<K;i++)
    polyvec_sigmam1_ntt(r->proj[i],a->proj[i],r->stride,a->stride,r->len);
}

/* FIXME: NTT rep! */
void polx_sigma5(polx r, const polx a) {
  size_t i;

  r->width = a->width;
  for(i=0;i<K;i++)
    poly_sigma5(r->proj[i],a->proj[i]);
}

void polxvec_sigma5(polxvec r, const polxvec a) {
  size_t i;

  assert(r->len == a->len);
  polxvec_copywidths(r,a);
  for(i=0;i<K;i++)
    polyvec_sigma5(r->proj[i],a->proj[i],r->stride,a->stride,r->len);
}

void polx_sigma5inv(polx r, const polx a) {
  size_t i;

  r->width = a->width;
  for(i=0;i<K;i++)
    poly_sigma5inv(r->proj[i],a->proj[i]);
}

void polxvec_sigma5inv(polxvec r, const polxvec a) {
  size_t i;

  assert(r->len == a->len);
  polxvec_copywidths(r,a);
  for(i=0;i<K;i++)
    polyvec_sigma5(r->proj[i],a->proj[i],r->stride,a->stride,r->len);
}

void polx_flip(polx r, const polx a) {
  size_t i;

  r->width = 1 + a->width;
  for(i=0;i<K;i++)
    poly_flip_ntt(r->proj[i],a->proj[i],primes[i]);
}

void polxvec_flip(polxvec r, const polxvec a) {
  size_t i;

  assert(r->len == a->len);
  assert(!r->widths || a->widths);
  if(r->widths) {
    for(i=0;i<r->len;i++)
      r->widths[r->stride*i] = 1 + a->widths[a->stride*i];
  }
  for(i=0;i<K;i++)
    polyvec_flip_ntt(r->proj[i],a->proj[i],r->stride,a->stride,r->len,primes[i]);
}
