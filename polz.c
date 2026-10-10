#include <stdint.h>
#include <immintrin.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "aesctr.h"
#include "data.h"
#include "fips202.h"
#include "polx.h"
#include "poly.h"
#include "polz.h"
#include "malloc.h"
#include "gaussian.h"
#include "revac_profile.h"

#define REVAC_GAUSSIAN_PAR_CHUNK_COEFFS 65536U

static int revac_gaussian_parallel_enabled(void) {
  const char *value = getenv("REVAC_LABRADOR_GAUSSIAN_PARALLEL");
  /* Each chunk owns its AES state, randbuf and output interval.  Earlier
   * failures attributed to this path were reproduced only when unrelated
   * experimental LNP sections shared mutable polynomial state.  With those
   * sections disabled, isolated stress tests (50 accumulator proofs plus 20
   * complete shows) verified every proof.  Keep an explicit serial rollback
   * for small hosts and diagnostics without changing proof parameters. */
  if (value != NULL && value[0] != '\0')
    return value[0] != '0';
  return 1;
}

static int revac_gaussian_threads(size_t work) {
#ifdef _OPENMP
  const char *value = getenv("REVAC_LABRADOR_GAUSSIAN_THREADS");
  char *end = NULL;
  long requested;
  if(value == NULL || value[0] == '\0')
    return revac_labrador_team_threads(work);
  requested = strtol(value, &end, 10);
  if(end == value || *end != '\0' || requested < 1)
    return 1;
  if(requested > omp_get_max_threads()) requested = omp_get_max_threads();
  if((size_t)requested > work) requested = (long)work;
  return (int)requested;
#else
  (void)work;
  return 1;
#endif
}

static void revac_store_u64_le(uint8_t out[8], uint64_t value) {
  for (size_t i = 0; i < 8; i++)
    out[i] = (uint8_t)(value >> (8 * i));
}

/* Sample fixed-size, independently domain-separated coefficient chunks.  The
 * chunk boundaries do not depend on the OpenMP thread count, so a proof run is
 * reproducible across machines with different core counts.  Each chunk is an
 * independent invocation of the unchanged discrete-Gaussian sampler; only
 * the PRG stream partitioning changes, not sigma or the response format. */
static void revac_gaussian_i32_parallel(
    int32_t *coeffs, size_t nelems, unsigned int log2sd,
    const uint8_t seed[16], uint64_t nonce, int threads) {
  const size_t chunks =
      (nelems + REVAC_GAUSSIAN_PAR_CHUNK_COEFFS - 1) /
      REVAC_GAUSSIAN_PAR_CHUNK_COEFFS;

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) if(chunks > 1 && threads > 1) schedule(static)
#endif
  for (size_t chunk = 0; chunk < chunks; chunk++) {
    uint8_t domain[40];
    uint8_t chunk_key[16];
    aes128ctr_ctx state;
    const size_t off = chunk * REVAC_GAUSSIAN_PAR_CHUNK_COEFFS;
    const size_t count = MIN(REVAC_GAUSSIAN_PAR_CHUNK_COEFFS,
                             nelems - off);

    memcpy(domain, seed, 16);
    revac_store_u64_le(domain + 16, nonce);
    revac_store_u64_le(domain + 24, chunk);
    memcpy(domain + 32, "LNP-GAUS", 8);
    shake128(chunk_key, sizeof(chunk_key), domain, sizeof(domain));
    aes128ctr_init(&state, chunk_key, 0);
    gaussian_i32(coeffs + off, (unsigned int)count, &state, log2sd);
  }
}

void polz_print(const polz a) {
  size_t i,j;

  for(i=0;i<N;i++) {
    printf("%2zu: ",i);
    printf("%5d",a->limbs[0]->c[i]);
    for(j=1;j<L;j++)
      printf(" + %5d * 2^%zu",a->limbs[j]->c[i],14*j);
    printf("\n");
  }
}

void polzvec_setzero(polz *r, size_t len) {
  polyvec_setzero(r[0]->limbs,1,L*len);
}

void polzvec_copy(polz *r, const polz *a, size_t len) {
  size_t i,j,k;
  __m512i f;

  for(i=0;i<len;i++) {
    for(j=0;j<L;j++) {
      for(k=0;k<N/32;k++) {
        f = _mm512_load_si512(&a[i]->limbs[j]->v[k]);
        _mm512_store_si512(&r[i]->limbs[j]->v[k],f);
      }
    }
  }
}

void polz_getcoeff(zz r, const polz a, int k) {
  size_t i;

  for(i=0;i<L;i++)
    r->limbs[i] = a->limbs[i]->c[k];
}

void polz_setcoeff(polz r, const zz a, int k) {
  size_t i;

  for(i=0;i<L;i++)
    r->limbs[i]->c[k] = a->limbs[i];
}

void polz_setcoeff_fromint64(polz r, int64_t a, int k) {
  size_t i;

  for(i=0;i<L-1;i++) {
    r->limbs[i]->c[k] = a & 0x3FFF;
    a >>= 14;
  }
 r->limbs[L-1]->c[k] = a;
}

void polzvec_fromint64vec(polz *r, size_t len, size_t deg, const int64_t v[len*deg*N]) {
  size_t i,j,k;

  for(i=0;i<len;i++)
    for(j=0;j<deg;j++)
      for(k=0;k<N;k++)
        polz_setcoeff_fromint64(r[i*deg+j],v[i*deg*N+k*deg+j],k);
}

int polz_iszero_constcoeff(const polz a) {
  size_t i;
  int64_t r;

  r = 0;
  for(i=0;i<L;i++)
    r |= (uint16_t)a->limbs[i]->c[0];

  r = -r >> 63;
  r += 1;
  return r;
}

// expects centered input
int polz_iszero(const polz a) {
  size_t i,j;
  int64_t r;

  r = 0;
  for(i=0;i<L;i++)
    for(j=0;j<N;j++)
      r |= (uint16_t)a->limbs[i]->c[j];

  r = -r >> 63;
  r += 1;
  return r;
}

int polzvec_iszero(const polz *a, size_t len) {
  size_t i;
  int64_t r;

  r = 0;
  for(i=0;i<len;i++)
    r += polz_iszero(a[i]);

  r -= len;
  r >>= 63;
  r += 1;
  return r;
}

static inline size_t uniform(zz r, const uint8_t *buf) {
  size_t i,j,s,bits;

  i = 0;
  s = 8;
  for(j=0;j<L;j++) {
    bits = (j < L-1) ? 14 : LOGQ - 14*(L-1);
    r->limbs[j] = (int16_t)buf[i] >> (8-s);
    while(s < bits) {
      r->limbs[j] |= (int16_t)buf[++i] << s;
      s += 8;
    }
    s = s - bits;
    r->limbs[j] &= (1 << bits) - 1;
  }

  if(zz_less_than(r,modulus->q))
    return 1;
  else
    return 0;
}

void polzvec_uniform(polz *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  size_t i,j,k;
  uint8_t *buf;
  const size_t chunk = (4096/AES128CTR_BLOCKBYTES+QBYTES-1)/QBYTES*AES128CTR_BLOCKBYTES/N;
  size_t nblocks = chunk*N*QBYTES/AES128CTR_BLOCKBYTES;
  aes128ctr_ctx aesctx;
  zz t;

  aes128ctr_init(&aesctx,seed,nonce);

  while(len >= chunk) {
    k = 0;
    buf = (uint8_t*)r + chunk*N*(2*L - QBYTES);
    aes128ctr_squeezeblocks(buf,nblocks,&aesctx);
    for(i=0;i<chunk;i++) {
      for(j=0;j<N;j++) {
        k += uniform(t,&buf[(i*N+j)*QBYTES]);
        polz_setcoeff(r[i],t,j);
      }
    }

    if(k == chunk*N) {
      r += chunk;
      len -= chunk;
    }
  }

  while(len) {
    k = 0;
    nblocks = (len*N*QBYTES+AES128CTR_BLOCKBYTES-1)/AES128CTR_BLOCKBYTES;
    __attribute__((aligned(64)))
    uint8_t buf2[nblocks*AES128CTR_BLOCKBYTES];
    aes128ctr_squeezeblocks(buf2,nblocks,&aesctx);
    for(i=0;i<len;i++) {
      for(j=0;j<N;j++) {
        k += uniform(t,&buf2[(i*N+j)*QBYTES]);
        polz_setcoeff(r[i],t,j);
      }
    }

    if(k == len*N)
      break;
  }
}

void polz_bitpack(uint8_t r[POLZBYTES], const polz a) {
  size_t i,j,k,bits;
  __m512i f,g,h;

  k = 0;
  h = _mm512_setzero_si512();
  for(i=0;i<L;i++) {
    bits = (i<L-1) ? 14 : LOGQ-14*(L-1);
    for(j=0;j<N/32;j++) {
      f = _mm512_load_si512(&a->limbs[i]->v[j]);
      g = _mm512_slli_epi16(f,k);
      h = _mm512_add_epi16(g,h);
      if(k >= 16-bits) {
        _mm512_storeu_si512(r,h);
        h = _mm512_srli_epi16(f,16-k);
        r += 64;
      }
      k = (k + bits) & 0xF;
    }
  }
  if(k)
    _mm512_storeu_si512(r,h);
}

void polzvec_bitpack(uint8_t *r, const polz *a, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_bitpack(&r[i*POLZBYTES],a[i]);
}

void polz_bitunpack(polz r, const uint8_t buf[POLZBYTES]) {
  int i,j,k,bits;
  __m512i f,g,h,mask;
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);
  const __m512i maskhi = _mm512_srli_epi16(mask14,14*L-LOGQ);

  k = 0;
  h = _mm512_setzero_si512();
  for(i=0;i<L;i++) {
    if(i < L-1) {
      bits = 14;
      mask = mask14;
    }
    else {
      bits = LOGQ - 14*(L-1);
      mask = maskhi;
    }

    for(j=0;j<N/32;j++) {
      if(k < bits) {
        f = _mm512_load_si512((__m512i*)buf);
        g = _mm512_slli_epi16(f,k);
        h = _mm512_add_epi16(g,h);
        g = _mm512_and_si512(h,mask);
        h = _mm512_srli_epi16(f,bits-k);
        buf += 64;
      }
      else {
        g = _mm512_and_si512(h,mask);
        h = _mm512_srli_epi16(h,bits);
      }
      _mm512_store_si512(&r->limbs[i]->v[j],g);
      k = (k - bits) & 0xF;
    }
  }
}

void polzvec_almostuniform(polz *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  size_t i;
  uint8_t *buf;
  size_t chunk,nblocks;
  aes128ctr_ctx aesctx;

  chunk = AES128CTR_BLOCKBYTES/64;
  nblocks = chunk*POLZBYTES/AES128CTR_BLOCKBYTES;
  aes128ctr_init(&aesctx,seed,nonce);

  while(len >= chunk) {
    buf = (uint8_t*)r + chunk*(N*2*L - POLZBYTES);
    aes128ctr_squeezeblocks(buf,nblocks,&aesctx);
    for(i=0;i<chunk;i++)
      polz_bitunpack(r[i],&buf[i*POLZBYTES]);
    len -= chunk;
    r += chunk;
  }

  if(len) {
    nblocks = (len*POLZBYTES+AES128CTR_BLOCKBYTES-1)/AES128CTR_BLOCKBYTES;
    __attribute__((aligned(64)))
    uint8_t buf2[nblocks*AES128CTR_BLOCKBYTES];
    aes128ctr_squeezeblocks(buf2,nblocks,&aesctx);
    for(i=0;i<len;i++)
      polz_bitunpack(r[i],&buf2[i*POLZBYTES]);
  }
}

double polzvec_norm(const polz *a, size_t len) {
  size_t i,j,k;
  long double r,t;
  long double q,hq;

  q = ldexpl(1,LOGQ) - QOFF;
  hq = (q-1)/2;

  r = 0;
  for(i=0;i<len;i++) {
    for(j=0;j<N;j++) {
      t = 0;
      for(k=0;k<L;k++)
        t += ldexpl(a[i]->limbs[k]->c[j],14*k);
      if(t > hq) t -= q;
      r += t*t;
    }
  }

  return sqrtl(r);
}

void polz_reduce(polz r) {
  size_t i,j;
  __m512i f,g,h;
  const __m512i qoff = _mm512_set1_epi16(QOFF);
  const __m512i qofft4 = _mm512_slli_epi16(qoff,2);
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);
  const __m512i maskhi = _mm512_srli_epi16(mask14,14*L-LOGQ);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&r->limbs[L-1]->v[i]);
    g = _mm512_srai_epi16(f,LOGQ-14*(L-1));
    f = _mm512_and_si512(f,maskhi);
    _mm512_store_si512(&r->limbs[L-1]->v[i],f);

    f = _mm512_load_si512(&r->limbs[0]->v[i]);
    h = _mm512_mullo_epi16(g,qoff);
    g = _mm512_mulhi_epi16(g,qofft4);
    h = _mm512_and_si512(h,mask14);
    f = _mm512_add_epi16(f,h);
    h = _mm512_srai_epi16(f,14);
    f = _mm512_and_si512(f,mask14);
    g = _mm512_add_epi16(g,h);
    _mm512_store_si512(&r->limbs[0]->v[i],f);

    for(j=1;j<L;j++) {
      f = _mm512_load_si512(&r->limbs[j]->v[i]);
      f = _mm512_add_epi16(f,g);
      if(j<L-1) {
        g = _mm512_srai_epi16(f,14);
        f = _mm512_and_si512(f,mask14);
      }
      _mm512_store_si512(&r->limbs[j]->v[i],f);
    }
  }
}

void polzvec_reduce(polz *r, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_reduce(r[i]);
}

void polz_caddq(polz r) {
  size_t i,j;
  __m512i f,c;
  __mmask32 mask;
  __m512i qq[L];
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<L;i++)
    qq[i] = _mm512_set1_epi16(modulus->q->limbs[i]);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&r->limbs[L-1]->v[i]);
    mask = _mm512_movepi16_mask(f);
    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&r->limbs[j]->v[i]);
      f = _mm512_mask_add_epi16(f,mask,f,qq[j]);
      if(j > 0) f = _mm512_add_epi16(f,c);
      if(j < L-1) {
        c = _mm512_srai_epi16(f,14);
        f = _mm512_and_si512(f,mask14);
      }
      _mm512_store_si512(&r->limbs[j]->v[i],f);
    }
  }
}

void polzvec_caddq(polz *r, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_caddq(r[i]);
}

void polz_center(polz r) {
  size_t i,j;
  __m512i f,c;
  __mmask32 mask;
  __m512i qq[L], hq[L];
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<L;i++)
    qq[i] = _mm512_set1_epi16(modulus->q->limbs[i]);
  for(i=0;i<L-1;i++) {
    hq[i] = _mm512_srli_epi16(qq[i],1);
    f = _mm512_slli_epi16(qq[i+1],13);
    hq[i] = _mm512_add_epi16(hq[i],f);
    hq[i] = _mm512_and_si512(hq[i],mask14);
  }
  hq[L-1] = _mm512_srli_epi16(qq[L-1],1);

  for(i=0;i<N/32;i++) {
    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&r->limbs[j]->v[i]);
      f = _mm512_sub_epi16(hq[j],f);
      if(j > 0) f = _mm512_add_epi16(f,c);
      if(j < L-1) c = _mm512_srai_epi16(f,14);
    }
    mask = _mm512_movepi16_mask(f);

    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&r->limbs[j]->v[i]);
      f = _mm512_mask_sub_epi16(f,mask,f,qq[j]);
      if(j > 0) f = _mm512_add_epi16(f,c);
      if(j < L-1) {
        c = _mm512_srai_epi16(f,14);
        f = _mm512_and_si512(f,mask14);
      }
      _mm512_store_si512(&r->limbs[j]->v[i],f);
    }
  }
}

void polzvec_center(polz *r, size_t len) {
  const int threads = revac_labrador_team_threads(len);

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) if(revac_labrador_parallel_work(len, 32)) schedule(static)
  #endif
  for(size_t i=0;i<len;i++)
    polz_center(r[i]);
}

void polz_topoly_montgomery(poly r, const polz a, const pdata prime) {
  size_t i,j;
  __m512i f,g,h;
  const __m512i pinvt4 = _mm512_slli_epi16(_mm512_set1_epi16(prime->pinv),2);
  const __m512i p = _mm512_set1_epi16(prime->p);

  for(i=0;i<N/32;i++) {
    g = _mm512_load_si512(&a->limbs[0]->v[i]);
    for(j=1;j<L;j++) {
      h = _mm512_mullo_epi16(g,pinvt4);
      f = _mm512_load_si512(&a->limbs[j]->v[i]);
      g = _mm512_srai_epi16(g,14);
      h = _mm512_mulhi_epu16(h,p);
      g = _mm512_add_epi16(f,g);
      g = _mm512_sub_epi16(g,h);
    }
    _mm512_store_si512(&r->v[i],g);
  }

  //Scaling factor multiplied during NTT
  //poly_scale(r,r,(1LL << (14*(L-1)+16)) % prime->p,prime);
}

void polzvec_topolyvec_montgomery(poly *r, const polz *a, ssize_t stride, size_t len, const pdata prime) {
  size_t i;

  for(i=0;i<len;i++)
    polz_topoly_montgomery(r[stride*i],a[i],prime);
}

void polz_topoly(poly r, const polz a) {
  size_t i;
  __m512i f;

  for(i=0;i<N/32;i++) {  //FIXME
    f = _mm512_load_si512(&a->limbs[0]->v[i]);
    _mm512_store_si512(&r->v[i],f);
  }
}

void polzvec_topolyvec(poly *r, const polz *a, ssize_t stride, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_topoly(r[stride*i],a[i]);
}

void polz_frompoly(polz r, const poly a) {
  size_t i,j;
  __m512i f;
  const __m512i zero = _mm512_setzero_si512();

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->v[i]);
    _mm512_store_si512(&r->limbs[0]->v[i],f);
  }

  for(i=1;i<L;i++)
    for(j=0;j<N/32;j++)
      _mm512_store_si512(&r->limbs[i]->v[j],zero);
}

void polzvec_frompolyvec(polz *r, const poly *a, ssize_t stride, size_t len) {
  const int threads = revac_labrador_team_threads(len);

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) if(revac_labrador_parallel_work(len, 32)) schedule(static)
  #endif
  for(size_t i=0;i<len;i++)
    polz_frompoly(r[i],a[stride*i]);
}

void polz_topolx(polx r, const polz a) {
  size_t i;

  r->width = ldexp(1,2*LOGQ)/12.0;
  for(i=0;i<K;i++)
    polz_topoly_montgomery(r->proj[i],a,primes[i]);

  polx_ntt(r,r);
}

void polzvec_topolxvec(polxvec r, const polz *a, size_t off, ssize_t stride, size_t len) {
  size_t i,k;
  const size_t nblocks = (len + 31) / 32;
  const int threads = revac_labrador_team_threads(nblocks);

  assert(off < r->len);
  assert(off+stride*(len-1) < r->len);
  polxvec_setwidths1(r,off,stride,len,ldexp(1,2*LOGQ)/12.0);
  off *= r->stride;
  stride *= r->stride;
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) \
      if(revac_labrador_conversion_parallel_enabled(\
             "REVAC_LABRADOR_TOPOLX_PARALLEL") && \
         revac_labrador_parallel_work(nblocks, 2)) \
      schedule(static) private(i,k)
  #endif
  for(size_t block=0;block<nblocks;block++) {
    const size_t j = block * 32;
    k = MIN(32,len-j);
    for(i=0;i<K;i++) {
      polzvec_topolyvec_montgomery(&r->proj[i][off+j*stride],&a[j],stride,k,primes[i]);
      polyvec_ntt(&r->proj[i][off+stride*j],&r->proj[i][off+stride*j],stride,stride,k,primes[i]);
    }
  }
}

static void zz_poly_mul(polz r, const zz a, const poly b) {
  size_t i,j;
  int bits;
  __m512i e[L],f,ff;
  __m512i g,h,k,l;
  __m512i mask,qoff;

  for(i=0;i<L;i++)
    e[i] =  _mm512_set1_epi16(a->limbs[i]);  // 2^bits

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&b->v[i]);  // 2^14
    ff = _mm512_slli_epi16(f,2);  // 2^16

    bits = LOGQ - 14*(L-1);
    mask = _mm512_set1_epi16((1 << bits) - 1);
    h = _mm512_mullo_epi16(f,e[L-1]);  // 2^16
    g = _mm512_mulhi_epu16(f,e[L-1]);  // 2^(bits-2)
    k = _mm512_srli_epi16(h,bits);  // 2^(16-bits)
    g = _mm512_slli_epi16(g,16-bits);  // 2^14
    h = _mm512_and_si512(h,mask);  // 2^bits
    k = _mm512_add_epi16(k,g);  // 2^14
    _mm512_store_si512(&r->limbs[L-1]->v[i],h);  // 2^bits

    // Assumes L > 1
    bits = 14;
    mask = _mm512_set1_epi16(0x3FFF);
    qoff = _mm512_set1_epi16(QOFF);  // 2^13
    h = _mm512_mullo_epi16(f,e[0]);  // 2^16
    l = _mm512_mullo_epi16(qoff,k);  // 2^16
    qoff = _mm512_slli_epi16(qoff,2);  // 2^15
    g = _mm512_mulhi_epu16(ff,e[0]);  // 2^14
    k = _mm512_mulhi_epu16(qoff,k);  // 2^13
    h = _mm512_and_si512(h,mask);  // 2^14
    l = _mm512_and_si512(l,mask);  // 2^14
    h = _mm512_add_epi16(h,l);  // 2^15-1
    l = _mm512_srli_epi16(h,bits);  // 2
    h = _mm512_and_si512(h,mask);  // 2^14
    k = _mm512_add_epi16(k,g);  // 2^14+2^13-1
    k = _mm512_add_epi16(k,l);  // 2^14+2^13
    _mm512_store_si512(&r->limbs[0]->v[i],h);  // 2^14

    for(j=1;j<L-1;j++) {
      h = _mm512_mullo_epi16(f,e[j]);  // 2^16
      g = _mm512_mulhi_epu16(ff,e[j]);  // 2^13
      h = _mm512_and_si512(h,mask);  // 2^14
      h = _mm512_add_epi16(h,k);  // 2^15+2^13-2
      k = _mm512_srli_epi16(h,bits);  // 2+1
      h = _mm512_and_si512(h,mask);  // 2^14
      k = _mm512_add_epi16(k,g);  // 2^13+2
      _mm512_store_si512(&r->limbs[j]->v[i],h);  // 2^14
    }

    g = _mm512_load_si512(&r->limbs[L-1]->v[i]);  // 2^bits
    g = _mm512_add_epi16(g,k);  // 2^bits+2^14+2^13-1
    _mm512_store_si512(&r->limbs[L-1]->v[i],g);  // 2^bits+2^14+2^13-1
  }
}

static void zz_poly_fma(polz r, const zz a, const poly b) {
  size_t i,j;
  int bits;
  __m512i e[L],f,ff;
  __m512i g,h,k,l;
  __m512i mask,qoff;

  for(i=0;i<L;i++)
    e[i] = _mm512_set1_epi16(a->limbs[i]);  // 2^bits

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&b->v[i]);  // 2^14 (p_i)
    ff = _mm512_slli_epi16(f,2);  // 2^16

    bits = LOGQ - 14*(L-1);
    mask = _mm512_set1_epi16((1 << bits) - 1);
    h = _mm512_mullo_epi16(f,e[L-1]);  // 2^16
    g = _mm512_mulhi_epu16(f,e[L-1]);  // 2^(bits-2)
    k = _mm512_srli_epi16(h,bits);  // 2^(16-bits)
    g = _mm512_slli_epi16(g,16-bits);  // 2^14
    h = _mm512_and_si512(h,mask);  // 2^bits
    g = _mm512_add_epi16(g,k);  // 2^14
    k = _mm512_load_si512(&r->limbs[L-1]->v[i]);  // 2^bits+2^15
    h = _mm512_add_epi16(h,k);  // 2^(bits+1)+2^15-1
    k = _mm512_srli_epi16(h,bits);  // 2^(15-bits)+2
    h = _mm512_and_si512(h,mask);  // 2^bits
    k = _mm512_add_epi16(k,g);  // 2^14+2^12+1
    _mm512_store_si512(&r->limbs[L-1]->v[i],h);  // 2^bits

    bits = 14;
    qoff = _mm512_set1_epi16(QOFF);  // 2^13
    mask = _mm512_set1_epi16(0x3FFF);

    h = _mm512_mullo_epi16(f,e[0]);  // 2^16
    l = _mm512_mullo_epi16(qoff,k);  // 2^16
    qoff = _mm512_slli_epi16(qoff,2);  // 2^15
    g = _mm512_mulhi_epu16(ff,e[0]);  // 2^14
    k = _mm512_mulhi_epu16(qoff,k);  // 2^14
    h = _mm512_and_si512(h,mask);  // 2^14
    l = _mm512_and_si512(l,mask);  // 2^14
    h = _mm512_add_epi16(h,l);  // 2^15-1
    l = _mm512_load_si512(&r->limbs[0]->v[i]);  // 2^14
    h = _mm512_add_epi16(h,l);  // 2^15+2^14-2
    l = _mm512_srli_epi16(h,bits);  // 2+1
    h = _mm512_and_si512(h,mask);  // 2^14
    k = _mm512_add_epi16(k,g);  // 2^15-1
    k = _mm512_add_epi16(k,l);  // 2^15+1
    _mm512_store_si512(&r->limbs[0]->v[i],h);  // 2^14

    for(j=1;j<L-1;j++) {
      h = _mm512_mullo_epi16(f,e[j]);  // 2^16
      g = _mm512_mulhi_epu16(ff,e[j]);  // 2^14
      h = _mm512_and_si512(h,mask);  // 2^14
      h = _mm512_add_epi16(h,k);  // 2^15+2^14
      k = _mm512_load_si512(&r->limbs[j]->v[i]);  // 2^14
      h = _mm512_add_epi16(h,k);  // 2^16-1
      k = _mm512_srli_epi16(h,bits);  // 2^2
      h = _mm512_and_si512(h,mask);  // 2^14
      k = _mm512_add_epi16(k,g);  // 2^14+2^2-1
      _mm512_store_si512(&r->limbs[j]->v[i],h);  // 2^14
    }

    g = _mm512_load_si512(&r->limbs[L-1]->v[i]);  // 2^bits
    g = _mm512_add_epi16(g,k);  // 2^bits+2^15
    _mm512_store_si512(&r->limbs[L-1]->v[i],g);  // 2^bits+2^15
  }
}

/* Explicit CRT mod q:
 * a_i = a mod p_i, |a| <= (P-1)/2
 * t_i = (P/p_i)^-1 mod p_i
 * alpha_i = a_it_i mod p_i
 * Explicit CRT: a = (\sum_i alpha_i/p_i - round(\sum_i alpha_i/p_i))P
 * Explicit CRT mod q: \sum_i alpha_i(P/p_i mod q) - round(\sum_i alpha_i/p_i)(P mod q)
 */
void polz_frompolx(polz r, const polx a) {
  size_t i;
  polx b;
  poly k;
  __m512i f;
  const __m512i shift = _mm512_set1_epi16(1024);

  polx_invntt(b,a);
  polyvec_setzero(&k,0,1);
  for(i=0;i<K;i++) {
    poly_scale(b->proj[i],b->proj[i],primes[i]->t,primes[i]);  // alpha_i
    poly_caddp(b->proj[i],primes[i]);  // needed for higher precision
    poly_quot_add(k,b->proj[i],primes[i]);  // TODO: Overflow possible for K > 8
  }

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&k->v[i]);
    f = _mm512_add_epi16(f,shift);
    f = _mm512_srli_epi16(f,11);  // round(\sum_i alpha_i/p_i)
    _mm512_store_si512(&k->v[i],f);
  }

  zz_poly_mul(r,modulus->pmq,k);
  for(i=0;i<K;i++)
    zz_poly_fma(r,modulus->xvec[i],b->proj[i]);

  polz_reduce(r);
}

static void frompolxvec(polz *r, const polxvec a, size_t off, ssize_t stride, size_t len) {
  size_t i,j;
  poly b[K][len];
  poly k[len];
  __m512i f;
  const __m512i shift = _mm512_set1_epi16(1024);

  for(i=0;i<K;i++)
    polyvec_invntt(b[i],&a->proj[i][a->stride*off],1,a->stride*stride,len,primes[i]);
  polyvec_setzero(k,1,len);
  for(i=0;i<K;i++) {
    polyvec_scale(b[i],b[i],1,1,len,primes[i]->t,primes[i]);  // alpha_i
    polyvec_caddp(b[i],1,len,primes[i]);  // needed for higher precision
    polyvec_quot_add(k,b[i],1,1,len,primes[i]);  // TODO: Overflow possible for K > 8
  }

  for(i=0;i<len;i++) {
    for(j=0;j<N/32;j++) {
      f = _mm512_load_si512(&k[i]->v[j]);
      f = _mm512_add_epi16(f,shift);
      f = _mm512_srli_epi16(f,11);  // round(\sum_i alpha_i/p_i)
      _mm512_store_si512(&k[i]->v[j],f);
    }
  }

  for(i=0;i<len;i++) {
    zz_poly_mul(r[i],modulus->pmq,k[i]);
    for(j=0;j<K;j++)
      zz_poly_fma(r[i],modulus->xvec[j],b[j][i]);
  }

  polzvec_reduce(r,len);
}

static int64_t exact_positive_mod(__int128 value, int64_t mod) {
  int64_t reduced = (int64_t)(value % mod);
  return reduced < 0 ? reduced + mod : reduced;
}

static int64_t exact_center_mod_q(__int128 value) {
  const int64_t q = ((int64_t)1 << LOGQ) - QOFF;

  value %= q;
  if(value > q/2)
    value -= q;
  else if(value < -q/2)
    value += q;
  return (int64_t)value;
}

static int exact_crt_initialized;
static int64_t exact_inv_scale[K];
static int64_t exact_crt_inv[K];
static __int128 exact_crt_factor[K];
static __int128 exact_crt_product;

static void exact_crt_init(void) {
  size_t k;

  if(exact_crt_initialized) return;
#ifdef _OPENMP
#pragma omp critical(polz_exact_crt_init)
#endif
  {
    if(!exact_crt_initialized){
      polx one, one_coeffs;

      polx_monomial(one, 0, 1);
      polx_invntt(one_coeffs, one);
      exact_crt_product = 1;
      for(k=0;k<K;k++)
        exact_crt_product *= primes[k]->p;
      for(k=0;k<K;k++){
        const int64_t scale = exact_positive_mod(
            one_coeffs->proj[k]->c[0], primes[k]->p);
        exact_inv_scale[k] = exact_positive_mod(
            pow_simple((int16_t)scale, primes[k]->p-2, primes[k]),
            primes[k]->p);
        exact_crt_factor[k] = exact_crt_product / primes[k]->p;
        exact_crt_inv[k] = exact_positive_mod(
            pow_simple((int16_t)exact_positive_mod(
                           exact_crt_factor[k], primes[k]->p),
                       primes[k]->p-2, primes[k]),
            primes[k]->p);
      }
      exact_crt_initialized = 1;
    }
  }
}

static void exact_frompolxvec(polz *r, const polxvec a, size_t off,
                              ssize_t stride, size_t len) {
  size_t coefficient, i, k;
  poly residues[K][len];

  for(k=0;k<K;k++)
    polyvec_invntt(residues[k], &a->proj[k][a->stride*off], 1,
                   a->stride*stride, len, primes[k]);

  for(i=0;i<len;i++){
    for(coefficient=0;coefficient<N;coefficient++){
      __int128 value = 0;

      for(k=0;k<K;k++){
        const int64_t prime = primes[k]->p;
        const int64_t residue = exact_positive_mod(
            (__int128)residues[k][i]->c[coefficient] * exact_inv_scale[k],
            prime);
        const int64_t alpha = exact_positive_mod(
            (__int128)residue * exact_crt_inv[k], prime);
        value += exact_crt_factor[k] * alpha;
      }
      value -= ((value + exact_crt_product/2) / exact_crt_product) *
               exact_crt_product;
      polz_setcoeff_fromint64(r[i], exact_center_mod_q(value), coefficient);
    }
  }
}

void polzvec_exact_frompolxvec(polz *r, const polxvec a, size_t off,
                               ssize_t stride, size_t len) {
  const size_t nblocks = (len + 31) / 32;
  const int threads = revac_labrador_team_threads(nblocks);

  assert(off < a->len);
  assert(off+stride*(len-1) < a->len);
  exact_crt_init();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) \
    if(revac_labrador_conversion_parallel_enabled(\
           "REVAC_LABRADOR_FROMPOLX_PARALLEL") && \
       revac_labrador_parallel_work(nblocks, 2)) schedule(static)
#endif
  for(size_t block=0;block<nblocks;block++){
    const size_t j = block * 32;
    const size_t k = MIN(32,len-j);
    exact_frompolxvec(&r[j], a, off+j*stride, stride, k);
  }
}

void polzvec_frompolxvec(polz *r, const polxvec a, size_t off, ssize_t stride, size_t len) {
  const size_t nblocks = (len + 31) / 32;
  const int threads = revac_labrador_team_threads(nblocks);

  assert(off < a->len);
  assert(off+stride*(len-1) < a->len);
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) \
      if(revac_labrador_conversion_parallel_enabled(\
             "REVAC_LABRADOR_FROMPOLX_PARALLEL") && \
         revac_labrador_parallel_work(nblocks, 2)) schedule(static)
  #endif
  for(size_t block=0;block<nblocks;block++) {
    const size_t j = block * 32;
    const size_t k = MIN(32,len-j);
    frompolxvec(&r[j],a,off+j*stride,stride,k);
  }
}

void polz_add(polz r, const polz a, const polz b) {
  size_t i,j;
  __m512i f,g,c;
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<N/32;i++) {
    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&a->limbs[j]->v[i]);
      g = _mm512_load_si512(&b->limbs[j]->v[i]);
      f = _mm512_add_epi16(f,g);
      if(j > 0) f = _mm512_add_epi16(f,c);
      if(j < L-1) {
        c = _mm512_srai_epi16(f,14);
        f = _mm512_and_si512(f,mask14);
      }
      _mm512_store_si512(&r->limbs[j]->v[i],f);
    }
  }
}

void polzvec_add(polz *r, const polz *a, const polz *b, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_add(r[i],a[i],b[i]);
}

void polz_sub(polz r, const polz a, const polz b) {
  size_t i,j;
  __m512i f,g,c;
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<N/32;i++) {
    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&a->limbs[j]->v[i]);
      g = _mm512_load_si512(&b->limbs[j]->v[i]);
      f = _mm512_sub_epi16(f,g);
      if(j > 0) f = _mm512_add_epi16(f,c);
      if(j < L-1) {
        c = _mm512_srai_epi16(f,14);
        f = _mm512_and_si512(f,mask14);
      }
      _mm512_store_si512(&r->limbs[j]->v[i],f);
    }
  }
}

void polzvec_sub(polz *r, const polz *a, const polz *b, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_sub(r[i],a[i],b[i]);
}

void polz_slli(polz r, const polz a, int s) {
  size_t i,j;
  __m512i f,g,h;
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<N/32;i++) {
    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&a->limbs[j]->v[i]);
      g = _mm512_slli_epi16(f,s);
      if(j) g = _mm512_add_epi16(g,h);
      if(j<L-1) {
        g = _mm512_and_si512(g,mask14);
        h = _mm512_srai_epi16(f,14-s);
      }
      _mm512_store_si512(&r->limbs[j]->v[i],g);
    }
  }
}

void polzvec_slli(polz *r, const polz *a, size_t len, int s) {
  size_t i;

  for(i=0;i<len;i++)
    polz_slli(r[i],a[i],s);
}

void polz_mul(polz r, const polz a, const polz b) {
  polx f,g;

  polz_topolx(f,a);
  polz_topolx(g,b);
  polx_mul(f,f,g);
  polz_frompolx(r,f);
}

void polz_poly_mul(polz r, const polz a, const poly b) {
  polx f,g;

  polz_topolx(f,a);
  polx_frompoly(g,b,ldexp(1,28)/12.0);
  polx_mul(f,f,g);
  polz_frompolx(r,f);
}

// expects centered input
void polz_split(poly lo, polz hi, const polz a, size_t d) {
  size_t i,j;
  __m512i f,g,h;
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->limbs[0]->v[i]);
    g = _mm512_slli_epi16(f,16-d);
    g = _mm512_srai_epi16(g,16-d);  // mod 2^d sign extended
    _mm512_store_si512(&lo->v[i],g);
    f = _mm512_sub_epi16(f,g);  // zero mod 2^d
    f = _mm512_srai_epi16(f,d);
    for(j=0;j<L-1;j++) {
      g = _mm512_load_si512(&a->limbs[j+1]->v[i]);
      h = _mm512_slli_epi16(g,14-d);
      h = _mm512_and_si512(h,mask14);
      f = _mm512_add_epi16(f,h);
      h = _mm512_srai_epi16(f,14);  // carry
      f = _mm512_and_si512(f,mask14);
      _mm512_store_si512(&hi->limbs[j]->v[i],f);
      f = _mm512_srai_epi16(g,d);
      f = _mm512_add_epi16(f,h);
    }
    _mm512_store_si512(&hi->limbs[L-1]->v[i],f);
  }
}

void polzvec_split(poly *lo, polz *hi, const polz *a, ssize_t stride, size_t len, size_t d) {
  size_t i;

  for(i=0;i<len;i++)
    polz_split(lo[stride*i],hi[i],a[i],d);
}

/* Expects centered input */
void polz_bindec(poly *r, const polz a, ssize_t stride, size_t t) {
  size_t i,j,k,s;
  __m512i f,g;
  const __m512i mask = _mm512_set1_epi16(1);

  for(i=0;i<N/32;i++) {
    k = s = 0;
    for(j=0;j<t;j++) {
      if(!s) {
        f = _mm512_load_si512(&a->limbs[k++]->v[i]);
        s = 14;
      }
      g = _mm512_and_si512(f,mask);
      f = _mm512_srai_epi16(f,1);
      _mm512_store_si512(&r[stride*j]->v[i],g);
      s -= 1;
    }
  }
}

void polzvec_bindec(poly *r, const polz *a, size_t len, ssize_t stride, size_t t) {
  size_t i;

  assert(stride >= (ssize_t)len);
  for(i=0;i<len;i++)
    polz_bindec(&r[i],a[i],stride,t);
}

void polz_bindec_topolxvec(polxvec r, const polz a, ssize_t stride, size_t t) {
  poly b[t];
  polxvec rs;

  polz_bindec(b,a,1,t);
  polxvec_init_subvec2(rs,r,0,stride,t);
  polxvec_frompolyvec(rs,b,1,t,1/2.0);
}

void polzvec_bindec_topolxvec(polxvec r, const polz *a, size_t len, ssize_t stride, size_t t) {
  size_t i,j,k;
  poly b[32*t];
  polxvec rs;

  assert(stride >= (ssize_t)len);
  assert(r->len == t*stride);
  for(j=0;(k=MIN(32,len-j));j+=k) {
    polzvec_bindec(b,&a[j],k,k,t);
    for(i=0;i<t;i++) {
      polxvec_init_subvec2(rs,r,stride*i+j,1,k);
      polxvec_frompolyvec(rs,&b[i*k],1,k,1/2.0);
    }
  }
  if((ssize_t)len < stride)
    for(i=0;i<t;i++)
      polxvec_setzero(r,stride*i+len,1,stride-len);
}

/* Expects centered input */
void polz_decompose(poly *r, const polz a, ssize_t stride, size_t t, size_t d) {
  size_t i,j,k,s;
  __m512i f,g,h;
  const __m512i mask = _mm512_set1_epi16((1<<d) - 1);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->limbs[0]->v[i]);
    k = 1;
    s = 14;
    for(j=0;j<t-1;j++) {
      if(s < d && k < L) {
        g = _mm512_load_si512(&a->limbs[k++]->v[i]);
        h = _mm512_slli_epi16(g,s);
        h = _mm512_and_si512(h,mask);
        f = _mm512_add_epi16(f,h);
        h = _mm512_slli_epi16(f,16-d);
        h = _mm512_srai_epi16(h,16-d);  // mod 2^d sign extended
        _mm512_store_si512(&r[stride*j]->v[i],h);
        f = _mm512_sub_epi16(f,h);  // zero mod 2^d
        f = _mm512_srai_epi16(f,d);
        g = _mm512_srai_epi16(g,d-s);
        f = _mm512_add_epi16(f,g);
        s += 14-d;
      }
      else {
        g = _mm512_slli_epi16(f,16-d);
        g = _mm512_srai_epi16(g,16-d);  // mod 2^d sign extended
        _mm512_store_si512(&r[stride*j]->v[i],g);
        f = _mm512_sub_epi16(f,g);  // zero mod 2^d
        f = _mm512_srai_epi16(f,d);
        s -= d;
      }
    }

    if(k < L) {
      g = _mm512_load_si512(&a->limbs[k++]->v[i]);
      g = _mm512_slli_epi16(g,s);
      f = _mm512_add_epi16(f,g);
      f = _mm512_slli_epi16(f,2);
      f = _mm512_srai_epi16(f,2);
    }
    _mm512_store_si512(&r[stride*(t-1)]->v[i],f);
  }
}

void polzvec_decompose(poly *r, const polz *a, size_t len, ssize_t stride, size_t t, size_t d) {
  size_t i;

  if(0 && t==1) {  // FIXME: why did it work before?
    polzvec_topolyvec(r,a,1,len);
    return;
  }

  if(stride >= (ssize_t)len) {
    for(i=0;i<len;i++)
      polz_decompose(&r[i],a[i],stride,t,d);
  }
  else {
    for(i=0;i<len;i++)
      polz_decompose(&r[t*stride*i],a[i],stride,t,d);
  }
}

void polz_decompose_topolxvec(polxvec r, const polz a, ssize_t stride, size_t t, size_t d) {
  poly b[t];
  polxvec rs;

  assert(r->len >= t*stride);
  polz_decompose(b,a,1,t,d);
  polxvec_init_subvec2(rs,r,0,stride,t);
  polxvec_frompolyvec(rs,b,1,t,ldexp(1,2*d)/12);
}

void polzvec_decompose_topolxvec(polxvec r, const polz *a, size_t len, ssize_t stride, size_t t, size_t d) {
  size_t i,j,k;
  double width;
  poly b[32*t];
  polxvec rs;

  assert(stride >= (ssize_t)len);
  width = ldexp(1,2*d)/12;
  for(j=0;(k=MIN(32,len-j));j+=k) {
    polzvec_decompose(b,&a[j],k,k,t,d);
    for(i=0;i<t;i++) {
      polxvec_init_subvec2(rs,r,stride*i+j,1,k);
      polxvec_frompolyvec(rs,&b[k*i],1,k,width);
    }
  }
  if((ssize_t)len < stride)
    for(i=0;i<t;i++)
      polxvec_setzero(r,stride*i+len,1,stride-len);
}

void polz_reconstruct(polz r, const poly *a, ssize_t stride, size_t t, size_t d) {
  size_t i,j,k,s;
  __m512i f,g,h;
  const __m512i mask = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<N/32;i++) {
    h = _mm512_setzero_si512();
    k = s = 0;
    for(j=0;j<L;j++) {
      while(k < t && (j == L-1 || s < 14-d)) {
        f = _mm512_load_si512(&a[stride*k++]->v[i]);
        f = _mm512_slli_epi16(f,s);
        h = _mm512_add_epi16(h,f);
        s += d;
      }  // k == t || (j < L-1 && 14-d <= s < 14)
      if(k < t) {
        f = _mm512_load_si512(&a[stride*k++]->v[i]);
        g = _mm512_slli_epi16(f,s);
        g = _mm512_and_si512(g,mask);
        h = _mm512_add_epi16(h,g);
        g = _mm512_and_si512(h,mask);
        _mm512_store_si512(&r->limbs[j]->v[i],g);
        h = _mm512_srai_epi16(h,14);
        g = _mm512_srai_epi16(f,14-s);
        h = _mm512_add_epi16(h,g);
        s -= 14-d;
      }  // k == t || 0 <= s < d
      else if(j < L-1) {
        g = _mm512_and_si512(h,mask);
        _mm512_store_si512(&r->limbs[j]->v[i],g);
        h = _mm512_srai_epi16(h,14);
      }
      else
        _mm512_store_si512(&r->limbs[j]->v[i],h);
    }
  }
}

void polzvec_reconstruct(polz *r, const poly *a, size_t len, ssize_t stride, size_t t, size_t d) {
  size_t i;

  if(t==1) {
    polzvec_frompolyvec(r,a,1,len);
    return;
  }

  for(i=0;i<len;i++)
    polz_reconstruct(r[i],&a[i],stride,t,d);
}

void polz_sigmam1(polz r, const polz a) {
  polyvec_sigmam1(r->limbs,a->limbs,1,1,L);
}

void polzvec_sigmam1(polz *r, const polz *a, size_t len) {
  polyvec_sigmam1(r[0]->limbs,a[0]->limbs,1,1,L*len);
}

void polz_printint64(const polz a) {
  size_t i;
  zz coeff;
  int64_t tmp;

  for(i=0;i<N;i++) {
    polz_getcoeff(coeff, a, i);
    tmp = int64_fromzz(coeff);
    printf("%2zu: %ld\n",i,tmp);
  }
}

/* The rejection checks only need centered signed 64-bit coefficients.  The
 * old path copied every coefficient into a temporary zz object before
 * reconstructing the same value.  Keep the limb order and sign convention
 * identical to int64_fromzz(), but read the limbs in place. */
static inline int64_t polz_coeff_i64_fast(const polz a, size_t k) {
  int64_t value = a->limbs[0]->c[k];
  for (size_t limb = 1; limb < L; limb++)
    value += (int64_t)a->limbs[limb]->c[k] << (14 * limb);
  return value;
}

int64_t polzvec_sprodz(const polz *a, const polz *b, size_t len) {
  int64_t sum = 0;
  const int threads = revac_labrador_lnp_norm_threads(len);

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) if(threads > 1 && len >= 16 && !omp_in_parallel()) schedule(static) reduction(+:sum)
  #endif
  for (size_t i = 0; i < len; i++) {
    int64_t ia, ib;
    for(size_t j = 0; j < N; j++) {
      ia = polz_coeff_i64_fast(a[i], j);
      ib = polz_coeff_i64_fast(b[i], j);

      sum += ia * ib;
    }
  }
  return sum;
}

/* Fused dot product and squared norm of an already centered response. */
void polzvec_response_products(const polz *z, const polz *v, size_t len,
                                int64_t *zv, int64_t *vv) {
  int64_t dot = 0, norm = 0;
  const int threads = revac_labrador_lnp_norm_threads(len);
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) if(threads > 1 && len >= 16 && !omp_in_parallel()) schedule(static) reduction(+:dot,norm)
  #endif
  for(size_t i=0;i<len;i++){
    for(size_t j=0;j<N;j++){
      const int64_t iz = polz_coeff_i64_fast(z[i], j);
      const int64_t iv = polz_coeff_i64_fast(v[i], j);
      dot += iz * iv;
      norm += iv * iv;
    }
  }
  *zv = dot;
  *vv = norm;
}

/* Compute <mask+shift,shift> and <shift,shift> in one pass.  Rejection
 * sampling needs both products; doing them separately decodes every shift
 * coefficient twice.  The caller keeps mask and shift in their centered
 * integer representatives, so this is exactly the protocol's z=mask+shift
 * product without constructing z before the candidate has been accepted. */
void polzvec_rejection_products(const polz *mask, const polz *shift,
                                size_t len, int64_t *zv, int64_t *vv) {
  int64_t mask_shift = 0;
  int64_t shift_sq = 0;
  const int threads = revac_labrador_lnp_norm_threads(len);

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) if(threads > 1 && len >= 16 && !omp_in_parallel()) schedule(static) reduction(+:mask_shift,shift_sq)
  #endif
  for (size_t i = 0; i < len; i++) {
    int64_t imask, ishift;
    for (size_t j = 0; j < N; j++) {
      imask = polz_coeff_i64_fast(mask[i], j);
      ishift = polz_coeff_i64_fast(shift[i], j);
      mask_shift += imask * ishift;
      shift_sq += ishift * ishift;
    }
  }

  *zv = mask_shift + shift_sq;
  *vv = shift_sq;
}

void polzvec_gaussian_scratch(polz *r, size_t len, unsigned int log2sd,
                              const uint8_t seed[16], uint64_t nonce,
                              int32_t *coeffs) {
  aes128ctr_ctx state;
  const int threads = revac_gaussian_threads(len);
  const size_t chunks = (N * len + REVAC_GAUSSIAN_PAR_CHUNK_COEFFS - 1) /
                        REVAC_GAUSSIAN_PAR_CHUNK_COEFFS;

  if (revac_gaussian_parallel_enabled() && threads > 1 && chunks > 1) {
    revac_gaussian_i32_parallel(coeffs, N * len, log2sd, seed, nonce,
                                MIN(threads, (int)chunks));
  } else {
    aes128ctr_init(&state, seed, nonce);
    gaussian_i32(coeffs, N * len, &state, log2sd);
  }

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) if(threads > 1 && len >= 16 && !omp_in_parallel()) schedule(static)
  #endif
  for (size_t i = 0; i < len; i++) {
    for (size_t j = 0; j < N; j++) {
      polz_setcoeff_fromint64(r[i], coeffs[i * N + j], j);
    }
  }
}


void polzvec_gaussian(polz *r, size_t len, unsigned int log2sd,
                        const uint8_t seed[16], uint64_t nonce) {
  int32_t *coeffs;

  coeffs = _malloc (N * len * sizeof(int32_t));
  polzvec_gaussian_scratch(r, len, log2sd, seed, nonce, coeffs);
  free (coeffs);
}
