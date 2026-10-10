#include "labrados_python.h"
#include "malloc.h"
#include "fips202.h"
#include "comkey.h"
#include "timing.h"
#include "revac_profile.h"
#include <stdlib.h>
#include <string.h>

#if defined(__GLIBC__)
extern int malloc_trim(size_t pad);
extern int mallopt(int param, int value);
#define PY_M_TRIM_THRESHOLD (-1)
#define PY_M_MMAP_THRESHOLD (-3)
#define PY_M_ARENA_MAX      (-8)
#endif

static const uint8_t PY_PROOF_MAGIC[] = {
  'R','E','V','A','C','-','L','A','B','R','A','D','O','R','-','P','R','O','O','F',1
};
static unsigned int py_proof_encoding(void) {
  const char *value = getenv("REVAC_LABRADOR_PROOF_ENCODING");
  if (value != NULL && strcmp(value, "v1") == 0) return 1U;
  if (value != NULL && strcmp(value, "v2") == 0) return 2U;
  if (value != NULL && strcmp(value, "v3") == 0) return 3U;
  return 4U;
}

static unsigned int py_proof_field_mask(void) {
  return py_proof_encoding() >= 2U ? 7U : 0U;
}

static int py_wr(uint8_t **out, size_t *rem, const void *src, size_t len);
static int py_rd(const uint8_t **in, size_t *rem, void *dst, size_t len);
static int py_wr_u64(uint8_t **out, size_t *rem, uint64_t x);
static int py_rd_u64(const uint8_t **in, size_t *rem, uint64_t *x);

static int py_statement_threads(size_t independent_work) {
#ifdef _OPENMP
  const char *value = getenv("REVAC_LABRADOR_STATEMENT_THREADS");
  char *end = NULL;
  long requested = 1;

  if (value != NULL && value[0] != '\0') {
    requested = strtol(value, &end, 10);
    if (end == value || *end != '\0' || requested < 1) requested = 1;
  }
  if (requested > omp_get_max_threads()) requested = omp_get_max_threads();
  if ((size_t)requested > independent_work) requested = (long)independent_work;
  if (omp_in_parallel()) requested = 1;
  return requested > 0 ? (int)requested : 1;
#else
  (void)independent_work;
  return 1;
#endif
}

static int py_q96_compact_hash_enabled(void) {
  const char *value = getenv("REVAC_LABRADOR_Q96_COMPACT_HASH");
  return value != NULL && value[0] != '\0' && value[0] != '0';
}

static unsigned int py_polyvec_signed_bits(const poly *vec, size_t len) {
  int32_t min_value = 0;
  int32_t max_value = 0;

  for (size_t i = 0; i < len; i++) {
    for (size_t j = 0; j < N; j++) {
      const int32_t value = vec[i]->c[j];
      if (value < min_value) min_value = value;
      if (value > max_value) max_value = value;
    }
  }
  for (unsigned int bits = 1; bits < 16; bits++) {
    const int32_t lower = -((int32_t)1 << (bits - 1));
    const int32_t upper = ((int32_t)1 << (bits - 1)) - 1;
    if (min_value >= lower && max_value <= upper) return bits;
  }
  return 16;
}

static unsigned int py_int32vec_signed_bits(const int32_t *vec, size_t len) {
  int64_t min_value = 0;
  int64_t max_value = 0;

  for (size_t i = 0; i < len; i++) {
    const int64_t value = vec[i];
    if (value < min_value) min_value = value;
    if (value > max_value) max_value = value;
  }
  for (unsigned int bits = 1; bits < 32; bits++) {
    const int64_t lower = -((int64_t)1 << (bits - 1));
    const int64_t upper = ((int64_t)1 << (bits - 1)) - 1;
    if (min_value >= lower && max_value <= upper) return bits;
  }
  return 32;
}

static size_t py_signed_polyvec_nbytes(size_t len, unsigned int bits) {
  if (bits == 0 || bits > 16 || len > SIZE_MAX / N ||
      len * N > (SIZE_MAX - 7) / bits) return SIZE_MAX;
  return (len * N * bits + 7) / 8;
}

static size_t py_signed_int32vec_nbytes(size_t len, unsigned int bits) {
  if (bits == 0 || bits > 32 || len > (SIZE_MAX - 7) / bits) {
    return SIZE_MAX;
  }
  return (len * bits + 7) / 8;
}

static int py_wr_signed_int32vec(uint8_t **out, size_t *rem,
                                 const int32_t *vec, size_t len) {
  const unsigned int bits = py_int32vec_signed_bits(vec, len);
  const size_t nbytes = py_signed_int32vec_nbytes(len, bits);
  const uint8_t width = (uint8_t)bits;
  size_t bitpos = 0;

  if (nbytes == SIZE_MAX || py_wr(out, rem, &width, 1) || *rem < nbytes) {
    return 1;
  }
  memset(*out, 0, nbytes);
  for (size_t i = 0; i < len; i++) {
    const uint32_t value = (uint32_t)vec[i];
    for (unsigned int bit = 0; bit < bits; bit++, bitpos++) {
      (*out)[bitpos / 8] |=
          (uint8_t)(((value >> bit) & 1U) << (bitpos % 8));
    }
  }
  *out += nbytes;
  *rem -= nbytes;
  return 0;
}

static int py_rd_signed_int32vec(const uint8_t **in, size_t *rem,
                                 int32_t *vec, size_t len) {
  uint8_t width;
  size_t nbytes;
  size_t bitpos = 0;

  if (py_rd(in, rem, &width, 1) || width == 0 || width > 32) return 1;
  nbytes = py_signed_int32vec_nbytes(len, width);
  if (nbytes == SIZE_MAX || *rem < nbytes) return 1;
  for (size_t i = 0; i < len; i++) {
    uint32_t value = 0;
    for (unsigned int bit = 0; bit < width; bit++, bitpos++) {
      value |= (uint32_t)(((*in)[bitpos / 8] >> (bitpos % 8)) & 1U)
               << bit;
    }
    int64_t signed_value = value;
    if ((value & ((uint32_t)1 << (width - 1))) != 0) {
      signed_value -= (int64_t)1 << width;
    }
    vec[i] = (int32_t)signed_value;
  }
  *in += nbytes;
  *rem -= nbytes;
  return py_int32vec_signed_bits(vec, len) == width ? 0 : 1;
}

static int py_wr_signed_polyvec_blob(uint8_t **out, size_t *rem,
                                     const poly *vec, size_t len) {
  const unsigned int bits = py_polyvec_signed_bits(vec, len);
  const size_t nbytes = py_signed_polyvec_nbytes(len, bits);
  const uint8_t width = (uint8_t)bits;
  const uint32_t mask = ((uint32_t)1 << bits) - 1U;
  size_t bitpos = 0;

  if (nbytes == SIZE_MAX || py_wr(out, rem, &width, 1) ||
      py_wr_u64(out, rem, nbytes) || *rem < nbytes) return 1;
  memset(*out, 0, nbytes);
  for (size_t i = 0; i < len; i++) {
    for (size_t j = 0; j < N; j++) {
      const uint32_t value = (uint16_t)vec[i]->c[j] & mask;
      for (unsigned int bit = 0; bit < bits; bit++, bitpos++) {
        (*out)[bitpos / 8] |=
            (uint8_t)(((value >> bit) & 1U) << (bitpos % 8));
      }
    }
  }
  *out += nbytes;
  *rem -= nbytes;
  return 0;
}

static int py_rd_signed_polyvec_blob(const uint8_t **in, size_t *rem,
                                     poly *vec, size_t len) {
  uint8_t width;
  uint64_t got;
  size_t nbytes;
  size_t bitpos = 0;

  if (py_rd(in, rem, &width, 1) || width == 0 || width > 16) return 1;
  nbytes = py_signed_polyvec_nbytes(len, width);
  if (nbytes == SIZE_MAX || py_rd_u64(in, rem, &got) || got != nbytes ||
      *rem < nbytes) return 1;
  for (size_t i = 0; i < len; i++) {
    for (size_t j = 0; j < N; j++) {
      uint32_t value = 0;
      for (unsigned int bit = 0; bit < width; bit++, bitpos++) {
        value |= (uint32_t)(((*in)[bitpos / 8] >> (bitpos % 8)) & 1U)
                 << bit;
      }
      int32_t signed_value = (int32_t)value;
      if ((value & ((uint32_t)1 << (width - 1))) != 0) {
        signed_value -= (int32_t)((uint32_t)1 << width);
      }
      vec[i]->c[j] = (int16_t)signed_value;
    }
  }
  *in += nbytes;
  *rem -= nbytes;
  return py_polyvec_signed_bits(vec, len) == width ? 0 : 1;
}

static size_t py_lab_proof_polz_len(const lab_params pp) {
  if (pp->compressed) {
    return 4 * pp->kappa[2];
  }
  if (!pp->tail) {
    return 2 * pp->kappa[1] + LIFTS;
  }
  return pp->len[LAB_INCOM] + pp->len[LAB_QUADG] +
         pp->len[LAB_LIFT] + pp->len[LAB_LING];
}

static size_t py_pack_lab_proof_polz_len(const pack_params pp, size_t round) {
  /* The ZK/LNP round reserves one entry in pp->p/pi->p, but that entry is a
     deliberately zeroed placeholder rather than a lab_proof.  Treat it as an
     empty blob; interpreting it as an ordinary non-tail Labrador round would
     invent LIFTS coefficients and dereference a NULL message pointer. */
  if (pp->zkp != NULL && round == pp->zkround) {
    return 0;
  }
  return py_lab_proof_polz_len(pp->p[round]);
}

static void py_lab_proof_set_ptrs(lab_proof pi, const lab_params pp, polz *arr) {
  pi->m[0] = arr;
  if (arr == NULL) {
    pi->m[1] = NULL;
    pi->m[2] = NULL;
    pi->m[3] = NULL;
  } else if (pp->compressed) {
    pi->m[1] = &pi->m[0][pp->kappa[2]];
    pi->m[2] = &pi->m[1][pp->kappa[2]];
    pi->m[3] = &pi->m[2][pp->kappa[2]];
  } else if (!pp->tail) {
    pi->m[1] = NULL;
    pi->m[2] = &pi->m[0][pp->kappa[1]];
    pi->m[3] = &pi->m[2][LIFTS];
  } else {
    pi->m[1] = NULL;
    pi->m[2] = &pi->m[0][pp->len[LAB_INCOM] + pp->len[LAB_QUADG]];
    pi->m[3] = &pi->m[2][pp->len[LAB_LIFT]];
  }
}

static size_t py_lnp_proof_polz_len(const lnp_params pp) {
  return 4 * pp->kappa_linfmsis + 256 / N;
}

static void py_lnp_proof_set_ptrs(lnp_proof pi, const lnp_params pp, polz *arr) {
  pi->m[0] = arr;
  pi->m[1] = &pi->m[0][pp->kappa_linfmsis];
  pi->m[2] = &pi->m[0][2 * pp->kappa_linfmsis];
  pi->m[3] = &pi->m[0][2 * pp->kappa_linfmsis + 256 / N];
  pi->m[4] = &pi->m[0][3 * pp->kappa_linfmsis + 256 / N];
}

static int py_wr(uint8_t **out, size_t *rem, const void *src, size_t len) {
  if (*rem < len) {
    return 1;
  }
  if (len != 0) {
    memcpy(*out, src, len);
  }
  *out += len;
  *rem -= len;
  return 0;
}

static int py_rd(const uint8_t **in, size_t *rem, void *dst, size_t len) {
  if (*rem < len) {
    return 1;
  }
  if (len != 0) {
    memcpy(dst, *in, len);
  }
  *in += len;
  *rem -= len;
  return 0;
}

static int py_wr_u64(uint8_t **out, size_t *rem, uint64_t x) {
  uint8_t buf[8];
  for (size_t i = 0; i < sizeof(buf); i++) {
    buf[i] = (uint8_t)(x >> (8 * i));
  }
  return py_wr(out, rem, buf, sizeof(buf));
}

static int py_rd_u64(const uint8_t **in, size_t *rem, uint64_t *x) {
  uint8_t buf[8];
  if (py_rd(in, rem, buf, sizeof(buf))) {
    return 1;
  }
  *x = 0;
  for (size_t i = 0; i < sizeof(buf); i++) {
    *x |= ((uint64_t)buf[i]) << (8 * i);
  }
  return 0;
}

static size_t py_blob_nbytes(size_t payload_len) {
  return sizeof(uint64_t) + payload_len;
}

static int py_wr_blob(uint8_t **out, size_t *rem, const void *src, size_t len) {
  if (py_wr_u64(out, rem, len)) {
    return 1;
  }
  return py_wr(out, rem, src, len);
}

static int py_rd_blob_exact(const uint8_t **in, size_t *rem, void *dst, size_t len) {
  uint64_t got;
  if (py_rd_u64(in, rem, &got) || got != len) {
    return 1;
  }
  return py_rd(in, rem, dst, len);
}

static int py_wr_polz_payload(uint8_t **out, size_t *rem,
                              const polz *src, size_t len) {
  __attribute__((aligned(64))) polz canonical;
  size_t nbytes;

  if (len > SIZE_MAX / POLZBYTES) return 1;
  nbytes = len * POLZBYTES;
  if (*rem < nbytes) return 1;
  for (size_t i = 0; i < len; i++) {
    memcpy(canonical, src[i], sizeof(polz));
    polz_reduce(canonical);
    polz_caddq(canonical);
    polz_bitpack(*out + i * POLZBYTES, canonical);
  }
  *out += nbytes;
  *rem -= nbytes;
  return 0;
}

static int py_rd_polz_payload(const uint8_t **in, size_t *rem,
                              polz *dst, size_t len) {
  __attribute__((aligned(64))) uint8_t packed[POLZBYTES];
  size_t nbytes;

  if (len > SIZE_MAX / POLZBYTES) return 1;
  nbytes = len * POLZBYTES;
  if (*rem < nbytes) return 1;
  for (size_t i = 0; i < len; i++) {
    /* polz_bitunpack uses aligned AVX-512 loads, while an inner proof blob can
       begin at any byte offset in its envelope/socket frame. */
    memcpy(packed, *in + i * POLZBYTES, POLZBYTES);
    polz_bitunpack(dst[i], packed);
  }
  *in += nbytes;
  *rem -= nbytes;
  return 0;
}

static int py_wr_polz_blob(uint8_t **out, size_t *rem,
                           const polz *src, size_t len) {
  if (len > SIZE_MAX / POLZBYTES ||
      py_wr_u64(out, rem, len * POLZBYTES)) return 1;
  return py_wr_polz_payload(out, rem, src, len);
}

static int py_rd_polz_blob(const uint8_t **in, size_t *rem,
                           polz *dst, size_t len) {
  uint64_t got;

  if (len > SIZE_MAX / POLZBYTES || py_rd_u64(in, rem, &got) ||
      got != len * POLZBYTES) return 1;
  return py_rd_polz_payload(in, rem, dst, len);
}

/* LNP m2 is a centered short response whose signed representative is checked
   directly by the verifier.  The four commitment vectors are ring elements
   and can use 38-bit residues; m2 deliberately remains in raw limb form. */
static size_t py_lnp_field_blob_nbytes(const lnp_params pp) {
  return 4 * pp->kappa_linfmsis * POLZBYTES +
         (256 / N) * sizeof(polz);
}

static int py_wr_lnp_field_blob(uint8_t **out, size_t *rem,
                                const lnp_params pp, const lnp_proof pi) {
  const size_t k = pp->kappa_linfmsis;
  const size_t centered = 256 / N;
  const size_t payload = py_lnp_field_blob_nbytes(pp);

  if (py_wr_u64(out, rem, payload) ||
      py_wr_polz_payload(out, rem, pi->m[0], 2 * k) ||
      py_wr(out, rem, pi->m[2], centered * sizeof(polz)) ||
      py_wr_polz_payload(out, rem, pi->m[3], 2 * k)) return 1;
  return 0;
}

static int py_rd_lnp_field_blob(const uint8_t **in, size_t *rem,
                                const lnp_params pp, lnp_proof pi) {
  const size_t k = pp->kappa_linfmsis;
  const size_t centered = 256 / N;
  uint64_t got;

  if (py_rd_u64(in, rem, &got) || got != py_lnp_field_blob_nbytes(pp) ||
      py_rd_polz_payload(in, rem, pi->m[0], 2 * k) ||
      py_rd(in, rem, pi->m[2], centered * sizeof(polz)) ||
      py_rd_polz_payload(in, rem, pi->m[3], 2 * k)) return 1;
  return 0;
}

void py_comkey_reset(void) {
  comkey_reset();
}

size_t py_comkey_nbytes(void) {
  size_t len = comkey_nonce ? comkey->len : 0;
  size_t nbytes = sizeof(uint64_t);
  nbytes += py_blob_nbytes(len * sizeof(double));
  for (size_t i = 0; i < K; i++) {
    nbytes += py_blob_nbytes(len * sizeof(poly));
  }
  return nbytes;
}

int py_comkey_serialize(uint8_t *out, size_t outlen) {
  uint8_t *cur = out;
  size_t rem = outlen;
  size_t len = comkey_nonce ? comkey->len : 0;

  if (outlen != py_comkey_nbytes()) {
    return 1;
  }
  if (py_wr_u64(&cur, &rem, len) ||
      py_wr_blob(&cur, &rem, len ? comkey->widths : NULL, len * sizeof(double))) {
    return 2;
  }
  for (size_t i = 0; i < K; i++) {
    if (py_wr_blob(&cur, &rem, len ? comkey->proj[i] : NULL, len * sizeof(poly))) {
      return 3;
    }
  }
  return rem == 0 ? 0 : 4;
}

int py_comkey_deserialize(const uint8_t *in, size_t inlen) {
  const uint8_t *cur = in;
  size_t rem = inlen;
  uint64_t len64;

  if (py_rd_u64(&cur, &rem, &len64)) {
    return 1;
  }
  comkey_reset();
  if (len64 == 0) {
    return rem == 0 ? 0 : 2;
  }
  polxvec_init(comkey, len64, 1);
  if (py_rd_blob_exact(&cur, &rem, comkey->widths, len64 * sizeof(double))) {
    return 3;
  }
  for (size_t i = 0; i < K; i++) {
    if (py_rd_blob_exact(&cur, &rem, comkey->proj[i], len64 * sizeof(poly))) {
      return 4;
    }
  }
  comkey_nonce = (len64 + 31) / 32;
  return rem == 0 ? 0 : 5;
}

void py_init_witness(witness wt, size_t r, size_t n[]){
  size_t i, nn;

  witness_init(wt, r, r);

  nn = 0;
  for(i=0;i<r;i++){
    wt->n[i] = n[i];
    nn += n[i];
  }

  wt->s[0] = _aligned_alloc(64, nn * sizeof(poly));
  for(i=1;i<wt->r;i++){
    wt->s[i] = &wt->s[i-1][wt->n[i-1]];
  }
}

int py_set_witness_vector(
  witness wt, 
  size_t idx, 
  size_t n, 
  size_t deg, 
  const int64_t s[]
)
{
  if(idx >= wt->r) {
    fprintf(stderr,"ERROR in py_set_witness_vector(): "
                   "Witness vector %zu does not exist\n", idx);
    return 1;
  }
  if(n*deg != wt->n[idx]) {
    fprintf(stderr,"ERROR in py_set_witness_vector(): "
                   "Mismatch of witness vector length\n");
    return 2;
  }

  polyvec_fromint64vec(wt->s[idx], s, 1, n, deg, NULL);
  return 0;
}

int py_print_witness_vector(witness wt, size_t idx){
  if(idx >= wt->r){
    fprintf(stderr,"ERROR in py_print_witness_vector(): "
                   "Witness vector %zu does not exist\n", idx);
    return 1;
  }
  polyvec_print(wt->s[idx], 1, wt->n[idx]);
  return 0;
}

void py_init_statement(
  statement st, 
  size_t r, 
  size_t n[], 
  uint64_t normsq[],
  uint64_t normsq_req[], 
  normtype normty[], 
  size_t num_rq_cnst,
  size_t num_zq_cnst,
  size_t num_int_cnst
)
{
  size_t i;

  statement_init(st, r, r);
  for(i=0;i<r;i++){
    st->n[i] = n[i];
    st->normsq[i] = normsq[i];
    st->normsq_req[i] = normsq_req[i];
    st->normty[i] = normty[i];
  }
  
  rqcnstset_init(st->rqcnst, num_rq_cnst, 0);
  st->rqcnst->nsparse = 0; // update when appending
  zqcnstset_init(st->zqcnst, 0, num_zq_cnst, 0, 0, num_int_cnst);
  st->zqcnst->nint = 0; // update when appending

  memset(st->h, 0, sizeof(st->h));

  comkey_init(1);
}

void py_statement_hash(uint8_t out[16], const statement st) {
  memcpy(out, st->h, 16);
}

int py_append_constraint(
  statement st,
  size_t nvec,
  const size_t idx[],
  const size_t n[],
  size_t deg,
  int64_t *phi,
  int64_t *b,
  int full
)
{
  size_t i, j, *off, off_tmp, *phi_off, total_phi;
  sparsecnst *cnst;
  double width_unif = ldexp(1, 2*LOGQ);
  __attribute__((aligned(16)))
  uint8_t hashbuf[deg * POLZBYTES];
  polz t[deg];
  shake128incctx shakectx;
  int statement_threads;

  if(!full && deg > 1){
    fprintf(stderr, "ERROR in py_append_constraint(): "
                    "Degree of integer constraint is too large\n");
    return 1;
  }

  for(i=0;i<nvec;i++){
    j = idx[i];
    if(j >= st->r){
      fprintf(stderr, "ERROR in py_append_constraint(): "
                      "Witness vector %zu does not exist\n", j);
      return 2;
    }
    if(n[i]*deg != st->n[j]){
      fprintf(stderr, "ERROR in py_append_constraint(): "
                      "Mismatch in witness vector length (vector %zu)\n", j);
      return 3;
    }
  }

  off = _malloc(st->r * sizeof(size_t));

  off[0] = 0;
  off_tmp = 0;
  for(i=1;i<st->r;i++){
    off_tmp += st->n[i-1];
    off[i] = off_tmp;
  }
  
  if(full){
    cnst = &st->rqcnst->sparse[st->rqcnst->nsparse];
    st->rqcnst->nsparse++;
    st->rqcnst->sparse_nchal += deg;
  }
  else{
    cnst = &st->zqcnst->sparse[st->zqcnst->nsparse];
    st->zqcnst->nsparse++;
    st->zqcnst->sparse_nchal += deg;
  }
  
  sparsecnst_init(*cnst, deg);
  linfunc_init((*cnst)->lin, deg, nvec, nvec);

  phi_off = _malloc((nvec + 1) * sizeof(size_t));
  phi_off[0] = 0;
  for(i=0;i<nvec;i++) phi_off[i+1] = phi_off[i] + n[i];
  total_phi = phi_off[nvec];
  if(total_phi > SIZE_MAX / deg ||
     total_phi * deg > SIZE_MAX / POLZBYTES){
    free(phi_off);
    free(off);
    return 4;
  }
  if(!py_q96_compact_hash_enabled()){
    (*cnst)->transcript_lin_len = total_phi * deg * POLZBYTES;
    if((*cnst)->transcript_lin_len != 0)
      (*cnst)->transcript_lin = _malloc((*cnst)->transcript_lin_len);
  }

  shake128_inc_init(&shakectx);
  shake128_inc_absorb(&shakectx,st->h,16);

  memset(t, 0, sizeof(polz) * deg);
  polzvec_fromint64vec(t, 1, deg, b);
  polzvec_topolxvec((*cnst)->b, t, 0, 1, deg);
  polxvec_setwidths1((*cnst)->b, 0, 1, deg, width_unif);
  polzvec_bitpack(hashbuf, t, deg);
  shake128_inc_absorb(&shakectx, hashbuf, deg * POLZBYTES);
  
  for(i=0;i<nvec;i++){
    j = idx[i];
    (*cnst)->lin->off[i] = off[j];
    polxvec_init((*cnst)->lin->phi[i], st->n[j], 1);
  }

  statement_threads = py_statement_threads(nvec);
#ifdef _OPENMP
#pragma omp parallel for num_threads(statement_threads) \
    if(statement_threads > 1 && nvec > 1) schedule(static)
#endif
  for(i=0;i<nvec;i++){
    const int64_t *part = phi + phi_off[i] * deg * N;
    polz values[deg];

    polxvec_fromint64vec((*cnst)->lin->phi[i], part, n[i], deg, width_unif);
    for(size_t k=0;k<n[i];k++){
      const int64_t *poly_values = part + k * deg * N;
      /* Preserve exact structural zeros in the width metadata.  The old
       * blanket width_unif assignment marked padding coordinates as dense,
       * hiding scalar extension-ring coefficients from the native
       * aggregation fast path.  Width zero is exact here because the input
       * coefficients are still available in the flat input. */
      for(size_t l=0;l<deg;l++){
        size_t coefficient;
        int nonzero = 0;
        for(coefficient=0;coefficient<N;coefficient++){
          nonzero |= poly_values[l*N + coefficient] != 0;
        }
        (*cnst)->lin->phi[i]->widths[deg*k + l] =
            nonzero ? width_unif : 0;
      }
      if(!py_q96_compact_hash_enabled()){
        memset(values, 0, sizeof(polz) * deg);
        polzvec_fromint64vec(values, 1, deg, poly_values);
        polzvec_bitpack(
            (*cnst)->transcript_lin + (phi_off[i] + k) * deg * POLZBYTES,
            values, deg);
      }
    }
  }
  if((*cnst)->transcript_lin_len != 0)
    shake128_inc_absorb(&shakectx, (*cnst)->transcript_lin,
                        (*cnst)->transcript_lin_len);

  shake128_inc_finalize(&shakectx);
  shake128_inc_squeeze(st->h, 16, &shakectx);

  free(phi_off);
  free(off);
  return 0;
}

static int py_append_quadratic_monomials_impl(
  statement st,
  size_t nlin,
  size_t nprod,
  const size_t idx_lin[],
  const size_t idx_prod1[],
  const size_t idx_prod2[],
  const size_t len_phi[],
  size_t deg,
  int64_t *a,
  const uint32_t *monomial_exponents,
  int64_t *phi,
  int64_t *b
)
{
  size_t i, j, k;
  sparsecnst *cnst;
  polz t;
  int64_t quadcoeffs[N] = {};
  shake128incctx shakectx;
  int ret;

  ret = py_append_constraint(st, nlin, idx_lin, len_phi, deg, phi, b, 1);

  if(ret){
    fprintf(stderr, "ERROR in py_append_quadratic(): "
                    "error while creating the linear part of the constraint\n");
    return 1;
  }

  for(i=0;i<nprod;i++){
    j = idx_prod1[i];
    k = idx_prod2[i];
    if(monomial_exponents != NULL && monomial_exponents[i] >= N){
      fprintf(stderr, "ERROR in py_append_quadratic(): "
                      "Quadratic monomial exponent is out of range\n");
      return 5;
    }
    if(j >= st->r){
      fprintf(stderr, "ERROR in py_append_quadratic(): "
                      "Witness vector %zu does not exist\n", j);
      return 2;
    }
    if(k >= st->r){
      fprintf(stderr, "ERROR in py_append_quadratic(): "
                      "Witness vector %zu does not exist\n", k);
      return 3;
    }
    if(st->n[j] != st->n[k]){
      fprintf(stderr, "ERROR in py_append_quadratic(): "
                      "Mismatch in witness lengths (vectors %zu, %zu)\n", j, k);
      return 4;
    }
  }

  cnst = &st->rqcnst->sparse[st->rqcnst->nsparse-1];

  quadfunc_init((*cnst)->quad, nprod, nprod);

  /* In the q88 deployment profile, st->h is reset to a domain-separated
   * digest of the PP, link-matrix context and compiler ABI before dynamic
   * constraints are appended.  That digest binds these deterministic
   * coefficients, so only the canonical dynamic RHS is chained here. */
  if(py_q96_compact_hash_enabled()){
    for(i=0;i<nprod;i++){
      j = idx_prod1[i];
      k = idx_prod2[i];
      if(st->normsq[j] <= st->normsq[k]){
        (*cnst)->quad->rows[i] = j;
        (*cnst)->quad->cols[i] = k;
      }
      else{
        (*cnst)->quad->rows[i] = k;
        (*cnst)->quad->cols[i] = j;
      }
      polx_monomial((*cnst)->quad->coeffs[i],
                    monomial_exponents == NULL ? 0 : monomial_exponents[i],
                    a[i]);
    }
    return 0;
  }

  (*cnst)->transcript_quad_len = nprod * POLZBYTES;
  if((*cnst)->transcript_quad_len != 0)
    (*cnst)->transcript_quad = _malloc((*cnst)->transcript_quad_len);

  shake128_inc_init(&shakectx);
  shake128_inc_absorb(&shakectx,st->h,16);
  
  for(i=0;i<nprod;i++){
    if(st->normsq[idx_prod1[i]] <= st->normsq[idx_prod2[i]]){
      (*cnst)->quad->rows[i] = idx_prod1[i];
      (*cnst)->quad->cols[i] = idx_prod2[i];
    }
    else{
      (*cnst)->quad->rows[i] = idx_prod2[i];
      (*cnst)->quad->cols[i] = idx_prod1[i];
    }
    
    polx_monomial((*cnst)->quad->coeffs[i],
                  monomial_exponents == NULL ? 0 : monomial_exponents[i],
                  a[i]);

    memset(quadcoeffs, 0, sizeof(quadcoeffs));
    quadcoeffs[monomial_exponents == NULL ? 0 : monomial_exponents[i]] = a[i];
    memset(&t, 0, sizeof(t));
    polzvec_fromint64vec(&t, 1, 1, quadcoeffs);
    polzvec_bitpack((*cnst)->transcript_quad + i * POLZBYTES, &t, 1);
  }

  if((*cnst)->transcript_quad_len != 0)
    shake128_inc_absorb(&shakectx, (*cnst)->transcript_quad,
                        (*cnst)->transcript_quad_len);

  shake128_inc_finalize(&shakectx);
  shake128_inc_squeeze(st->h, 16, &shakectx);

  return 0;
}

int py_append_quadratic(
  statement st, size_t nlin, size_t nprod, const size_t idx_lin[],
  const size_t idx_prod1[], const size_t idx_prod2[],
  const size_t len_phi[], size_t deg, int64_t *a, int64_t *phi, int64_t *b
)
{
  return py_append_quadratic_monomials_impl(
      st, nlin, nprod, idx_lin, idx_prod1, idx_prod2, len_phi, deg, a,
      NULL, phi, b);
}

int py_append_quadratic_monomials(
  statement st, size_t nlin, size_t nprod, const size_t idx_lin[],
  const size_t idx_prod1[], const size_t idx_prod2[],
  const size_t len_phi[], size_t deg, int64_t *a,
  const uint32_t *monomial_exponents, int64_t *phi, int64_t *b
)
{
  if(nprod != 0 && monomial_exponents == NULL) return 6;
  return py_append_quadratic_monomials_impl(
      st, nlin, nprod, idx_lin, idx_prod1, idx_prod2, len_phi, deg, a,
      monomial_exponents, phi, b);
}

int py_append_deg0_constraint(
  statement st,
  size_t idx,
  size_t deg
)
{
  size_t i, off;
  intcnst *cnst;

  if(idx >= st->r){
    fprintf(stderr, "ERROR in py_append_deg0_constraint(): "
                    "Witness vector %zu does not exist\n", idx);
    return 1;
  }
  if(st->n[idx] != deg){
    fprintf(stderr, "ERROR in py_append_deg0_constraint(): "
                    "Mismatch in witness degree (vector %zu\n)", idx);
    return 2;
  }

  cnst = &st->zqcnst->intc[st->zqcnst->nint];
  st->zqcnst->nint++;
  st->zqcnst->int_nchal += deg*N-1;

  off = 0;
  for(i=0;i<idx;i++){
    off += st->n[i];
  }
  (*cnst)->off = off;
  (*cnst)->rank = deg;

  return 0;
}

/* Add one degree-zero predicate for every polynomial in a witness vector.
 * A single rank>1 intcnst means "one integer followed by zero polynomials",
 * so it is not equivalent to this batched operation. */
int py_append_deg0_vector_constraints(
  statement st,
  size_t idx,
  size_t count
)
{
  size_t i, off;
  intcnst *cnst;

  if(idx >= st->r || count == 0 || st->n[idx] != count){
    fprintf(stderr, "ERROR in py_append_deg0_vector_constraints(): "
                    "invalid witness vector or count\n");
    return 1;
  }

  off = 0;
  for(i=0;i<idx;i++) off += st->n[i];
  for(i=0;i<count;i++){
    cnst = &st->zqcnst->intc[st->zqcnst->nint++];
    (*cnst)->off = off + i;
    (*cnst)->rank = 1;
    st->zqcnst->int_nchal += N - 1;
  }
  return 0;
}

/* Add one degree-zero predicate for every logical extension-ring polynomial.
 * Each logical polynomial occupies `deg` consecutive base-ring polynomials in
 * the witness.  This is the vector analogue of py_append_deg0_constraint(). */
int py_append_deg0_logical_vector_constraints(
  statement st,
  size_t idx,
  size_t count,
  size_t deg
)
{
  size_t i, off;
  intcnst *cnst;

  if(idx >= st->r || count == 0 || deg == 0 ||
     count > SIZE_MAX / deg || st->n[idx] != count * deg){
    fprintf(stderr, "ERROR in py_append_deg0_logical_vector_constraints(): "
                    "invalid witness vector, count, or degree\n");
    return 1;
  }

  off = 0;
  for(i=0;i<idx;i++) off += st->n[i];
  for(i=0;i<count;i++){
    cnst = &st->zqcnst->intc[st->zqcnst->nint++];
    (*cnst)->off = off + i * deg;
    (*cnst)->rank = deg;
    st->zqcnst->int_nchal += deg * N - 1;
  }
  return 0;
}

/* Deterministically expand a seed-defined matrix bank into an int32 buffer.
 * This is a bridge primitive: it mirrors the application's SHAKE256 block
 * derivation and SHAKE128 q-uniform rejection sampler, while avoiding tens of
 * thousands of Python hashlib calls during worker cold start. */
int py_seeded_uniform_matrix_bank(
  int32_t *out,
  const uint8_t seed[32],
  size_t matrices,
  size_t blocks,
  size_t rows,
  size_t cols,
  size_t degree,
  uint32_t q,
  uint32_t qbits,
  uint32_t domain_separator
)
{
  size_t mb;
  const uint64_t mask = (((uint64_t)1) << qbits) - 1;

  if(out == NULL || seed == NULL || matrices == 0 || blocks == 0 ||
     rows == 0 || cols == 0 || degree == 0 || qbits == 0 || qbits > 20 ||
     q == 0) return 1;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(mb=0;mb<matrices*blocks;mb++){
    const size_t matrix_index = mb / blocks;
    const size_t block = mb % blocks;
    uint8_t block_seed[32];
    shake256incctx block_ctx;
    size_t row, col;

    shake256_inc_init(&block_ctx);
    shake256_inc_absorb(&block_ctx, seed, 32);
    shake256_inc_absorb(&block_ctx, (const uint8_t *)&matrix_index,
                        sizeof(matrix_index));
    shake256_inc_absorb(&block_ctx, (const uint8_t *)&block, sizeof(block));
    shake256_inc_finalize(&block_ctx);
    shake256_inc_squeeze(block_seed, sizeof(block_seed), &block_ctx);

    for(row=0;row<rows;row++){
      for(col=0;col<cols;col++){
        uint8_t buf[2*SHAKE128_RATE];
        shake128incctx ctx;
        size_t produced = 0, off = 0, available = sizeof(buf);
        int32_t *poly_out = out +
          ((((matrix_index*blocks + block)*rows + row)*cols + col)*degree);

        shake128_inc_init(&ctx);
        shake128_inc_absorb(&ctx, block_seed, sizeof(block_seed));
        shake128_inc_absorb(&ctx, (const uint8_t *)&domain_separator,
                            sizeof(domain_separator));
        shake128_inc_absorb(&ctx, (const uint8_t *)&row, sizeof(row));
        shake128_inc_absorb(&ctx, (const uint8_t *)&col, sizeof(col));
        shake128_inc_finalize(&ctx);
        shake128_inc_squeeze(buf, sizeof(buf), &ctx);

        while(produced < degree){
          uint64_t word, candidate;
          size_t k;
          if(available < 5){
            for(k=0;k<available;k++) buf[k] = buf[off+k];
            shake128_inc_squeeze(&buf[available], SHAKE128_RATE, &ctx);
            off = 0;
            available += SHAKE128_RATE;
          }
          word = (uint64_t)buf[off] |
                 ((uint64_t)buf[off+1] << 8) |
                 ((uint64_t)buf[off+2] << 16) |
                 ((uint64_t)buf[off+3] << 24) |
                 ((uint64_t)buf[off+4] << 32);
          candidate = word & mask;
          if(candidate < q) poly_out[produced++] = (int32_t)candidate;
          if(produced < degree){
            candidate = (word >> 20) & mask;
            if(candidate < q) poly_out[produced++] = (int32_t)candidate;
          }
          off += 5;
          available -= 5;
        }
      }
    }
  }
  return 0;
}

/* Compose A-hat=A_0+sum r_i A_i over the integers, derive the exact q-lift,
 * signed-decompose the centered value, and write all resulting vectors
 * directly into an already initialized resident witness.  This removes the
 * largest remaining per-request Python integer/list construction path. */
int py_refresh_set_witness(
  witness wt,
  const size_t *chunk_idx,
  const size_t *lift_idx,
  const int32_t *bank,
  const uint8_t *bits,
  size_t bit_count,
  size_t blocks,
  size_t rows,
  size_t cols,
  size_t degree,
  int64_t q,
  int64_t base,
  size_t chunks
)
{
  size_t row, entry, coeff, bit, chunk;
  const size_t entries = blocks * cols;
  const size_t vector_coeffs = entries * degree;
  int64_t *out;
  int error = 0;

  if(wt == NULL || chunk_idx == NULL || lift_idx == NULL || bank == NULL ||
     bits == NULL || bit_count == 0 || blocks == 0 || rows == 0 ||
     cols == 0 || degree != N || q <= 1 || base <= 1 || chunks == 0)
    return 1;
  out = (int64_t *)calloc((chunks + 1) * vector_coeffs, sizeof(int64_t));
  if(out == NULL) return 2;

  for(row=0;row<rows;row++){
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) private(bit, chunk)
#endif
    for(entry=0;entry<entries;entry++){
      for(coeff=0;coeff<degree;coeff++){
        const size_t block = entry / cols;
        const size_t col = entry % cols;
        const size_t bank_stride = blocks * rows * cols * degree;
        const size_t bank_offset =
          (((block * rows + row) * cols + col) * degree + coeff);
        const size_t out_offset = entry * degree + coeff;
        int64_t raw = bank[bank_offset];
        int64_t centered, value;

        for(bit=0;bit<bit_count;bit++){
          if(bits[bit]) raw += bank[(bit + 1) * bank_stride + bank_offset];
        }
        centered = raw % q;
        if(centered > q / 2) centered -= q;
        out[chunks * vector_coeffs + out_offset] =
          (raw - centered) / q;
        value = centered;
        for(chunk=0;chunk<chunks;chunk++){
          int64_t digit = (value + base / 2) % base;
          if(digit < 0) digit += base;
          digit -= base / 2;
          out[chunk * vector_coeffs + out_offset] = digit;
          value = (value - digit) / base;
        }
        if(value != 0){
#ifdef _OPENMP
#pragma omp atomic write
#endif
          error = 3;
        }
      }
    }
    if(error != 0) break;
    for(chunk=0;chunk<chunks;chunk++){
      if(py_set_witness_vector(wt, chunk_idx[row * chunks + chunk],
                               entries, 1,
                               out + chunk * vector_coeffs) != 0){
        error = 4;
        break;
      }
    }
    if(error != 0) break;
    if(py_set_witness_vector(wt, lift_idx[row], entries, 1,
                             out + chunks * vector_coeffs) != 0){
      error = 5;
      break;
    }
  }
  free(out);
  return error;
}

int py_statement_set_norms(statement st, size_t count,
                           const uint64_t *normsq)
{
  size_t i;
  if(st == NULL || normsq == NULL || count != st->r) return 1;
  for(i=0;i<count;i++) st->normsq[i] = normsq[i];
  return 0;
}

/* Append A-hat binding constraints using true sparse selector parts.  The
 * generic Python bridge represented each one-hot selector as a full 48-poly
 * vector.  Here every selector is a single monomial at the corresponding
 * witness offset, preserving the relation while removing zero storage and
 * zero multiplications from both Labrador folds. */
int py_append_refresh_constraints(
  statement st,
  size_t r_idx,
  const size_t *chunk_idx,
  const size_t *lift_idx,
  const int32_t *bank,
  size_t bit_count,
  size_t blocks,
  size_t rows,
  size_t cols,
  size_t degree,
  int64_t q,
  int64_t base,
  size_t chunks
)
{
  size_t *woff;
  size_t row, entry, bit, chunk, coeff, i;
  const size_t entries = blocks * cols;
  const size_t bank_stride = blocks * rows * cols * degree;
  const double width_unif = ldexp(1, 2*LOGQ);
  int64_t coeff64[N];
  uint8_t hashbuf[POLZBYTES];
  polz integer;

  if(st == NULL || chunk_idx == NULL || lift_idx == NULL || bank == NULL ||
     r_idx >= st->r || st->n[r_idx] != bit_count || degree != N ||
     bit_count == 0 || blocks == 0 || rows == 0 || cols == 0 ||
     chunks == 0 || q <= 1 || base <= 1) return 1;
  for(row=0;row<rows;row++){
    if(lift_idx[row] >= st->r || st->n[lift_idx[row]] != entries) return 2;
    for(chunk=0;chunk<chunks;chunk++){
      const size_t idx = chunk_idx[row * chunks + chunk];
      if(idx >= st->r || st->n[idx] != entries) return 3;
    }
  }

  woff = _malloc(st->r * sizeof(size_t));
  woff[0] = 0;
  for(i=1;i<st->r;i++) woff[i] = woff[i-1] + st->n[i-1];

  for(row=0;row<rows;row++){
    for(entry=0;entry<entries;entry++){
      sparsecnst *slot = &st->rqcnst->sparse[st->rqcnst->nsparse++];
      shake128incctx hashctx;
      const size_t block = entry / cols;
      const size_t col = entry % cols;
      const size_t bank_offset =
        (((block * rows + row) * cols + col) * degree);
      int64_t scale = 1;

      st->rqcnst->sparse_nchal += 1;
      sparsecnst_init(*slot, 1);
      linfunc_init((*slot)->lin, 1, chunks + 2, chunks + 2);
      shake128_inc_init(&hashctx);
      shake128_inc_absorb(&hashctx, st->h, sizeof(st->h));

      for(coeff=0;coeff<degree;coeff++)
        coeff64[coeff] = -(int64_t)bank[bank_offset + coeff];
      polzvec_fromint64vec(&integer, 1, 1, coeff64);
      polzvec_topolxvec((*slot)->b, &integer, 0, 1, 1);
      polxvec_setwidths1((*slot)->b, 0, 1, 1, width_unif);
      polzvec_bitpack(hashbuf, &integer, 1);
      shake128_inc_absorb(&hashctx, hashbuf, sizeof(hashbuf));

      (*slot)->lin->off[0] = woff[r_idx];
      polxvec_init((*slot)->lin->phi[0], bit_count, 1);
      for(bit=0;bit<bit_count;bit++){
        for(coeff=0;coeff<degree;coeff++)
          coeff64[coeff] = bank[(bit + 1) * bank_stride +
                                bank_offset + coeff];
        polzvec_fromint64vec(&integer, 1, 1, coeff64);
        polzvec_topolxvec((*slot)->lin->phi[0], &integer, bit, 1, 1);
        polzvec_bitpack(hashbuf, &integer, 1);
        shake128_inc_absorb(&hashctx, hashbuf, sizeof(hashbuf));
      }
      polxvec_setwidths1((*slot)->lin->phi[0], 0, 1, bit_count, width_unif);

      memset(coeff64, 0, sizeof(coeff64));
      for(chunk=0;chunk<chunks;chunk++){
        (*slot)->lin->off[chunk + 1] =
          woff[chunk_idx[row * chunks + chunk]] + entry;
        polxvec_init((*slot)->lin->phi[chunk + 1], 1, 1);
        coeff64[0] = -scale;
        polzvec_fromint64vec(&integer, 1, 1, coeff64);
        polzvec_topolxvec((*slot)->lin->phi[chunk + 1], &integer, 0, 1, 1);
        polxvec_setwidths1((*slot)->lin->phi[chunk + 1], 0, 1, 1,
                          width_unif);
        polzvec_bitpack(hashbuf, &integer, 1);
        shake128_inc_absorb(&hashctx, hashbuf, sizeof(hashbuf));
        coeff64[0] = 0;
        scale *= base;
      }
      (*slot)->lin->off[chunks + 1] = woff[lift_idx[row]] + entry;
      polxvec_init((*slot)->lin->phi[chunks + 1], 1, 1);
      coeff64[0] = -q;
      polzvec_fromint64vec(&integer, 1, 1, coeff64);
      polzvec_topolxvec((*slot)->lin->phi[chunks + 1], &integer, 0, 1, 1);
      polxvec_setwidths1((*slot)->lin->phi[chunks + 1], 0, 1, 1, width_unif);
      polzvec_bitpack(hashbuf, &integer, 1);
      shake128_inc_absorb(&hashctx, hashbuf, sizeof(hashbuf));
      shake128_inc_finalize(&hashctx);
      shake128_inc_squeeze(st->h, sizeof(st->h), &hashctx);
    }
  }
  free(woff);
  return 0;
}

/* q88 stores every A-hat radix chunk either as a blocks-wide vector for one
 * CRT-lane/phase pair, or as a blocks*phases-wide vector containing complete
 * degree-(phases*N) polynomials.  Append the same coefficient-wise
 * equalities as the generic bridge above, but address the selected base-ring
 * polynomial directly instead of materializing a one-hot coefficient
 * vector.  The packed form lets the full-ring compiler keep the hidden
 * Remark-4.3 selector bits in the base ring (bit_count backend polynomials)
 * rather than repeating every bit in all extension phases. */
int py_append_refresh_phase_constraints(
  statement st,
  size_t r_idx,
  const size_t *chunk_idx,
  const size_t *lift_idx,
  const int32_t *bank,
  size_t bit_count,
  size_t blocks,
  size_t lanes,
  size_t phases,
  size_t degree,
  const uint32_t *moduli,
  int64_t base,
  size_t chunks
)
{
  size_t *woff;
  size_t lane, phase, block, bit, chunk, coeff, i;
  const size_t bank_matrix_stride = lanes * phases * blocks * degree;
  const double width_unif = ldexp(1, 2*LOGQ);
  int64_t coeff64[N];
  uint8_t hashbuf[POLZBYTES];
  polz integer;

  if(st == NULL || chunk_idx == NULL || lift_idx == NULL || bank == NULL ||
     moduli == NULL || r_idx >= st->r || st->n[r_idx] != bit_count ||
     degree != N || bit_count == 0 || blocks == 0 || lanes == 0 ||
     phases == 0 || chunks == 0 || base <= 1) return 1;
  for(lane=0;lane<lanes;lane++){
    if(moduli[lane] <= 1) return 2;
    for(phase=0;phase<phases;phase++){
      const size_t pair = lane * phases + phase;
      if(lift_idx[pair] >= st->r ||
         (st->n[lift_idx[pair]] != blocks &&
          st->n[lift_idx[pair]] != blocks * phases))
        return 3;
      for(chunk=0;chunk<chunks;chunk++){
        const size_t idx = chunk_idx[pair * chunks + chunk];
        if(idx >= st->r ||
           (st->n[idx] != blocks && st->n[idx] != blocks * phases))
          return 4;
      }
    }
  }

  woff = _malloc(st->r * sizeof(size_t));
  woff[0] = 0;
  for(i=1;i<st->r;i++) woff[i] = woff[i-1] + st->n[i-1];

  for(lane=0;lane<lanes;lane++){
    for(phase=0;phase<phases;phase++){
      const size_t pair = lane * phases + phase;
      for(block=0;block<blocks;block++){
        sparsecnst *slot = &st->rqcnst->sparse[st->rqcnst->nsparse++];
        shake128incctx hashctx;
        const size_t bank_offset =
          (((lane * phases + phase) * blocks + block) * degree);
        int64_t scale = 1;

        st->rqcnst->sparse_nchal += 1;
        sparsecnst_init(*slot, 1);
        linfunc_init((*slot)->lin, 1, chunks + 2, chunks + 2);
        shake128_inc_init(&hashctx);
        shake128_inc_absorb(&hashctx, st->h, sizeof(st->h));

        for(coeff=0;coeff<degree;coeff++)
          coeff64[coeff] = -(int64_t)bank[bank_offset + coeff];
        polzvec_fromint64vec(&integer, 1, 1, coeff64);
        polzvec_topolxvec((*slot)->b, &integer, 0, 1, 1);
        polxvec_setwidths1((*slot)->b, 0, 1, 1, width_unif);
        polzvec_bitpack(hashbuf, &integer, 1);
        shake128_inc_absorb(&hashctx, hashbuf, sizeof(hashbuf));

        (*slot)->lin->off[0] = woff[r_idx];
        polxvec_init((*slot)->lin->phi[0], bit_count, 1);
        for(bit=0;bit<bit_count;bit++){
          for(coeff=0;coeff<degree;coeff++)
            coeff64[coeff] = bank[(bit + 1) * bank_matrix_stride +
                                  bank_offset + coeff];
          polzvec_fromint64vec(&integer, 1, 1, coeff64);
          polzvec_topolxvec((*slot)->lin->phi[0], &integer, bit, 1, 1);
          polzvec_bitpack(hashbuf, &integer, 1);
          shake128_inc_absorb(&hashctx, hashbuf, sizeof(hashbuf));
        }
        polxvec_setwidths1((*slot)->lin->phi[0], 0, 1, bit_count,
                          width_unif);

        memset(coeff64, 0, sizeof(coeff64));
        for(chunk=0;chunk<chunks;chunk++){
          const size_t idx = chunk_idx[pair * chunks + chunk];
          const size_t selected = st->n[idx] == blocks
                                  ? block : block * phases + phase;
          (*slot)->lin->off[chunk + 1] = woff[idx] + selected;
          polxvec_init((*slot)->lin->phi[chunk + 1], 1, 1);
          coeff64[0] = -scale;
          polzvec_fromint64vec(&integer, 1, 1, coeff64);
          polzvec_topolxvec((*slot)->lin->phi[chunk + 1], &integer,
                           0, 1, 1);
          polxvec_setwidths1((*slot)->lin->phi[chunk + 1], 0, 1, 1,
                            width_unif);
          polzvec_bitpack(hashbuf, &integer, 1);
          shake128_inc_absorb(&hashctx, hashbuf, sizeof(hashbuf));
          coeff64[0] = 0;
          scale *= base;
        }
        {
          const size_t idx = lift_idx[pair];
          const size_t selected = st->n[idx] == blocks
                                  ? block : block * phases + phase;
          (*slot)->lin->off[chunks + 1] = woff[idx] + selected;
        }
        polxvec_init((*slot)->lin->phi[chunks + 1], 1, 1);
        coeff64[0] = -(int64_t)moduli[lane];
        polzvec_fromint64vec(&integer, 1, 1, coeff64);
        polzvec_topolxvec((*slot)->lin->phi[chunks + 1], &integer,
                         0, 1, 1);
        polxvec_setwidths1((*slot)->lin->phi[chunks + 1], 0, 1, 1,
                          width_unif);
        polzvec_bitpack(hashbuf, &integer, 1);
        shake128_inc_absorb(&hashctx, hashbuf, sizeof(hashbuf));
        shake128_inc_finalize(&hashctx);
        shake128_inc_squeeze(st->h, sizeof(st->h), &hashctx);
      }
    }
  }
  free(woff);
  return 0;
}

/* Signed-decompose exact C-generated lift polynomials and place every chunk
 * directly in its resident witness slot. */
int py_decompose_set_witness(witness wt, const size_t *indices,
                             const int64_t *source, size_t poly_count,
                             size_t degree, int64_t base, size_t chunks)
{
  size_t poly_index, coeff, chunk;
  int64_t *out;
  int error = 0;
  if(wt == NULL || indices == NULL || source == NULL || poly_count == 0 ||
     degree != N || base <= 1 || chunks == 0) return 1;
  out = _aligned_alloc(64, chunks * degree * sizeof(int64_t));
  for(poly_index=0;poly_index<poly_count;poly_index++){
    for(coeff=0;coeff<degree;coeff++){
      int64_t value = source[poly_index * degree + coeff];
      for(chunk=0;chunk<chunks;chunk++){
        int64_t digit = (value + base / 2) % base;
        if(digit < 0) digit += base;
        digit -= base / 2;
        out[chunk * degree + coeff] = digit;
        value = (value - digit) / base;
      }
      if(value != 0) error = 2;
    }
    if(error != 0) break;
    for(chunk=0;chunk<chunks;chunk++){
      if(py_set_witness_vector(
           wt, indices[poly_index * chunks + chunk], 1, 1,
           out + chunk * degree) != 0){
        error = 3;
        break;
      }
    }
    if(error != 0) break;
  }
  free(out);
  return error;
}

int py_malloc_trim(void)
{
#if defined(__GLIBC__)
  return malloc_trim(0);
#else
  return 0;
#endif
}

int py_malloc_configure(size_t mmap_threshold, size_t trim_threshold,
                        int arena_max)
{
#if defined(__GLIBC__)
  int ok = 1;
  if(mmap_threshold > 0)
    ok &= mallopt(PY_M_MMAP_THRESHOLD, (int)mmap_threshold);
  if(trim_threshold > 0)
    ok &= mallopt(PY_M_TRIM_THRESHOLD, (int)trim_threshold);
  if(arena_max > 0) ok &= mallopt(PY_M_ARENA_MAX, arena_max);
  return ok ? 0 : 1;
#else
  (void)mmap_threshold;
  (void)trim_threshold;
  (void)arena_max;
  return 0;
#endif
}

static void py_hash_existing_sparse(statement st, const sparsecnst cnst)
{
  size_t i, j;
  uint8_t hashbuf[POLZBYTES];
  polx value;
  polz integer;
  shake128incctx ctx;

  shake128_inc_init(&ctx);
  shake128_inc_absorb(&ctx, st->h, sizeof(st->h));
  for(j=0;j<cnst->b->len;j++){
    polx_frompolxvec(value, cnst->b, j);
    polz_frompolx(integer, value);
    polzvec_bitpack(hashbuf, &integer, 1);
    shake128_inc_absorb(&ctx, hashbuf, sizeof(hashbuf));
  }
  if(py_q96_compact_hash_enabled()){
    shake128_inc_finalize(&ctx);
    shake128_inc_squeeze(st->h, sizeof(st->h), &ctx);
    return;
  }
  if(cnst->transcript_lin != NULL){
    shake128_inc_absorb(&ctx, cnst->transcript_lin,
                        cnst->transcript_lin_len);
  }
  else{
    for(i=0;i<cnst->lin->nparts;i++){
      for(j=0;j<cnst->lin->phi[i]->len;j++){
        polx_frompolxvec(value, cnst->lin->phi[i], j);
        polz_frompolx(integer, value);
        polzvec_bitpack(hashbuf, &integer, 1);
        shake128_inc_absorb(&ctx, hashbuf, sizeof(hashbuf));
      }
    }
  }
  shake128_inc_finalize(&ctx);
  shake128_inc_squeeze(st->h, sizeof(st->h), &ctx);

  if(cnst->quad->len != 0){
    shake128_inc_init(&ctx);
    shake128_inc_absorb(&ctx, st->h, sizeof(st->h));
    if(cnst->transcript_quad != NULL){
      shake128_inc_absorb(&ctx, cnst->transcript_quad,
                          cnst->transcript_quad_len);
    }
    else{
      for(i=0;i<cnst->quad->len;i++){
        polz_frompolx(integer, cnst->quad->coeffs[i]);
        polzvec_bitpack(hashbuf, &integer, 1);
        shake128_inc_absorb(&ctx, hashbuf, sizeof(hashbuf));
      }
    }
    shake128_inc_finalize(&ctx);
    shake128_inc_squeeze(st->h, sizeof(st->h), &ctx);
  }
}

/* Rebind the request-specific RHS of a static statement tail and reconstruct
 * exactly the same transcript hash that py_append_constraint/quadratic would
 * have produced.  The worker is single-threaded, so the cached statement is
 * never observed while it is being rebound. */
int py_statement_rebind_tail(
  statement st,
  size_t first,
  size_t count,
  const int64_t *rhs,
  const uint8_t static_hash[16]
)
{
  size_t i, rhs_off = 0;
  double width_unif = ldexp(1, 2*LOGQ);

  /* ``tail`` historically required the rebound interval to end at the final
   * sparse constraint.  Also accepting a proper prefix is useful to cache an
   * intermediate transcript hash: callers can rebind a large, infrequently
   * changing prefix once, save st->h, and then rebind a small dynamic tail on
   * every request.  The hash operation itself is unchanged. */
  if(st == NULL || rhs == NULL || static_hash == NULL || count == 0 ||
     first >= st->rqcnst->nsparse ||
     count > st->rqcnst->nsparse - first) return 1;
  memcpy(st->h, static_hash, sizeof(st->h));
  for(i=0;i<count;i++){
    sparsecnst *cnst = &st->rqcnst->sparse[first+i];
    const size_t rank = (*cnst)->b->len;
    polz values[rank];
    if(rank == 0 || rhs_off > SIZE_MAX - rank * N) return 2;
    polzvec_fromint64vec(values, 1, rank, rhs + rhs_off);
    polzvec_topolxvec((*cnst)->b, values, 0, 1, rank);
    polxvec_setwidths1((*cnst)->b, 0, 1, rank, width_unif);
    py_hash_existing_sparse(st, *cnst);
    rhs_off += rank * N;
  }
  return 0;
}

void py_statement_hash_reset(statement st, const uint8_t hash[16])
{
  memcpy(st->h, hash, sizeof(st->h));
}

int py_statement_rebind_constraint_raw(
  statement st,
  size_t index,
  size_t nvec,
  const size_t n[],
  size_t deg,
  const int64_t *phi,
  const int64_t *b,
  size_t nprod,
  const int64_t *a
)
{
  size_t i, k;
  sparsecnst *cnst;
  double width_unif = ldexp(1, 2*LOGQ);
  uint8_t hashbuf[deg * POLZBYTES];
  polz values[deg];
  shake128incctx ctx;

  if(st == NULL || index >= st->rqcnst->nsparse || deg == 0 ||
     phi == NULL || b == NULL || (nprod != 0 && a == NULL)) return 1;
  cnst = &st->rqcnst->sparse[index];
  if((*cnst)->b->len != deg || (*cnst)->lin->nparts != nvec ||
     (*cnst)->quad->len != nprod) return 2;

  /* The raw rebind API can replace the coefficient transcript without
   * replacing the resident coefficient graph.  Drop any append-time cache so
   * a later RHS-only rebind takes the established graph-conversion fallback
   * instead of using bytes for the former coefficients. */
  free((*cnst)->transcript_lin);
  free((*cnst)->transcript_quad);
  (*cnst)->transcript_lin = NULL;
  (*cnst)->transcript_lin_len = 0;
  (*cnst)->transcript_quad = NULL;
  (*cnst)->transcript_quad_len = 0;

  polzvec_fromint64vec(values, 1, deg, b);
  polzvec_topolxvec((*cnst)->b, values, 0, 1, deg);
  polxvec_setwidths1((*cnst)->b, 0, 1, deg, width_unif);

  shake128_inc_init(&ctx);
  shake128_inc_absorb(&ctx, st->h, sizeof(st->h));
  polzvec_bitpack(hashbuf, values, deg);
  shake128_inc_absorb(&ctx, hashbuf, sizeof(hashbuf));
  for(i=0;i<nvec;i++){
    for(k=0;k<n[i];k++){
      memset(values, 0, sizeof(values));
      polzvec_fromint64vec(values, 1, deg, phi);
      polzvec_bitpack(hashbuf, values, deg);
      shake128_inc_absorb(&ctx, hashbuf, sizeof(hashbuf));
      phi += deg*N;
    }
  }
  shake128_inc_finalize(&ctx);
  shake128_inc_squeeze(st->h, sizeof(st->h), &ctx);

  if(nprod != 0){
    shake128_inc_init(&ctx);
    shake128_inc_absorb(&ctx, st->h, sizeof(st->h));
    for(i=0;i<nprod;i++){
      int64_t coeffs[N] = {0};
      coeffs[0] = a[i];
      memset(values, 0, sizeof(values));
      polzvec_fromint64vec(values, 1, 1, coeffs);
      polzvec_bitpack(hashbuf, values, 1);
      shake128_inc_absorb(&ctx, hashbuf, POLZBYTES);
    }
    shake128_inc_finalize(&ctx);
    shake128_inc_squeeze(st->h, sizeof(st->h), &ctx);
  }
  return 0;
}

int py_gen_params(dch_pack_params pp, const statement ist, int zk, int debug){
  size_t pibits;
  timing time;

  timing_start(&time, "Parameter generation");

  if(dch_pack_params_gen(pp, &pibits, ist, zk)){
    return 1;
  }
  if(debug){
    dch_pack_params_print(pp);
  }

  timing_end(&time);
  timing_print(&time, 0);

  printf("预估 proof 大小：%.2f KB\n", ((double) pibits)/8192);

  return 0;
}

int py_simple_verify(const statement st, const witness wt){
  return verify(st, wt);
}

void py_prove(
  dch_pack_proof pi,
  const statement ist, 
  const witness iwt,
  const dch_pack_params pp
)
{
  timing time;

  timing_buffer_init();

  timing_start(&time, "Prover");

  proof_hash_canonical_set(py_proof_field_mask() != 0);
  dch_pack_prove(pi, ist, iwt, pp);

  timing_end(&time);
  timing_print(&time, 0);

  timing_buffer_flush();
}

int py_last_prove_tail_valid(void){
  return dch_pack_last_prove_relation_valid();
}

int py_verify(
  const statement ist, 
  const dch_pack_params pp, 
  const dch_pack_proof pi
)
{
  int ret;
  timing time;

  timing_buffer_init();

  timing_start(&time, "Verifier");

  ret = dch_pack_verify(ist, pp, pi);

  timing_end(&time);
  timing_print(&time, 0);

  timing_buffer_flush();

  return ret;
}

size_t py_proof_nbytes(const dch_pack_params pp, const dch_pack_proof pi) {
  const unsigned int encoding = py_proof_encoding();
  const unsigned int fields = py_proof_field_mask();
  size_t nbytes = sizeof(PY_PROOF_MAGIC);

  nbytes += py_blob_nbytes(pp->pp_dch->kappa_outer *
      ((fields & 1U) ? POLZBYTES : sizeof(polz)));
  nbytes += sizeof(uint64_t);
  for (size_t i = 0; i < pp->pp_pack->np; i++) {
    nbytes += py_blob_nbytes(
        py_pack_lab_proof_polz_len(pp->pp_pack, i) *
        ((fields & 2U) ? POLZBYTES : sizeof(polz)));
    if (encoding >= 4U) {
      const unsigned int bits = py_int32vec_signed_bits(
          pi->pi_pack->p[i]->p, 256);
      const size_t packed = py_signed_int32vec_nbytes(256, bits);
      if (packed == SIZE_MAX) return SIZE_MAX;
      nbytes += 1 + packed;
    } else {
      nbytes += 256 * sizeof(int32_t);
    }
  }
  nbytes += sizeof(uint64_t);
  if (pp->pp_pack->zkp != NULL) {
    nbytes += py_blob_nbytes((fields & 4U) ?
        py_lnp_field_blob_nbytes(pp->pp_pack->zkp[0]) :
        py_lnp_proof_polz_len(pp->pp_pack->zkp[0]) * sizeof(polz));
  }

  nbytes += sizeof(uint64_t);
  for (size_t i = 0; i < pi->pi_pack->owt->r; i++) {
    nbytes += sizeof(uint64_t);
    if (encoding >= 3U) {
      const unsigned int bits = py_polyvec_signed_bits(
          pi->pi_pack->owt->s[i], pi->pi_pack->owt->n[i]);
      const size_t packed = py_signed_polyvec_nbytes(
          pi->pi_pack->owt->n[i], bits);
      if (packed == SIZE_MAX) return SIZE_MAX;
      nbytes += 1 + py_blob_nbytes(packed);
    } else {
      nbytes += py_blob_nbytes(pi->pi_pack->owt->n[i] * sizeof(poly));
    }
  }
  return nbytes;
}

static size_t py_dch_params_nbytes(const dch_params pp) {
  size_t nbytes = sizeof(*pp);
  nbytes += py_blob_nbytes(pp->r * sizeof(size_t));
  nbytes += py_blob_nbytes(pp->r * sizeof(uint64_t));
  nbytes += py_blob_nbytes(pp->r * sizeof(normtype));
  nbytes += py_blob_nbytes(pp->nexact * sizeof(size_t));
  nbytes += py_blob_nbytes(pp->nincom * sizeof(size_t));
  nbytes += py_blob_nbytes(pp->nincom * sizeof(size_t));
  nbytes += py_blob_nbytes(pp->nincom * sizeof(size_t));
  nbytes += py_blob_nbytes(pp->nexact * sizeof(size_t));
  nbytes += py_blob_nbytes(pp->nexact * sizeof(size_t));
  nbytes += py_blob_nbytes(pp->nincom * sizeof(size_t));
  nbytes += py_blob_nbytes(pp->nexact * sizeof(size_t));
  nbytes += py_blob_nbytes(pp->nexact * sizeof(size_t));
  nbytes += py_blob_nbytes(pp->nexact * sizeof(size_t));
  nbytes += py_blob_nbytes(pp->nexact * sizeof(size_t));
  nbytes += py_blob_nbytes(pp->nexact * sizeof(size_t));
  return nbytes;
}

static size_t py_lab_params_nbytes(const lab_params pp) {
  size_t nbytes = sizeof(*pp);
  nbytes += py_blob_nbytes(pp->r_old * sizeof(size_t));
  nbytes += py_blob_nbytes(pp->r * sizeof(size_t));
  nbytes += py_blob_nbytes(pp->r * sizeof(uint64_t));
  return nbytes;
}

size_t py_params_nbytes(const dch_pack_params pp) {
  size_t nbytes = sizeof(PY_PROOF_MAGIC);
  nbytes += py_dch_params_nbytes(pp->pp_dch);
  nbytes += sizeof(*pp->pp_pack);
  nbytes += sizeof(uint64_t);
  for (size_t i = 0; i < pp->pp_pack->np; i++) {
    nbytes += py_lab_params_nbytes(pp->pp_pack->p[i]);
  }
  nbytes += sizeof(uint64_t);
  if (pp->pp_pack->zkp != NULL) {
    nbytes += sizeof(*pp->pp_pack->zkp[0]);
  }
  return nbytes;
}

static int py_wr_dch_params(uint8_t **cur, size_t *rem, const dch_params pp) {
  if (py_wr(cur, rem, pp, sizeof(*pp)) ||
      py_wr_blob(cur, rem, pp->n, pp->r * sizeof(size_t)) ||
      py_wr_blob(cur, rem, pp->normsq, pp->r * sizeof(uint64_t)) ||
      py_wr_blob(cur, rem, pp->normty, pp->r * sizeof(normtype)) ||
      py_wr_blob(cur, rem, pp->exact_map, pp->nexact * sizeof(size_t)) ||
      py_wr_blob(cur, rem, pp->incom_offw, pp->nincom * sizeof(size_t)) ||
      py_wr_blob(cur, rem, pp->incom_lenw, pp->nincom * sizeof(size_t)) ||
      py_wr_blob(cur, rem, pp->kappa_inner, pp->nincom * sizeof(size_t)) ||
      py_wr_blob(cur, rem, pp->off_exact, pp->nexact * sizeof(size_t)) ||
      py_wr_blob(cur, rem, pp->off_sigma, pp->nexact * sizeof(size_t)) ||
      py_wr_blob(cur, rem, pp->off_incom, pp->nincom * sizeof(size_t)) ||
      py_wr_blob(cur, rem, pp->len_exact_merge, pp->nexact * sizeof(size_t)) ||
      py_wr_blob(cur, rem, pp->idx_exact, pp->nexact * sizeof(size_t)) ||
      py_wr_blob(cur, rem, pp->idx_exact_merge, pp->nexact * sizeof(size_t)) ||
      py_wr_blob(cur, rem, pp->idx_sigma, pp->nexact * sizeof(size_t)) ||
      py_wr_blob(cur, rem, pp->nparts_exact_merge, pp->nexact * sizeof(size_t))) {
    return 1;
  }
  return 0;
}

static int py_wr_lab_params(uint8_t **cur, size_t *rem, const lab_params pp) {
  if (py_wr(cur, rem, pp, sizeof(*pp)) ||
      py_wr_blob(cur, rem, pp->rr, pp->r_old * sizeof(size_t)) ||
      py_wr_blob(cur, rem, pp->n, pp->r * sizeof(size_t)) ||
      py_wr_blob(cur, rem, pp->normsq, pp->r * sizeof(uint64_t))) {
    return 1;
  }
  return 0;
}

int py_params_serialize(uint8_t *out, size_t outlen, const dch_pack_params pp) {
  uint8_t *cur = out;
  size_t rem = outlen;

  if (outlen != py_params_nbytes(pp)) {
    return 1;
  }
  if (py_wr(&cur, &rem, PY_PROOF_MAGIC, sizeof(PY_PROOF_MAGIC)) ||
      py_wr_dch_params(&cur, &rem, pp->pp_dch) ||
      py_wr(&cur, &rem, pp->pp_pack, sizeof(*pp->pp_pack)) ||
      py_wr_u64(&cur, &rem, pp->pp_pack->np)) {
    return 2;
  }
  for (size_t i = 0; i < pp->pp_pack->np; i++) {
    if (py_wr_lab_params(&cur, &rem, pp->pp_pack->p[i])) {
      return 3;
    }
  }
  if (py_wr_u64(&cur, &rem, pp->pp_pack->zkp != NULL ? 1 : 0)) {
    return 4;
  }
  if (pp->pp_pack->zkp != NULL &&
      py_wr(&cur, &rem, pp->pp_pack->zkp[0], sizeof(*pp->pp_pack->zkp[0]))) {
    return 5;
  }
  return rem == 0 ? 0 : 6;
}

static int py_rd_alloc_blob(const uint8_t **cur, size_t *rem, void **out, size_t len) {
  *out = len ? _malloc(len) : NULL;
  if (len && *out == NULL) {
    return 1;
  }
  return py_rd_blob_exact(cur, rem, *out, len);
}

static int py_rd_dch_params(const uint8_t **cur, size_t *rem, dch_params pp) {
  if (py_rd(cur, rem, pp, sizeof(*pp))) {
    return 1;
  }
  if (py_rd_alloc_blob(cur, rem, (void **)&pp->n, pp->r * sizeof(size_t)) ||
      py_rd_alloc_blob(cur, rem, (void **)&pp->normsq, pp->r * sizeof(uint64_t)) ||
      py_rd_alloc_blob(cur, rem, (void **)&pp->normty, pp->r * sizeof(normtype)) ||
      py_rd_alloc_blob(cur, rem, (void **)&pp->exact_map, pp->nexact * sizeof(size_t)) ||
      py_rd_alloc_blob(cur, rem, (void **)&pp->incom_offw, pp->nincom * sizeof(size_t)) ||
      py_rd_alloc_blob(cur, rem, (void **)&pp->incom_lenw, pp->nincom * sizeof(size_t)) ||
      py_rd_alloc_blob(cur, rem, (void **)&pp->kappa_inner, pp->nincom * sizeof(size_t)) ||
      py_rd_alloc_blob(cur, rem, (void **)&pp->off_exact, pp->nexact * sizeof(size_t)) ||
      py_rd_alloc_blob(cur, rem, (void **)&pp->off_sigma, pp->nexact * sizeof(size_t)) ||
      py_rd_alloc_blob(cur, rem, (void **)&pp->off_incom, pp->nincom * sizeof(size_t)) ||
      py_rd_alloc_blob(cur, rem, (void **)&pp->len_exact_merge, pp->nexact * sizeof(size_t)) ||
      py_rd_alloc_blob(cur, rem, (void **)&pp->idx_exact, pp->nexact * sizeof(size_t)) ||
      py_rd_alloc_blob(cur, rem, (void **)&pp->idx_exact_merge, pp->nexact * sizeof(size_t)) ||
      py_rd_alloc_blob(cur, rem, (void **)&pp->idx_sigma, pp->nexact * sizeof(size_t)) ||
      py_rd_alloc_blob(cur, rem, (void **)&pp->nparts_exact_merge, pp->nexact * sizeof(size_t))) {
    return 2;
  }
  return 0;
}

static int py_rd_lab_params(const uint8_t **cur, size_t *rem, lab_params pp) {
  size_t *rr = NULL, *n = NULL;
  uint64_t *normsq = NULL;

  if (py_rd(cur, rem, pp, sizeof(*pp)) ||
      py_rd_alloc_blob(cur, rem, (void **)&rr, pp->r_old * sizeof(size_t)) ||
      py_rd_alloc_blob(cur, rem, (void **)&n, pp->r * sizeof(size_t)) ||
      py_rd_alloc_blob(cur, rem, (void **)&normsq, pp->r * sizeof(uint64_t))) {
    return 1;
  }
  pp->rr = _malloc((pp->r_old + pp->r) * sizeof(size_t) + pp->r * sizeof(uint64_t));
  if (pp->rr == NULL) {
    return 2;
  }
  pp->n = &pp->rr[pp->r_old];
  pp->normsq = (uint64_t *)&pp->n[pp->r];
  memcpy(pp->rr, rr, pp->r_old * sizeof(size_t));
  memcpy(pp->n, n, pp->r * sizeof(size_t));
  memcpy(pp->normsq, normsq, pp->r * sizeof(uint64_t));
  free(rr);
  free(n);
  free(normsq);
  return 0;
}

int py_params_deserialize(dch_pack_params pp, const uint8_t *in, size_t inlen) {
  const uint8_t *cur = in;
  size_t rem = inlen;
  uint8_t magic[sizeof(PY_PROOF_MAGIC)];
  uint64_t got;

  memset(pp, 0, sizeof(*pp));
  if (py_rd(&cur, &rem, magic, sizeof(magic)) ||
      memcmp(magic, PY_PROOF_MAGIC, sizeof(PY_PROOF_MAGIC)) != 0) {
    return 1;
  }
  if (py_rd_dch_params(&cur, &rem, pp->pp_dch) ||
      py_rd(&cur, &rem, pp->pp_pack, sizeof(*pp->pp_pack)) ||
      py_rd_u64(&cur, &rem, &got)) {
    return 2;
  }
  if (got != pp->pp_pack->np) {
    return 3;
  }
  pp->pp_pack->p = pp->pp_pack->np ? _malloc(pp->pp_pack->np * sizeof(lab_params)) : NULL;
  for (size_t i = 0; i < pp->pp_pack->np; i++) {
    if (py_rd_lab_params(&cur, &rem, pp->pp_pack->p[i])) {
      return 4;
    }
  }
  if (py_rd_u64(&cur, &rem, &got)) {
    return 5;
  }
  if (got) {
    pp->pp_pack->zkp = _malloc(sizeof(lnp_params));
    if (py_rd(&cur, &rem, pp->pp_pack->zkp[0], sizeof(*pp->pp_pack->zkp[0]))) {
      return 6;
    }
  } else {
    pp->pp_pack->zkp = NULL;
  }
  return rem == 0 ? 0 : 7;
}

int py_params_is_zk(const dch_pack_params pp) {
  return pp->pp_pack->zkp != NULL;
}

int py_proof_serialize(uint8_t *out, size_t outlen, const dch_pack_params pp,
                       const dch_pack_proof pi) {
  const unsigned int encoding = py_proof_encoding();
  const unsigned int fields = py_proof_field_mask();
  const char *profile = getenv("REVAC_LABRADOR_PROFILE");
  uint8_t *cur = out;
  size_t rem = outlen;
  uint8_t magic[sizeof(PY_PROOF_MAGIC)];

  if (outlen != py_proof_nbytes(pp, pi)) {
    return 1;
  }
  if (profile != NULL && strcmp(profile, "0") != 0) {
    size_t raw_owt = 0;
    fprintf(stderr, "[Labrador proof owt] r=%zu", pi->pi_pack->owt->r);
    for (size_t i = 0; i < pi->pi_pack->owt->r; i++) {
      const unsigned int bits = py_polyvec_signed_bits(
          pi->pi_pack->owt->s[i], pi->pi_pack->owt->n[i]);
      raw_owt += pi->pi_pack->owt->n[i] * sizeof(poly);
      fprintf(stderr, "; v%zu:n=%zu,bits=%u", i,
              pi->pi_pack->owt->n[i], bits);
    }
    fprintf(stderr, "; raw=%zu 字节\n", raw_owt);
    fprintf(stderr, "[Labrador proof layout] dch=%zu",
            pp->pp_dch->kappa_outer * POLZBYTES);
    for (size_t i = 0; i < pp->pp_pack->np; i++) {
      const size_t lab_len = py_pack_lab_proof_polz_len(pp->pp_pack, i);
      fprintf(stderr, "; lab%zu=%zu,p%zu=%zu(%u bits)", i,
              lab_len * POLZBYTES, i, 256 * sizeof(int32_t),
              py_int32vec_signed_bits(pi->pi_pack->p[i]->p, 256));
    }
    if (pp->pp_pack->zkp != NULL) {
      fprintf(stderr, "; lnp=%zu",
              py_lnp_field_blob_nbytes(pp->pp_pack->zkp[0]));
    }
    fprintf(stderr, "\n");
  }
  memcpy(magic, PY_PROOF_MAGIC, sizeof(magic));
  if (encoding == 2U) magic[sizeof(magic) - 1] = 0x17U;
  if (encoding == 3U) magic[sizeof(magic) - 1] = 0x27U;
  if (encoding >= 4U) magic[sizeof(magic) - 1] = 0x37U;
  if (py_wr(&cur, &rem, magic, sizeof(magic))) {
    return 2;
  }

  if (((fields & 1U) && py_wr_polz_blob(&cur, &rem, pi->pi_dch->com,
                                        pp->pp_dch->kappa_outer)) ||
      (!(fields & 1U) && py_wr_blob(
          &cur, &rem, pi->pi_dch->com,
          pp->pp_dch->kappa_outer * sizeof(polz)))) {
    return 3;
  }

  if (py_wr_u64(&cur, &rem, pp->pp_pack->np)) {
    return 4;
  }
  for (size_t i = 0; i < pp->pp_pack->np; i++) {
    size_t lab_len = py_pack_lab_proof_polz_len(pp->pp_pack, i);
    if (((fields & 2U) && py_wr_polz_blob(
            &cur, &rem, pi->pi_pack->p[i]->m[0], lab_len)) ||
        (!(fields & 2U) && py_wr_blob(
            &cur, &rem, pi->pi_pack->p[i]->m[0],
            lab_len * sizeof(polz)))) {
      return 5;
    }
    if ((encoding >= 4U && py_wr_signed_int32vec(
            &cur, &rem, pi->pi_pack->p[i]->p, 256)) ||
        (encoding < 4U && py_wr(
            &cur, &rem, pi->pi_pack->p[i]->p,
            256 * sizeof(int32_t)))) {
      return 6;
    }
  }

  if (py_wr_u64(&cur, &rem, pp->pp_pack->zkp != NULL ? 1 : 0)) {
    return 7;
  }
  if (pp->pp_pack->zkp != NULL) {
    if (((fields & 4U) && py_wr_lnp_field_blob(
             &cur, &rem, pp->pp_pack->zkp[0],
             pi->pi_pack->zkp[0])) ||
        (!(fields & 4U) && py_wr_blob(
            &cur, &rem, pi->pi_pack->zkp[0]->m[0],
            py_lnp_proof_polz_len(pp->pp_pack->zkp[0]) * sizeof(polz)))) {
      return 8;
    }
  }

  if (py_wr_u64(&cur, &rem, pi->pi_pack->owt->r)) {
    return 9;
  }
  for (size_t i = 0; i < pi->pi_pack->owt->r; i++) {
    if (py_wr_u64(&cur, &rem, pi->pi_pack->owt->n[i])) {
      return 10;
    }
  }
  for (size_t i = 0; i < pi->pi_pack->owt->r; i++) {
    if ((encoding >= 3U && py_wr_signed_polyvec_blob(
             &cur, &rem, pi->pi_pack->owt->s[i],
             pi->pi_pack->owt->n[i])) ||
        (encoding < 3U && py_wr_blob(
             &cur, &rem, pi->pi_pack->owt->s[i],
             pi->pi_pack->owt->n[i] * sizeof(poly)))) {
      return 10;
    }
  }
  return rem == 0 ? 0 : 11;
}

int py_proof_deserialize(dch_pack_proof pi, const dch_pack_params pp,
                         const uint8_t *in, size_t inlen) {
  const uint8_t *cur = in;
  size_t rem = inlen;
  uint8_t magic[sizeof(PY_PROOF_MAGIC)];
  uint64_t got;
  unsigned int encoding;
  unsigned int fields;

  memset(pi, 0, sizeof(*pi));

  if (py_rd(&cur, &rem, magic, sizeof(magic))) {
    return 1;
  }
  if (memcmp(magic, PY_PROOF_MAGIC, sizeof(PY_PROOF_MAGIC)) == 0) {
    encoding = 1U;
    fields = 0;
  } else if (memcmp(magic, PY_PROOF_MAGIC,
                    sizeof(PY_PROOF_MAGIC) - 1) == 0 &&
             magic[sizeof(magic) - 1] == 0x17U) {
    encoding = 2U;
    fields = 7U;
  } else if (memcmp(magic, PY_PROOF_MAGIC,
                    sizeof(PY_PROOF_MAGIC) - 1) == 0 &&
             magic[sizeof(magic) - 1] == 0x27U) {
    encoding = 3U;
    fields = 7U;
  } else if (memcmp(magic, PY_PROOF_MAGIC,
                    sizeof(PY_PROOF_MAGIC) - 1) == 0 &&
             magic[sizeof(magic) - 1] == 0x37U) {
    encoding = 4U;
    fields = 7U;
  } else {
    return 1;
  }
  proof_hash_canonical_set(fields != 0);

  size_t dch_len = pp->pp_dch->kappa_outer * sizeof(polz);
  pi->pi_dch->com = dch_len ? _aligned_alloc(64, dch_len) : NULL;
  if (dch_len && pi->pi_dch->com == NULL) {
    return 2;
  }
  if (((fields & 1U) && py_rd_polz_blob(
          &cur, &rem, pi->pi_dch->com, pp->pp_dch->kappa_outer)) ||
      (!(fields & 1U) &&
       py_rd_blob_exact(&cur, &rem, pi->pi_dch->com, dch_len))) {
    return 3;
  }

  if (py_rd_u64(&cur, &rem, &got) || got != pp->pp_pack->np) {
    return 4;
  }
  pi->pi_pack->np = pp->pp_pack->np;
  pi->pi_pack->p = pi->pi_pack->np ? _malloc(pi->pi_pack->np * sizeof(lab_proof)) : NULL;
  for (size_t i = 0; i < pi->pi_pack->np; i++) {
    size_t lab_len = py_pack_lab_proof_polz_len(pp->pp_pack, i);
    polz *arr = lab_len ? _aligned_alloc(64, lab_len * sizeof(polz)) : NULL;
    if (lab_len && arr == NULL) {
      return 5;
    }
    if (((fields & 2U) && py_rd_polz_blob(&cur, &rem, arr, lab_len)) ||
        (!(fields & 2U) &&
         py_rd_blob_exact(&cur, &rem, arr, lab_len * sizeof(polz)))) {
      return 6;
    }
    py_lab_proof_set_ptrs(pi->pi_pack->p[i], pp->pp_pack->p[i], arr);
    if ((encoding >= 4U && py_rd_signed_int32vec(
            &cur, &rem, pi->pi_pack->p[i]->p, 256)) ||
        (encoding < 4U && py_rd(
            &cur, &rem, pi->pi_pack->p[i]->p,
            256 * sizeof(int32_t)))) {
      return 7;
    }
  }

  if (py_rd_u64(&cur, &rem, &got)) {
    return 8;
  }
  if ((pp->pp_pack->zkp != NULL) != (got != 0)) {
    return 9;
  }
  if (pp->pp_pack->zkp != NULL) {
    size_t zkp_len = py_lnp_proof_polz_len(pp->pp_pack->zkp[0]);
    polz *arr = _aligned_alloc(64, zkp_len * sizeof(polz));
    if (arr == NULL) {
      return 10;
    }
    pi->pi_pack->zkp = _malloc(sizeof(lnp_proof));
    py_lnp_proof_set_ptrs(pi->pi_pack->zkp[0], pp->pp_pack->zkp[0], arr);
    if (((fields & 4U) && py_rd_lnp_field_blob(
            &cur, &rem, pp->pp_pack->zkp[0],
            pi->pi_pack->zkp[0])) ||
        (!(fields & 4U) &&
         py_rd_blob_exact(&cur, &rem, arr, zkp_len * sizeof(polz)))) {
      return 11;
    }
  } else {
    pi->pi_pack->zkp = NULL;
  }

  if (py_rd_u64(&cur, &rem, &got)) {
    return 12;
  }
  witness_init(pi->pi_pack->owt, got, got);
  size_t nn = 0;
  for (size_t i = 0; i < pi->pi_pack->owt->r; i++) {
    if (py_rd_u64(&cur, &rem, &got)) {
      return 13;
    }
    pi->pi_pack->owt->n[i] = got;
    nn += got;
  }
  pi->pi_pack->owt->s[0] = nn ? _aligned_alloc(64, nn * sizeof(poly)) : NULL;
  if (nn && pi->pi_pack->owt->s[0] == NULL) {
    return 14;
  }
  for (size_t i = 1; i < pi->pi_pack->owt->r; i++) {
    pi->pi_pack->owt->s[i] = &pi->pi_pack->owt->s[i - 1][pi->pi_pack->owt->n[i - 1]];
  }
  for (size_t i = 0; i < pi->pi_pack->owt->r; i++) {
    if ((encoding >= 3U && py_rd_signed_polyvec_blob(
             &cur, &rem, pi->pi_pack->owt->s[i],
             pi->pi_pack->owt->n[i])) ||
        (encoding < 3U && py_rd_blob_exact(
             &cur, &rem, pi->pi_pack->owt->s[i],
             pi->pi_pack->owt->n[i] * sizeof(poly)))) {
      return 15;
    }
  }
  return rem == 0 ? 0 : 16;
}

void py_free_witness(witness wt){
  witness_free(wt);
}

void py_free_statement(statement st){
  statement_free(st);
  comkey_free();
}

void py_free_params(dch_pack_params pp){
  dch_pack_params_free(pp);
}

void py_free_proof(dch_pack_proof pi){
  dch_pack_proof_free(pi);
}
