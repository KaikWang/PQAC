/*
 * Isolated ablation benchmark for the accumulator improvements.
 *
 * The private base-2 reference below is built with the current compiler,
 * polynomial arithmetic, OpenMP policy, streamed-tail layout and redundant
 * composite removal.  Therefore it changes only the gadget radix/dimension,
 * rather than comparing the current code with a historical checkout.
 */
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "arith.h"
#include "fips202.h"
#include "poly_q_sampling.h"
#include "random.h"
#include "randombytes.h"
#include "revac_context.h"
#include "revac_labrador_bridge.h"
#include "revac_sync.h"
#include "revac_ta.h"
#include "revac_uav.h"

#define REF_B2_DIM PARAM_Q_BITLEN
#define AB_REPETITIONS_DEFAULT 3u
#define AB_WARMUPS_DEFAULT 1u
#define AB_MAX_SAMPLES 1000u
/* Must match LINEAR_DECOMP_CHUNKS in revac_acc_memver_import_labrador.py. */
#define AB_LABRADOR_LINEAR_DECOMP_CHUNKS 4u
#define AB_ZK_REPETITIONS_DEFAULT 1u
#define AB_ZK_WARMUPS_DEFAULT 1u
#define AB_RANDOMIZED_REFRESH_BITS_DEFAULT 128u

typedef struct { poly_q_mat_d_d blocks[REF_B2_DIM]; } ref_b2_gmat_t;
typedef struct { ref_b2_gmat_t bit[ACC_GENPP_HANDLE_BITS]; } ref_b2_Ba_t;
typedef struct {
  poly_q_mat_d_d blocks[REF_B2_DIM][REF_B2_DIM];
} ref_b2_ginv_mat_t;
typedef struct {
  poly_q_vec_d left;
  poly_q_vec_d right[ACC_GENPP_HANDLE_BITS][REF_B2_DIM];
} ref_b2_delta_t;
typedef struct {
  double min, max, total, sum_squares;
  double samples[AB_MAX_SAMPLES];
  size_t count;
} timing_stats_t;

/*
 * Keep a data dependency on every measured result without XOR cancellation:
 * an even number of identical benchmark iterations used to print 0x0, which
 * looked like a failed check despite being expected.  This is not a
 * cryptographic digest and is never used by the scheme.
 */
static volatile uint64_t benchmark_digest = UINT64_C(0x6a09e667f3bcc909);

static void benchmark_mix(uint64_t value) {
  uint64_t mixed = benchmark_digest ^
      (value + UINT64_C(0x9e3779b97f4a7c15) +
       (benchmark_digest << 6) + (benchmark_digest >> 2));

  mixed = (mixed << 27) | (mixed >> (64 - 27));
  benchmark_digest = mixed * UINT64_C(0x94d049bb133111eb) +
                     UINT64_C(0x2545f4914f6cdd1d);
}

static double monotonic_seconds(void) {
  struct timespec now;
  (void)clock_gettime(CLOCK_MONOTONIC, &now);
  return (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
}

static int expect(int condition, const char *message) {
  if (!condition) { fprintf(stderr, "失败：%s\n", message); return 0; }
  return 1;
}

static size_t env_count(const char *name, size_t fallback) {
  const char *value = getenv(name);
  char *end = NULL;
  unsigned long parsed;
  if (value == NULL || value[0] == '\0') return fallback;
  parsed = strtoul(value, &end, 10);
  if (end == value || *end != '\0' || parsed == 0 || parsed > 1000) {
    fprintf(stderr, "忽略无效的 %s=%s，使用 %zu。\n", name, value, fallback);
    return fallback;
  }
  return (size_t)parsed;
}
static size_t env_nonnegative_count(const char *name, size_t fallback) {
  const char *value = getenv(name);
  char *end = NULL;
  unsigned long parsed;
  if (value == NULL || value[0] == '\0') return fallback;
  parsed = strtoul(value, &end, 10);
  if (end == value || *end != '\0' || parsed > 1000) {
    fprintf(stderr, "忽略无效的 %s=%s，使用 %zu。\n", name, value, fallback);
    return fallback;
  }
  return (size_t)parsed;
}

static int env_enabled(const char *name) {
  const char *value = getenv(name);
  return value != NULL && value[0] != '\0' && strcmp(value, "0") != 0 &&
         strcmp(value, "false") != 0 && strcmp(value, "FALSE") != 0;
}

static void stats_init(timing_stats_t *s) {
  s->min = s->max = s->total = s->sum_squares = 0.0;
  s->count = 0;
}
static void stats_add(timing_stats_t *s, double v) {
  if (s->count >= AB_MAX_SAMPLES) abort();
  if (s->count == 0 || v < s->min) s->min = v;
  if (s->count == 0 || v > s->max) s->max = v;
  s->total += v; s->sum_squares += v * v; s->samples[s->count++] = v;
}
static double stats_mean(const timing_stats_t *s) { return s->count == 0 ? 0.0 : s->total / (double)s->count; }
static double stats_stddev(const timing_stats_t *s) {
  double variance;
  if (s->count < 2) return 0.0;
  variance = (s->sum_squares - s->total * s->total / (double)s->count) /
             (double)(s->count - 1);
  return sqrt(variance > 0.0 ? variance : 0.0);
}
static void print_stats(const char *label, const timing_stats_t *s) {
  const double mean = stats_mean(s), stddev = stats_stddev(s);
  printf("  %-40s 平均 %8.3f ms，标准差 %7.3f ms，最小 %8.3f ms，最大 %8.3f ms，CV %.2f%%（%zu 次）\n",
         label, mean * 1000.0, stddev * 1000.0, s->min * 1000.0,
         s->max * 1000.0, mean == 0.0 ? 0.0 : stddev / mean * 100.0,
         s->count);
  printf("    原始样本（ms）：");
  for (size_t i = 0; i < s->count; i++) {
    printf("%s%.3f", i == 0 ? "" : ", ", s->samples[i] * 1000.0);
  }
  printf("\n");
}
static int handle_bit(uint32_t value, size_t bit) { return (int)((value >> bit) & 1u); }

static size_t reference_b2_parallel_threads(void) {
  size_t workers = 1;
#ifdef _OPENMP
  const char *configured = getenv("REVAC_ACC_THREADS");
  int runtime_max = omp_get_max_threads();
  if (runtime_max > 1) workers = (size_t)runtime_max;
  if (configured != NULL && configured[0] != '\0') {
    char *end = NULL;
    unsigned long requested = strtoul(configured, &end, 10);
    if (end != configured && *end == '\0' && requested > 0) {
      workers = (size_t)requested;
      if (runtime_max > 0 && workers > (size_t)runtime_max) workers = (size_t)runtime_max;
    }
  }
#endif
  return workers > REF_B2_DIM ? REF_B2_DIM : workers;
}

static void ref_b2_gmat_init(ref_b2_gmat_t *m) {
  for (size_t i = 0; i < REF_B2_DIM; i++) poly_q_mat_d_d_init(m->blocks[i]);
}
static void ref_b2_gmat_clear(ref_b2_gmat_t *m) {
  for (size_t i = 0; i < REF_B2_DIM; i++) poly_q_mat_d_d_clear(m->blocks[i]);
}
static void ref_b2_Ba_init(ref_b2_Ba_t *B) {
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) ref_b2_gmat_init(&B->bit[i]);
}
static void ref_b2_Ba_clear(ref_b2_Ba_t *B) {
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) ref_b2_gmat_clear(&B->bit[i]);
}
static void ref_b2_ginv_init(ref_b2_ginv_mat_t *m) {
  for (size_t i = 0; i < REF_B2_DIM; i++) for (size_t j = 0; j < REF_B2_DIM; j++) poly_q_mat_d_d_init(m->blocks[i][j]);
}
static void ref_b2_ginv_clear(ref_b2_ginv_mat_t *m) {
  for (size_t i = 0; i < REF_B2_DIM; i++) for (size_t j = 0; j < REF_B2_DIM; j++) poly_q_mat_d_d_clear(m->blocks[i][j]);
}
static void ref_b2_delta_init(ref_b2_delta_t *d) {
  poly_q_vec_d_init(d->left);
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) for (size_t j = 0; j < REF_B2_DIM; j++) poly_q_vec_d_init(d->right[i][j]);
}
static void ref_b2_delta_clear(ref_b2_delta_t *d) {
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) for (size_t j = 0; j < REF_B2_DIM; j++) poly_q_vec_d_clear(d->right[i][j]);
  poly_q_vec_d_clear(d->left);
}

static void sample_uniform_Ba_b2(ref_b2_Ba_t *B, const uint8_t seed[SEED_BYTES]) {
  uint8_t block_seed[SEED_BYTES]; keccak_state state;
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) for (size_t j = 0; j < REF_B2_DIM; j++) {
    shake256_init(&state); shake256_absorb(&state, seed, SEED_BYTES);
    shake256_absorb(&state, (const uint8_t *)&i, sizeof(i));
    shake256_absorb(&state, (const uint8_t *)&j, sizeof(j));
    shake256_finalize(&state); shake256_squeeze(block_seed, SEED_BYTES, &state);
    poly_q_mat_d_d_uniform(B->bit[i].blocks[j], block_seed, DOMAIN_SEPARATOR_DS, 0);
  }
}
/* Test-local copy of parameter expansion: identical seed derivation, but no
   production helper is made public just for this benchmark. */
static void sample_uniform_Ba_b3(acc_genpp_Ba_t *B, const uint8_t seed[SEED_BYTES]) {
  uint8_t block_seed[SEED_BYTES]; keccak_state state;
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
    shake256_init(&state); shake256_absorb(&state, seed, SEED_BYTES);
    shake256_absorb(&state, (const uint8_t *)&i, sizeof(i));
    shake256_absorb(&state, (const uint8_t *)&j, sizeof(j));
    shake256_finalize(&state); shake256_squeeze(block_seed, SEED_BYTES, &state);
    poly_q_mat_d_d_uniform(B->bit[i].blocks[j], block_seed, DOMAIN_SEPARATOR_DS, 0);
  }
}

static void ref_b2_gadget_block(poly_q_mat_d_d out, size_t digit) {
  poly_q_mat_d_d_zero(out);
  for (size_t row = 0; row < PARAM_D; row++) poly_q_set_coeff(out->rows[row]->entries[row], 0, (coeff_q)(UINT64_C(1) << digit));
}
static void ref_b2_Di(ref_b2_gmat_t *out, const ref_b2_Ba_t *B, size_t bit, int value) {
  poly_q_mat_d_d gadget; poly_q_mat_d_d_init(gadget);
  for (size_t j = 0; j < REF_B2_DIM; j++) {
    if (value) poly_q_mat_d_d_set(out->blocks[j], B->bit[bit].blocks[j]);
    else { ref_b2_gadget_block(gadget, j); poly_q_mat_d_d_sub(out->blocks[j], gadget, B->bit[bit].blocks[j]); }
  }
  poly_q_mat_d_d_clear(gadget);
}

/* Exact unsigned binary G^{-1}; q < 2^PARAM_Q_BITLEN.  The reconstruction
   guard makes the reference decomposition self-checking. */
static void ref_b2_ginv_vec(poly_q_vec_d out[REF_B2_DIM], const poly_q_vec_d in) {
  for (size_t d = 0; d < REF_B2_DIM; d++) poly_q_vec_d_zero(out[d]);
  for (size_t row = 0; row < PARAM_D; row++) for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
    uint64_t value = (uint64_t)poly_q_get_coeff(in->entries[row], coeff), recovered = 0;
    for (size_t d = 0; d < REF_B2_DIM; d++) { uint64_t bit = (value >> d) & 1u; recovered |= bit << d; if (bit) poly_q_set_coeff(out[d]->entries[row], coeff, 1); }
    if (recovered != value) abort();
  }
}
static void ref_b2_ginv_mat(ref_b2_ginv_mat_t *out, const ref_b2_gmat_t *in) {
  for (size_t r = 0; r < REF_B2_DIM; r++) for (size_t c = 0; c < REF_B2_DIM; c++) poly_q_mat_d_d_zero(out->blocks[r][c]);
  for (size_t block = 0; block < REF_B2_DIM; block++) for (size_t row = 0; row < PARAM_D; row++) for (size_t col = 0; col < PARAM_D; col++) for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
    uint64_t value = (uint64_t)poly_q_get_coeff(in->blocks[block]->rows[row]->entries[col], coeff), recovered = 0;
    for (size_t d = 0; d < REF_B2_DIM; d++) { uint64_t bit = (value >> d) & 1u; recovered |= bit << d; if (bit) poly_q_set_coeff(out->blocks[d][block]->rows[row]->entries[col], coeff, 1); }
    if (recovered != value) abort();
  }
}

static void ref_b2_mul_ginv_column(ref_b2_gmat_t *out, const ref_b2_gmat_t *lhs, const ref_b2_ginv_mat_t *rhs, size_t output) {
  poly_q acc, product; poly_q_init(acc); poly_q_init(product);
  for (size_t row = 0; row < PARAM_D; row++) for (size_t col = 0; col < PARAM_D; col++) {
    poly_q_zero(acc);
    for (size_t gadget = 0; gadget < REF_B2_DIM; gadget++) for (size_t inner = 0; inner < PARAM_D; inner++) {
      poly_q_mul(product, lhs->blocks[gadget]->rows[row]->entries[inner], rhs->blocks[gadget][output]->rows[inner]->entries[col]);
      poly_q_add(acc, acc, product);
    }
    poly_q_set(out->blocks[output]->rows[row]->entries[col], acc);
  }
  poly_q_clear(product); poly_q_clear(acc);
}
static void ref_b2_mul_ginv(ref_b2_gmat_t *out, const ref_b2_gmat_t *lhs, const ref_b2_ginv_mat_t *rhs) {
  const size_t workers = reference_b2_parallel_threads();
  if (workers <= 1) { for (size_t c = 0; c < REF_B2_DIM; c++) ref_b2_mul_ginv_column(out, lhs, rhs, c); return; }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(workers)
  for (int c = 0; c < REF_B2_DIM; c++) ref_b2_mul_ginv_column(out, lhs, rhs, (size_t)c);
#else
  for (size_t c = 0; c < REF_B2_DIM; c++) ref_b2_mul_ginv_column(out, lhs, rhs, c);
#endif
}

/* Same modern streamed EvalF (31 composites) as the production base-3 code. */
static void ref_b2_evalf(ref_b2_gmat_t *out, const ref_b2_Ba_t *B, uint32_t selector) {
  ref_b2_gmat_t current, D; ref_b2_ginv_mat_t inverse;
  ref_b2_gmat_init(&current); ref_b2_gmat_init(&D); ref_b2_ginv_init(&inverse);
  ref_b2_Di(&current, B, ACC_GENPP_HANDLE_BITS - 1, handle_bit(selector, ACC_GENPP_HANDLE_BITS - 1));
  for (size_t offset = 1; offset < ACC_GENPP_HANDLE_BITS; offset++) {
    size_t bit = ACC_GENPP_HANDLE_BITS - 1 - offset;
    ref_b2_ginv_mat(&inverse, &current); ref_b2_Di(&D, B, bit, handle_bit(selector, bit)); ref_b2_mul_ginv(&current, &D, &inverse);
  }
  for (size_t j = 0; j < REF_B2_DIM; j++) poly_q_mat_d_d_set(out->blocks[j], current.blocks[j]);
  ref_b2_ginv_clear(&inverse); ref_b2_gmat_clear(&D); ref_b2_gmat_clear(&current);
}

/* Point differs at its low bit.  As in production LOCAL EvalFX it performs
   the full 30-composite tail recurrence; only T1 is then consumed. */
static void ref_b2_evalfx(ref_b2_delta_t *delta, const ref_b2_Ba_t *B, uint32_t selector, const poly_q_vec_d u) {
  ref_b2_gmat_t current, D; ref_b2_ginv_mat_t inverse;
  poly_q_vec_d ginv_u[REF_B2_DIM], tmp_in[REF_B2_DIM], tmp_out[REF_B2_DIM];
  ref_b2_gmat_init(&current); ref_b2_gmat_init(&D); ref_b2_ginv_init(&inverse);
  for (size_t j = 0; j < REF_B2_DIM; j++) { poly_q_vec_d_init(ginv_u[j]); poly_q_vec_d_init(tmp_in[j]); poly_q_vec_d_init(tmp_out[j]); }
  ref_b2_Di(&current, B, ACC_GENPP_HANDLE_BITS - 1, handle_bit(selector, ACC_GENPP_HANDLE_BITS - 1));
  for (size_t bit = ACC_GENPP_HANDLE_BITS - 2;; bit--) {
    ref_b2_ginv_mat(&inverse, &current); ref_b2_Di(&D, B, bit, handle_bit(selector, bit)); ref_b2_mul_ginv(&current, &D, &inverse);
    if (bit == 1) break;
  }
  ref_b2_ginv_vec(ginv_u, u); poly_q_vec_d_zero(delta->left);
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) for (size_t j = 0; j < REF_B2_DIM; j++) poly_q_vec_d_zero(delta->right[i][j]);
  for (size_t j = 0; j < REF_B2_DIM; j++) poly_q_vec_d_set(tmp_in[j], ginv_u[j]);
  ref_b2_ginv_mat(&inverse, &current);
  for (size_t row = 0; row < REF_B2_DIM; row++) {
    poly_q_vec_d_zero(tmp_out[row]);
    for (size_t col = 0; col < REF_B2_DIM; col++) poly_q_mat_d_d_muladd_vec_d(tmp_out[row], inverse.blocks[row][col], tmp_in[col]);
    if (handle_bit(selector, 0)) poly_q_vec_d_set(delta->right[0][row], tmp_out[row]);
    else poly_q_vec_d_neg(delta->right[0][row], tmp_out[row]);
  }
  for (size_t j = 0; j < REF_B2_DIM; j++) { poly_q_vec_d_clear(tmp_out[j]); poly_q_vec_d_clear(tmp_in[j]); poly_q_vec_d_clear(ginv_u[j]); }
  ref_b2_ginv_clear(&inverse); ref_b2_gmat_clear(&D); ref_b2_gmat_clear(&current);
}

static void witness_set(acc_genpp_wit_t *dst, const acc_genpp_wit_t *src) {
  poly_q_vec_d_set(dst->left, src->left); poly_q_vec_d_set(dst->v12, src->v12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) poly_q_vec_d_set(dst->v2[i], src->v2[i]);
  poly_q_set(dst->tag, src->tag);
  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) poly_q_vec_d_set(dst->right[bit][j], src->right[bit][j]);
  memcpy(dst->refresh_r, src->refresh_r, sizeof(dst->refresh_r));
  for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) poly_q_vec_d_set(dst->refresh_z[j], src->refresh_z[j]);
}
static void uav_clone(revac_uav_t *dst, const revac_uav_t *src) {
  dst->id = src->id; dst->handle = src->handle; dst->last_seq = src->last_seq;
  dst->registered = src->registered; dst->issued = src->issued; witness_set(&dst->witness, &src->witness);
}

static size_t witness_full_bytes_for_dim(size_t gadget_dim) {
  const size_t vec_bits = (size_t)PARAM_D * PARAM_N * PARAM_Q_BITLEN;
  const size_t vec_count = 2 + ACC_GENPP_TRUNC_K_L +
                           ACC_GENPP_HANDLE_BITS * gadget_dim;
  return (vec_count * vec_bits + PARAM_N + 7) / 8;
}

static size_t accumulator_Ba_bytes_for_dim(size_t gadget_dim) {
  const size_t coefficients = (size_t)ACC_GENPP_HANDLE_BITS * gadget_dim *
                              PARAM_D * PARAM_D * PARAM_N;
  return (coefficients * PARAM_Q_BITLEN + 7) / 8;
}

static size_t randomized_matrix_block_bytes(void) {
  const size_t coefficients = (size_t)ACC_GENPP_GADGET_DIM * PARAM_D *
                              PARAM_D * PARAM_N;
  return (coefficients * PARAM_Q_BITLEN + 7) / 8;
}

static size_t randomized_witness_block_bytes(void) {
  const size_t coefficients = (size_t)ACC_GENPP_GADGET_DIM * PARAM_D *
                              PARAM_N;
  return (coefficients * PARAM_Q_BITLEN + 7) / 8;
}

static void randomized_matrix_expand(acc_genpp_gmat_t *out,
                                     const uint8_t master[SEED_BYTES],
                                     size_t matrix_index) {
  uint8_t block_seed[SEED_BYTES];
  keccak_state state;

  for (size_t block = 0; block < ACC_GENPP_GADGET_DIM; block++) {
    shake256_init(&state);
    shake256_absorb(&state, master, SEED_BYTES);
    shake256_absorb(&state, (const uint8_t *)&matrix_index,
                    sizeof(matrix_index));
    shake256_absorb(&state, (const uint8_t *)&block, sizeof(block));
    shake256_finalize(&state);
    shake256_squeeze(block_seed, sizeof(block_seed), &state);
    poly_q_mat_d_d_uniform(out->blocks[block], block_seed,
                           DOMAIN_SEPARATOR_DS, 0);
  }
}

static void randomized_matrix_set(acc_genpp_gmat_t *dst,
                                  const acc_genpp_gmat_t *src) {
  for (size_t block = 0; block < ACC_GENPP_GADGET_DIM; block++) {
    poly_q_mat_d_d_set(dst->blocks[block], src->blocks[block]);
  }
}

static void randomized_matrix_add(acc_genpp_gmat_t *dst,
                                  const acc_genpp_gmat_t *src) {
  for (size_t block = 0; block < ACC_GENPP_GADGET_DIM; block++) {
    poly_q_mat_d_d_add(dst->blocks[block], dst->blocks[block],
                       src->blocks[block]);
  }
}

static void randomized_bits(uint8_t *bits, size_t bit_count,
                            const uint8_t seed[SEED_BYTES], size_t run) {
  keccak_state state;
  uint8_t byte = 0;

  shake256_init(&state);
  shake256_absorb(&state, seed, SEED_BYTES);
  shake256_absorb(&state, (const uint8_t *)&run, sizeof(run));
  shake256_finalize(&state);
  for (size_t i = 0; i < bit_count; i++) {
    if ((i & 7u) == 0) shake256_squeeze(&byte, 1, &state);
    bits[i] = (byte >> (i & 7u)) & 1u;
  }
}

static void randomized_matrix_compose_cached(
    acc_genpp_gmat_t *out, const acc_genpp_gmat_t *bank,
    const uint8_t *bits, size_t bit_count) {
  randomized_matrix_set(out, &bank[0]);
  for (size_t i = 0; i < bit_count; i++) {
    if (bits[i] != 0) randomized_matrix_add(out, &bank[i + 1]);
  }
}

static void randomized_matrix_compose_seeded(
    acc_genpp_gmat_t *out, const uint8_t master[SEED_BYTES],
    const uint8_t *bits, size_t bit_count) {
  acc_genpp_gmat_t expanded;

  acc_genpp_gmat_init(&expanded);
  randomized_matrix_expand(out, master, 0);
  for (size_t i = 0; i < bit_count; i++) {
    if (bits[i] == 0) continue;
    randomized_matrix_expand(&expanded, master, i + 1);
    randomized_matrix_add(out, &expanded);
  }
  acc_genpp_gmat_clear(&expanded);
}

static void randomized_small_vec(poly_q_vec_d out,
                                 const uint8_t seed[SEED_BYTES],
                                 size_t block) {
  keccak_state state;
  uint8_t buf[SHAKE256_RATE];
  size_t pos = sizeof(buf);

  poly_q_vec_d_zero(out);
  shake256_init(&state);
  shake256_absorb(&state, seed, SEED_BYTES);
  shake256_absorb(&state, (const uint8_t *)&block, sizeof(block));
  shake256_finalize(&state);
  for (size_t row = 0; row < PARAM_D; row++) {
    for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
      uint8_t sample;
      do {
        if (pos == sizeof(buf)) {
          shake256_squeeze(buf, sizeof(buf), &state);
          pos = 0;
        }
        sample = buf[pos++] & 3u;
      } while (sample == 3u);
      poly_q_set_coeff(out->entries[row], coeff, (coeff_q)sample - 1);
    }
  }
}

static void randomized_matrix_mul_vec(
    poly_q_vec_d out, const acc_genpp_gmat_t *matrix,
    const poly_q_vec_d rhs[ACC_GENPP_GADGET_DIM]) {
  poly_q_vec_d_zero(out);
  for (size_t block = 0; block < ACC_GENPP_GADGET_DIM; block++) {
    poly_q_mat_d_d_muladd_vec_d(out, matrix->blocks[block], rhs[block]);
  }
}

static int randomized_vec_equal(const poly_q_vec_d lhs,
                                const poly_q_vec_d rhs) {
  for (size_t row = 0; row < PARAM_D; row++) {
    for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
      if (poly_q_get_coeff(lhs->entries[row], coeff) !=
          poly_q_get_coeff(rhs->entries[row], coeff)) return 0;
    }
  }
  return 1;
}

/* Exact counts for the current grouped-right Labrador layout.  These are
   statement/witness-layout counts, not a claim that packed proof bytes scale
   exactly linearly with them. */
static size_t labrador_vector_count_for_dim(size_t gadget_dim) {
  const size_t same_handle_opening =
      (size_t)2 * PARAM_D +
      (size_t)PARAM_D * AB_LABRADOR_LINEAR_DECOMP_CHUNKS;
  return ACC_GENPP_HANDLE_BITS + 1 + 1 +
         (size_t)PARAM_D * AB_LABRADOR_LINEAR_DECOMP_CHUNKS * 3 +
         (size_t)ACC_GENPP_TRUNC_K_L * PARAM_D +
         gadget_dim * PARAM_D + same_handle_opening;
}

static size_t labrador_poly_count_for_dim(size_t gadget_dim) {
  const size_t same_handle_opening =
      (size_t)2 * PARAM_D +
      (size_t)PARAM_D * AB_LABRADOR_LINEAR_DECOMP_CHUNKS;
  return ACC_GENPP_HANDLE_BITS + ACC_GENPP_HANDLE_BITS + 1 +
         (size_t)PARAM_D * AB_LABRADOR_LINEAR_DECOMP_CHUNKS * 3 +
         (size_t)ACC_GENPP_TRUNC_K_L * PARAM_D +
         (size_t)ACC_GENPP_HANDLE_BITS * gadget_dim * PARAM_D +
         same_handle_opening;
}

static void print_kib_comparison(const char *label, size_t base2,
                                 size_t base3) {
  const double saved = 100.0 * (1.0 - (double)base3 / (double)base2);
  printf("  %-30s base-2 %9.2f KiB；base-3 %9.2f KiB；减少 %6.2f%%（%.2fx）\n",
         label, base2 / 1024.0, base3 / 1024.0, saved,
         (double)base2 / (double)base3);
}

static void print_structural_ablation(void) {
  const size_t b2_full = witness_full_bytes_for_dim(REF_B2_DIM);
  const size_t b3_full = acc_genpp_wit_full_size_bytes();
  const size_t compact = acc_genpp_wit_compact_size_bytes();
  const size_t b2_ba = accumulator_Ba_bytes_for_dim(REF_B2_DIM);
  const size_t b3_ba = accumulator_Ba_bytes_for_dim(ACC_GENPP_GADGET_DIM);
  const size_t b2_vectors = labrador_vector_count_for_dim(REF_B2_DIM);
  const size_t b3_vectors = labrador_vector_count_for_dim(ACC_GENPP_GADGET_DIM);
  const size_t b2_polys = labrador_poly_count_for_dim(REF_B2_DIM);
  const size_t b3_polys = labrador_poly_count_for_dim(ACC_GENPP_GADGET_DIM);

  printf("\n=== 消融 0：平衡 base-3 的结构性收益 ===\n");
  printf("完整见证、Ba 公参按 q 系数的无压缩紧凑表示计算；Labrador 项为当前 grouped-right 布局的精确关系规模。\n");
  print_kib_comparison("累加器完整见证", b2_full, b3_full);
  printf("  %-30s %9.2f KiB（与 gadget 维度无关）\n",
         "累加器紧凑见证（不含 right）", compact / 1024.0);
  print_kib_comparison("累加器 Ba 公共矩阵", b2_ba, b3_ba);
  printf("  Labrador 见证向量：base-2 %zu；base-3 %zu；减少 %.2f%%。\n",
         b2_vectors, b3_vectors,
         100.0 * (1.0 - (double)b3_vectors / (double)b2_vectors));
  printf("  Labrador 见证多项式：base-2 %zu；base-3 %zu；减少 %.2f%%。\n",
         b2_polys, b3_polys,
         100.0 * (1.0 - (double)b3_polys / (double)b2_polys));
  printf("  注：上述两种布局都包含同一句柄承诺的 24 个固定见证向量/多项式；该固定开销不归因于 base-3。\n");
  printf("  注：Labrador 打包 proof 的实际字节数取决于参数和随机性；下面的在线测试会报告 base-3 实测值。\n");
}

static int run_radix_ablation(size_t warmups, size_t repetitions) {
  uint8_t seed[SEED_BYTES]; ref_b2_Ba_t b2; acc_genpp_Ba_t b3;
  ref_b2_gmat_t b2_f_out; acc_genpp_gmat_t b3_f_out;
  ref_b2_delta_t b2_fx_out; acc_genpp_wit_t b3_fx_out; poly_q_vec_d u;
  timing_stats_t b2_f, b3_f, b2_fx, b3_fx;
  const uint32_t selector = UINT32_C(0x7a35c69d), point = selector ^ UINT32_C(1);
  memset(seed, 0, sizeof(seed)); for (size_t i = 0; i < sizeof(seed); i++) seed[i] = (uint8_t)(0xa5u ^ (uint8_t)(17u * i));
  ref_b2_Ba_init(&b2); acc_genpp_Ba_init(&b3); ref_b2_gmat_init(&b2_f_out); acc_genpp_gmat_init(&b3_f_out);
  ref_b2_delta_init(&b2_fx_out); acc_genpp_wit_init(&b3_fx_out); poly_q_vec_d_init(u);
  sample_uniform_Ba_b2(&b2, seed); sample_uniform_Ba_b3(&b3, seed); poly_q_vec_d_uniform(u, seed, DOMAIN_SEPARATOR_U);
  stats_init(&b2_f); stats_init(&b3_f); stats_init(&b2_fx); stats_init(&b3_fx);
  printf("\n=== 消融 A：平衡 base-3 对 EvalF / EvalFX 的影响 ===\n");
  printf("对照条件：同一可执行文件、同一 -march=native/OpenMP、固定公参种子；两者均使用当前流式 tail/去冗余复合。\n");
  printf("base-2 维度：%u；平衡 base-3 维度：%u；工作线程：base-2=%zu，base-3=%zu。\n", (unsigned)REF_B2_DIM, (unsigned)ACC_GENPP_GADGET_DIM, reference_b2_parallel_threads(), acc_genpp_parallel_threads());
  for (size_t run = 0; run < warmups + repetitions; run++) {
    double start, elapsed;
    start = monotonic_seconds(); ref_b2_evalf(&b2_f_out, &b2, selector); elapsed = monotonic_seconds() - start;
    benchmark_mix((uint64_t)poly_q_get_coeff(b2_f_out.blocks[0]->rows[0]->entries[0], 0)); if (run >= warmups) stats_add(&b2_f, elapsed);
    start = monotonic_seconds(); acc_genpp_evalf(&b3_f_out, &b3, selector); elapsed = monotonic_seconds() - start;
    benchmark_mix((uint64_t)poly_q_get_coeff(b3_f_out.blocks[0]->rows[0]->entries[0], 0)); if (run >= warmups) stats_add(&b3_f, elapsed);
    start = monotonic_seconds(); ref_b2_evalfx(&b2_fx_out, &b2, selector, u); elapsed = monotonic_seconds() - start;
    benchmark_mix((uint64_t)poly_q_get_coeff(b2_fx_out.right[0][0]->entries[0], 0)); if (run >= warmups) stats_add(&b2_fx, elapsed);
    start = monotonic_seconds(); acc_genpp_evalfx_apply(&b3_fx_out, &b3, selector, point, u); elapsed = monotonic_seconds() - start;
    benchmark_mix((uint64_t)poly_q_get_coeff(b3_fx_out.right[0][0]->entries[0], 0)); if (run >= warmups) stats_add(&b3_fx, elapsed);
  }
  print_stats("base-2 EvalF（当前优化参考）", &b2_f); print_stats("平衡 base-3 EvalF（生产实现）", &b3_f);
  print_stats("base-2 EvalFX（当前优化参考）", &b2_fx); print_stats("平衡 base-3 EvalFX（生产实现）", &b3_fx);
  printf("  EvalF 加速比（base-2 / base-3）：%.2fx\n", stats_mean(&b2_f) / stats_mean(&b3_f));
  printf("  EvalFX 加速比（base-2 / base-3）：%.2fx\n", stats_mean(&b2_fx) / stats_mean(&b3_fx));
  poly_q_vec_d_clear(u); acc_genpp_wit_clear(&b3_fx_out); ref_b2_delta_clear(&b2_fx_out);
  acc_genpp_gmat_clear(&b3_f_out); ref_b2_gmat_clear(&b2_f_out); acc_genpp_Ba_clear(&b3); ref_b2_Ba_clear(&b2);
  return 1;
}

/*
 * Test-only algebraic prototype of Appendix-E randomized refresh.
 *
 * It constructs A_r exactly and samples an extra short m' block z, then uses
 * the existing trapdoor sampler for the adjusted target c-A_r z.  Hence the
 * resulting pair satisfies [M_x | A_r](s,z)=c and is suitable for measuring
 * the structural and arithmetic cost.  It deliberately does not claim that
 * this independently sampled distribution is the production SampleLeft
 * distribution for the enlarged matrix, nor does it alter the production
 * MemVer/Labrador relation.
 */
static int run_randomized_refresh_ablation(
    const revac_ta_t *ta, const revac_uav_t *member,
    size_t random_bits, size_t warmups, size_t repetitions) {
  uint8_t master[SEED_BYTES], refresh_seed[SEED_BYTES];
  uint8_t *bits = NULL;
  acc_genpp_gmat_t *bank = NULL;
  acc_genpp_gmat_t cached_ar, seeded_ar;
  acc_genpp_wit_t cached_witness, seeded_witness;
  poly_q_vec_d cached_z[ACC_GENPP_GADGET_DIM];
  poly_q_vec_d seeded_z[ACC_GENPP_GADGET_DIM];
  poly_q_vec_d contribution, target, recombined;
  timing_stats_t cached_stats, seeded_stats;
  double setup_started, setup_elapsed;
  uint64_t cached_norm2 = 0, seeded_norm2 = 0;
  const size_t matrix_bytes = randomized_matrix_block_bytes();
  const size_t extra_witness_bytes = randomized_witness_block_bytes();
  const size_t randomizer_bytes = (random_bits + 7) / 8;
  const size_t current_witness_bytes = acc_genpp_wit_full_size_bytes();
  const size_t current_wire_bytes = REVAC_SYNC_WIRE_HEADER_BYTES +
                                    REVAC_CONTEXT_SIGNATURE_BYTES +
                                    REVAC_SYNC_WITNESS_BYTES;
  const size_t current_polys =
      labrador_poly_count_for_dim(ACC_GENPP_GADGET_DIM);
  const size_t structural_polys = current_polys +
      (size_t)ACC_GENPP_GADGET_DIM * PARAM_D +
      (random_bits + PARAM_N - 1) / PARAM_N;
  int ok = 0;

  printf("\n=== 消融 D：随机化矩阵异步 FreshSync（测试内原型）===\n");
  printf("随机量长度：%zu bit；刷新窗口示例：每个见证最多 400 次更新，FreshSync 后窗口计数归零。\n",
         random_bits);
  printf("说明：400 是策略阈值，不是本测试给出的安全定理；这里测量 A_r、扩展见证等式及采样成本。\n");

  bits = calloc(random_bits, 1);
  bank = calloc(random_bits + 1, sizeof(*bank));
  if (bits == NULL || bank == NULL) goto cleanup;
  memset(master, 0, sizeof(master));
  for (size_t i = 0; i < sizeof(master); i++) {
    master[i] = (uint8_t)(0x6du ^ (uint8_t)(i * 29u));
  }

  setup_started = monotonic_seconds();
  for (size_t i = 0; i <= random_bits; i++) {
    acc_genpp_gmat_init(&bank[i]);
    randomized_matrix_expand(&bank[i], master, i);
  }
  setup_elapsed = monotonic_seconds() - setup_started;

  acc_genpp_gmat_init(&cached_ar);
  acc_genpp_gmat_init(&seeded_ar);
  acc_genpp_wit_init(&cached_witness);
  acc_genpp_wit_init(&seeded_witness);
  poly_q_vec_d_init(contribution);
  poly_q_vec_d_init(target);
  poly_q_vec_d_init(recombined);
  for (size_t block = 0; block < ACC_GENPP_GADGET_DIM; block++) {
    poly_q_vec_d_init(cached_z[block]);
    poly_q_vec_d_init(seeded_z[block]);
  }
  stats_init(&cached_stats);
  stats_init(&seeded_stats);

  for (size_t run = 0; run < warmups + repetitions; run++) {
    double started, elapsed;

    shake256(refresh_seed, sizeof(refresh_seed),
             (const uint8_t *)&run, sizeof(run));
    randomized_bits(bits, random_bits, master, run);

    started = monotonic_seconds();
    randomized_matrix_compose_cached(&cached_ar, bank, bits, random_bits);
    for (size_t block = 0; block < ACC_GENPP_GADGET_DIM; block++) {
      randomized_small_vec(cached_z[block], refresh_seed, block);
    }
    randomized_matrix_mul_vec(contribution, &cached_ar, cached_z);
    poly_q_vec_d_sub(target, ta->acc, contribution);
    acc_genpp_sampleleft(&cached_witness, &ta->acc_pp, &ta->acc_td,
                         target, member->handle, refresh_seed);
    elapsed = monotonic_seconds() - started;
    if (!expect(acc_genpp_memver(&ta->acc_pp, target, member->handle,
                                 &cached_witness),
                "缓存矩阵 FreshSync 的基础见证关系") ) goto initialized_cleanup;
    poly_q_vec_d_add(recombined, target, contribution);
    if (!expect(randomized_vec_equal(recombined, ta->acc),
                "缓存矩阵 FreshSync 的扩展关系")) goto initialized_cleanup;
    cached_norm2 = acc_genpp_wit_norm2(&cached_witness);
    for (size_t block = 0; block < ACC_GENPP_GADGET_DIM; block++) {
      cached_norm2 += poly_q_vec_d_norm2(cached_z[block]);
    }
    benchmark_mix(cached_norm2);
    if (run >= warmups) stats_add(&cached_stats, elapsed);

    started = monotonic_seconds();
    randomized_matrix_compose_seeded(&seeded_ar, master, bits, random_bits);
    for (size_t block = 0; block < ACC_GENPP_GADGET_DIM; block++) {
      randomized_small_vec(seeded_z[block], refresh_seed, block);
    }
    randomized_matrix_mul_vec(contribution, &seeded_ar, seeded_z);
    poly_q_vec_d_sub(target, ta->acc, contribution);
    acc_genpp_sampleleft(&seeded_witness, &ta->acc_pp, &ta->acc_td,
                         target, member->handle, refresh_seed);
    elapsed = monotonic_seconds() - started;
    if (!expect(acc_genpp_memver(&ta->acc_pp, target, member->handle,
                                 &seeded_witness),
                "种子展开 FreshSync 的基础见证关系")) goto initialized_cleanup;
    poly_q_vec_d_add(recombined, target, contribution);
    if (!expect(randomized_vec_equal(recombined, ta->acc),
                "种子展开 FreshSync 的扩展关系")) goto initialized_cleanup;
    seeded_norm2 = acc_genpp_wit_norm2(&seeded_witness);
    for (size_t block = 0; block < ACC_GENPP_GADGET_DIM; block++) {
      seeded_norm2 += poly_q_vec_d_norm2(seeded_z[block]);
    }
    benchmark_mix(seeded_norm2);
    if (run >= warmups) stats_add(&seeded_stats, elapsed);
  }

  printf("  一次性显式矩阵库：%zu 个 A_i，共 %.2f MiB；生成 %.3f 秒。\n",
         random_bits + 1,
         (double)(random_bits + 1) * matrix_bytes / (1024.0 * 1024.0),
         setup_elapsed);
  printf("  种子展开公参：%u 字节主种子（不含安全证明所需的域分离元数据）。\n",
         (unsigned)SEED_BYTES);
  printf("  当前完整累加器见证：%.2f KiB。\n",
         current_witness_bytes / 1024.0);
  printf("  随机化见证结构投影：%.2f KiB = 当前见证 + %.2f KiB 的 m' 块 + %zu 字节 r；增加 %.2f%%。\n",
         (current_witness_bytes + extra_witness_bytes + randomizer_bytes) /
             1024.0,
         extra_witness_bytes / 1024.0, randomizer_bytes,
         100.0 * (double)(extra_witness_bytes + randomizer_bytes) /
             (double)current_witness_bytes);
  printf("  当前 MemWitSync wire：%.2f KiB；随机化未压缩投影：%.2f KiB。\n",
         current_wire_bytes / 1024.0,
         (current_wire_bytes + extra_witness_bytes + randomizer_bytes) /
             1024.0);
  printf("  Labrador witness 多项式结构下界：%zu -> %zu（增加 %.2f%%）；尚未计入 A_r 构造的二次约束。\n",
         current_polys, structural_polys,
         100.0 * ((double)structural_polys / (double)current_polys - 1.0));
  print_stats("FreshSync：显式 A_i 缓存 + 扩展关系采样", &cached_stats);
  print_stats("FreshSync：32-byte 种子按需展开 + 扩展关系采样",
              &seeded_stats);
  printf("  fresh 扩展见证 norm^2（最后一轮）：缓存=%" PRIu64
         "，种子展开=%" PRIu64 "；两条路径代数关系均通过。\n",
         cached_norm2, seeded_norm2);
  printf("  限制：这是正确等式的测试内采样原型；未修改生产 MemVer，也未声称等同于附录 E 的联合 SampleLeft 分布。\n");
  printf("  限制：随机化 Labrador proof 的大小/prove/verify 必须在扩展 ZK 关系完成后实测，不能由上述多项式比例直接推断。\n");
  ok = 1;

initialized_cleanup:
  for (size_t block = 0; block < ACC_GENPP_GADGET_DIM; block++) {
    poly_q_vec_d_clear(seeded_z[block]);
    poly_q_vec_d_clear(cached_z[block]);
  }
  poly_q_vec_d_clear(recombined);
  poly_q_vec_d_clear(target);
  poly_q_vec_d_clear(contribution);
  acc_genpp_wit_clear(&seeded_witness);
  acc_genpp_wit_clear(&cached_witness);
  acc_genpp_gmat_clear(&seeded_ar);
  acc_genpp_gmat_clear(&cached_ar);
cleanup:
  if (bank != NULL) {
    for (size_t i = 0; i <= random_bits; i++) {
      /* Only initialized after allocation in the normal path. */
      if (bits != NULL) acc_genpp_gmat_clear(&bank[i]);
    }
  }
  free(bank);
  free(bits);
  return ok;
}

/* This is intentionally the accumulator proof only: no anonymous-signature
   proving/verifying is included in either timer.  The bridge invokes the
   production Labrador script in zk=True/prove mode and its normal verifier. */
static int run_labrador_zk_benchmark(const revac_public_state_t *public_state,
                                     revac_uav_t *uav,
                                     size_t warmups,
                                     size_t repetitions) {
  static const uint8_t compact_magic[] =
      "REVAC-LABRADOR38-MEMVER\x03";
  revac_show_context_t ctx;
  timing_stats_t prove_stats, verify_stats;
  size_t proof_min = 0, proof_max = 0;
  const double layout_ratio =
      (double)labrador_poly_count_for_dim(REF_B2_DIM) /
      (double)labrador_poly_count_for_dim(ACC_GENPP_GADGET_DIM);
  int compact_profile_tamper_checked = 0;
  int ok = 1;

  printf("\n=== 消融 E：累加器 Labrador 零知识证明（仅服务器在线执行）===\n");
  if (!revac_labrador_online_enabled()) {
    printf("已跳过：REVAC_LABRADOR_ONLINE 未启用。服务器设置为 1 后，本测试会运行 zk=True 的真实 prove/verify。\n");
    return 1;
  }
  revac_show_context_init(&ctx);
  revac_show_context_from_public(&ctx, public_state, uav->pk.seed);
  stats_init(&prove_stats);
  stats_init(&verify_stats);
  for (size_t run = 0; run < warmups + repetitions; run++) {
    revac_show_proof_t proof;
    poly_q_vec_d link_randomness[2];
    double started, elapsed;

    revac_show_proof_init(&proof);
    poly_q_vec_d_init(link_randomness[0]);
    poly_q_vec_d_init(link_randomness[1]);
    revac_same_handle_commitment_generate(proof.handle_commitment,
                                          link_randomness, public_state,
                                          uav->handle);
    started = monotonic_seconds();
    if (!revac_labrador_memver_prove(&proof, public_state, uav,
                                     link_randomness)) {
      poly_q_vec_d_clear(link_randomness[1]);
      poly_q_vec_d_clear(link_randomness[0]);
      revac_show_proof_clear(&proof);
      ok = 0;
      break;
    }
    elapsed = monotonic_seconds() - started;
    if (!expect(proof.acc_zk_kind == REVAC_ACC_ZK_LABRADOR_MEMVER &&
                    proof.acc_zk_proof_len > 0,
                "Labrador prover 必须产生 zk=True 的累加器证明")) {
      poly_q_vec_d_clear(link_randomness[1]);
      poly_q_vec_d_clear(link_randomness[0]);
      revac_show_proof_clear(&proof);
      ok = 0;
      break;
    }
    if (run == 0 && warmups != 0) {
      printf("  冷启动（worker 启动/公参加载/参数生成，不计入稳定统计）：prove %.3f ms", elapsed * 1000.0);
    }
    if (run >= warmups) {
      stats_add(&prove_stats, elapsed);
      if (prove_stats.count == 1 || proof.acc_zk_proof_len < proof_min) proof_min = proof.acc_zk_proof_len;
      if (prove_stats.count == 1 || proof.acc_zk_proof_len > proof_max) proof_max = proof.acc_zk_proof_len;
    }

    started = monotonic_seconds();
    if (!revac_labrador_memver_verify(public_state, &ctx, &proof)) {
      poly_q_vec_d_clear(link_randomness[1]);
      poly_q_vec_d_clear(link_randomness[0]);
      revac_show_proof_clear(&proof);
      ok = 0;
      break;
    }
    elapsed = monotonic_seconds() - started;
    if (run == 0 && warmups != 0) {
      printf("；verify %.3f ms。\n", elapsed * 1000.0);
    }
    if (run >= warmups) stats_add(&verify_stats, elapsed);

    /* Compact worker proofs bind a canonical public-parameter profile instead
       of retransmitting 18.5 KiB of parameters.  Exercise the GCS rejection
       path once, outside the reported prove/verify samples. */
    if (!compact_profile_tamper_checked &&
        proof.acc_zk_proof_len > sizeof(compact_magic) - 1 + 32 &&
        memcmp(proof.acc_zk_proof, compact_magic,
               sizeof(compact_magic) - 1) == 0) {
      const size_t profile_offset = sizeof(compact_magic) - 1;
      const uint8_t saved = proof.acc_zk_proof[profile_offset];

      proof.acc_zk_proof[profile_offset] ^= UINT8_C(1);
      if (revac_labrador_memver_verify(public_state, &ctx, &proof)) {
        fprintf(stderr, "失败：GCS 接受了被篡改的 compact profile ID。\n");
        ok = 0;
      }
      proof.acc_zk_proof[profile_offset] = saved;
      if (!ok) {
        poly_q_vec_d_clear(link_randomness[1]);
        poly_q_vec_d_clear(link_randomness[0]);
        revac_show_proof_clear(&proof);
        break;
      }
      compact_profile_tamper_checked = 1;
      printf("  compact profile ID 篡改拒绝检查：通过。\n");
    }
    poly_q_vec_d_clear(link_randomness[1]);
    poly_q_vec_d_clear(link_randomness[0]);
    revac_show_proof_clear(&proof);
  }
  if (ok) {
    print_stats("累加器 ZK prove：预热后 C bridge + zk=True Labrador", &prove_stats);
    print_stats("累加器 ZK verify：预热后 C bridge + Labrador", &verify_stats);
    printf("  base-3 实测累加器 ZK payload：%.2f–%.2f KiB。\n",
           proof_min / 1024.0, proof_max / 1024.0);
    printf("  base-2 线性布局投影（非实测）：proof payload 约 %.2f–%.2f KiB；prove 约 %.3f 秒；verify 约 %.3f 秒。\n",
           proof_min * layout_ratio / 1024.0, proof_max * layout_ratio / 1024.0,
           stats_mean(&prove_stats) * layout_ratio,
           stats_mean(&verify_stats) * layout_ratio);
    printf("  该投影仅按 Labrador 见证多项式数量 %zu/%zu=%.3fx 缩放；最终论文应以可运行的 base-2 证明器实测替代。\n",
           labrador_poly_count_for_dim(REF_B2_DIM),
           labrador_poly_count_for_dim(ACC_GENPP_GADGET_DIM), layout_ratio);
  }
  revac_show_context_clear(&ctx);
  return ok;
}

static int run_sync_and_update_ablation(size_t warmups, size_t repetitions,
                                        size_t randomized_refresh_bits,
                                        size_t zk_warmups,
                                        size_t zk_repetitions,
                                        int zk_only) {
  revac_ta_t ta; revac_public_state_t public_state; revac_uav_t revoked, member, local_uav, delta_uav;
  acc_genpp_wit_t cached_sync, fallback_sync; acc_genpp_update_t update;
  uint8_t msg0[PARAM_M * PARAM_N / 8], msg1[PARAM_M * PARAM_N / 8], *wire = NULL;
  size_t wire_len = 0, cached_composites = 0, fallback_composites = 0, local_composites = 0, delta_composites = 0, cached_updates = 0, fallback_updates = 0;
  uint64_t member_id = UINT64_C(0xab10000000000002);
  int have_delta4_member = 0;
  timing_stats_t cached_stats, fallback_stats, local_stats, delta_stats; int ok = 0;
  revac_ta_init(&ta); revac_uav_init(&revoked); revac_uav_init(&member); revac_uav_init(&local_uav); revac_uav_init(&delta_uav);
  acc_genpp_wit_init(&cached_sync); acc_genpp_wit_init(&fallback_sync);
  stats_init(&cached_stats); stats_init(&fallback_stats); stats_init(&local_stats); stats_init(&delta_stats);
  printf("\n=== 消融 B：base-3 下 TA 复用 tail/DELTA6 的 MemWitSync ===\n");
  revac_ta_setup_q19(&ta); revac_ta_public_state(&public_state, &ta);
  revac_uav_keygen(&revoked, UINT64_C(0xab10000000000001), ta.ta_pk.seed);
  /* DELTA4 legitimately falls back for a handle sharing its low four bits
     with the revoked handle.  Select a deterministic candidate that does
     not share that prefix so this is an acceleration, not fallback, test. */
  for (size_t attempt = 0; attempt < 128; attempt++, member_id++) {
    revac_uav_keygen(&member, member_id, ta.ta_pk.seed);
    if (((revac_handle_from_identity(&member.pk, member.id) ^
          revac_handle_from_identity(&revoked.pk, revoked.id)) &
         UINT32_C(0x0f)) != 0) {
      have_delta4_member = 1;
      break;
    }
  }
  randombytes(msg0, sizeof(msg0)); randombytes(msg1, sizeof(msg1));
  if (!expect(have_delta4_member, "找到可用的 DELTA4 成员句柄") ||
      !expect(revac_ta_register_uav(&ta, &revoked), "撤销对象注册") || !expect(revac_ta_register_uav(&ta, &member), "成员注册") || !expect(revac_ta_issue_uav(&ta, &revoked, msg0), "撤销对象签发") || !expect(revac_ta_issue_uav(&ta, &member, msg1), "成员签发") || !expect(revac_ta_revoke_uav_context(&ta, revoked.handle, &update), "生成 DELTA4/DELTA6 撤销更新")) goto cleanup;
  revac_ta_public_state(&public_state, &ta);
  if (!expect(ta.sync_delta_cache[update.seq - 1].used, "TA 私有 DELTA6 缓存存在") || !expect(revac_ta_export_context(&ta, update.seq, REVAC_CONTEXT_PROFILE_DELTA4, &wire, &wire_len), "导出 BML DELTA4")) goto cleanup;
  if (zk_only) {
    size_t composites = 0;
    uav_clone(&delta_uav, &member);
    if (!expect(revac_uav_apply_update_context(
                    &delta_uav, &public_state, &update, wire, wire_len,
                    0, NULL, &composites),
                "ZK-only 成员应用 DELTA4") ||
        !expect(composites == 0, "ZK-only DELTA4 不应执行 tail 复合")) {
      goto cleanup;
    }
    ok = run_labrador_zk_benchmark(&public_state, &delta_uav, zk_warmups,
                                   zk_repetitions);
    goto cleanup;
  }
  for (size_t run = 0; run < warmups + repetitions; run++) {
    const size_t slot = update.seq - 1; const int saved_used = ta.sync_delta_cache[slot].used; double start, elapsed; size_t hits = 0, composites = 0;
    /* acc_genpp_memwitsync constructs a fresh witness.  Recreate the output
       containers too, so repetitions do not accidentally retain a sampled
       truncation component from an earlier measurement. */
    acc_genpp_wit_clear(&cached_sync); acc_genpp_wit_init(&cached_sync);
    acc_genpp_wit_clear(&fallback_sync); acc_genpp_wit_init(&fallback_sync);
    start = monotonic_seconds(); if (!revac_ta_memwitsync_cached(&cached_sync, &ta, member.handle, &hits, &composites)) goto cleanup; elapsed = monotonic_seconds() - start;
    if (!expect(hits == 1 && composites == 0, "DELTA6 命中时不应重算 tail") || !expect(acc_genpp_memver(&ta.acc_pp, ta.acc, member.handle, &cached_sync), "DELTA6 同步见证正确")) goto cleanup;
    if (run >= warmups) { stats_add(&cached_stats, elapsed); cached_updates = hits; cached_composites = composites; }
    ta.sync_delta_cache[slot].used = 0; hits = 0; composites = 0; start = monotonic_seconds();
    if (!revac_ta_memwitsync_cached(&fallback_sync, &ta, member.handle, &hits, &composites)) { ta.sync_delta_cache[slot].used = saved_used; goto cleanup; }
    elapsed = monotonic_seconds() - start; ta.sync_delta_cache[slot].used = saved_used;
    if (!expect(hits == 0 && composites > 0, "强制未命中应走本地 tail 链") || !expect(acc_genpp_memver(&ta.acc_pp, ta.acc, member.handle, &fallback_sync), "本地回退同步见证正确")) goto cleanup;
    if (run >= warmups) { stats_add(&fallback_stats, elapsed); fallback_updates = hits; fallback_composites = composites; }
    uav_clone(&local_uav, &member); composites = 0; start = monotonic_seconds();
    if (!revac_uav_apply_update_local(&local_uav, &public_state, &update, &composites)) goto cleanup; elapsed = monotonic_seconds() - start;
    if (!expect(composites > 0, "UAV LOCAL 应计算 tail") || !expect(acc_genpp_memver(&ta.acc_pp, ta.acc, local_uav.handle, &local_uav.witness), "UAV LOCAL 更新正确")) goto cleanup;
    if (run >= warmups) { stats_add(&local_stats, elapsed); local_composites = composites; }
    uav_clone(&delta_uav, &member); composites = 0; start = monotonic_seconds();
    if (!revac_uav_apply_update_context(&delta_uav, &public_state, &update, wire, wire_len, 0, NULL, &composites)) goto cleanup; elapsed = monotonic_seconds() - start;
    if (!expect(composites == 0, "DELTA4 应避免 UAV tail 复合") || !expect(acc_genpp_memver(&ta.acc_pp, ta.acc, delta_uav.handle, &delta_uav.witness), "UAV DELTA4 更新正确")) goto cleanup;
    if (run >= warmups) { stats_add(&delta_stats, elapsed); delta_composites = composites; }
  }
  printf("TA 私有 DELTA6 大小：%.2f KiB；公共 BML DELTA4 大小：%.2f KiB。\n", ta.sync_delta_cache[update.seq - 1].delta_data_len / 1024.0, wire_len / 1024.0);
  print_stats("TA MemWitSync：DELTA6 命中", &cached_stats); print_stats("TA MemWitSync：强制本地 tail 回退", &fallback_stats);
  printf("  DELTA6 加速比（回退 / 命中）：%.2fx；命中：%zu 条缓存更新、%zu 次复合；回退：%zu 条缓存更新、%zu 次复合。\n", stats_mean(&fallback_stats) / stats_mean(&cached_stats), cached_updates, cached_composites, fallback_updates, fallback_composites);
  printf("\n=== 消融 C：base-3 下 Δ 增量对 UAV 见证更新的影响 ===\n");
  print_stats("UAV LOCAL：认证广播后本地 EvalFX", &local_stats); print_stats("UAV DELTA4：认证广播后直接叠加 Δ", &delta_stats);
  printf("  Δ 更新加速比（LOCAL / DELTA4）：%.2fx；LOCAL=%zu 次复合，DELTA4=%zu 次复合。\n", stats_mean(&local_stats) / stats_mean(&delta_stats), local_composites, delta_composites);
  printf("  两路径均已对最终见证执行 MemVer 正确性检查。\n"); ok = 1;
  if (!run_randomized_refresh_ablation(&ta, &delta_uav,
                                       randomized_refresh_bits, warmups,
                                       repetitions)) {
    ok = 0;
  }
  if (!run_labrador_zk_benchmark(&public_state, &delta_uav, zk_warmups,
                                 zk_repetitions)) {
    ok = 0;
  }
cleanup:
  free(wire); acc_genpp_wit_clear(&fallback_sync); acc_genpp_wit_clear(&cached_sync); revac_uav_clear(&delta_uav); revac_uav_clear(&local_uav); revac_uav_clear(&member); revac_uav_clear(&revoked); revac_ta_clear(&ta); return ok;
}

int main(void) {
  const size_t warmups = env_nonnegative_count("REVAC_ABLATION_WARMUPS", AB_WARMUPS_DEFAULT);
  const size_t repetitions = env_count("REVAC_ABLATION_REPETITIONS", AB_REPETITIONS_DEFAULT);
  const size_t randomized_refresh_bits =
      env_count("REVAC_RANDOMIZED_REFRESH_BITS",
                AB_RANDOMIZED_REFRESH_BITS_DEFAULT);
  const size_t zk_warmups = env_nonnegative_count("REVAC_ABLATION_ZK_WARMUPS", AB_ZK_WARMUPS_DEFAULT);
  const size_t zk_repetitions = env_count("REVAC_ABLATION_ZK_REPETITIONS", AB_ZK_REPETITIONS_DEFAULT);
  const int zk_only = env_enabled("REVAC_ABLATION_ZK_ONLY");
  int ok;
  arith_setup(); random_init();
  printf("可撤销累加器改进消融测试（不修改生产方案）\n");
  printf("参数 profile：%s\n", REVAC_PARAM_PROFILE_NAME);
  printf("参数 profile ID：%s；发布状态：%s\n",
         REVAC_PARAM_PROFILE_ID_HEX,
         REVAC_PARAM_RELEASE_READY ? "frozen" : "candidate");
  printf("Eval/同步预热 %zu 次、计时 %zu 次；Labrador ZK 预热 %zu 次、稳定计时 %zu 次。\n", warmups, repetitions, zk_warmups, zk_repetitions);
  printf("随机化 FreshSync 原型使用 %zu-bit r。\n",
         randomized_refresh_bits);
  printf("可用 REVAC_ABLATION_WARMUPS、REVAC_ABLATION_REPETITIONS、REVAC_RANDOMIZED_REFRESH_BITS、REVAC_ABLATION_ZK_WARMUPS、REVAC_ABLATION_ZK_REPETITIONS 调整。\n");
  if (!zk_only) print_structural_ablation();
  ok = (zk_only || run_radix_ablation(warmups, repetitions)) &&
       run_sync_and_update_ablation(warmups, repetitions,
                                    randomized_refresh_bits, zk_warmups,
                                    zk_repetitions, zk_only);
  arith_teardown(); printf("基准滚动校验哈希（非密码学，仅防止省略计算）：0x%" PRIx64 "\n", benchmark_digest);
  if (!ok) return 1;
  printf("改进消融测试通过。\n"); return 0;
}
