#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "acc_genpp.h"
#include "arith.h"
#include "random.h"

#ifndef ACC_GENPP_GADGET_BASE
#define ACC_GENPP_GADGET_BASE 2
#endif

#define DEFAULT_REPETITIONS 5
#define DEFAULT_WARMUPS 1

static volatile uint64_t benchmark_sink;

static double now_seconds(void) {
  struct timespec ts;

  if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) {
    perror("clock_gettime");
    exit(1);
  }
  return (double)ts.tv_sec + (double)ts.tv_nsec / 1000000000.0;
}

static size_t env_count(const char *name, size_t fallback) {
  const char *value = getenv(name);
  char *end = NULL;
  unsigned long parsed;

  if (value == NULL || value[0] == '\0') {
    return fallback;
  }
  parsed = strtoul(value, &end, 10);
  if (end == value || *end != '\0' || parsed == 0 || parsed > 1000) {
    fprintf(stderr, "%s must be an integer in [1, 1000]\n", name);
    exit(2);
  }
  return (size_t)parsed;
}

static int compare_double(const void *lhs, const void *rhs) {
  const double a = *(const double *)lhs;
  const double b = *(const double *)rhs;
  return (a > b) - (a < b);
}

static void print_summary(const char *name, const double *samples, size_t count) {
  double *sorted = malloc(count * sizeof(*sorted));
  double sum = 0.0;
  double variance = 0.0;
  double median;

  if (sorted == NULL) {
    fprintf(stderr, "benchmark summary allocation failed\n");
    exit(2);
  }
  for (size_t i = 0; i < count; i++) {
    sorted[i] = samples[i];
    sum += samples[i];
  }
  qsort(sorted, count, sizeof(*sorted), compare_double);
  median = count % 2 == 0
               ? (sorted[count / 2 - 1] + sorted[count / 2]) / 2.0
               : sorted[count / 2];
  for (size_t i = 0; i < count; i++) {
    const double delta = samples[i] - sum / (double)count;
    variance += delta * delta;
  }
  variance /= (double)count;
  printf("%s mean_s=%.9f median_s=%.9f min_s=%.9f max_s=%.9f stddev_s=%.9f\n",
         name, sum / (double)count, median, sorted[0], sorted[count - 1],
         sqrt(variance));
  free(sorted);
}

static uint64_t consume_evalf(const acc_genpp_gmat_t *value) {
  uint64_t checksum = 0;

  for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
    checksum ^= (uint64_t)poly_q_get_coeff(value->blocks[j]->rows[0]->entries[0],
                                           j % PARAM_N);
  }
  return checksum;
}

static uint64_t consume_evalfx(const acc_genpp_wit_t *value) {
  uint64_t checksum = 0;

  checksum ^= (uint64_t)poly_q_get_coeff(value->left->entries[0], 0);
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    const size_t j = i % ACC_GENPP_GADGET_DIM;
    checksum ^= (uint64_t)poly_q_get_coeff(value->right[i][j]->entries[0],
                                           i % PARAM_N);
  }
  return checksum;
}

int main(void) {
  static const uint32_t selectors[] = {
      UINT32_C(0x10203040), UINT32_C(0x89abcdef), UINT32_C(0x31415926),
      UINT32_C(0xdeadbeef), UINT32_C(0x76543210), UINT32_C(0xa5a5a5a5),
      UINT32_C(0x0f1e2d3c), UINT32_C(0xc001d00d)};
  static const uint32_t points[] = {
      UINT32_C(0x55667788), UINT32_C(0x01234567), UINT32_C(0x27182818),
      UINT32_C(0xcafebabe), UINT32_C(0xfedcba98), UINT32_C(0x5a5a5a5a),
      UINT32_C(0x4b5a6978), UINT32_C(0x13579bdf)};
  const size_t repetitions = env_count("REVAC_BENCH_REPS", DEFAULT_REPETITIONS);
  const size_t warmups = env_count("REVAC_BENCH_WARMUPS", DEFAULT_WARMUPS);
  double *evalf_samples = malloc(repetitions * sizeof(*evalf_samples));
  double *tail30_samples = malloc(repetitions * sizeof(*tail30_samples));
  double *evalfx_samples = malloc(repetitions * sizeof(*evalfx_samples));
  acc_genpp_pp_t pp;
  acc_genpp_td_t td;
  acc_genpp_state_t st;
  acc_genpp_gmat_t evalf_out;
  acc_genpp_wit_t evalfx_out;
  poly_q_vec_d accumulator;

  if (evalf_samples == NULL || tail30_samples == NULL ||
      evalfx_samples == NULL) {
    fprintf(stderr, "benchmark sample allocation failed\n");
    return 2;
  }

  arith_setup();
  random_init();
  acc_genpp_pp_init(&pp);
  acc_genpp_td_init(&td);
  acc_genpp_state_init(&st);
  acc_genpp_gmat_init(&evalf_out);
  acc_genpp_wit_init(&evalfx_out);
  poly_q_vec_d_init(accumulator);
  acc_genpp_gen(&pp, &td, &st, accumulator);

  printf("config gadget_base=%d gadget_dim=%d handle_bits=%d acc_threads=%zu repetitions=%zu warmups=%zu\n",
         ACC_GENPP_GADGET_BASE, ACC_GENPP_GADGET_DIM,
         ACC_GENPP_HANDLE_BITS, acc_genpp_parallel_threads(),
         repetitions, warmups);

  for (size_t i = 0; i < warmups; i++) {
    const size_t input = i % (sizeof(selectors) / sizeof(selectors[0]));
    acc_genpp_evalf(&evalf_out, &pp.Ba, selectors[input]);
    benchmark_sink ^= consume_evalf(&evalf_out);
    acc_genpp_evalf_tail1(&evalf_out, &pp.Ba, selectors[input]);
    benchmark_sink ^= consume_evalf(&evalf_out);
    acc_genpp_evalfx_apply(&evalfx_out, &pp.Ba, selectors[input],
                           points[input], pp.ua);
    benchmark_sink ^= consume_evalfx(&evalfx_out);
  }

  for (size_t i = 0; i < repetitions; i++) {
    const size_t input = (i + warmups) %
                         (sizeof(selectors) / sizeof(selectors[0]));
    double start = now_seconds();
    acc_genpp_evalf(&evalf_out, &pp.Ba, selectors[input]);
    evalf_samples[i] = now_seconds() - start;
    benchmark_sink ^= consume_evalf(&evalf_out);
    printf("sample iteration=%zu operation=EvalF seconds=%.9f\n",
           i + 1, evalf_samples[i]);

    start = now_seconds();
    acc_genpp_evalf_tail1(&evalf_out, &pp.Ba, selectors[input]);
    tail30_samples[i] = now_seconds() - start;
    benchmark_sink ^= consume_evalf(&evalf_out);
    printf("sample iteration=%zu operation=Tail30_T31_to_T1 seconds=%.9f\n",
           i + 1, tail30_samples[i]);

    start = now_seconds();
    acc_genpp_evalfx_apply(&evalfx_out, &pp.Ba, selectors[input],
                           points[input], pp.ua);
    evalfx_samples[i] = now_seconds() - start;
    benchmark_sink ^= consume_evalfx(&evalfx_out);
    printf("sample iteration=%zu operation=EvalFX_apply seconds=%.9f\n",
           i + 1, evalfx_samples[i]);
  }

  print_summary("EvalF", evalf_samples, repetitions);
  print_summary("Tail30_T31_to_T1", tail30_samples, repetitions);
  print_summary("EvalFX_apply", evalfx_samples, repetitions);
  printf("benchmark_sink=%llu\n", (unsigned long long)benchmark_sink);

  poly_q_vec_d_clear(accumulator);
  acc_genpp_wit_clear(&evalfx_out);
  acc_genpp_gmat_clear(&evalf_out);
  acc_genpp_state_clear(&st);
  acc_genpp_td_clear(&td);
  acc_genpp_pp_clear(&pp);
  arith_teardown();
  free(evalfx_samples);
  free(tail30_samples);
  free(evalf_samples);
  return 0;
}
