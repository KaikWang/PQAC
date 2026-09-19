#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "acc_genpp.h"
#include "arith.h"
#include "random.h"

static double monotonic_seconds(void) {
  struct timespec now;

  (void)clock_gettime(CLOCK_MONOTONIC, &now);
  return (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
}

static int gmat_equal(const acc_genpp_gmat_t *lhs,
                      const acc_genpp_gmat_t *rhs) {
  for (size_t block = 0; block < ACC_GENPP_GADGET_DIM; block++) {
    if (!poly_q_mat_d_d_equal(lhs->blocks[block], rhs->blocks[block])) {
      return 0;
    }
  }
  return 1;
}

static int witness_delta_equal(const acc_genpp_wit_t *lhs,
                               const acc_genpp_wit_t *rhs) {
  if (!poly_q_vec_d_equal(lhs->left, rhs->left)) {
    return 0;
  }
  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
    for (size_t gadget = 0; gadget < ACC_GENPP_GADGET_DIM; gadget++) {
      if (!poly_q_vec_d_equal(lhs->right[bit][gadget],
                              rhs->right[bit][gadget])) {
        return 0;
      }
    }
  }
  return 1;
}

int main(void) {
  const uint32_t selector = UINT32_C(0x89abcdef);
  const uint32_t point = UINT32_C(0x01234567);
  acc_genpp_pp_t pp;
  acc_genpp_td_t td;
  acc_genpp_state_t state;
  acc_genpp_gmat_t serial_evalf, parallel_evalf;
  acc_genpp_wit_t serial_evalfx, parallel_evalfx;
  poly_q_vec_d accumulator;
  double serial_evalf_seconds;
  double parallel_evalf_seconds;
  double serial_evalfx_seconds;
  double parallel_evalfx_seconds;
  double started;
  size_t parallel_threads;
  int ok = 0;

  arith_setup();
  random_init();
  acc_genpp_pp_init(&pp);
  acc_genpp_td_init(&td);
  acc_genpp_state_init(&state);
  acc_genpp_gmat_init(&serial_evalf);
  acc_genpp_gmat_init(&parallel_evalf);
  acc_genpp_wit_init(&serial_evalfx);
  acc_genpp_wit_init(&parallel_evalfx);
  poly_q_vec_d_init(accumulator);
  acc_genpp_gen(&pp, &td, &state, accumulator);

  (void)setenv("REVAC_ACC_THREADS", "1", 1);
  started = monotonic_seconds();
  acc_genpp_evalf(&serial_evalf, &pp.Ba, selector);
  serial_evalf_seconds = monotonic_seconds() - started;
  started = monotonic_seconds();
  acc_genpp_evalfx_apply(&serial_evalfx, &pp.Ba, selector, point, pp.ua);
  serial_evalfx_seconds = monotonic_seconds() - started;

  (void)setenv("REVAC_ACC_THREADS", "12", 1);
  parallel_threads = acc_genpp_parallel_threads();
  started = monotonic_seconds();
  acc_genpp_evalf(&parallel_evalf, &pp.Ba, selector);
  parallel_evalf_seconds = monotonic_seconds() - started;
  started = monotonic_seconds();
  acc_genpp_evalfx_apply(&parallel_evalfx, &pp.Ba, selector, point, pp.ua);
  parallel_evalfx_seconds = monotonic_seconds() - started;

  if (!gmat_equal(&serial_evalf, &parallel_evalf)) {
    fprintf(stderr, "serial/parallel EvalF mismatch\n");
    goto cleanup;
  }
  if (!witness_delta_equal(&serial_evalfx, &parallel_evalfx)) {
    fprintf(stderr, "serial/parallel EvalFX mismatch\n");
    goto cleanup;
  }

  printf("parallel workers: %zu\n", parallel_threads);
  printf("EvalF serial %.6f s, parallel %.6f s, speedup %.2fx\n",
         serial_evalf_seconds, parallel_evalf_seconds,
         serial_evalf_seconds / parallel_evalf_seconds);
  printf("EvalFX serial %.6f s, parallel %.6f s, speedup %.2fx\n",
         serial_evalfx_seconds, parallel_evalfx_seconds,
         serial_evalfx_seconds / parallel_evalfx_seconds);
  printf("serial/parallel coefficient equality: passed\n");
  ok = 1;

cleanup:
  poly_q_vec_d_clear(accumulator);
  acc_genpp_wit_clear(&parallel_evalfx);
  acc_genpp_wit_clear(&serial_evalfx);
  acc_genpp_gmat_clear(&parallel_evalf);
  acc_genpp_gmat_clear(&serial_evalf);
  acc_genpp_state_clear(&state);
  acc_genpp_td_clear(&td);
  acc_genpp_pp_clear(&pp);
  arith_teardown();
  return ok ? 0 : 1;
}
