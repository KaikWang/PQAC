#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "acc96_sampleleft.h"

static double seconds(void) {
  struct timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);
  return (double)now.tv_sec + (double)now.tv_nsec / 1e9;
}

int main(void) {
  const uint8_t public_seed[ACC96_SEED_BYTES] = {0x96, 0x51, 0x10};
  const uint8_t trapdoor_seed[ACC96_SEED_BYTES] = {0x96, 0x52, 0x10};
  const uint8_t accumulator_seed[ACC96_SEED_BYTES] = {0x96, 0x53, 0x10};
  const uint32_t handle = UINT32_C(0x00345678);
  const long double test_bound = 1.0e40L;
  acc96_ctx_t ctx;
  acc96_accumulator_t *state = calloc(1, sizeof(*state));
  acc96_relation_witness_t *witness = calloc(1, sizeof(*witness));
  acc96_sampleleft_ctx_t *sampler = NULL;
  double started;
  double construction_seconds;
  double issue_seconds;
  long double norm_square;
  int ok = 0;
  if (state == NULL || witness == NULL || !acc96_ctx_init(&ctx)) goto cleanup;
  if (!acc96_accgen(state, public_seed, trapdoor_seed, accumulator_seed,
                    &ctx))
    goto cleanup;
  started = seconds();
  sampler = acc96_sampleleft_ctx_new(&state->trapdoor, &ctx, 384, 0);
  construction_seconds = seconds() - started;
  if (sampler == NULL) goto cleanup;
  started = seconds();
  if (!acc96_add(state, witness, handle, test_bound,
                 acc96_sampleleft_secure, sampler, &ctx))
    goto cleanup;
  issue_seconds = seconds() - started;
  norm_square = acc96_relation_witness_norm2(witness, &ctx);
  if (!acc96_memver(state, handle, witness, test_bound, &ctx)) goto cleanup;
  printf("q88 secure SampleLeft passed: covariance cache %.3f s, Issue %.3f "
         "s, threads %zu, witness norm^2 %.6Le, max Arb radius %.3e.\n",
         construction_seconds, issue_seconds,
         acc96_sampleleft_thread_count(sampler), norm_square,
         acc96_sampleleft_maximum_arb_radius(sampler));
  ok = 1;
cleanup:
  acc96_sampleleft_ctx_free(sampler);
  free(witness);
  free(state);
  return ok ? 0 : 1;
}
