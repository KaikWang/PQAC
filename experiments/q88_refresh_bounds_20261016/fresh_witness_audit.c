#define _POSIX_C_SOURCE 200809L
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "acc96_sampleleft.h"
#include "acc96_truncated_params_generated.h"

/* Clear-relation mechanism check. No synthetic large deletion counter,
   production cap changes, secret dumps, or proof/security certification. */
static int report(const char *label, const acc96_relation_witness_t *w,
                  const acc96_ctx_t *ctx) {
  const acc96_poly_t *head[2 + ACC96_TRUNC_RIGHT];
  head[0] = &w->left;
  head[1] = &w->v12;
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++) head[2+i] = &w->v2[i];
  acc96_u128 maximum = 0;
  for (size_t i = 0; i < 2 + ACC96_TRUNC_RIGHT; i++) {
    for (size_t j = 0; j < ACC96_N; j++) {
      acc96_i128 c = acc96_coeff_to_centered(&head[i]->coeff[j], ctx);
      acc96_u128 a = c < 0 ? (acc96_u128)(-(c+1))+1 : (acc96_u128)c;
      if (a > maximum) maximum = a;
    }
  }
  int ternary = 1;
  for (size_t bit = 0; bit < ACC96_HANDLE_BITS; bit++)
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++)
      for (size_t j = 0; j < ACC96_N; j++) {
        acc96_i128 c = acc96_coeff_to_centered(&w->right.right[bit][digit].coeff[j], ctx);
        ternary &= c >= -1 && c <= 1;
      }
  for (size_t block = 0; block < ACC96_REFRESH_BLOCKS; block++)
    for (size_t j = 0; j < ACC96_N; j++) {
      acc96_i128 c = acc96_coeff_to_centered(&w->refresh_z[block].coeff[j], ctx);
      ternary &= c >= -1 && c <= 1;
    }
  const int capped = maximum <= ((UINT64_C(1)<<30)-1);
  printf("{\"phase\":\"%s\",\"head_max_abs\":%" PRIu64
         ",\"head_cap_pass\":%s,\"right_z_ternary\":%s,\"norm_squared\":%.15Le}\n",
         label, maximum > UINT64_MAX ? UINT64_MAX : (uint64_t)maximum,
         capped ? "true" : "false", ternary ? "true" : "false",
         acc96_relation_witness_norm2(w, ctx));
  fflush(stdout);
  return capped && ternary;
}

int main(void) {
  const uint8_t ps[ACC96_SEED_BYTES] = {0x16,0x10,0x26,0xa1};
  const uint8_t ts[ACC96_SEED_BYTES] = {0x16,0x10,0x26,0xa2};
  const uint8_t as[ACC96_SEED_BYTES] = {0x16,0x10,0x26,0xa3};
  const uint32_t member = UINT32_C(0x345678);
  acc96_ctx_t ctx;
  acc96_accumulator_t *state = calloc(1,sizeof(*state));
  acc96_relation_witness_t *current = calloc(1,sizeof(*current));
  acc96_relation_witness_t *next = calloc(1,sizeof(*next));
  acc96_relation_witness_t *other = calloc(1,sizeof(*other));
  acc96_sampleleft_ctx_t *sampler = NULL;
  int ok = 0;
  if (!state || !current || !next || !other || !acc96_ctx_init(&ctx)) goto cleanup;
  if (!acc96_accgen(state,ps,ts,as,&ctx)) goto cleanup;
  sampler = acc96_sampleleft_ctx_new(&state->trapdoor,&ctx,384,0);
  if (!sampler || !acc96_add(state,current,member,
          ACC96_PUBLIC_WITNESS_NORM_SQ_BOUND,acc96_sampleleft_secure,sampler,&ctx) ||
      !report("initial_fresh",current,&ctx)) goto cleanup;
  for (unsigned round = 0; round < 2; round++) {
    const uint32_t revoked = UINT32_C(0x456780) + round;
    acc96_relation_update_t update;
    if (!acc96_add(state,other,revoked,ACC96_PUBLIC_WITNESS_NORM_SQ_BOUND,
                   acc96_sampleleft_secure,sampler,&ctx) ||
        !acc96_delete(state,&update,revoked,&ctx) ||
        !acc96_memwitup(current,state,member,&update,&ctx) ||
        !acc96_memver(state,member,current,ACC96_PUBLIC_WITNESS_NORM_SQ_BOUND,&ctx))
      goto cleanup;
    printf("{\"phase\":\"updated_old_witness\",\"sequence\":%u,\"norm_squared\":%.15Le}\n",
           state->sequence,acc96_relation_witness_norm2(current,&ctx));
    if (!acc96_sampleleft_secure(next,&state->pp,&state->trapdoor,
                                &state->accumulator,member,sampler,&ctx) ||
        !acc96_memver(state,member,next,ACC96_PUBLIC_WITNESS_NORM_SQ_BOUND,&ctx) ||
        acc96_relation_memver(&state->pp,&state->accumulator,member,next,1.0L,&ctx) ||
        !report("refreshed_same_pp_current_accumulator",next,&ctx)) goto cleanup;
    acc96_relation_witness_t *swap = current; current = next; next = swap;
  }
  puts("Q88_REFRESH_CLEAR_MECHANISM_PASS");
  ok = 1;
cleanup:
  acc96_sampleleft_ctx_free(sampler);
  free(other); free(next); free(current); free(state);
  return ok ? 0 : 1;
}
