#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arith.h"
#include "random.h"
#include "revac_sync.h"
#include "revac_ta.h"
#include "revac_uav.h"

#define CHECK(c) do { if (!(c)) { \
  fprintf(stderr,"fresh sync failure at line %d: %s\n",__LINE__,#c); \
  goto cleanup; } } while (0)

typedef struct { const acc96_relation_witness_t *w; unsigned calls; } bad_sample_t;

static int bad_sample(acc96_relation_witness_t *out,
    const acc96_relation_pp_t *pp, const acc96_trapdoor_t *trapdoor,
    const acc96_poly_t *acc, uint32_t handle, void *opaque,
    const acc96_ctx_t *ctx) {
  (void)pp; (void)trapdoor; (void)acc; (void)handle; (void)ctx;
  bad_sample_t *bad = opaque;
  bad->calls++;
  memcpy(out,bad->w,sizeof(*out));
  return 1;
}

static void set_signed(acc96_coeff_t *out, acc96_i128 v) {
  acc96_coeff_from_u128(out,v < 0 ? acc96_modulus()-(acc96_u128)(-v) : (acc96_u128)v);
}

static int fetch(void *opaque, const revac_sync_request_t *request,
                 uint8_t **wire, size_t *wire_len) {
  return revac_ta_export_memwitsync(opaque,request->handle,request,wire,wire_len);
}

int main(void) {
  revac_ta_t *ta = calloc(1,sizeof(*ta));
  revac_uav_t *uav = calloc(1,sizeof(*uav));
  acc96_relation_witness_t *bad = calloc(1,sizeof(*bad));
  acc96_relation_witness_t *decoded = calloc(1,sizeof(*decoded));
  revac_public_state_t public_state;
  revac_sync_request_t request;
  uint8_t message[PARAM_M*PARAM_N/8] = {0};
  uint8_t *wire = NULL;
  size_t wire_len = 0;
  uint32_t seq = 0;
  int initialized = 0, ok = 0;
  arith_setup(); random_init();
  CHECK(ta && uav && bad && decoded);
  revac_ta_init(ta); revac_uav_init(uav); initialized = 1;
  CHECK(revac_ta_setup(ta));
  revac_uav_keygen(uav,UINT64_C(161001),ta->ta_pk.seed);
  CHECK(revac_ta_register_uav(ta,uav));
  CHECK(revac_ta_issue_uav(ta,uav,message));
  CHECK(acc96_relation_fresh_witness_check(uav->acc96_witness,&ta->acc96_ctx));
  revac_ta_public_state(&public_state,ta);
  CHECK(revac_sync_request_create(&request,uav->handle,uav->last_seq,&public_state));

  /* Add an exact kernel vector to a valid witness, using the real TrapGen
     relation trunc_b[0]=R0+A*R1. The result exceeds F but still satisfies
     the actual membership equation and existing public L2 bound. */
  memcpy(bad,uav->acc96_witness,sizeof(*bad));
  const acc96_i128 a = (acc96_i128)ACC96_FRESH_COEFF_BOUND+1;
  for (size_t i=0;i<ACC96_N;i++) {
    acc96_i128 left = acc96_coeff_to_centered(&bad->left.coeff[i],&ta->acc96_ctx);
    acc96_i128 second = acc96_coeff_to_centered(&bad->v12.coeff[i],&ta->acc96_ctx);
    left += a*acc96_coeff_to_centered(&ta->acc96_state->trapdoor.r[0][0].coeff[i],&ta->acc96_ctx);
    second += a*acc96_coeff_to_centered(&ta->acc96_state->trapdoor.r[1][0].coeff[i],&ta->acc96_ctx);
    if (i==0) left -= a*ACC96_TRUNC_BASE*ACC96_TRUNC_BASE;
    set_signed(&bad->left.coeff[i],left);
    set_signed(&bad->v12.coeff[i],second);
  }
  set_signed(&bad->v2[0].coeff[0],a+acc96_coeff_to_centered(&bad->v2[0].coeff[0],&ta->acc96_ctx));
  CHECK(acc96_relation_memver(public_state.acc96_pp,public_state.acc96_acc,
        uav->handle,bad,public_state.acc96_witness_norm_sq_bound,&ta->acc96_ctx));
  CHECK(!acc96_relation_fresh_witness_check(bad,&ta->acc96_ctx));
  CHECK(!revac_q96_sync_response_pack(&wire,&wire_len,&request,public_state.acc96_seq,
        public_state.acc96_acc,bad,&ta->acc96_ctx,&ta->ta_sk,&ta->ta_pk));
  CHECK(wire==NULL && wire_len==0);
  bad_sample_t injection = {bad,0};
  acc96_sampleleft_fn saved_fn = ta->acc96_sampler;
  void *saved_opaque = ta->acc96_sampler_opaque;
  ta->acc96_sampler = bad_sample; ta->acc96_sampler_opaque = &injection;
  int rejected = !revac_ta_export_memwitsync(ta,uav->handle,&request,&wire,&wire_len);
  ta->acc96_sampler = saved_fn; ta->acc96_sampler_opaque = saved_opaque;
  CHECK(rejected && injection.calls==1 && wire==NULL && wire_len==0);
  CHECK(acc96_relation_fresh_witness_check(uav->acc96_witness,&ta->acc96_ctx));

  CHECK(revac_ta_export_memwitsync(ta,uav->handle,&request,&wire,&wire_len));
  CHECK(revac_q96_sync_response_unpack(decoded,&seq,wire,wire_len,&request,&public_state));
  CHECK(acc96_relation_fresh_witness_check(decoded,&ta->acc96_ctx));
  wire[REVAC_Q96_SYNC_WIRE_HEADER_BYTES] ^= 1;
  CHECK(!revac_q96_sync_response_unpack(decoded,&seq,wire,wire_len,&request,&public_state));
  free(wire); wire = NULL;
  CHECK(revac_uav_recover_updates(uav,&public_state,fetch,ta));
  CHECK(acc96_relation_fresh_witness_check(uav->acc96_witness,&ta->acc96_ctx));
  CHECK(revac_uav_memver_q96(uav,&public_state));
  puts("Q88_FRESH_GATE_SIGNED_SYNC_AND_UAV_INSTALL_PASS");
  ok = 1;
cleanup:
  free(wire);
  if (initialized) { revac_uav_clear(uav); revac_ta_clear(ta); }
  free(decoded); free(bad); free(uav); free(ta);
  arith_teardown();
  return ok ? 0 : 1;
}
