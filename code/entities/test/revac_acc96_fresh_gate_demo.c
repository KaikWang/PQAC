#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "acc96_accumulator.h"
#include "acc96_proof_layout.h"

#define CHECK(condition) do { if (!(condition)) { \
  fprintf(stderr, "fresh gate failure at line %d: %s\n", __LINE__, #condition); \
  goto cleanup; } } while (0)

static void set_signed(acc96_coeff_t *out, int64_t value) {
  acc96_u128 encoded = value < 0
      ? acc96_modulus() - (acc96_u128)(-value) : (acc96_u128)value;
  acc96_coeff_from_u128(out, encoded);
}

static int copied_sample(acc96_relation_witness_t *out,
    const acc96_relation_pp_t *pp, const acc96_trapdoor_t *trapdoor,
    const acc96_poly_t *acc, uint32_t handle, void *opaque,
    const acc96_ctx_t *ctx) {
  (void)pp; (void)trapdoor; (void)acc; (void)handle; (void)ctx;
  memcpy(out, opaque, sizeof(*out));
  return 1;
}

int main(void) {
  acc96_ctx_t ctx, uninitialized = {0};
  acc96_relation_witness_t *w = calloc(1, sizeof(*w));
  acc96_relation_witness_t *out = calloc(1, sizeof(*out));
  acc96_accumulator_t *state = calloc(1, sizeof(*state));
  acc96_poly_t *head[2+ACC96_TRUNC_RIGHT];
  const uint32_t handle = UINT32_C(0x123456);
  const int64_t f = ACC96_FRESH_COEFF_BOUND;
  int ok = 0;
  CHECK(w && out && state && acc96_ctx_init(&ctx));
  CHECK(ACC96_WINDOW_FRESH_COEFF_BOUND == f);
  acc96_coeff_from_u128(&w->tag.coeff[0], 1);
  CHECK(acc96_relation_fresh_witness_check(w, &ctx));
  CHECK(!acc96_relation_fresh_witness_check(NULL, &ctx));
  CHECK(!acc96_relation_fresh_witness_check(w, NULL));
  CHECK(!acc96_relation_fresh_witness_check(w, &uninitialized));
  head[0] = &w->left; head[1] = &w->v12;
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++) head[2+i] = &w->v2[i];
  for (size_t i = 0; i < 2+ACC96_TRUNC_RIGHT; i++) {
    acc96_coeff_t *c = &head[i]->coeff[ACC96_N-1];
    set_signed(c, f); CHECK(acc96_relation_fresh_witness_check(w, &ctx));
    set_signed(c, -f); CHECK(acc96_relation_fresh_witness_check(w, &ctx));
    set_signed(c, f+1); CHECK(!acc96_relation_fresh_witness_check(w, &ctx));
    set_signed(c, -f-1); CHECK(!acc96_relation_fresh_witness_check(w, &ctx));
    set_signed(c, 0);
  }
  /* Noncanonical CRT residues must not become accepted by reconstruction. */
  w->left.coeff[0].limb[0] = ACC96_P0;
  CHECK(!acc96_relation_fresh_witness_check(w, &ctx));
  set_signed(&w->left.coeff[0], 0);
  acc96_coeff_t *ternary[] = {
    &w->right.right[0][0].coeff[0],
    &w->right.right[ACC96_HANDLE_BITS-1][ACC96_GADGET_DIGITS-1].coeff[ACC96_N-1],
    &w->refresh_z[0].coeff[0],
    &w->refresh_z[ACC96_REFRESH_BLOCKS-1].coeff[ACC96_N-1],
  };
  for (size_t i = 0; i < sizeof(ternary)/sizeof(*ternary); i++) {
    set_signed(ternary[i], -1); CHECK(acc96_relation_fresh_witness_check(w, &ctx));
    set_signed(ternary[i], 1); CHECK(acc96_relation_fresh_witness_check(w, &ctx));
    set_signed(ternary[i], 2); CHECK(!acc96_relation_fresh_witness_check(w, &ctx));
    set_signed(ternary[i], -2); CHECK(!acc96_relation_fresh_witness_check(w, &ctx));
    set_signed(ternary[i], 0);
  }
  set_signed(&w->tag.coeff[0], 0);
  CHECK(!acc96_relation_fresh_witness_check(w, &ctx));
  set_signed(&w->tag.coeff[0], 1);
  set_signed(&w->tag.coeff[ACC96_N-1], 1);
  CHECK(!acc96_relation_fresh_witness_check(w, &ctx));
  set_signed(&w->tag.coeff[ACC96_N-1], 0);

  /* A genuinely equation-valid, low-L2 candidate supplied by a callback must
     still fail Add if it exceeds F. No production trapdoor/sampler is mocked
     as secure: this synthetic state isolates the admission gate. */
  set_signed(&w->left.coeff[0], f+1);
  acc96_poly_set(&state->accumulator, &w->left);
  CHECK(acc96_relation_memver(&state->pp, &state->accumulator, handle,
                             w, 1.0e30L, &ctx));
  CHECK(!acc96_add(state, out, handle, 1.0e30L, copied_sample, w, &ctx));
  CHECK(state->sequence == 0);
  CHECK(acc96_poly_equal(&state->accumulator, &w->left));
  for (size_t i = 0; i < ACC96_MAX_MEMBERS; i++)
    CHECK(state->members[i].status == ACC96_MEMBER_EMPTY);
  const unsigned char *bytes = (const unsigned char *)out;
  for (size_t i = 0; i < sizeof(*out); i++) CHECK(bytes[i] == 0);

  /* The larger post-update encoding check must remain distinct from F. */
  CHECK(acc96_window_proof_witness_check(w, &ctx, NULL));
  set_signed(&w->left.coeff[0], f);
  acc96_poly_set(&state->accumulator, &w->left);
  CHECK(acc96_add(state, out, handle, 1.0e30L, copied_sample, w, &ctx));
  CHECK(state->members[0].status == ACC96_MEMBER_ACTIVE);
  puts("Q88_FRESH_GATE_BOUNDARIES_AND_ADMISSION_PASS");
  ok = 1;
cleanup:
  free(state); free(out); free(w);
  return ok ? 0 : 1;
}
