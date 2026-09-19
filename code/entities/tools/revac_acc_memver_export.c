#include <stdio.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdlib.h>

#include "arith.h"
#include "poly_q_sampling.h"
#include "random.h"
#include "randombytes.h"
#include "revac_ta.h"
#include "revac_uav.h"

static void json_poly(const poly_q p) {
  printf("[");
  for (size_t i = 0; i < PARAM_N; i++) {
    if (i) {
      printf(",");
    }
    printf("%ld", (long)poly_q_get_coeff_centered(p, i));
  }
  printf("]");
}

static void json_vec_d(const poly_q_vec_d v) {
  printf("[");
  for (size_t i = 0; i < PARAM_D; i++) {
    if (i) {
      printf(",");
    }
    json_poly(v->entries[i]);
  }
  printf("]");
}

static void json_lift(int64_t lift[PARAM_D][PARAM_N]) {
  printf("[");
  for (size_t row = 0; row < PARAM_D; row++) {
    if (row != 0) {
      printf(",");
    }
    printf("[");
    for (size_t k = 0; k < PARAM_N; k++) {
      if (k != 0) {
        printf(",");
      }
      printf("%" PRId64, lift[row][k]);
    }
    printf("]");
  }
  printf("]");
}

static void json_mat_d_d(const poly_q_mat_d_d m) {
  printf("[");
  for (size_t r = 0; r < PARAM_D; r++) {
    if (r) {
      printf(",");
    }
    printf("[");
    for (size_t c = 0; c < PARAM_D; c++) {
      if (c) {
        printf(",");
      }
      json_poly(m->rows[r]->entries[c]);
    }
    printf("]");
  }
  printf("]");
}

int main(void) {
  revac_ta_t ta;
  revac_public_state_t public_state;
  revac_uav_t revoked_uav, uav;
  acc_genpp_update_t revoke_update;
  poly_q_mat_d_d Atrap, Gblock, Mblock;
  poly_q_vec_d lhs;
  poly_q tag_times_bexpi;
  int64_t lift[PARAM_D][PARAM_N];
  uint8_t revoked_msg[PARAM_M * PARAM_N / 8];
  uint8_t msg[PARAM_M * PARAM_N / 8];
  const char *fresh_env = getenv("REVAC_EXPORT_FRESH");
  int export_fresh = fresh_env != NULL && fresh_env[0] != '\0' &&
                     fresh_env[0] != '0';

  arith_setup();
  random_init();
  revac_ta_init(&ta);
  revac_uav_init(&revoked_uav);
  revac_uav_init(&uav);
  poly_q_mat_d_d_init(Atrap);
  poly_q_mat_d_d_init(Gblock);
  poly_q_mat_d_d_init(Mblock);
  poly_q_vec_d_init(lhs);
  poly_q_init(tag_times_bexpi);

  revac_ta_setup_q19(&ta);
  revac_ta_public_state(&public_state, &ta);
  revac_uav_keygen(&revoked_uav, UINT64_C(0x1000000000000041),
                   ta.ta_pk.seed);
  revac_uav_keygen(&uav, UINT64_C(0x1000000000000042), ta.ta_pk.seed);
  randombytes(revoked_msg, sizeof(revoked_msg));
  randombytes(msg, sizeof(msg));
  if (!revac_ta_register_uav(&ta, &revoked_uav) ||
      !revac_ta_register_uav(&ta, &uav) ||
      !revac_ta_issue_uav(&ta, &revoked_uav, revoked_msg) ||
      !revac_ta_issue_uav(&ta, &uav, msg) ||
      (!export_fresh &&
       !revac_ta_revoke_uav(&ta, revoked_uav.handle, &revoke_update))) {
    fprintf(stderr, "failed to create export instance\n");
    return 1;
  }
  if (!export_fresh) {
    revac_ta_public_state(&public_state, &ta);
    revac_uav_apply_update(&uav, &public_state, &revoke_update);
  }
  if (!acc_genpp_memver(&ta.acc_pp, ta.acc, uav.handle, &uav.witness)) {
    fprintf(stderr, "updated survivor witness failed before export\n");
    return 1;
  }
  if (!acc_genpp_memver_exact_lift(lift, &ta.acc_pp, ta.acc, uav.handle,
                                   &uav.witness, ACC_GENPP_LIFT_AUTO, NULL)) {
    fprintf(stderr, "failed to compute exact lift before export\n");
    return 1;
  }
  poly_q_mat_d_d_uniform(Atrap, ta.acc_pp.A.seed, DOMAIN_SEPARATOR_A, 0);

  poly_q_vec_d_set(lhs, uav.witness.left);
  poly_q_mat_d_d_muladd_vec_d(lhs, Atrap, uav.witness.v12);
  coeff_q bexpi = ACC_GENPP_TRUNC_B_POW_L;
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_mat_d_d_set(Mblock, ta.acc_pp.trunc_B[i]);
    poly_q_mul_scalar(tag_times_bexpi, uav.witness.tag, bexpi);
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_q_sub(Mblock->rows[j]->entries[j], Mblock->rows[j]->entries[j], tag_times_bexpi);
    }
    poly_q_mat_d_d_mulsub_vec_d(lhs, Mblock, uav.witness.v2[i]);
    bexpi *= ACC_GENPP_TRUNC_B;
  }
  for (size_t h = 0; h < ACC_GENPP_HANDLE_BITS; h++) {
    int xbit = (int)((uav.handle >> h) & 1u);
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_mat_d_d_set(Mblock, ta.acc_pp.Ba.bit[h].blocks[j]);
      if (xbit) {
        poly_q_mat_d_d_zero(Gblock);
        for (size_t r = 0; r < PARAM_D; r++) {
          poly_q_set_coeff(Gblock->rows[r]->entries[r], 0,
                           (coeff_q)acc_genpp_gadget_power(j));
        }
        poly_q_mat_d_d_sub(Mblock, Mblock, Gblock);
      }
      poly_q_mat_d_d_muladd_vec_d(lhs, Mblock, uav.witness.right[h][j]);
    }
  }

  printf("{\n");
  printf("\"params\":{\"n\":%d,\"d\":%d,\"q\":%ld,\"handle_bits\":%d,"
         "\"gadget_base\":%d,\"gadget_dim\":%d,\"gadget_signed\":%d,"
         "\"trunc_k_l\":%d,\"deletion_count\":%u},\n",
         PARAM_N, PARAM_D, (long)PARAM_Q, ACC_GENPP_HANDLE_BITS,
         ACC_GENPP_GADGET_BASE, ACC_GENPP_GADGET_DIM,
         ACC_GENPP_GADGET_SIGNED, ACC_GENPP_TRUNC_K_L,
         (unsigned)ta.acc_st.ctr);
  printf("\"x\":%u,\n", uav.handle);
  printf("\"x_bits\":[");
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    if (i) {
      printf(",");
    }
    printf("%u", (unsigned)((uav.handle >> i) & 1u));
  }
  printf("],\n");
  printf("\"Atrap\":");
  json_mat_d_d(Atrap);
  printf(",\n\"trunc_B\":[");
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    if (i) {
      printf(",");
    }
    json_mat_d_d(ta.acc_pp.trunc_B[i]);
  }
  printf("],\n\"Ba\":[");
  for (size_t h = 0; h < ACC_GENPP_HANDLE_BITS; h++) {
    if (h) {
      printf(",");
    }
    printf("[");
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      if (j) {
        printf(",");
      }
      json_mat_d_d(ta.acc_pp.Ba.bit[h].blocks[j]);
    }
    printf("]");
  }
  printf("],\n\"acc\":");
  json_vec_d(ta.acc);
  printf(",\n\"memver_lhs\":");
  json_vec_d(lhs);
  printf(",\n\"lift\":");
  json_lift(lift);
  printf(",\n\"witness\":{");
  printf("\"left\":");
  json_vec_d(uav.witness.left);
  printf(",\"v12\":");
  json_vec_d(uav.witness.v12);
  printf(",\"v2\":[");
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    if (i) {
      printf(",");
    }
    json_vec_d(uav.witness.v2[i]);
  }
  printf("],\"tag\":");
  json_poly(uav.witness.tag);
  printf(",\"right\":[");
  for (size_t h = 0; h < ACC_GENPP_HANDLE_BITS; h++) {
    if (h) {
      printf(",");
    }
    printf("[");
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      if (j) {
        printf(",");
      }
      json_vec_d(uav.witness.right[h][j]);
    }
    printf("]");
  }
  printf("]}}\n");
  printf("\n");

  poly_q_clear(tag_times_bexpi);
  poly_q_vec_d_clear(lhs);
  poly_q_mat_d_d_clear(Mblock);
  poly_q_mat_d_d_clear(Gblock);
  poly_q_mat_d_d_clear(Atrap);
  revac_uav_clear(&uav);
  revac_uav_clear(&revoked_uav);
  revac_ta_clear(&ta);
  arith_teardown();
  return 0;
}
