#include "revocable_ac.h"

#include "fips202.h"
#include "poly_q_sampling.h"
#include "randombytes.h"
#include "revac_ta.h"
#include <stdlib.h>
#include <string.h>

static void revac_absorb_u64(keccak_state *state, uint64_t x) {
  uint8_t buf[8];

  for (size_t i = 0; i < sizeof(buf); i++) {
    buf[i] = (uint8_t)(x >> (8 * i));
  }
  shake256_absorb(state, buf, sizeof(buf));
}

static void revac_absorb_poly_q(keccak_state *state, const poly_q p) {
  for (size_t i = 0; i < PARAM_N; i++) {
    revac_absorb_u64(state, (uint64_t)poly_q_get_coeff(p, i));
  }
}

static void revac_absorb_u32(keccak_state *state, uint32_t x) {
  uint8_t buf[4];

  for (size_t i = 0; i < sizeof(buf); i++) {
    buf[i] = (uint8_t)(x >> (8 * i));
  }
  shake256_absorb(state, buf, sizeof(buf));
}

int revac_public_accumulator_digest(
    uint8_t digest[REVAC_CONTEXT_HASH_BYTES],
    const revac_public_state_t *public_state) {
  static const uint8_t domain[] = "REVAC-ACCUMULATOR-STATE-v1";
  keccak_state state;

  if (digest == NULL || public_state == NULL) return 0;
  shake256_init(&state);
  shake256_absorb(&state, domain, sizeof(domain) - 1);
  revac_absorb_u32(&state, (uint32_t)public_state->acc_profile);
  if (public_state->acc_profile == REVAC_ACC_PROFILE_Q19) {
    if (public_state->acc == NULL) return 0;
    for (size_t i = 0; i < PARAM_D; i++) {
      revac_absorb_poly_q(&state, public_state->acc->entries[i]);
    }
  } else if (public_state->acc_profile == REVAC_ACC_PROFILE_Q96) {
    uint8_t *packed;
    size_t packed_len;

    if (public_state->acc96_acc == NULL) return 0;
    packed = malloc(ACC96_POLY_PACKED_BYTES);
    if (packed == NULL) return 0;
    packed_len = acc96_poly_pack(packed, public_state->acc96_acc);
    if (packed_len != ACC96_POLY_PACKED_BYTES) {
      free(packed);
      return 0;
    }
    shake256_absorb(&state, packed, packed_len);
    memset(packed, 0, packed_len);
    free(packed);
  } else {
    return 0;
  }
  shake256_finalize(&state);
  shake256_squeeze(digest, REVAC_CONTEXT_HASH_BYTES, &state);
  return 1;
}

int revac_show_context_matches_public(
    const revac_show_context_t *ctx,
    const revac_public_state_t *public_state) {
  uint8_t expected[REVAC_CONTEXT_HASH_BYTES];

  if (ctx == NULL || public_state == NULL ||
      ctx->acc_profile != public_state->acc_profile ||
      ctx->acc_seq != public_state->acc_seq ||
      !revac_public_accumulator_digest(expected, public_state)) {
    return 0;
  }
  return memcmp(ctx->acc_digest, expected, sizeof(expected)) == 0;
}

uint32_t revac_handle_from_identity(const user_pk_t *pk, uint64_t id) {
  keccak_state state;
  uint8_t out[4];

  shake256_init(&state);
  revac_absorb_u64(&state, id);
  shake256_absorb(&state, pk->seed, sizeof(pk->seed));
  for (size_t i = 0; i < PARAM_D; i++) {
    revac_absorb_poly_q(&state, pk->t->entries[i]);
  }
  shake256_finalize(&state);
  shake256_squeeze(out, sizeof(out), &state);
  const uint32_t handle = ((uint32_t)out[0]) |
                          ((uint32_t)out[1] << 8) |
                          ((uint32_t)out[2] << 16) |
                          ((uint32_t)out[3] << 24);
  return handle & ((UINT32_C(1) << ACC96_HANDLE_BITS) - 1U);
}

static void revac_same_handle_to_vec(poly_q_vec_d out, uint32_t handle) {
  poly_q_vec_d_zero(out);
  for (size_t bit = 0; bit < ACC96_HANDLE_BITS; bit++) {
    poly_q_set_coeff(out->entries[0], bit,
                     (coeff_q)((handle >> bit) & 1U));
  }
}

void revac_same_handle_public_matrices(
    poly_q_mat_d_d A_link,
    poly_q_mat_d_d D_link,
    const revac_public_state_t *public_state) {
  poly_q_mat_d_d_uniform(A_link, public_state->context_pp_hash,
                         DOMAIN_SEPARATOR_REVAC_LINK_A, 0);
  poly_q_mat_d_d_uniform(D_link, public_state->context_pp_hash,
                         DOMAIN_SEPARATOR_REVAC_LINK_D, 0);
}

void revac_same_handle_commitment_generate(
    poly_q_vec_d commitment,
    poly_q_vec_d randomness[2],
    const revac_public_state_t *public_state,
    uint32_t handle) {
  poly_q_mat_d_d A_link, D_link;
  poly_q_vec_d handle_vec, tmp;
  uint8_t seed[SEED_BYTES];

  poly_q_mat_d_d_init(A_link);
  poly_q_mat_d_d_init(D_link);
  poly_q_vec_d_init(handle_vec);
  poly_q_vec_d_init(tmp);
  randombytes(seed, sizeof(seed));
  poly_q_vec_d_bin_uniform(randomness[0], seed,
                           DOMAIN_SEPARATOR_REVAC_LINK_R, 0);
  poly_q_vec_d_bin_uniform(randomness[1], seed,
                           DOMAIN_SEPARATOR_REVAC_LINK_R, PARAM_D);
  revac_same_handle_public_matrices(A_link, D_link, public_state);
  revac_same_handle_to_vec(handle_vec, handle);
  poly_q_vec_d_set(commitment, randomness[0]);
  poly_q_mat_d_d_muladd_vec_d(commitment, A_link, randomness[1]);
  poly_q_mat_d_d_mul_vec_d(tmp, D_link, handle_vec);
  poly_q_vec_d_add(commitment, commitment, tmp);

  memset(seed, 0, sizeof(seed));
  poly_q_vec_d_clear(tmp);
  poly_q_vec_d_clear(handle_vec);
  poly_q_mat_d_d_clear(D_link);
  poly_q_mat_d_d_clear(A_link);
}

static void revac_same_handle_muladd_exact(
    int64_t out[PARAM_N], const poly_q lhs, const poly_q rhs) {
  for (size_t i = 0; i < PARAM_N; i++) {
    const int64_t a = (int64_t)poly_q_get_coeff_centered(lhs, i);
    for (size_t j = 0; j < PARAM_N; j++) {
      const int64_t b = (int64_t)poly_q_get_coeff_centered(rhs, j);
      const size_t degree = i + j;
      if (degree < PARAM_N) {
        out[degree] += a * b;
      } else {
        out[degree - PARAM_N] -= a * b;
      }
    }
  }
}

int revac_same_handle_exact_lift(
    int64_t lift[PARAM_D][PARAM_N],
    const revac_public_state_t *public_state,
    uint32_t handle,
    const poly_q_vec_d randomness[2],
    const poly_q_vec_d commitment) {
  poly_q_mat_d_d A_link, D_link;
  poly_q_vec_d handle_vec;
  int row_ok[PARAM_D];
  const char *parallel_env = getenv("REVAC_LIFT_PARALLEL");
  const int parallel_rows =
      parallel_env == NULL || parallel_env[0] == '\0' ||
      parallel_env[0] != '0';

  if (public_state == NULL || randomness == NULL || commitment == NULL) {
    return 0;
  }
  poly_q_mat_d_d_init(A_link);
  poly_q_mat_d_d_init(D_link);
  poly_q_vec_d_init(handle_vec);
  revac_same_handle_public_matrices(A_link, D_link, public_state);
  revac_same_handle_to_vec(handle_vec, handle);
#ifdef _OPENMP
#pragma omp parallel for num_threads(PARAM_D) if(parallel_rows) schedule(static)
#endif
  for (size_t row = 0; row < PARAM_D; row++) {
    int64_t expression[PARAM_N] = {0};
    row_ok[row] = 1;
    for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
      expression[coeff] = (int64_t)poly_q_get_coeff_centered(
          randomness[0]->entries[row], coeff);
    }
    for (size_t col = 0; col < PARAM_D; col++) {
      revac_same_handle_muladd_exact(expression,
                                    A_link->rows[row]->entries[col],
                                    randomness[1]->entries[col]);
      revac_same_handle_muladd_exact(expression,
                                    D_link->rows[row]->entries[col],
                                    handle_vec->entries[col]);
    }
    for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
      const int64_t target = (int64_t)poly_q_get_coeff_centered(
          commitment->entries[row], coeff);
      const int64_t difference = target - expression[coeff];
      if ((difference % PARAM_Q) != 0) {
        row_ok[row] = 0;
        break;
      }
      lift[row][coeff] = difference / PARAM_Q;
    }
  }
  poly_q_vec_d_clear(handle_vec);
  poly_q_mat_d_d_clear(D_link);
  poly_q_mat_d_d_clear(A_link);
  for (size_t row = 0; row < PARAM_D; row++) {
    if (!row_ok[row]) {
      return 0;
    }
  }
  return 1;
}

void revac_show_proof_init(revac_show_proof_t *proof) {
  show_proof_init(&proof->sig_show);
  memset(proof->nonce, 0, sizeof(proof->nonce));
  poly_q_vec_d_init(proof->handle_commitment);
  proof->acc_zk_kind = REVAC_ACC_ZK_NONE;
  proof->acc_zk_proof = NULL;
  proof->acc_zk_proof_len = 0;
}

void revac_show_proof_clear(revac_show_proof_t *proof) {
  free(proof->acc_zk_proof);
  proof->acc_zk_proof = NULL;
  proof->acc_zk_proof_len = 0;
  proof->acc_zk_kind = REVAC_ACC_ZK_NONE;
  poly_q_vec_d_clear(proof->handle_commitment);
  show_proof_clear(&proof->sig_show);
}

int revac_show_proof_set_acc_zk(revac_show_proof_t *proof,
                                revac_acc_zk_kind_t kind,
                                const uint8_t *bytes,
                                size_t len) {
  uint8_t *copy = NULL;

  if ((bytes == NULL && len != 0) ||
      (len == 0 && kind != REVAC_ACC_ZK_NONE) ||
      (len != 0 &&
       kind != REVAC_ACC_ZK_LABRADOR_MEMVER &&
       kind != REVAC_ACC_ZK_LABRADOR_MEMVER_DEFLATE &&
       kind != REVAC_ACC_ZK_Q96_LABRADOR_MEMVER &&
       kind != REVAC_ACC_ZK_Q96_LABRADOR_MEMVER_DEFLATE)) {
    return 0;
  }
  if (len != 0) {
    copy = malloc(len);
    if (copy == NULL) {
      return 0;
    }
    memcpy(copy, bytes, len);
  }
  free(proof->acc_zk_proof);
  proof->acc_zk_proof = copy;
  proof->acc_zk_proof_len = len;
  proof->acc_zk_kind = kind;
  return 1;
}

void revac_show_context_init(revac_show_context_t *ctx) {
  memset(ctx->user_seed, 0, sizeof(ctx->user_seed));
  ctx->acc_profile = REVAC_ACC_PROFILE_Q19;
  ctx->acc_seq = 0;
  memset(ctx->acc_digest, 0, sizeof(ctx->acc_digest));
  poly_q_vec_d_init(ctx->acc);
}

void revac_show_context_clear(revac_show_context_t *ctx) {
  poly_q_vec_d_clear(ctx->acc);
}

void revac_show_context_from_ta(revac_show_context_t *ctx,
                                const revac_ta_t *ta,
                                const uint8_t user_seed[SEED_BYTES]) {
  memcpy(ctx->user_seed, user_seed, SEED_BYTES);
  ctx->acc_profile = ta->acc_profile;
  if (ta->acc_profile == REVAC_ACC_PROFILE_Q96) {
    ctx->acc_seq = ta->acc96_state != NULL ? ta->acc96_state->sequence : 0;
    poly_q_vec_d_zero(ctx->acc);
  } else {
    ctx->acc_seq = ta->acc_st.ctr;
    poly_q_vec_d_set(ctx->acc, ta->acc);
  }
  {
    revac_public_state_t public_state;
    revac_ta_public_state(&public_state, ta);
    if (!revac_public_accumulator_digest(ctx->acc_digest, &public_state)) {
      memset(ctx->acc_digest, 0, sizeof(ctx->acc_digest));
    }
  }
}

void revac_show_context_from_public(
    revac_show_context_t *ctx,
    const revac_public_state_t *public_state,
    const uint8_t user_seed[SEED_BYTES]) {
  memcpy(ctx->user_seed, user_seed, SEED_BYTES);
  ctx->acc_profile = public_state->acc_profile;
  ctx->acc_seq = public_state->acc_seq;
  if (public_state->acc_profile == REVAC_ACC_PROFILE_Q96 ||
      public_state->acc == NULL) {
    poly_q_vec_d_zero(ctx->acc);
  } else {
    poly_q_vec_d_set(ctx->acc, public_state->acc);
  }
  if (!revac_public_accumulator_digest(ctx->acc_digest, public_state)) {
    memset(ctx->acc_digest, 0, sizeof(ctx->acc_digest));
  }
}

void revac_show_lang_init(poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
                          poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K],
                          poly_qshow_mat_k_k A3_embed[PARAM_D],
                          poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D],
                          poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
                          poly_qshow_mat_k_k Dx_embed[PARAM_D],
                          poly_qshow_vec_k u_embed[PARAM_D]) {
  for (size_t i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_init(u_embed[i]);
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_init(A_embed[i][j]);
      poly_qshow_mat_k_k_init(Ds_embed[i][j]);
      poly_qshow_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
    }
    for (size_t j = 0; j < PARAM_D * PARAM_K; j++) {
      poly_qshow_mat_k_k_init(B_embed[i][j]);
    }
    poly_qshow_mat_k_k_init(A3_embed[i]);
    for (size_t j = 0; j < PARAM_M; j++) {
      poly_qshow_mat_k_k_init(D_embed[i][j]);
    }
    poly_qshow_mat_k_k_init(Dx_embed[i]);
  }
}

void revac_show_lang_clear(poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
                           poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K],
                           poly_qshow_mat_k_k A3_embed[PARAM_D],
                           poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D],
                           poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
                           poly_qshow_mat_k_k Dx_embed[PARAM_D],
                           poly_qshow_vec_k u_embed[PARAM_D]) {
  for (size_t i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_clear(u_embed[i]);
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_clear(A_embed[i][j]);
      poly_qshow_mat_k_k_clear(Ds_embed[i][j]);
      poly_qshow_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
    }
    for (size_t j = 0; j < PARAM_D * PARAM_K; j++) {
      poly_qshow_mat_k_k_clear(B_embed[i][j]);
    }
    poly_qshow_mat_k_k_clear(A3_embed[i]);
    for (size_t j = 0; j < PARAM_M; j++) {
      poly_qshow_mat_k_k_clear(D_embed[i][j]);
    }
    poly_qshow_mat_k_k_clear(Dx_embed[i]);
  }
}

void revac_show_public_embed(poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
                             poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K],
                             poly_qshow_mat_k_k A3_embed[PARAM_D],
                             poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D],
                             poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
                             poly_qshow_mat_k_k Dx_embed[PARAM_D],
                             poly_qshow_vec_k u_embed[PARAM_D],
                             const sep_pk_t *ta_pk,
                             const uint8_t user_seed[SEED_BYTES]) {
  poly_q_mat_d_d A, Ds[2], Dx;
  poly_q_mat_d_m D;
  poly_q_vec_d a3, u;

  poly_q_mat_d_d_init(A);
  poly_q_mat_d_d_init(Ds[0]);
  poly_q_mat_d_d_init(Ds[1]);
  poly_q_mat_d_d_init(Dx);
  poly_q_mat_d_m_init(D);
  poly_q_vec_d_init(a3);
  poly_q_vec_d_init(u);

  poly_q_mat_d_d_uniform(A, ta_pk->seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_vec_d_uniform(a3, ta_pk->seed, DOMAIN_SEPARATOR_A3);
  poly_q_mat_d_d_uniform(Ds[0], ta_pk->seed, DOMAIN_SEPARATOR_DS, 0);
  poly_q_mat_d_d_uniform(Ds[1], ta_pk->seed, DOMAIN_SEPARATOR_DS, PARAM_D);
  poly_q_mat_d_m_uniform(D, ta_pk->seed, DOMAIN_SEPARATOR_D);
  poly_q_mat_d_d_uniform(Dx, user_seed, DOMAIN_SEPARATOR_DX, 0);
  poly_q_vec_d_uniform(u, ta_pk->seed, DOMAIN_SEPARATOR_U);

  for (size_t i = 0; i < PARAM_D; i++) {
    poly_qshow_subring_embed_vec_k(u_embed[i], u->entries[i], PARAM_Q1_SHOW);
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qshow_subring_embed_mat_k_k(A_embed[i][j], A->rows[i]->entries[j], PARAM_Q1_SHOW);
      poly_qshow_subring_embed_mat_k_k(Ds_embed[i][j], Ds[0]->rows[i]->entries[j], PARAM_Q1_SHOW);
      poly_qshow_subring_embed_mat_k_k(Ds_embed[i][j + PARAM_D], Ds[1]->rows[i]->entries[j], PARAM_Q1_SHOW);
    }
    for (size_t j = 0; j < PARAM_D * PARAM_K; j++) {
      poly_qshow_subring_embed_mat_k_k(B_embed[i][j], ta_pk->B[j / PARAM_D]->rows[i]->entries[j % PARAM_D], PARAM_Q1_SHOW);
    }
    poly_qshow_subring_embed_mat_k_k(A3_embed[i], a3->entries[i], PARAM_Q1_SHOW);
    for (size_t j = 0; j < PARAM_M; j++) {
      poly_qshow_subring_embed_mat_k_k(D_embed[i][j], D->rows[i]->entries[j], PARAM_Q1_SHOW);
    }
    poly_qshow_subring_embed_mat_k_k(Dx_embed[i], Dx->rows[i]->entries[0], PARAM_Q1_SHOW);
  }

  poly_q_vec_d_clear(u);
  poly_q_vec_d_clear(a3);
  poly_q_mat_d_m_clear(D);
  poly_q_mat_d_d_clear(Dx);
  poly_q_mat_d_d_clear(Ds[1]);
  poly_q_mat_d_d_clear(Ds[0]);
  poly_q_mat_d_d_clear(A);
}

void revac_show_derive_crs_seed(uint8_t crs_seed[CRS_SEED_BYTES],
                                const sep_pk_t *ta_pk,
                                const uint8_t user_seed[SEED_BYTES],
                                const uint8_t nonce[REVAC_SHOW_NONCE_BYTES],
                                revac_acc_profile_t acc_profile,
                                uint32_t acc_seq,
                                const uint8_t acc_digest[REVAC_CONTEXT_HASH_BYTES]) {
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, ta_pk->seed, SEED_BYTES);
  shake256_absorb(&state, user_seed, SEED_BYTES);
  shake256_absorb(&state, nonce, REVAC_SHOW_NONCE_BYTES);
  revac_absorb_u32(&state, (uint32_t)acc_profile);
  revac_absorb_u32(&state, acc_seq);
  shake256_absorb(&state, acc_digest, REVAC_CONTEXT_HASH_BYTES);
  shake256_finalize(&state);
  shake256_squeeze(crs_seed, CRS_SEED_BYTES, &state);
}
