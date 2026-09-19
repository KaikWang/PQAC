#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arith.h"
#include "huffman_compress.h"
#include "random.h"
#include "randombytes.h"
#include "revac_gcs.h"
#include "revac_proof_codec.h"
#include "revac_ta.h"
#include "revac_uav.h"

static int expect_true(int cond, const char *msg) {
  if (!cond) {
    fprintf(stderr, "FAIL: %s\n", msg);
    return 0;
  }
  return 1;
}

static uint32_t load_u32_le(const uint8_t in[4]) {
  return ((uint32_t)in[0]) |
         ((uint32_t)in[1] << 8) |
         ((uint32_t)in[2] << 16) |
         ((uint32_t)in[3] << 24);
}

static void print_size(const char *label, size_t bytes) {
  printf("%s: %zu bytes (%.2f KiB)\n", label, bytes, bytes / 1024.0);
}

static void print_signature_breakdown(const revac_show_signature_size_breakdown_t *sizes) {
  print_size("  signature magic", sizes->magic);
  print_size("  tA", sizes->tA);
  print_size("  compression hints", sizes->hints);
  print_size("  tB", sizes->tB);
  print_size("  z3", sizes->z3);
  print_size("  h", sizes->h);
  print_size("  t1", sizes->t1);
  print_size("  challenge and counter", sizes->challenge_and_counter);
  print_size("  z1", sizes->z1);
  print_size("  z2", sizes->z2);
  print_size("  signature show proof total", sizes->total);
}

static int test_show_compression_identities(void) {
  static const int64_t values[] = {
      0, 1, ((int64_t)1 << (PARAM_D_ROUND_SHOW - 1)) - 1,
      (int64_t)1 << (PARAM_D_ROUND_SHOW - 1),
      ((int64_t)1 << PARAM_D_ROUND_SHOW) - 1,
      PARAM_Q_SHOW - 2, PARAM_Q_SHOW - 1};
  static const int64_t r_values[] = {
      0, PARAM_GAMMA_SHOW / 2 - 1, PARAM_GAMMA_SHOW / 2,
      PARAM_GAMMA_SHOW - 1, PARAM_Q_SHOW - PARAM_GAMMA_SHOW / 2,
      PARAM_Q_SHOW - 2, PARAM_Q_SHOW - 1};
  static const int64_t z_values[] = {0, 2, -2, 3, -3, -5, 5};
  poly_qshow_vec_d value, high, low, recombined;
  poly_qshow_vec_d r, z, sum, hints, recovered, expected;
  int ok;

  poly_qshow_vec_d_init(value);
  poly_qshow_vec_d_init(high);
  poly_qshow_vec_d_init(low);
  poly_qshow_vec_d_init(recombined);
  poly_qshow_vec_d_init(r);
  poly_qshow_vec_d_init(z);
  poly_qshow_vec_d_init(sum);
  poly_qshow_vec_d_init(hints);
  poly_qshow_vec_d_init(recovered);
  poly_qshow_vec_d_init(expected);
  poly_qshow_vec_d_zero(value);
  poly_qshow_vec_d_zero(r);
  poly_qshow_vec_d_zero(z);

  for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); i++) {
    poly_qshow_set_coeff(value->entries[0], i, values[i]);
    poly_qshow_set_coeff(r->entries[0], i, r_values[i]);
    poly_qshow_set_coeff(z->entries[0], i, z_values[i]);
  }

  show_power2round_vec_d(high, low, value);
  poly_qshow_vec_d_mul_scalar(recombined, high,
                              (int64_t)1 << PARAM_D_ROUND_SHOW);
  poly_qshow_vec_d_add(recombined, recombined, low);
  ok = poly_qshow_vec_d_equal(recombined, value);

  poly_qshow_vec_d_add(sum, r, z);
  show_make_g_hints_vec_d(hints, z, r);
  show_use_g_hints_vec_d(recovered, hints, r);
  show_highbits_vec_d(expected, sum);
  ok = ok && poly_qshow_vec_d_equal(recovered, expected);

  poly_qshow_vec_d_clear(expected);
  poly_qshow_vec_d_clear(recovered);
  poly_qshow_vec_d_clear(hints);
  poly_qshow_vec_d_clear(sum);
  poly_qshow_vec_d_clear(z);
  poly_qshow_vec_d_clear(r);
  poly_qshow_vec_d_clear(recombined);
  poly_qshow_vec_d_clear(low);
  poly_qshow_vec_d_clear(high);
  poly_qshow_vec_d_clear(value);
  return ok;
}

static int test_same_handle_exact_lift(
    const revac_public_state_t *public_state, uint32_t handle) {
  poly_q_vec_d commitment, randomness[2];
  int64_t lift[PARAM_D][PARAM_N];
  coeff_q original;
  int ok;

  poly_q_vec_d_init(commitment);
  poly_q_vec_d_init(randomness[0]);
  poly_q_vec_d_init(randomness[1]);
  revac_same_handle_commitment_generate(commitment, randomness,
                                        public_state, handle);
  ok = revac_same_handle_exact_lift(lift, public_state, handle,
                                    randomness, commitment);

  /* A one-coefficient public-commitment modification cannot retain an exact
   * integer q-lift for the unchanged short opening. */
  original = poly_q_get_coeff(commitment->entries[0], 0);
  poly_q_set_coeff(commitment->entries[0], 0,
                   (coeff_q)((original + 1) % PARAM_Q));
  ok = ok && !revac_same_handle_exact_lift(lift, public_state, handle,
                                           randomness, commitment);

  poly_q_vec_d_clear(randomness[1]);
  poly_q_vec_d_clear(randomness[0]);
  poly_q_vec_d_clear(commitment);
  return ok;
}

int main(void) {
  revac_ta_t ta;
  revac_public_state_t public_state;
  revac_uav_t revoked_uav, uav, other_uav;
  revac_gcs_t gcs;
  acc_genpp_update_t revoke_update;
  revac_show_proof_t proof;
  revac_show_proof_t proof2;
  revac_show_proof_t unpacked;
  revac_show_context_t ctx;
  revac_show_context_t other_ctx;
  uint8_t revoked_msg[PARAM_M * PARAM_N / 8];
  uint8_t msg[PARAM_M * PARAM_N / 8];
  uint8_t other_msg[PARAM_M * PARAM_N / 8];
  uint8_t nonce[REVAC_SHOW_NONCE_BYTES];
  uint8_t nonce2[REVAC_SHOW_NONCE_BYTES];
  uint8_t *wire = NULL;
  uint8_t *tampered_wire = NULL;
  size_t wire_len = 0;
  revac_show_signature_size_breakdown_t signature_sizes;
  int online;
  int ok = 1;

  arith_setup();
  random_init();

  ok = expect_true(test_show_compression_identities(),
                   "Power2Round/generalized-hint identity failed");
  if (!ok) {
    goto teardown;
  }

  revac_ta_init(&ta);
  revac_uav_init(&revoked_uav);
  revac_uav_init(&uav);
  revac_uav_init(&other_uav);
  revac_gcs_init(&gcs);
  revac_show_proof_init(&proof);
  revac_show_proof_init(&proof2);
  revac_show_proof_init(&unpacked);
  revac_show_context_init(&ctx);
  revac_show_context_init(&other_ctx);

  revac_ta_setup_q19(&ta);
  revac_ta_public_state(&public_state, &ta);
  revac_uav_keygen(&revoked_uav, UINT64_C(0x1000000000000122),
                   ta.ta_pk.seed);
  revac_uav_keygen(&uav, UINT64_C(0x1000000000000123), ta.ta_pk.seed);
  revac_uav_keygen(&other_uav, UINT64_C(0x1000000000000124),
                   ta.ta_pk.seed);
  revac_gcs_keygen(&gcs, UINT64_C(0x2000000000000123),
                   ta.ta_pk.seed, UINT64_C(123456));
  randombytes(revoked_msg, sizeof(revoked_msg));
  randombytes(msg, sizeof(msg));
  randombytes(other_msg, sizeof(other_msg));
  randombytes(nonce, sizeof(nonce));

  ok = ok && expect_true(revac_ta_register_uav(&ta, &revoked_uav),
                         "revocation-target UAV registration failed");
  ok = ok && expect_true(revac_ta_register_uav(&ta, &uav),
                         "surviving UAV registration failed");
  ok = ok && expect_true(revac_ta_register_uav(&ta, &other_uav),
                         "second surviving UAV registration failed");
  ok = ok && expect_true(revac_ta_register_gcs(&ta, &gcs),
                         "GCS registration failed");
  ok = ok && expect_true(revac_ta_issue_uav(&ta, &revoked_uav, revoked_msg),
                         "revocation-target UAV issuance failed");
  ok = ok && expect_true(revac_ta_issue_uav(&ta, &uav, msg),
                         "surviving UAV issuance failed");
  ok = ok && expect_true(revac_ta_issue_uav(&ta, &other_uav, other_msg),
                         "second surviving UAV issuance failed");
  ok = ok && expect_true(
      revac_ta_revoke_uav(&ta, revoked_uav.handle, &revoke_update),
      "pre-show revocation failed");
  if (!ok) {
    goto cleanup;
  }
  revac_ta_public_state(&public_state, &ta);
  ok = expect_true(test_same_handle_exact_lift(&public_state, uav.handle),
                   "same-handle commitment exact-lift test failed");
  if (!ok) {
    goto cleanup;
  }
  revac_uav_apply_update(&uav, &public_state, &revoke_update);
  revac_uav_apply_update(&other_uav, &public_state, &revoke_update);
  ok = expect_true(acc_genpp_memver(&ta.acc_pp, ta.acc, uav.handle,
                                    &uav.witness),
                   "surviving UAV witness failed after MemWitUp");
  ok = ok && expect_true(acc_genpp_memver(&ta.acc_pp, ta.acc,
                                          other_uav.handle,
                                          &other_uav.witness),
                         "second surviving UAV witness failed after MemWitUp");
  if (!ok) {
    goto cleanup;
  }

  revac_show_context_from_public(&ctx, &public_state, uav.pk.seed);
  revac_show_context_from_public(&other_ctx, &public_state,
                                 other_uav.pk.seed);
  ok = expect_true(revac_uav_show_prove(&proof, &public_state, &uav, msg, nonce),
                   "show proof generation failed");
  if (!ok) {
    goto cleanup;
  }
  {
    const coeff_q saved =
        poly_q_get_coeff(uav.witness.left->entries[0], 0);
    poly_q_set_coeff(uav.witness.left->entries[0], 0,
                     (saved + 1) % PARAM_Q);
    randombytes(nonce2, sizeof(nonce2));
    ok = expect_true(!revac_uav_show_prove(&proof2, &public_state, &uav,
                                            msg, nonce2),
                     "MemVer cache accepted a directly modified witness");
    poly_q_set_coeff(uav.witness.left->entries[0], 0, saved);
    if (!ok) {
      goto cleanup;
    }
    printf("MemVer 摘要缓存篡改回归：已拒绝。\n");
  }
  ok = expect_true(revac_gcs_verify_show_signature(&public_state, &gcs, &ctx,
                                                   &proof),
                   "in-memory signature show proof did not verify");
  if (!ok) {
    goto cleanup;
  }

  print_size("signature show proof raw packed estimate",
             revac_show_signature_raw_size_bytes());
  print_size("signature show proof compact bound",
             revac_show_signature_compact_bound_bytes());
  ok = expect_true(revac_show_signature_compact_size_breakdown(&proof.sig_show,
                                                                &signature_sizes),
                   "could not measure compact signature fields");
  if (!ok) {
    goto cleanup;
  }
  print_signature_breakdown(&signature_sizes);
  /* The shared-handle opening adds two binary R_q^d response vectors to z1.
   * Keep a separate post-link budget so this intentional security cost is not
   * mistaken for a compact-codec regression. */
  ok = expect_true(signature_sizes.total <= 70 * 1024,
                   "same-handle signature proof exceeded 70 KiB budget");
  if (!ok) {
    goto cleanup;
  }
  print_size("in-memory/wire-equivalent proof before codec",
             revac_show_proof_wire_size_bytes(&proof));
  print_size("accumulator zk proof payload before codec",
             proof.acc_zk_proof_len);

  ok = expect_true(revac_show_proof_pack_compact(&wire, &wire_len, &proof, 1),
                   "compact proof packing failed");
  if (!ok) {
    goto cleanup;
  }
  print_size("compact revocable show proof wire bytes", wire_len);

  tampered_wire = malloc(wire_len + 1);
  ok = expect_true(tampered_wire != NULL, "tamper-test allocation failed");
  if (!ok) {
    goto cleanup;
  }
  memcpy(tampered_wire, wire, wire_len);
  tampered_wire[6] ^= 1; /* first nonce byte after RVACP2 */
  ok = expect_true(!revac_gcs_verify_show_signature_wire(&public_state, &gcs, &ctx,
                                                          tampered_wire,
                                                          wire_len),
                   "tampered show nonce unexpectedly verified");
  if (!ok) {
    goto cleanup;
  }
  memcpy(tampered_wire, wire, wire_len);
  tampered_wire[6 + REVAC_SHOW_NONCE_BYTES] ^= 1;
  ok = expect_true(!revac_gcs_verify_show_signature_wire(
                       &public_state, &gcs, &ctx, tampered_wire, wire_len),
                   "tampered same-handle commitment unexpectedly verified");
  if (!ok) {
    goto cleanup;
  }
  memcpy(tampered_wire, wire, wire_len);
  tampered_wire[6 + REVAC_SHOW_NONCE_BYTES +
                REVAC_SAME_HANDLE_COMMITMENT_BYTES + 4 + 6] ^= 1;
  /* first tA byte */
  ok = expect_true(!revac_gcs_verify_show_signature_wire(&public_state, &gcs, &ctx,
                                                          tampered_wire,
                                                          wire_len),
                   "tampered show commitment unexpectedly verified");
  if (!ok) {
    goto cleanup;
  }
  memcpy(tampered_wire, wire, wire_len);
  memset(tampered_wire + 6 + REVAC_SHOW_NONCE_BYTES +
                         REVAC_SAME_HANDLE_COMMITMENT_BYTES + 4 + 6,
         0xff, 5);
  ok = expect_true(!revac_gcs_verify_show_signature_wire(&public_state, &gcs, &ctx,
                                                          tampered_wire,
                                                          wire_len),
                   "non-canonical rounded tA unexpectedly decoded");
  if (!ok) {
    goto cleanup;
  }
  memcpy(tampered_wire, wire, wire_len);
  tampered_wire[6 + REVAC_SHOW_NONCE_BYTES +
                REVAC_SAME_HANDLE_COMMITMENT_BYTES + 4 + 6 +
                COMPACT_TA_BYTES + 4] ^= 1; /* first generalized hint */
  ok = expect_true(!revac_gcs_verify_show_signature_wire(&public_state, &gcs, &ctx,
                                                          tampered_wire,
                                                          wire_len),
                   "tampered compression hint unexpectedly verified");
  if (!ok) {
    goto cleanup;
  }
  ok = expect_true(!revac_gcs_verify_show_signature_wire(&public_state, &gcs, &ctx,
                                                          wire,
                                                          wire_len - 1),
                   "truncated compact proof unexpectedly verified");
  if (!ok) {
    goto cleanup;
  }
  memcpy(tampered_wire, wire, wire_len);
  tampered_wire[wire_len] = 0;
  ok = expect_true(!revac_gcs_verify_show_signature_wire(&public_state, &gcs, &ctx,
                                                          tampered_wire,
                                                          wire_len + 1),
                   "compact proof with trailing byte unexpectedly decoded");
  if (!ok) {
    goto cleanup;
  }
  memcpy(tampered_wire, wire, wire_len);
  {
    size_t sig_len_off = 6 + REVAC_SHOW_NONCE_BYTES +
                         REVAC_SAME_HANDLE_COMMITMENT_BYTES;
    size_t kind_off = sig_len_off + 4 +
                      load_u32_le(tampered_wire + sig_len_off);
    ok = expect_true(kind_off <= wire_len && 4 <= wire_len - kind_off,
                     "could not locate accumulator proof kind");
    if (!ok) {
      goto cleanup;
    }
    memset(tampered_wire + kind_off, 0xff, 4);
  }
  ok = expect_true(!revac_gcs_verify_show_signature_wire(&public_state, &gcs, &ctx,
                                                          tampered_wire,
                                                          wire_len),
                   "unknown accumulator proof kind unexpectedly decoded");
  if (!ok) {
    goto cleanup;
  }

  ok = expect_true(revac_show_proof_unpack_compact(&unpacked, wire, wire_len),
                   "compact proof unpacking failed");
  if (!ok) {
    goto cleanup;
  }
  print_size("accumulator zk proof payload after codec",
             unpacked.acc_zk_proof_len);

  ok = expect_true(revac_gcs_verify_show_signature(&public_state, &gcs, &ctx, &unpacked),
                   "unpacked signature show proof did not verify");
  if (!ok) {
    goto cleanup;
  }

  online = unpacked.acc_zk_kind == REVAC_ACC_ZK_LABRADOR_MEMVER ||
           unpacked.acc_zk_kind == REVAC_ACC_ZK_LABRADOR_MEMVER_DEFLATE;
  if (online) {
    ok = expect_true(revac_gcs_verify_show(&public_state, &gcs, &ctx, &unpacked),
                     "unpacked revocable show proof did not verify");
    if (!ok) {
      goto cleanup;
    }

    /* Generate a valid proof for a second UAV with a different handle, then
     * splice its accumulator component behind the first UAV's signature
     * proof.  Before the shared commitment was wired into both relations,
     * this x_A/x_B component-mixing shape had no common-x equation. */
    randombytes(nonce2, sizeof(nonce2));
    ok = expect_true(uav.handle != other_uav.handle,
                     "mixing-test UAV handles unexpectedly collided");
    ok = ok && expect_true(revac_uav_show_prove(&proof2, &public_state,
                                                &other_uav, other_msg,
                                                nonce2),
                     "second linked show proof generation failed");
    ok = ok && expect_true(!poly_q_vec_d_equal(proof.handle_commitment,
                                                proof2.handle_commitment),
                           "independent same-handle commitments collided");
    ok = ok && expect_true(revac_gcs_verify_show(&public_state, &gcs,
                                                 &other_ctx,
                                                 &proof2),
                           "second linked show proof did not verify");
    if (!ok) {
      goto cleanup;
    }
    ok = expect_true(revac_show_proof_set_acc_zk(
                         &unpacked, proof2.acc_zk_kind,
                         proof2.acc_zk_proof, proof2.acc_zk_proof_len),
                     "could not construct mixed-component proof");
    ok = ok && expect_true(!revac_gcs_verify_show(&public_state, &gcs, &ctx,
                                                  &unpacked),
                           "mixed signature/accumulator proof unexpectedly verified");
    if (!ok) {
      goto cleanup;
    }
    printf("same-handle component-mixing regression: rejected.\n");
  }

  printf("codec demo passed%s.\n",
         online ? " with accumulator proof verification" :
                  " with signature proof verification");

cleanup:
  free(tampered_wire);
  free(wire);
  revac_show_context_clear(&other_ctx);
  revac_show_context_clear(&ctx);
  revac_show_proof_clear(&unpacked);
  revac_show_proof_clear(&proof2);
  revac_show_proof_clear(&proof);
  revac_gcs_clear(&gcs);
  revac_uav_clear(&other_uav);
  revac_uav_clear(&uav);
  revac_uav_clear(&revoked_uav);
  revac_ta_clear(&ta);
teardown:
  arith_teardown();
  return ok ? 0 : 1;
}
