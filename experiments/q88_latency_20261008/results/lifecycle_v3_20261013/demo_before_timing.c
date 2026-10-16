#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "randombytes.h"
#include "acc96_exact_lift.h"
#include "acc96_labrador_wire.h"
#include "acc96_proof_layout.h"
#include "revac_gcs.h"
#include "revac_labrador_bridge.h"
#include "revac_proof_codec.h"
#include "revac_q96_link_wire.h"
#include "revac_sync.h"
#include "revac_ta.h"
#include "revac_uav.h"

static int require(int condition, const char *message) {
  if (!condition) {
    fprintf(stderr, "[失败] %s\n", message);
    return 0;
  }
  return 1;
}

static double monotonic_seconds(void) {
  struct timespec now;
  (void)clock_gettime(CLOCK_MONOTONIC, &now);
  return (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
}

static size_t q96_show_repetitions(void) {
  const char *value = getenv("REVAC_Q96_SHOW_REPETITIONS");
  char *end = NULL;
  unsigned long parsed;

  if (value == NULL || value[0] == '\0') return 1;
  parsed = strtoul(value, &end, 10);
  if (end == value || *end != '\0' || parsed < 1 || parsed > 20) return 1;
  return (size_t)parsed;
}

static int q96_prewarm_enabled(void) {
  const char *value = getenv("REVAC_LABRADOR_Q96_PREWARM");
  return value == NULL || value[0] == '\0' || value[0] != '0';
}

static size_t first_differing_bit(uint32_t left, uint32_t right) {
  uint32_t difference = left ^ right;
  size_t bit = 0;

  if (difference == 0) return ACC96_HANDLE_BITS;
  while ((difference & 1U) == 0U) {
    difference >>= 1;
    bit++;
  }
  return bit;
}

static int select_test_identity(const user_pk_t *pk, uint32_t revoked_handle,
                                uint64_t first_id, size_t minimum_bit,
                                size_t maximum_bit, uint64_t *identity) {
  const uint64_t attempts = UINT64_C(1) << 20;

  if (pk == NULL || identity == NULL || minimum_bit > maximum_bit ||
      maximum_bit >= ACC96_HANDLE_BITS || first_id > UINT64_MAX - attempts)
    return 0;
  for (uint64_t offset = 0; offset < attempts; offset++) {
    const uint64_t candidate = first_id + offset;
    const uint32_t handle = revac_handle_from_identity(pk, candidate);
    const size_t bit = first_differing_bit(handle, revoked_handle);
    if (bit >= minimum_bit && bit <= maximum_bit) {
      *identity = candidate;
      return 1;
    }
  }
  return 0;
}

static int write_optional_artifact(const char *environment_name,
                                   const uint8_t *data, size_t length) {
  const char *path = getenv(environment_name);
  FILE *file;
  int ok;
  if (path == NULL || path[0] == '\0') return 1;
  file = fopen(path, "wb");
  if (file == NULL) return 0;
  ok = fwrite(data, 1, length, file) == length;
  if (fclose(file) != 0) ok = 0;
  return ok;
}

typedef struct {
  const revac_ta_t *ta;
  uint32_t authenticated_handle;
} q96_sync_endpoint_t;

static int fetch_q96_sync_from_ta(void *opaque,
                                  const revac_sync_request_t *request,
                                  uint8_t **wire, size_t *wire_len) {
  const q96_sync_endpoint_t *endpoint = opaque;
  return endpoint != NULL &&
         revac_ta_export_memwitsync(endpoint->ta,
                                    endpoint->authenticated_handle,
                                    request, wire, wire_len);
}

int main(void) {
  uint8_t msg[PARAM_M * PARAM_N / 8] = {0};
  uint8_t nonce[REVAC_SHOW_NONCE_BYTES] = {0x96};
  revac_ta_t ta;
  revac_uav_t uav1;
  revac_uav_t uav2;
  revac_uav_t uav_local;
  revac_gcs_t gcs;
  revac_gcs_verify_guard_t verify_guard;
  revac_public_state_t public_state;
  revac_show_proof_t show;
  revac_show_context_t show_context;
  uint8_t first_acc_digest[REVAC_CONTEXT_HASH_BYTES];
  uint64_t fresh_witness_max_abs = 0;
  uint64_t updated_witness_max_abs = 0;
  acc96_relation_update_t update;
  acc96_relation_pp_t *expanded_pp = NULL;
  acc96_exact_lift_t *exact_lift = NULL;
  uint8_t *public_wire = NULL;
  uint8_t *prover_wire = NULL;
  uint8_t *link_wire = NULL;
  uint8_t *instance_wire = NULL;
  uint8_t *show_wire = NULL;
  size_t show_wire_len = 0;
  uint8_t *update_wire = NULL;
  size_t update_wire_len = 0;
  uint8_t *prefetch_wire = NULL;
  size_t prefetch_wire_len = 0;
  uint8_t *manifest_wire = NULL;
  size_t manifest_wire_len = 0;
  uint8_t *delta6_manifest_wire = NULL;
  size_t delta6_manifest_wire_len = 0;
  uint8_t *delta6_layer_wire[REVAC_Q96_DELTA_LAYERS] = {0};
  size_t delta6_layer_wire_len[REVAC_Q96_DELTA_LAYERS] = {0};
  size_t delta6_needed_layers = 0;
  size_t delta6_download_bytes = 0;
  size_t delta6_applied_layers = 0;
  int delta6_used_local = 0;
  size_t delta6_local_applied_layers = 0;
  int delta6_local_used_local = 0;
  size_t sync_wire_len = 0;
  double q96_update_started = 0.0;
  double q96_update_seconds = 0.0;
  double q96_delete_seconds = 0.0;
  double q96_prefetch_started = 0.0;
  double q96_prefetch_seconds = 0.0;
  int q96_update_used_tail4 = 0;
  poly_q_vec_d link_randomness[2];
  poly_q_vec_d link_commitment;
  uint8_t decoded_public_seed[ACC96_SEED_BYTES];
  unsigned exact_lift_bits = 0;
  int ok = 0;
  int verify_guard_initialized = 0;

  arith_setup();
  random_init();
  revac_ta_init(&ta);
  revac_uav_init(&uav1);
  revac_uav_init(&uav2);
  revac_uav_init(&uav_local);
  revac_gcs_init(&gcs);
  memset(&verify_guard, 0, sizeof(verify_guard));
  revac_show_proof_init(&show);
  revac_show_context_init(&show_context);
  poly_q_vec_d_init(link_commitment);
  for (size_t i = 0; i < 2; i++) poly_q_vec_d_init(link_randomness[i]);

  if (!require(revac_gcs_verify_guard_init(&verify_guard, 32.0, 8.0),
               "GCS 并行验证限流器初始化失败")) {
    goto cleanup;
  }
  verify_guard_initialized = 1;

  puts("累加器配置：q88 候选 D（兼容 acc96 API），24-bit handle，"
       "balanced base-26007 x 6，50000-delete window。");

  if (!require(revac_ta_setup(&ta), "默认 q88 Setup/ACB trapdoor 认证失败")) {
    goto cleanup;
  }
  revac_uav_keygen(&uav1, UINT64_C(1001), ta.ta_pk.seed);
  revac_uav_keygen(&uav2, UINT64_C(1002), ta.ta_pk.seed);
  revac_uav_keygen(&uav_local, UINT64_C(2002), ta.ta_pk.seed);
  {
    const uint32_t revoked_handle =
        revac_handle_from_identity(&uav1.pk, uav1.id);
    uint64_t delta6_identity = 0;
    uint64_t local_identity = 0;
    if (!require(select_test_identity(&uav2.pk, revoked_handle,
                                      UINT64_C(1002), 0, 0,
                                      &delta6_identity),
                 "无法构造确定性 DELTA6 命中测试句柄") ||
        !require(select_test_identity(&uav_local.pk, revoked_handle,
                                      UINT64_C(2002),
                                      REVAC_Q96_DELTA_LAYERS,
                                      ACC96_HANDLE_BITS - 1,
                                      &local_identity),
                 "无法构造确定性 DELTA6 LOCAL 回退测试句柄")) {
      goto cleanup;
    }
    uav2.id = delta6_identity;
    uav_local.id = local_identity;
  }
  revac_gcs_keygen(&gcs, UINT64_C(2001), ta.ta_pk.seed, UINT64_C(1));
  randombytes(msg, sizeof(msg));
  if (!require(revac_ta_register_uav(&ta, &uav1), "UAV1 注册失败") ||
      !require(revac_ta_register_uav(&ta, &uav2), "UAV2 注册失败") ||
      !require(revac_ta_register_uav(&ta, &uav_local),
               "LOCAL 回退 UAV 注册失败") ||
      !require(revac_ta_register_gcs(&ta, &gcs), "GCS 注册失败") ||
      !require(revac_ta_issue_uav(&ta, &uav1, msg), "q96 UAV1 Issue 失败") ||
      !require(revac_ta_issue_uav(&ta, &uav2, msg), "q96 UAV2 Issue 失败") ||
      !require(revac_ta_issue_uav(&ta, &uav_local, msg),
               "q96 LOCAL 回退 UAV Issue 失败")) {
    goto cleanup;
  }

  revac_ta_public_state(&public_state, &ta);
  expanded_pp = malloc(sizeof(*expanded_pp));
  public_wire = malloc(ACC96_LABRADOR_PUBLIC_BYTES);
  if (!require(expanded_pp != NULL && public_wire != NULL &&
                   acc96_labrador_public_encode(
                       public_wire, ACC96_LABRADOR_PUBLIC_BYTES,
                       public_state.acc96_public_seed,
                       public_state.acc96_pp) &&
                   acc96_labrador_public_decode(
                       expanded_pp, decoded_public_seed, public_wire,
                       ACC96_LABRADOR_PUBLIC_BYTES,
                       public_state.acc96_ctx) &&
                   memcmp(decoded_public_seed,
                          public_state.acc96_public_seed,
                          ACC96_SEED_BYTES) == 0 &&
                   memcmp(expanded_pp, public_state.acc96_pp,
                          sizeof(*expanded_pp)) == 0,
               "q96 二进制公参不能从 seed + trunc_b 确定性重建")) {
    goto cleanup;
  }
  free(expanded_pp);
  expanded_pp = NULL;
  revac_show_context_from_public(&show_context, &public_state, uav2.pk.seed);
  if (!require(public_state.acc_profile == REVAC_ACC_PROFILE_Q96,
               "public state 没有选择 q96") ||
      !require(public_state.window_id == 0 &&
                   show_context.window_id == public_state.window_id,
               "q96 Setup 没有从 window_id=0 开始") ||
      !require(memcmp(public_state.parameter_digest,
                      public_state.context_pp_hash,
                      REVAC_CONTEXT_HASH_BYTES) == 0 &&
                   memcmp(show_context.parameter_digest,
                          public_state.parameter_digest,
                          REVAC_CONTEXT_HASH_BYTES) == 0,
               "q96 show context 没有显式绑定参数摘要") ||
      !require(show_context.acc_profile == REVAC_ACC_PROFILE_Q96 &&
                   revac_show_context_matches_public(&show_context,
                                                     &public_state),
               "q96 show context 没有绑定真实累加器摘要") ||
      !require(acc96_window_proof_witness_check(
                   uav2.acc96_witness, public_state.acc96_ctx,
                   &fresh_witness_max_abs),
               "新签发 q96 witness 超出 radix-8191 window 编码门") ||
      !require(revac_gcs_verify_uav(&public_state, &gcs, &uav1, msg),
               "GCS q96 UAV1 明文门控失败") ||
      !require(revac_gcs_verify_uav(&public_state, &gcs, &uav2, msg),
               "GCS q96 UAV2 明文门控失败")) {
      goto cleanup;
  }
  show_context.window_id ^= UINT32_C(1);
  if (!require(!revac_show_context_matches_public(&show_context, &public_state),
               "篡改 q96 window_id 后 context 仍被接受")) {
    goto cleanup;
  }
  show_context.window_id ^= UINT32_C(1);
  show_context.parameter_digest[0] ^= UINT8_C(1);
  if (!require(!revac_show_context_matches_public(&show_context, &public_state),
               "篡改 q96 参数摘要后 context 仍被接受")) {
    goto cleanup;
  }
  show_context.parameter_digest[0] ^= UINT8_C(1);
  if (revac_labrador_online_enabled() && q96_prewarm_enabled()) {
    double started = monotonic_seconds();
    double prover_seconds;
    double verifier_seconds;

    if (!require(revac_uav_memver_cache_prewarm(&uav2, &public_state),
                 "q96 UAV MemVer 摘要缓存预热失败") ||
        !require(revac_labrador_q96_prover_prewarm(&public_state, &uav2),
                 "q96 UAV prover 离线预热失败")) {
      goto cleanup;
    }
    prover_seconds = monotonic_seconds() - started;
    started = monotonic_seconds();
    if (!require(revac_labrador_q96_verifier_prewarm(&public_state),
                 "q96 GCS verifier 离线预热失败")) {
      goto cleanup;
    }
    verifier_seconds = monotonic_seconds() - started;
    printf("q96 离线预热（不计入在线 show）：UAV %.3f 秒，GCS %.3f 秒。"
           "设置 REVAC_LABRADOR_Q96_PREWARM=0 可测冷启动。\n",
           prover_seconds, verifier_seconds);
  }
  if (revac_labrador_online_enabled()) {
    const size_t repetitions = q96_show_repetitions();
    for (size_t run = 0; run < repetitions; run++) {
      const double prove_started = monotonic_seconds();
      double prove_seconds;
      double verify_seconds;
      double pack_started, pack_seconds;

      if (run != 0) randombytes(nonce, sizeof(nonce));
      free(show_wire);
      show_wire = NULL;
      show_wire_len = 0;
      if (!require(revac_uav_show_prove(
                       &show, &public_state, &uav2, msg, nonce),
                   "q96 在线完整 show 生成失败")) {
        goto cleanup;
      }
      prove_seconds = monotonic_seconds() - prove_started;
      pack_started = monotonic_seconds();
      if (!require(show.acc_zk_kind == REVAC_ACC_ZK_Q96_LABRADOR_MEMVER,
                   "q96 show 没有使用独立 proof kind") ||
          !require(revac_show_proof_pack_compact(
                       &show_wire, &show_wire_len, &show, 1),
                   "q96 show compact/deflate 编码失败")) {
        goto cleanup;
      }
      pack_seconds = monotonic_seconds() - pack_started;
      prove_seconds += pack_seconds; /* Complete packet generation, including deflate. */
      {
        const double verify_started = monotonic_seconds();
        if (!require(revac_gcs_verify_show_wire_guarded_parallel(
                         &public_state, &gcs, &show_context, &verify_guard,
                         show_wire, show_wire_len),
                     "GCS 从 wire 验证 q96 匿名签名/累加器/"
                     "同句柄失败")) {
          goto cleanup;
        }
        verify_seconds = monotonic_seconds() - verify_started;
      }
      printf("q96 show 第 %zu/%zu 次：UAV 生成 %.3f 秒，"
             "GCS wire 验证 %.3f 秒；累加器 proof envelope %.2f KiB，"
             "压缩后完整 show packet %.2f KiB。\n",
             run + 1, repetitions, prove_seconds, verify_seconds,
             show.acc_zk_proof_len / 1024.0, show_wire_len / 1024.0);
      printf("[q88 full packet measurement] {\"run\":%zu,"
             "\"generation_seconds\":%.9f,\"verification_seconds\":%.9f,"
             "\"serialization_seconds\":%.9f,\"accumulator_bytes\":%zu,"
             "\"packet_bytes\":%zu}\n",
             run + 1, prove_seconds, verify_seconds, pack_seconds,
             show.acc_zk_proof_len, show_wire_len);
      if (!require(write_optional_artifact("REVAC_Q88_SHOW_PACKET_DUMP",
                                          show_wire, show_wire_len),
                   "写出完整 q88 show packet 失败") ||
          !require(write_optional_artifact("REVAC_Q88_ACCUMULATOR_DUMP",
                                          show.acc_zk_proof,
                                          show.acc_zk_proof_len),
                   "写出 q88 accumulator envelope 失败")) goto cleanup;
      {
        revac_public_state_t next_window_state = public_state;
        next_window_state.window_id++;
        if (!require(!revac_gcs_verify_show_accumulator(
                         &next_window_state, &gcs, &show_context, &show),
                     "旧窗口 q96 show 被新窗口状态接受")) {
          goto cleanup;
        }
      }
    }
  } else if (!require(!revac_uav_show_prove(
                           &show, &public_state, &uav2, msg, nonce),
                       "关闭在线证明时 q96 show 没有失败关闭")) {
    goto cleanup;
  }
  q96_prefetch_started = monotonic_seconds();
  if (!require(revac_ta_export_q96_tail4_prefetch(
                   &ta, uav1.handle, &prefetch_wire, &prefetch_wire_len),
               "TA 在撤销前导出 q96 T4 prefetch 失败")) {
    goto cleanup;
  }
  q96_prefetch_seconds = monotonic_seconds() - q96_prefetch_started;
  prefetch_wire[8] ^= UINT8_C(1);
  if (!require(!revac_uav_cache_q96_tail4_prefetch(
                   &uav2, &public_state, prefetch_wire, prefetch_wire_len),
               "篡改候选句柄的 q96 prefetch 未被拒绝")) {
    goto cleanup;
  }
  prefetch_wire[8] ^= UINT8_C(1);
  prefetch_wire[REVAC_Q96_PREFETCH_WIRE_HEADER_BYTES] ^= UINT8_C(1);
  if (!require(!revac_uav_cache_q96_tail4_prefetch(
                   &uav2, &public_state, prefetch_wire, prefetch_wire_len),
               "篡改签名的 q96 prefetch 未被拒绝")) {
    goto cleanup;
  }
  prefetch_wire[REVAC_Q96_PREFETCH_WIRE_HEADER_BYTES] ^= UINT8_C(1);
  prefetch_wire[prefetch_wire_len - 1] ^= UINT8_C(1);
  if (!require(!revac_uav_cache_q96_tail4_prefetch(
                   &uav2, &public_state, prefetch_wire, prefetch_wire_len),
               "篡改 T4 payload 的 q96 prefetch 未被拒绝")) {
    goto cleanup;
  }
  prefetch_wire[prefetch_wire_len - 1] ^= UINT8_C(1);
  if (!require(revac_uav_cache_q96_tail4_prefetch(
                   &uav2, &public_state, prefetch_wire, prefetch_wire_len),
               "UAV 在撤销前验签并缓存 q96 T4 prefetch 失败")) {
    goto cleanup;
  }
  q96_update_started = monotonic_seconds();
  if (!require(revac_ta_revoke_uav_q96(&ta, uav1.handle, &update),
               "q96 Delete 失败")) {
    goto cleanup;
  }
  q96_delete_seconds = monotonic_seconds() - q96_update_started;
  memcpy(first_acc_digest, show_context.acc_digest, sizeof(first_acc_digest));
  show_context.acc_digest[0] ^= 1U;
  if (!require(!revac_show_context_matches_public(&show_context, &public_state),
               "篡改 q96 累加器摘要后 context 仍被接受")) {
    goto cleanup;
  }
  show_context.acc_digest[0] ^= 1U;

  revac_ta_public_state(&public_state, &ta);
  if (!require(!revac_show_context_matches_public(&show_context, &public_state),
               "撤销后旧 q96 show context 仍匹配新状态")) {
    goto cleanup;
  }
  revac_show_context_from_public(&show_context, &public_state, uav2.pk.seed);
  if (!require(revac_ta_export_q96_update(
                   &ta, update.sequence, &update_wire, &update_wire_len),
               "TA 导出 q96 签名 BML 更新失败") ||
      !require(revac_ta_export_q96_update_manifest(
                   &ta, update.sequence, &manifest_wire, &manifest_wire_len),
               "TA 导出 q96 实时撤销 manifest 失败") ||
      !require(revac_ta_export_q96_delta6_manifest(
                   &ta, update.sequence, &delta6_manifest_wire,
                   &delta6_manifest_wire_len),
               "TA 导出 q88 DELTA6 manifest 失败") ||
      !require(write_optional_artifact("REVAC_Q96_UPDATE_ARTIFACT",
                                       update_wire, update_wire_len),
               "写出 q96 T4 packet artifact 失败") ||
      !require(write_optional_artifact("REVAC_Q96_MANIFEST_ARTIFACT",
                                       manifest_wire, manifest_wire_len),
               "写出 q96 manifest artifact 失败")) {
    goto cleanup;
  }
  {
    uint32_t difference = update.revoked_handle ^ uav2.handle;
    while ((difference & 1U) == 0U) {
      difference >>= 1;
      delta6_needed_layers++;
    }
    delta6_needed_layers++;
    if (delta6_needed_layers <= REVAC_Q96_DELTA_LAYERS) {
      for (size_t i = 0; i < delta6_needed_layers; i++) {
        if (!require(revac_ta_export_q96_delta6_layer(
                         &ta, update.sequence, i, &delta6_layer_wire[i],
                         &delta6_layer_wire_len[i]),
                     "TA 按需导出 q88 DELTA6 层失败"))
          goto cleanup;
        delta6_download_bytes += delta6_layer_wire_len[i];
      }
    } else {
      delta6_needed_layers = 0;
    }
  }
  if (!require(first_differing_bit(update.revoked_handle, uav2.handle) == 0,
               "确定性 DELTA6 命中句柄未落在第 0 层") ||
      !require(first_differing_bit(update.revoked_handle, uav_local.handle) >=
                   REVAC_Q96_DELTA_LAYERS,
               "确定性 LOCAL 句柄仍落在 DELTA6 覆盖范围") ||
      !require(revac_uav_apply_update_q96_delta6(
                   &uav_local, &public_state, delta6_manifest_wire,
                   delta6_manifest_wire_len, NULL, NULL, 0,
                   &delta6_local_applied_layers, &delta6_local_used_local),
               "完整签名 DELTA6 API 的 LOCAL 回退失败") ||
      !require(delta6_local_used_local && delta6_local_applied_layers == 0,
               "DELTA6 覆盖外句柄未报告零层 LOCAL 回退") ||
      !require(revac_uav_memver_q96(&uav_local, &public_state),
               "DELTA6 LOCAL 回退后的 witness 未通过 MemVer")) {
    goto cleanup;
  }
  update_wire[12] ^= UINT8_C(1);
  if (!require(!revac_uav_apply_update_q96_wire(
                   &uav2, &public_state, update_wire, update_wire_len),
               "篡改撤销句柄的 q96 广播未被拒绝")) {
    goto cleanup;
  }
  update_wire[12] ^= UINT8_C(1);
  update_wire[125] ^= UINT8_C(1);
  if (!require(!revac_uav_apply_update_q96_wire(
                   &uav2, &public_state, update_wire, update_wire_len),
               "篡改 q96 T4 codec 的广播未被拒绝")) {
    goto cleanup;
  }
  update_wire[125] ^= UINT8_C(1);
  update_wire[update_wire_len - 1] ^= UINT8_C(1);
  if (!require(!revac_uav_apply_update_q96_wire(
                   &uav2, &public_state, update_wire, update_wire_len),
               "篡改 q96 T4 payload 的广播未被拒绝")) {
    goto cleanup;
  }
  update_wire[update_wire_len - 1] ^= UINT8_C(1);

  /* A valid manifest alone is sufficient for an exact LOCAL fallback. */
  {
    revac_uav_t uncached;
    revac_uav_init(&uncached);
    uncached.handle = uav2.handle;
    uncached.last_seq = uav2.last_seq;
    uncached.issued = 1;
    uncached.acc96_witness = malloc(sizeof(*uncached.acc96_witness));
    if (!require(uncached.acc96_witness != NULL,
                 "分配 q96 manifest LOCAL 回退 witness 失败")) {
      revac_uav_clear(&uncached);
      goto cleanup;
    }
    memcpy(uncached.acc96_witness, uav2.acc96_witness,
           sizeof(*uncached.acc96_witness));
    if (!require(revac_uav_apply_update_q96_manifest(
                     &uncached, &public_state, manifest_wire,
                     manifest_wire_len),
                 "没有预分发 T4 时 manifest LOCAL 回退失败") ||
        !require(revac_uav_memver_q96(&uncached, &public_state),
                 "manifest LOCAL 回退后的 witness 未通过 MemVer")) {
      revac_uav_clear(&uncached);
      goto cleanup;
    }
    revac_uav_clear(&uncached);
  }
  manifest_wire[12] ^= UINT8_C(1);
  if (!require(!revac_uav_apply_update_q96_manifest(
                   &uav2, &public_state, manifest_wire, manifest_wire_len),
               "篡改撤销句柄的 q96 manifest 未被拒绝")) {
    goto cleanup;
  }
  manifest_wire[12] ^= UINT8_C(1);
  manifest_wire[88] ^= UINT8_C(1);
  if (!require(!revac_uav_apply_update_q96_manifest(
                   &uav2, &public_state, manifest_wire, manifest_wire_len),
               "篡改 T4 hash 的 q96 manifest 未被拒绝")) {
    goto cleanup;
  }
  manifest_wire[88] ^= UINT8_C(1);
  manifest_wire[REVAC_Q96_UPDATE_WIRE_HEADER_BYTES] ^= UINT8_C(1);
  if (!require(!revac_uav_apply_update_q96_manifest(
                   &uav2, &public_state, manifest_wire, manifest_wire_len),
               "篡改签名的 q96 manifest 未被拒绝")) {
    goto cleanup;
  }
  manifest_wire[REVAC_Q96_UPDATE_WIRE_HEADER_BYTES] ^= UINT8_C(1);
  delta6_manifest_wire[12] ^= UINT8_C(1);
  if (!require(!revac_uav_apply_update_q96_delta6(
                   &uav2, &public_state, delta6_manifest_wire,
                   delta6_manifest_wire_len,
                   (const uint8_t **)delta6_layer_wire,
                   delta6_layer_wire_len, delta6_needed_layers,
                   NULL, NULL),
               "篡改撤销句柄的 q88 DELTA6 manifest 未被拒绝"))
    goto cleanup;
  delta6_manifest_wire[12] ^= UINT8_C(1);
  if (delta6_needed_layers > 0) {
    const uint32_t last_seq_before = uav2.last_seq;
    acc96_relation_witness_t *witness_before =
        malloc(sizeof(*witness_before));
    if (!require(witness_before != NULL,
                 "分配 DELTA6 原子拒绝快照失败"))
      goto cleanup;
    memcpy(witness_before, uav2.acc96_witness, sizeof(*witness_before));
    delta6_layer_wire[0][delta6_layer_wire_len[0] - 1] ^= UINT8_C(1);
    if (!require(!revac_uav_apply_update_q96_delta6(
                     &uav2, &public_state, delta6_manifest_wire,
                     delta6_manifest_wire_len,
                     (const uint8_t **)delta6_layer_wire,
                     delta6_layer_wire_len, delta6_needed_layers,
                     NULL, NULL),
                 "篡改 q88 DELTA6 层未被拒绝") ||
        !require(uav2.last_seq == last_seq_before &&
                     memcmp(witness_before, uav2.acc96_witness,
                            sizeof(*witness_before)) == 0,
                 "篡改 DELTA6 层后 witness/sequence 发生变化")) {
      free(witness_before);
      goto cleanup;
    }
    delta6_layer_wire[0][delta6_layer_wire_len[0] - 1] ^= UINT8_C(1);
    if (!require(!revac_uav_apply_update_q96_delta6(
                     &uav2, &public_state, delta6_manifest_wire,
                     delta6_manifest_wire_len,
                     (const uint8_t **)delta6_layer_wire,
                     delta6_layer_wire_len, delta6_needed_layers - 1,
                     NULL, NULL),
                 "缺少 q88 DELTA6 必需层时未拒绝") ||
        !require(uav2.last_seq == last_seq_before &&
                     memcmp(witness_before, uav2.acc96_witness,
                            sizeof(*witness_before)) == 0,
                 "缺少 DELTA6 层后 witness/sequence 发生变化")) {
      free(witness_before);
      goto cleanup;
    }
    delta6_layer_wire[0][8] ^= UINT8_C(1);
    if (!require(!revac_uav_apply_update_q96_delta6(
                     &uav2, &public_state, delta6_manifest_wire,
                     delta6_manifest_wire_len,
                     (const uint8_t **)delta6_layer_wire,
                     delta6_layer_wire_len, delta6_needed_layers,
                     NULL, NULL),
                 "跨序列重放 q88 DELTA6 层未被拒绝")) {
      free(witness_before);
      goto cleanup;
    }
    delta6_layer_wire[0][8] ^= UINT8_C(1);
    delta6_layer_wire[0][12] = UINT8_C(1);
    if (!require(!revac_uav_apply_update_q96_delta6(
                     &uav2, &public_state, delta6_manifest_wire,
                     delta6_manifest_wire_len,
                     (const uint8_t **)delta6_layer_wire,
                     delta6_layer_wire_len, delta6_needed_layers,
                     NULL, NULL),
                 "错序 q88 DELTA6 层未被拒绝")) {
      free(witness_before);
      goto cleanup;
    }
    delta6_layer_wire[0][12] = UINT8_C(0);
    free(witness_before);
  }
  q96_update_used_tail4 =
      ((update.revoked_handle ^ uav2.handle) & UINT32_C(0xf)) != 0;
  if (!require(!revac_gcs_verify_uav(&public_state, &gcs, &uav1, msg),
               "已撤销 UAV1 仍通过 q96 门控")) {
    goto cleanup;
  }
  q96_update_started = monotonic_seconds();
  if (!require(revac_uav_apply_update_q96_delta6(
                   &uav2, &public_state, delta6_manifest_wire,
                   delta6_manifest_wire_len,
                   (const uint8_t **)delta6_layer_wire,
                   delta6_layer_wire_len, delta6_needed_layers,
                   &delta6_applied_layers, &delta6_used_local),
               "UAV2 验签 manifest 并选择性应用 q88 DELTA6 失败")) {
    goto cleanup;
  }
  q96_update_seconds = monotonic_seconds() - q96_update_started;
  if (!require(revac_gcs_verify_uav(&public_state, &gcs, &uav2, msg),
               "更新后 UAV2 q96 MemVer 失败") ||
      !require(!revac_ta_revoke_uav_q96(&ta, uav1.handle, &update),
               "重复 q96 Delete 没有被拒绝")) {
    goto cleanup;
  }
  if (!acc96_window_proof_witness_check(
          uav2.acc96_witness, public_state.acc96_ctx,
          &updated_witness_max_abs)) {
    fprintf(stderr, "[诊断] 更新后 q96 witness 最大绝对系数=%llu，"
                    "超过 radix-8191 gate=%lld。\n",
            (unsigned long long)updated_witness_max_abs,
            (long long)ACC96_WINDOW_COEFF_BOUND);
    goto cleanup;
  }
  exact_lift = malloc(sizeof(*exact_lift));
  if (!require(exact_lift != NULL &&
                   acc96_memver_exact_lift(
                       exact_lift, &exact_lift_bits,
                       public_state.acc96_pp, public_state.acc96_acc,
                       uav2.handle, uav2.acc96_witness,
                       public_state.acc96_ctx),
               "三条 q96 CRT MemVer 精确整数 lift 检查失败")) {
    goto cleanup;
  }
  for (size_t lane = 0; lane < ACC96_LIMBS; lane++) {
    for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
      int32_t limbs[ACC96_LIFT_PROOF_LIMBS];
      if (!require(acc96_lift_proof_limb_decompose(
                       limbs, exact_lift->coeff[lane][coefficient]) &&
                       acc96_lift_proof_limb_reconstruct(limbs) ==
                           exact_lift->coeff[lane][coefficient],
                   "q88 exact lift 超出 6xradix-8191 编码门")) {
        goto cleanup;
      }
    }
  }
  prover_wire = malloc(ACC96_LABRADOR_PROVER_BYTES);
  link_wire = malloc(REVAC_Q96_LINK_WIRE_BYTES);
  instance_wire = malloc(REVAC_Q96_INSTANCE_WIRE_BYTES);
  revac_same_handle_commitment_generate(link_commitment, link_randomness,
                                        &public_state, uav2.handle);
  if (!require(prover_wire != NULL &&
                   acc96_labrador_prover_encode(
                       prover_wire, ACC96_LABRADOR_PROVER_BYTES,
                       public_state.acc96_seq, uav2.handle,
                       public_state.acc96_acc, uav2.acc96_witness,
                       exact_lift, public_state.acc96_ctx) &&
                   memcmp(prover_wire, ACC96_LABRADOR_PROVER_MAGIC,
                          ACC96_LABRADOR_MAGIC_BYTES) == 0,
               "q96 C 端预分解 Labrador prover payload 生成失败")) {
    goto cleanup;
  }
  if (!require(link_wire != NULL &&
                   revac_q96_link_wire_encode(
                       link_wire, REVAC_Q96_LINK_WIRE_BYTES, &public_state,
                       uav2.handle, link_randomness, link_commitment),
               "q96/q19 同句柄 link prover payload 生成失败")) {
    goto cleanup;
  }
  if (!require(instance_wire != NULL &&
                   revac_q96_instance_wire_encode(
                       instance_wire, REVAC_Q96_INSTANCE_WIRE_BYTES,
                       &public_state, link_commitment),
               "q96 GCS 公开验证实例生成失败")) {
    goto cleanup;
  }
  if (!require(write_optional_artifact("REVAC_Q96_PUBLIC_OUT", public_wire,
                                       ACC96_LABRADOR_PUBLIC_BYTES) &&
                   write_optional_artifact("REVAC_Q96_PROVER_OUT", prover_wire,
                                           ACC96_LABRADOR_PROVER_BYTES) &&
                   write_optional_artifact("REVAC_Q96_LINK_OUT", link_wire,
                                           REVAC_Q96_LINK_WIRE_BYTES) &&
                   write_optional_artifact("REVAC_Q96_INSTANCE_OUT",
                                           instance_wire,
                                           REVAC_Q96_INSTANCE_WIRE_BYTES),
               "写出 q96 Labrador 测试 artifact 失败")) {
    goto cleanup;
  }
  {
    const acc96_coeff_t saved = uav2.acc96_witness->left.coeff[0];
    const uint64_t lane_moduli[ACC96_LIMBS] = {
        ACC96_P0, ACC96_P1, ACC96_P2};
    for (size_t lane = 0; lane < ACC96_LIMBS; lane++)
      uav2.acc96_witness->left.coeff[0].limb[lane] =
          (saved.limb[lane] + 1U) % lane_moduli[lane];
    if (!require(!acc96_memver_exact_lift(
                     exact_lift, NULL, public_state.acc96_pp,
                     public_state.acc96_acc, uav2.handle,
                     uav2.acc96_witness, public_state.acc96_ctx),
                 "篡改 q96 witness 后精确整数 lift 仍被接受")) {
      uav2.acc96_witness->left.coeff[0] = saved;
      goto cleanup;
    }
    uav2.acc96_witness->left.coeff[0] = saved;
  }
  if (!require(memcmp(first_acc_digest, show_context.acc_digest,
                      sizeof(first_acc_digest)) != 0,
               "Delete 后 q96 累加器摘要没有变化")) {
    goto cleanup;
  }
  {
    q96_sync_endpoint_t endpoint = {&ta, uav2.handle};
    revac_sync_request_t request;
    acc96_relation_witness_t *scratch = calloc(1, sizeof(*scratch));
    uint8_t *sync_wire = NULL;
    uint32_t response_seq = 0;

    if (!require(scratch != NULL &&
                     revac_sync_request_create(
                         &request, uav2.handle, uav2.last_seq,
                         &public_state) &&
                     !revac_ta_export_memwitsync(
                         &ta, uav1.handle, &request,
                         &sync_wire, &sync_wire_len) &&
                     revac_ta_export_memwitsync(
                         &ta, uav2.handle, &request,
                         &sync_wire, &sync_wire_len),
                 "TA q96 MemWitSync 生成失败")) {
      free(sync_wire);
      free(scratch);
      goto cleanup;
    }
    sync_wire[REVAC_Q96_SYNC_WIRE_HEADER_BYTES] ^= UINT8_C(1);
    if (!require(!revac_q96_sync_response_unpack(
                     scratch, &response_seq, sync_wire, sync_wire_len,
                     &request, &public_state),
                 "篡改 q96 MemWitSync 签名未被拒绝")) {
      free(sync_wire);
      free(scratch);
      goto cleanup;
    }
    sync_wire[REVAC_Q96_SYNC_WIRE_HEADER_BYTES] ^= UINT8_C(1);
    {
      revac_sync_request_t replayed = request;
      replayed.nonce[0] ^= UINT8_C(1);
      if (!require(!revac_q96_sync_response_unpack(
                       scratch, &response_seq, sync_wire, sync_wire_len,
                       &replayed, &public_state),
                   "不同 nonce 重放 q96 MemWitSync 未被拒绝")) {
        free(sync_wire);
        free(scratch);
        goto cleanup;
      }
    }
    free(sync_wire);
    free(scratch);

    /* Corrupt the local copy, then recover only through the authenticated
       TA endpoint.  The decoder verifies into scratch space and installs the
       witness atomically after clear q96 MemVer succeeds. */
    uav2.acc96_witness->left.coeff[0].limb[0] ^= UINT64_C(1);
    if (!require(!revac_uav_memver_q96(&uav2, &public_state),
                 "q96 同步前的损坏 witness 意外通过") ||
        !require(revac_uav_recover_updates(
                     &uav2, &public_state, fetch_q96_sync_from_ta,
                     &endpoint),
                 "UAV q96 远程 MemWitSync 恢复失败") ||
        !require(revac_uav_memver_q96(&uav2, &public_state),
                 "q96 MemWitSync 安装后 MemVer 失败")) {
      goto cleanup;
    }
  }
  ok = 1;
  printf("q88 候选 D 实体接入测试通过（兼容 acc96 API）："
         "同一 TA/UAV/GCS API 已覆盖 Setup、注册、"
         "Issue、明文 MemVer、Delete 与 MemWitUp。\n");
  printf("q88 Issue 已使用区间认证扰动 + 精确有理 Klein SampleLeft；"
         "q88 LaBRADOR 关闭时 fail-closed、开启时执行完整关系。\n");
  printf("真实 witness 已通过 radix-8191 window 编码门：新签发最大系数=%llu，"
         "一次更新后最大系数=%llu，上限=%lld。\n",
         (unsigned long long)fresh_witness_max_abs,
         (unsigned long long)updated_witness_max_abs,
         (long long)ACC96_WINDOW_COEFF_BOUND);
  printf("三条 CRT lane 共用同一份有符号 witness；精确整数 lift 最大为 %u bit。\n",
         exact_lift_bits);
  printf("q88 Labrador 本地桥接：公开静态数据 %.2f MiB，"
         "预分解 prover payload %.2f MiB（仅限 UAV 本机 Unix socket），"
         "GCS 每次公开实例 %.2f KiB。\n",
         ACC96_LABRADOR_PUBLIC_BYTES / (1024.0 * 1024.0),
         ACC96_LABRADOR_PROVER_BYTES / (1024.0 * 1024.0),
         REVAC_Q96_INSTANCE_WIRE_BYTES / 1024.0);
  printf("q88 跨设备撤销广播：撤销前生成并验签 T4 prefetch %.2f KiB"
         "（TA %.3f 秒）；"
         "实时签名 manifest %.2f KiB；%s 缓存命中更新 %.3f 秒；"
         "已验证缓存缺失时 LOCAL 回退及篡改拒绝。\n",
         prefetch_wire_len / 1024.0, q96_prefetch_seconds,
         manifest_wire_len / 1024.0,
         q96_update_used_tail4 ? "TAIL4" : "LOCAL 回退",
         q96_update_seconds);
  printf("q88 选择性 DELTA6：TA 融合 Delete/六层构造并保留兼容 T4 "
         "%.3f 秒；认证 manifest %.2f KiB；该 UAV %s，"
         "接收 %zu 层/%.2f KiB，应用 %zu 层，更新 %.3f 秒；"
         "缺层、错序、跨序列、篡改层和篡改句柄均已拒绝。\n",
         q96_delete_seconds, delta6_manifest_wire_len / 1024.0,
         delta6_used_local ? "LOCAL 回退" : "DELTA6 命中",
         delta6_needed_layers, delta6_download_bytes / 1024.0,
         delta6_applied_layers, q96_update_seconds);
  puts("q88 DELTA6 完整协议分支：已确定性验证第 0 层命中，以及首个差异"
       "位 >= 6 时零层下载并 LOCAL 更新。");
  printf("q88 远程 MemWitSync：TA 重新采样当前 witness，签名恢复 packet "
         "%.2f MiB；篡改拒绝、明文 MemVer 与原子安装均通过。\n",
         sync_wire_len / (1024.0 * 1024.0));

cleanup:
  free(link_wire);
  free(instance_wire);
  free(show_wire);
  free(update_wire);
  free(prefetch_wire);
  free(manifest_wire);
  free(delta6_manifest_wire);
  for (size_t i = 0; i < REVAC_Q96_DELTA_LAYERS; i++)
    free(delta6_layer_wire[i]);
  free(prover_wire);
  free(public_wire);
  free(exact_lift);
  free(expanded_pp);
  for (size_t i = 0; i < 2; i++) poly_q_vec_d_clear(link_randomness[i]);
  poly_q_vec_d_clear(link_commitment);
  revac_show_context_clear(&show_context);
  revac_show_proof_clear(&show);
  if (verify_guard_initialized) revac_gcs_verify_guard_clear(&verify_guard);
  revac_gcs_clear(&gcs);
  revac_uav_clear(&uav_local);
  revac_uav_clear(&uav2);
  revac_uav_clear(&uav1);
  revac_ta_clear(&ta);
  arith_teardown();
  return ok ? 0 : 1;
}
