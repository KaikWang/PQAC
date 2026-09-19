#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "arith.h"
#include "random.h"
#include "randombytes.h"
#include "revac_context.h"
#include "revac_gcs.h"
#include "revac_labrador_bridge.h"
#include "revac_sync.h"
#include "revac_ta.h"
#include "revac_uav.h"

#define REVAC_DEMO_SHOW_STEADY_REPETITIONS_DEFAULT 3u
#define REVAC_DEMO_SHOW_STEADY_REPETITIONS_MAX 50u

typedef struct {
  double min;
  double max;
  double total;
  double sum_squares;
  double samples[REVAC_DEMO_SHOW_STEADY_REPETITIONS_MAX];
  size_t count;
} revac_demo_timing_stats_t;

static double monotonic_seconds(void) {
  struct timespec now;

  (void)clock_gettime(CLOCK_MONOTONIC, &now);
  return (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
}

static double elapsed_seconds(double started) {
  return monotonic_seconds() - started;
}

static size_t show_steady_repetitions(void) {
  const char *value = getenv("REVAC_DEMO_SHOW_REPETITIONS");
  char *end = NULL;
  unsigned long parsed;

  if (value == NULL || value[0] == '\0') {
    return REVAC_DEMO_SHOW_STEADY_REPETITIONS_DEFAULT;
  }
  parsed = strtoul(value, &end, 10);
  if (end == value || *end != '\0' || parsed < 1 ||
      parsed > REVAC_DEMO_SHOW_STEADY_REPETITIONS_MAX) {
    return REVAC_DEMO_SHOW_STEADY_REPETITIONS_DEFAULT;
  }
  return (size_t)parsed;
}

static void timing_stats_add(revac_demo_timing_stats_t *stats, double value) {
  if (stats->count >= REVAC_DEMO_SHOW_STEADY_REPETITIONS_MAX) {
    abort();
  }
  if (stats->count == 0 || value < stats->min) stats->min = value;
  if (stats->count == 0 || value > stats->max) stats->max = value;
  stats->total += value;
  stats->sum_squares += value * value;
  stats->samples[stats->count++] = value;
}

static double timing_stats_mean(const revac_demo_timing_stats_t *stats) {
  return stats->count == 0 ? 0.0 : stats->total / (double)stats->count;
}

static double timing_stats_stddev(const revac_demo_timing_stats_t *stats) {
  double variance;

  if (stats->count < 2) return 0.0;
  variance = (stats->sum_squares -
              stats->total * stats->total / (double)stats->count) /
             (double)(stats->count - 1);
  return variance > 0.0 ? sqrt(variance) : 0.0;
}

static void print_timing_stats(const char *label,
                               const revac_demo_timing_stats_t *stats) {
  const double mean = timing_stats_mean(stats);
  const double stddev = timing_stats_stddev(stats);

  printf("    %-44s 平均 %8.3f ms，标准差 %7.3f ms，最小 %8.3f ms，最大 %8.3f ms（%zu 次）\n",
         label, mean * 1000.0, stddev * 1000.0, stats->min * 1000.0,
         stats->max * 1000.0, stats->count);
  printf("      原始样本（ms）：");
  for (size_t i = 0; i < stats->count; i++) {
    printf("%s%.3f", i == 0 ? "" : ", ", stats->samples[i] * 1000.0);
  }
  printf("\n");
}

static void print_section(const char *title) {
  printf("\n============================================================\n");
  printf("%s\n", title);
  printf("============================================================\n");
}

static int detailed_logs_enabled(void) {
  const char *value = getenv("REVAC_DEMO_VERBOSE");

  return value != NULL && value[0] != '\0' && strcmp(value, "0") != 0;
}

/* Keep the integration-test output short by default while preserving the
   component-level diagnostics behind REVAC_DEMO_VERBOSE=1.  stderr remains
   visible, and callers restore stdout before reporting a failed operation. */
static int quiet_stdout_begin(void) {
  int saved_stdout;
  int null_fd;

  (void)fflush(stdout);
  saved_stdout = dup(STDOUT_FILENO);
  if (saved_stdout < 0) {
    return -1;
  }
  null_fd = open("/dev/null", O_WRONLY);
  if (null_fd < 0 || dup2(null_fd, STDOUT_FILENO) < 0) {
    if (null_fd >= 0) {
      (void)close(null_fd);
    }
    (void)close(saved_stdout);
    return -1;
  }
  (void)close(null_fd);
  return saved_stdout;
}

static void quiet_stdout_end(int saved_stdout) {
  if (saved_stdout < 0) {
    return;
  }
  (void)fflush(stdout);
  (void)dup2(saved_stdout, STDOUT_FILENO);
  (void)close(saved_stdout);
}

static void diagnose_show_verify_failure(
    const revac_public_state_t *public_state,
    const revac_gcs_t *gcs,
    const revac_show_context_t *show_ctx,
    const uint8_t *wire,
    size_t wire_len) {
  const int signature_ok = revac_gcs_verify_show_signature_wire(
      public_state, gcs, show_ctx, wire, wire_len);
  const int accumulator_ok = revac_gcs_verify_show_accumulator_wire(
      public_state, gcs, show_ctx, wire, wire_len);

  /* Failure-only diagnostics: the normal benchmark path and its timings are
     unchanged.  Keeping the two component results visible prevents an
     occasional native-proof failure from being misdiagnosed as a scheduling
     or packet-codec regression. */
  printf("  [失败定位] 匿名签名组件：%s；累加器 Labrador 组件：%s\n",
         signature_ok ? "通过" : "失败",
         accumulator_ok ? "通过" : "失败");
}

static int expect_true(int cond, const char *msg) {
  if (!cond) {
    printf("%s\n", msg);
    return 0;
  }
  return 1;
}

static int expect_false(int cond, const char *msg) {
  return expect_true(!cond, msg);
}

/* In a cross-device deployment this callback is replaced by the network
   receiver.  The demo maps it directly to TA storage. */
static int fetch_context_from_ta(void *opaque,
                                 uint32_t seq,
                                 revac_context_profile_t profile,
                                 uint8_t **wire,
                                 size_t *wire_len) {
  return revac_ta_export_context((const revac_ta_t *)opaque, seq, profile,
                                 wire, wire_len);
}

typedef struct {
  const revac_ta_t *ta;
  const revac_public_state_t *public_state;
  uint32_t authenticated_handle;
  int tamper_checks_passed;
} revac_sync_demo_endpoint_t;

typedef struct {
  double setup_and_issuance_total;
  double phase1_setup;
  double phase2_registration;
  double phase3_issuance;
  double phase3b_uav_prove;
  double phase3b_guarded_parallel_verify;
  double phase3b_total;
  revac_demo_timing_stats_t show_steady_uav_prove;
  revac_demo_timing_stats_t show_steady_serial_verify;
  revac_demo_timing_stats_t show_steady_parallel_verify;
  double show_benchmark_extra;
  double phase4_revoke1;
  double phase4_online_update;
  double phase4_state_checks;
  double phase4_offline_recovery;
  double phase4_revoke2;
  double phase4_remote_sync;
  double phase4_total;
  double workflow_total;
  size_t show_packet_bytes;
  size_t online_context_bytes;
  size_t offline_context_bytes;
  size_t second_context_bytes;
  size_t sync_response_bytes;
  size_t online_local_composites;
  size_t offline_local_composites;
  revac_context_profile_t online_profile;
  revac_context_profile_t offline_profile;
  int labrador_online;
} revac_demo_report_t;

/* The correctness run immediately before this helper is the cold-start
 * sample.  Repeating with fresh nonces exercises the same production API
 * after the local worker/public-parameter caches are warm. */
static int benchmark_show_steady_state(
    revac_demo_report_t *report,
    const revac_public_state_t *public_state,
    const revac_gcs_t *gcs,
    const revac_show_context_t *show_ctx,
    revac_gcs_verify_guard_t *verify_guard,
    revac_uav_t *uav,
    const uint8_t msg[PARAM_M * PARAM_N / 8],
    int verbose) {
  const double benchmark_started = monotonic_seconds();
  const size_t repetitions = show_steady_repetitions();

  for (size_t run = 0; run < repetitions; run++) {
    uint8_t nonce[REVAC_SHOW_NONCE_BYTES];
    uint8_t *wire = NULL;
    size_t wire_len = 0;
    double started;
    int saved_stdout;
    int ok;

    randombytes(nonce, sizeof(nonce));
    started = monotonic_seconds();
    saved_stdout = verbose ? -1 : quiet_stdout_begin();
    ok = revac_uav_show_prove_wire(&wire, &wire_len, public_state, uav,
                                   msg, nonce);
    quiet_stdout_end(saved_stdout);
    if (!ok) {
      free(wire);
      printf("稳定态 show 基准第 %zu 次生成失败。\n", run + 1);
      return 0;
    }
    timing_stats_add(&report->show_steady_uav_prove,
                     elapsed_seconds(started));

    started = monotonic_seconds();
    saved_stdout = verbose ? -1 : quiet_stdout_begin();
    ok = revac_gcs_verify_show_wire(public_state, gcs, show_ctx,
                                    wire, wire_len);
    quiet_stdout_end(saved_stdout);
    if (report->labrador_online ? !ok : ok) {
      printf("稳定态 show 基准第 %zu 次串行完整验证失败。\n", run + 1);
      diagnose_show_verify_failure(public_state, gcs, show_ctx,
                                   wire, wire_len);
      free(wire);
      return 0;
    }
    timing_stats_add(&report->show_steady_serial_verify,
                     elapsed_seconds(started));

    started = monotonic_seconds();
    saved_stdout = verbose ? -1 : quiet_stdout_begin();
    ok = revac_gcs_verify_show_wire_guarded_parallel(
        public_state, gcs, show_ctx, verify_guard, wire, wire_len);
    quiet_stdout_end(saved_stdout);
    if (report->labrador_online ? !ok : ok) {
      printf("稳定态 show 基准第 %zu 次受保护并行验证结果不符合预期。\n",
             run + 1);
      diagnose_show_verify_failure(public_state, gcs, show_ctx,
                                   wire, wire_len);
      free(wire);
      return 0;
    }
    free(wire);
    timing_stats_add(&report->show_steady_parallel_verify,
                     elapsed_seconds(started));
  }
  report->show_benchmark_extra = elapsed_seconds(benchmark_started);
  return 1;
}

static void print_byte_count(const char *label, size_t bytes) {
  printf("  %-34s %10zu 字节  (%8.2f KiB)\n",
         label, bytes, bytes / 1024.0);
}

static void print_report(const revac_demo_report_t *report,
                         double runtime_initialization,
                         double runtime_teardown,
                         double program_total) {
  print_section("第二部分：大小统计（序列化证明与网络载荷）");
  printf("  show 模式                         %s\n",
         report->labrador_online
             ? "匿名签名 + zk=True Labrador 累加器证明"
             : "仅匿名签名（离线门控测试）");
  print_byte_count("压缩后的可撤销 show packet", report->show_packet_bytes);
  print_byte_count("在线 UAV2 下载的 DELTA4 context",
                   report->online_context_bytes);
  print_byte_count("离线 UAV3 恢复用 DELTA4 context",
                   report->offline_context_bytes);
  print_byte_count("第二次撤销缓存的 DELTA4 context",
                   report->second_context_bytes);
  print_byte_count("远程 MemWitSync 响应", report->sync_response_bytes);
  printf("  在线 context profile              %s，本地复合 %zu 次\n",
         revac_context_profile_name(report->online_profile),
         report->online_local_composites);
  printf("  离线恢复 context profile           %s，本地复合 %zu 次\n",
         revac_context_profile_name(report->offline_profile),
         report->offline_local_composites);
  printf("  注：这里统计的是 wire/payload 大小，不是进程峰值内存。\n");

  print_section("第三部分：时间基准汇总");
  printf("  %-48s %10zu 个\n", "累加器 EvalF/EvalFX 并行 worker",
         acc_genpp_parallel_threads());
  printf("  %-48s %10.6f 秒\n", "运行时初始化", runtime_initialization);
  printf("  %-48s %10.6f 秒\n", "业务阶段 1：初始化、注册与凭证签发",
         report->setup_and_issuance_total);
  printf("    %-46s %10.6f 秒\n", "TA 初始化 + KeyGen + AccGen + BML",
         report->phase1_setup);
  printf("    %-46s %10.6f 秒\n", "UAV/GCS 注册",
         report->phase2_registration);
  printf("    %-46s %10.6f 秒\n", "签发 3 个凭证 + Add + 验证",
         report->phase3_issuance);
  printf("  %-48s %10.6f 秒\n", "业务阶段 2：匿名 show 生成与验证",
         report->phase3b_total);
  printf("    %-46s %10.6f 秒\n", "UAV 完整 show（冷启动/正确性样本）",
         report->phase3b_uav_prove);
  printf("    %-46s %10.6f 秒\n", "GCS 预过滤/限流后并行完整验证（冷启动样本）",
         report->phase3b_guarded_parallel_verify);
  printf("    注：以下 %zu 次稳定态样本不计入业务流程总时间。\n",
         report->show_steady_uav_prove.count);
  print_timing_stats("UAV 完整 show（冷启动后稳定态）",
                     &report->show_steady_uav_prove);
  print_timing_stats("GCS 串行完整验证（稳定态对照）",
                     &report->show_steady_serial_verify);
  print_timing_stats("GCS 预过滤/限流后并行完整验证（稳定态）",
                     &report->show_steady_parallel_verify);
  if (timing_stats_mean(&report->show_steady_parallel_verify) > 0.0) {
    printf("    GCS 并行验证加速比                         %.2fx\n",
           timing_stats_mean(&report->show_steady_serial_verify) /
               timing_stats_mean(&report->show_steady_parallel_verify));
  }
  printf("  %-48s %10.6f 秒\n", "业务阶段 3：撤销、更新与同步",
         report->phase4_total);
  printf("    %-46s %10.6f 秒\n", "第一次撤销 + 捕获 context",
         report->phase4_revoke1);
  printf("    %-46s %10.6f 秒\n", "在线 UAV2 context 更新",
         report->phase4_online_update);
  printf("    %-46s %10.6f 秒\n", "撤销/最新/过期状态检查",
         report->phase4_state_checks);
  printf("    %-46s %10.6f 秒\n", "离线 UAV3 BML/context 恢复",
         report->phase4_offline_recovery);
  printf("    %-46s %10.6f 秒\n", "第二次撤销 + context 缓存检查",
         report->phase4_revoke2);
  printf("    %-46s %10.6f 秒\n", "远程 MemWitSync + 篡改检查",
         report->phase4_remote_sync);
  printf("  %-48s %10.6f 秒\n", "完整业务流程（三个业务阶段）",
         report->workflow_total);
  printf("  %-48s %10.6f 秒\n", "运行时清理", runtime_teardown);
  printf("  %-48s %10.6f 秒\n", "程序总计", program_total);
}

/* Cross-device MemWitSync endpoint.  Only this TA-side callback can access
   the accumulator trapdoor.  authenticated_handle models the identity bound
   to an mTLS/QUIC session; the UAV receives only a signed wire response. */
static int fetch_sync_from_ta(void *opaque,
                              const revac_sync_request_t *request,
                              uint8_t **wire,
                              size_t *wire_len) {
  revac_sync_demo_endpoint_t *endpoint = opaque;
  acc_genpp_wit_t scratch;
  revac_sync_request_t wrong_request;
  uint8_t *tampered = NULL;
  uint32_t response_seq = 0;
  int rejected;

  if (!revac_ta_export_memwitsync(endpoint->ta,
                                  endpoint->authenticated_handle,
                                  request, wire, wire_len)) {
    return 0;
  }
  tampered = malloc(*wire_len);
  if (tampered == NULL) {
    free(*wire);
    *wire = NULL;
    *wire_len = 0;
    return 0;
  }
  acc_genpp_wit_init(&scratch);

  memcpy(tampered, *wire, *wire_len);
  tampered[REVAC_SYNC_WIRE_HEADER_BYTES] ^= UINT8_C(1);
  rejected = !revac_sync_response_unpack(
      &scratch, &response_seq, tampered, *wire_len, request,
      endpoint->public_state);

  memcpy(tampered, *wire, *wire_len);
  tampered[*wire_len - 1] ^= UINT8_C(1);
  rejected = rejected && !revac_sync_response_unpack(
      &scratch, &response_seq, tampered, *wire_len, request,
      endpoint->public_state);

  wrong_request = *request;
  wrong_request.nonce[0] ^= UINT8_C(1);
  rejected = rejected && !revac_sync_response_unpack(
      &scratch, &response_seq, *wire, *wire_len, &wrong_request,
      endpoint->public_state);

  acc_genpp_wit_clear(&scratch);
  free(tampered);
  if (!rejected) {
    free(*wire);
    *wire = NULL;
    *wire_len = 0;
    return 0;
  }
  endpoint->tamper_checks_passed = 1;
  return 1;
}

static int revac_plan_process_demo(revac_demo_report_t *report) {
  int ok = 1;
  int operation_ok;
  int saved_stdout;
  const int verbose = detailed_logs_enabled();
  double phase_started;
  double major_stage_started;
  double step_started;
  double workflow_started;
  revac_ta_t ta;
  revac_public_state_t public_state;
  revac_sync_demo_endpoint_t sync_endpoint;
  revac_uav_t uav1, uav2, uav3;
  revac_gcs_t gcs;
  revac_gcs_verify_guard_t verify_guard;
  revac_show_context_t show1_ctx;
  acc_genpp_update_t revoke_uav1, revoke_uav2;
  uint8_t *show1_wire = NULL;
  uint8_t *context_probe = NULL;
  size_t show1_wire_len = 0;
  size_t context_probe_len = 0;
  uint8_t nonce1[REVAC_SHOW_NONCE_BYTES];
  uint8_t msg1[PARAM_M * PARAM_N / 8];
  uint8_t msg2[PARAM_M * PARAM_N / 8];
  uint8_t msg3[PARAM_M * PARAM_N / 8];

  memset(report, 0, sizeof(*report));
  memset(&verify_guard, 0, sizeof(verify_guard));
  report->online_profile = REVAC_CONTEXT_PROFILE_LOCAL;
  report->offline_profile = REVAC_CONTEXT_PROFILE_LOCAL;
  report->labrador_online = revac_labrador_online_enabled();

  print_section("第一部分：正确性测试");
  printf("测试对象：TA、3 个 UAV、GCS、匿名 show、撤销更新和远程同步\n");
  workflow_started = monotonic_seconds();

  revac_ta_init(&ta);
  revac_uav_init(&uav1);
  revac_uav_init(&uav2);
  revac_uav_init(&uav3);
  revac_gcs_init(&gcs);
  if (!revac_gcs_verify_guard_init(&verify_guard, 32.0, 32.0)) {
    printf("GCS 验证限流器初始化失败。\n");
    ok = 0;
    goto cleanup;
  }
  revac_show_context_init(&show1_ctx);

  printf("\n[业务阶段 1/3] 初始化、注册与凭证签发\n");
  major_stage_started = monotonic_seconds();
  printf("  - TA 初始化、匿名凭证密钥生成、AccGen、BML 初始化\n");
  phase_started = monotonic_seconds();
  revac_ta_setup_q19(&ta);
  revac_ta_public_state(&public_state, &ta);
  memset(&sync_endpoint, 0, sizeof(sync_endpoint));
  sync_endpoint.ta = &ta;
  sync_endpoint.public_state = &public_state;
  if (!expect_true(ta.acc_st.ctr == 0, "TA 累加器计数器没有正确初始化。") ||
      !expect_true(ta.acc_st.bml_len == 0, "TA 的 BML 初始状态不是空的。") ||
      !expect_true(gcs.cert_valid == 0, "注册前 GCS 证书应处于无效状态。")) {
    ok = 0;
    goto cleanup;
  }
  report->phase1_setup = elapsed_seconds(phase_started);
  printf("  [通过] TA、累加器和 BML 初始状态正确\n");

  printf("  - UAV/GCS 密钥生成、注册和句柄派生\n");
  phase_started = monotonic_seconds();
  revac_uav_keygen(&uav1, UINT64_C(0x1000000000000001), ta.ta_pk.seed);
  revac_uav_keygen(&uav2, UINT64_C(0x1000000000000002), ta.ta_pk.seed);
  revac_uav_keygen(&uav3, UINT64_C(0x1000000000000003), ta.ta_pk.seed);
  revac_gcs_keygen(&gcs, UINT64_C(0x2000000000000001),
                   ta.ta_pk.seed, UINT64_C(2026071401));
  randombytes(msg1, sizeof(msg1));
  randombytes(msg2, sizeof(msg2));
  randombytes(msg3, sizeof(msg3));
  randombytes(nonce1, sizeof(nonce1));

  if (!expect_true(revac_ta_register_uav(&ta, &uav1), "UAV1 注册失败。") ||
      !expect_true(revac_ta_register_uav(&ta, &uav2), "UAV2 注册失败。") ||
      !expect_true(revac_ta_register_uav(&ta, &uav3), "UAV3 注册失败。") ||
      !expect_true(revac_ta_register_gcs(&ta, &gcs), "GCS 注册失败。") ||
      !expect_true(uav1.handle == revac_handle_from_identity(&uav1.pk, uav1.id),
                   "UAV1 撤销句柄不等于 H(PK_UAV, ID_UAV)。") ||
      !expect_true(uav2.handle == revac_handle_from_identity(&uav2.pk, uav2.id),
                   "UAV2 撤销句柄不等于 H(PK_UAV, ID_UAV)。") ||
      !expect_true(uav3.handle == revac_handle_from_identity(&uav3.pk, uav3.id),
                   "UAV3 撤销句柄不等于 H(PK_UAV, ID_UAV)。") ||
      !expect_true(gcs.registered && gcs.cert_valid == REVAC_GCS_CERT_VALID,
                   "TA 没有正确安装 GCS 证书状态。")) {
    ok = 0;
    goto cleanup;
  }
  report->phase2_registration = elapsed_seconds(phase_started);
  printf("  [通过] 注册、GCS 证书和撤销句柄派生正确\n");

  printf("  - 匿名凭证签发、累加器 Add 和凭证验证\n");
  phase_started = monotonic_seconds();
  if (!expect_true(revac_ta_issue_uav(&ta, &uav1, msg1), "UAV1 凭证签发失败。") ||
      !expect_true(revac_ta_issue_uav(&ta, &uav2, msg2), "UAV2 凭证签发失败。") ||
      !expect_true(revac_ta_issue_uav(&ta, &uav3, msg3), "UAV3 凭证签发失败。") ||
      !expect_true(revac_gcs_verify_uav(&public_state, &gcs, &uav1, msg1),
                   "撤销前 GCS 验证 UAV1 失败。") ||
      !expect_true(revac_gcs_verify_uav(&public_state, &gcs, &uav2, msg2),
                   "撤销前 GCS 验证 UAV2 失败。") ||
      !expect_true(revac_gcs_verify_uav(&public_state, &gcs, &uav3, msg3),
                   "撤销前 GCS 验证 UAV3 失败。")) {
    ok = 0;
    goto cleanup;
  }
  gcs.cert_valid = 0;
  if (!expect_false(revac_gcs_verify_uav(&public_state, &gcs, &uav1, msg1),
                    "GCS 错误接受了无效 GCS 证书下的 UAV 凭证。")) {
    ok = 0;
    goto cleanup;
  }
  gcs.cert_valid = REVAC_GCS_CERT_VALID;
  report->phase3_issuance = elapsed_seconds(phase_started);
  report->setup_and_issuance_total = elapsed_seconds(major_stage_started);
  printf("  [通过] 签发、签名、MemVer 和 GCS 证书门控正确\n");

  printf("\n[业务阶段 2/3] 匿名 show proof 与累加器 ZK 门控\n");
  if (verbose) {
    printf("---------------- show proof 内部详细日志开始 ----------------\n");
  } else {
    printf("  - 正在生成并验证 show proof（内部逐轮日志已隐藏）\n");
  }
  phase_started = monotonic_seconds();
  revac_show_context_from_public(&show1_ctx, &public_state, uav1.pk.seed);
  step_started = monotonic_seconds();
  saved_stdout = verbose ? -1 : quiet_stdout_begin();
  operation_ok = revac_uav_show_prove_wire(&show1_wire, &show1_wire_len,
                                           &public_state, &uav1, msg1, nonce1);
  quiet_stdout_end(saved_stdout);
  if (!expect_true(operation_ok,
                   "UAV1 生成压缩的可撤销 show packet 失败。")) {
    ok = 0;
    goto cleanup;
  }
  report->phase3b_uav_prove = elapsed_seconds(step_started);
  report->show_packet_bytes = show1_wire_len;
  step_started = monotonic_seconds();
  saved_stdout = verbose ? -1 : quiet_stdout_begin();
  operation_ok = revac_gcs_verify_show_wire_guarded_parallel(
      &public_state, &gcs, &show1_ctx, &verify_guard,
      show1_wire, show1_wire_len);
  quiet_stdout_end(saved_stdout);
  if (report->labrador_online ? !operation_ok : operation_ok) {
    printf("GCS 受保护并行完整验证结果不符合预期。\n");
    diagnose_show_verify_failure(&public_state, &gcs, &show1_ctx,
                                 show1_wire, show1_wire_len);
    ok = 0;
    goto cleanup;
  }
  report->phase3b_guarded_parallel_verify = elapsed_seconds(step_started);
  /* The same nonce must be rejected before either expensive verifier starts. */
  if (!expect_false(revac_gcs_verify_show_wire_guarded_parallel(
                        &public_state, &gcs, &show1_ctx, &verify_guard,
                        show1_wire, show1_wire_len),
                    "GCS 防重放门控错误接受了重复 show nonce。")) {
    ok = 0;
    goto cleanup;
  }
  report->phase3b_total = elapsed_seconds(phase_started);
  if (!benchmark_show_steady_state(report, &public_state, &gcs, &show1_ctx,
                                   &verify_guard, &uav1, msg1, verbose)) {
    ok = 0;
    goto cleanup;
  }
  if (verbose) {
    printf("---------------- show proof 内部详细日志结束 ----------------\n");
  }
  if (report->labrador_online) {
    printf("  [通过] 匿名签名 proof、zk=True Labrador MemVer proof 及同一句柄共享承诺均验证成功\n");
  } else {
    printf("  [通过] 匿名签名 proof 有效，缺少累加器 proof 的完整请求被正确拒绝\n");
  }
  printf("  [通过] 冷启动 show 与冷启动后 %zu 次稳定态 show 均完成验证\n",
         report->show_steady_uav_prove.count);

  printf("\n[业务阶段 3/3] 撤销、DELTA4 广播、LOCAL 回退与远程同步\n");
  phase_started = monotonic_seconds();
  step_started = monotonic_seconds();
  if (!expect_true(revac_ta_revoke_uav_context(&ta, uav1.handle,
                                               &revoke_uav1),
                   "撤销 UAV1 失败。") ||
      !expect_true(revoke_uav1.seq == 1 && revoke_uav1.handle == uav1.handle,
                   "UAV1 撤销广播中的 seq/handle 不符合预期。") ||
      !expect_true(ta.acc_st.bml_len == 1, "BML 没有保存最新的撤销广播。")) {
    ok = 0;
    goto cleanup;
  }
  revac_ta_public_state(&public_state, &ta);
  sync_endpoint.authenticated_handle = uav3.handle;
  report->phase4_revoke1 = elapsed_seconds(step_started);

  step_started = monotonic_seconds();
  {
    const double goodput_mbps = 20.0;
    const double recurrence_ms = 158.6093566;
    revac_context_profile_t expected_profile =
        revac_context_profile_select(goodput_mbps, recurrence_ms);
    revac_context_profile_t used_profile = REVAC_CONTEXT_PROFILE_LOCAL;
    size_t download_bytes = 0;
    size_t local_composites = 0;

    if (!expect_true(revac_uav_apply_update_auto(
                         &uav2, &public_state, &revoke_uav1, goodput_mbps,
                         recurrence_ms, fetch_context_from_ta, &ta,
                         &used_profile, &local_composites, &download_bytes),
                     "在线 UAV2 使用 context 更新失败。") ||
        !expect_true(used_profile == expected_profile ||
                         (expected_profile == REVAC_CONTEXT_PROFILE_DELTA4 &&
                          used_profile == REVAC_CONTEXT_PROFILE_LOCAL),
                     "在线 UAV2 没有使用 DELTA4 或安全 LOCAL 回退。") ||
        !expect_true(download_bytes >=
                         revac_context_profile_wire_bytes(
                             REVAC_CONTEXT_PROFILE_DELTA4),
                     "在线 UAV2 下载的 DELTA4 字节数不正确。")) {
      ok = 0;
      goto cleanup;
    }
    report->online_profile = used_profile;
    report->online_context_bytes = download_bytes;
    report->online_local_composites = local_composites;
  }
  report->phase4_online_update = elapsed_seconds(step_started);
  step_started = monotonic_seconds();
  if (!expect_false(revac_gcs_verify_uav(&public_state, &gcs, &uav1, msg1),
                    "已撤销的 UAV1 仍然通过验证。") ||
      !expect_true(revac_gcs_verify_uav(&public_state, &gcs, &uav2, msg2),
                   "在线 UAV2 执行 MemWitUp 后验证失败。") ||
      !expect_false(revac_gcs_verify_uav(&public_state, &gcs, &uav3, msg3),
                    "持有过期见证的离线 UAV3 被错误验证通过。")) {
    ok = 0;
    goto cleanup;
  }
  report->phase4_state_checks = elapsed_seconds(step_started);

  step_started = monotonic_seconds();
  {
    const double goodput_mbps = 5.0;
    const double recurrence_ms = 158.6093566;
    revac_context_profile_t expected_profile =
        revac_context_profile_select(goodput_mbps, recurrence_ms);
    size_t download_bytes = 0;
    size_t local_composites = 0;

    if (!expect_true(revac_uav_recover_updates_auto(
                         &uav3, &public_state, goodput_mbps, recurrence_ms,
                         fetch_context_from_ta, &ta,
                         fetch_sync_from_ta, &sync_endpoint, &download_bytes,
                         &local_composites),
                     "离线 UAV3 使用 context 支持的 BML 恢复失败。") ||
        !expect_true(download_bytes >=
                         revac_context_profile_wire_bytes(
                             REVAC_CONTEXT_PROFILE_DELTA4),
                     "UAV3 的 BML DELTA4 字节数不正确。") ||
        !expect_true(revac_gcs_verify_uav(&public_state, &gcs, &uav3, msg3),
                     "UAV3 完成 BML 恢复后，GCS 验证失败。") ||
        !expect_true(uav3.last_seq == revoke_uav1.seq,
                     "BML 恢复后 UAV3 的 last_seq 没有更新。")) {
      ok = 0;
      goto cleanup;
    }
    report->offline_profile =
        download_bytes > 0 && local_composites == 0
            ? REVAC_CONTEXT_PROFILE_DELTA4
            : REVAC_CONTEXT_PROFILE_LOCAL;
    report->offline_context_bytes = download_bytes;
    report->offline_local_composites = local_composites;
  }
  report->phase4_offline_recovery = elapsed_seconds(step_started);
  printf("  [通过] 已撤销/最新/过期状态、在线更新和离线恢复均正确\n");

  step_started = monotonic_seconds();
  if (!expect_true(revac_ta_revoke_uav_context(&ta, uav2.handle,
                                               &revoke_uav2),
                   "撤销 UAV2 失败。") ||
      !expect_true(revac_ta_export_context(
                       &ta, revoke_uav2.seq, REVAC_CONTEXT_PROFILE_DELTA4,
                       &context_probe, &context_probe_len),
                   "TA 没有缓存第二次撤销对应的 context。") ||
      !expect_true(context_probe_len >= revac_context_profile_wire_bytes(
                                             REVAC_CONTEXT_PROFILE_DELTA4),
                   "第二次撤销的 context 大小不正确。") ||
      !expect_false(revac_gcs_verify_uav(&public_state, &gcs, &uav2, msg2),
                    "已撤销的 UAV2 仍然通过验证。") ||
      !expect_false(revac_gcs_verify_uav(&public_state, &gcs, &uav3, msg3),
                    "第二次撤销后，UAV3 的过期见证被错误验证通过。")) {
    ok = 0;
    goto cleanup;
  }
  revac_ta_public_state(&public_state, &ta);
  report->second_context_bytes = context_probe_len;
  free(context_probe);
  context_probe = NULL;
  context_probe_len = 0;
  report->phase4_revoke2 = elapsed_seconds(step_started);

  ta.acc_st.bml_len = 0;
  ta.acc_st.bml_start = 0;
  memset(ta.acc_st.bml, 0, sizeof(ta.acc_st.bml));
  revac_ta_public_state(&public_state, &ta);
  step_started = monotonic_seconds();
  {
    revac_sync_request_t unauthorized_request;
    uint8_t *unauthorized_wire = NULL;
    size_t unauthorized_len = 0;
    size_t sync_bytes = 0;

    sync_endpoint.tamper_checks_passed = 0;
    if (!expect_true(revac_sync_request_create(
                         &unauthorized_request, uav3.handle, uav3.last_seq,
                         &public_state),
                     "无法构造绑定认证身份的同步请求。") ||
        !expect_false(revac_ta_export_memwitsync(
                          &ta, uav2.handle, &unauthorized_request,
                          &unauthorized_wire, &unauthorized_len),
                      "TA 错误接受了冒用其他已认证 UAV 的同步请求。") ||
        !expect_true(revac_uav_recover_updates_auto(
                         &uav3, &public_state, 20.0, 158.6093566,
                         fetch_context_from_ta, &ta,
                         fetch_sync_from_ta, &sync_endpoint,
                         &sync_bytes, NULL),
                     "BML 不可用时，远程 MemWitSync 失败。") ||
        !expect_true(sync_endpoint.tamper_checks_passed,
                     "没有执行已签名同步响应的篡改检查。") ||
        !expect_true(sync_bytes == REVAC_SYNC_WIRE_HEADER_BYTES +
                                       REVAC_CONTEXT_SIGNATURE_BYTES +
                                       REVAC_SYNC_WITNESS_BYTES,
                     "远程 MemWitSync 的字节数不正确。") ||
        !expect_true(revac_gcs_verify_uav(&public_state, &gcs, &uav3, msg3),
                     "远程 MemWitSync 后 GCS 验证失败。") ||
        !expect_true(uav3.last_seq == revoke_uav2.seq,
                     "远程 MemWitSync 后 UAV3 的 last_seq 没有更新。")) {
      free(unauthorized_wire);
      ok = 0;
      goto cleanup;
    }
    free(unauthorized_wire);
    report->sync_response_bytes = sync_bytes;
  }
  report->phase4_remote_sync = elapsed_seconds(step_started);

  if (!expect_false(revac_ta_register_uav(&ta, &uav1),
                    "TA 错误接受了重复或已撤销的 UAV 句柄。")) {
    ok = 0;
    goto cleanup;
  }
  printf("  - 远程 MemWitSync、安全绑定和篡改拒绝\n");
  printf("  [通过] 签名、payload、nonce、对端身份绑定和篡改检查均正确\n");
  printf("  [通过] 重复/已撤销句柄被拒绝\n");
  report->phase4_total = elapsed_seconds(phase_started);
  report->workflow_total = elapsed_seconds(workflow_started) -
                           report->show_benchmark_extra;

  printf("\n正确性结论：全部 3 个业务阶段测试通过。\n");

cleanup:
  free(context_probe);
  free(show1_wire);
  revac_show_context_clear(&show1_ctx);
  revac_gcs_verify_guard_clear(&verify_guard);
  revac_gcs_clear(&gcs);
  revac_uav_clear(&uav3);
  revac_uav_clear(&uav2);
  revac_uav_clear(&uav1);
  revac_ta_clear(&ta);
  return ok;
}

int main(void) {
  int ok;
  revac_demo_report_t report;
  double initialization_started;
  double program_started;
  double teardown_started;
  double initialization_elapsed;
  double teardown_elapsed;
  double program_elapsed;

  (void)setvbuf(stdout, NULL, _IOLBF, 0);
  printf("[参数] profile=%s\n", REVAC_PARAM_PROFILE_NAME);
  printf("[参数] profile_id=%s；发布状态=%s\n",
         REVAC_PARAM_PROFILE_ID_HEX,
         REVAC_PARAM_RELEASE_READY ? "frozen" : "candidate（不可作为最终安全声明）");
  program_started = monotonic_seconds();
  initialization_started = monotonic_seconds();
  arith_setup();
  random_init();
  initialization_elapsed = elapsed_seconds(initialization_started);
  ok = revac_plan_process_demo(&report);
  teardown_started = monotonic_seconds();
  arith_teardown();
  teardown_elapsed = elapsed_seconds(teardown_started);
  program_elapsed = elapsed_seconds(program_started);

  if (ok) {
    print_report(&report, initialization_elapsed, teardown_elapsed,
                 program_elapsed);
    print_section("最终结果");
    printf("全部正确性测试通过；大小与时间数据已分类汇总。\n");
    return 0;
  }
  print_section("最终结果");
  printf("测试失败；请从第一部分最后一个失败信息开始排查。\n");
  return 1;
}
