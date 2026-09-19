#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "arith.h"
#include "random.h"
#include "randombytes.h"
#include "revac_ta.h"
#include "revac_uav.h"

static double monotonic_seconds(void) {
  struct timespec now;
  (void)clock_gettime(CLOCK_MONOTONIC, &now);
  return (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
}

static int expect(int condition, const char *message) {
  if (!condition) {
    fprintf(stderr, "失败：%s\n", message);
    return 0;
  }
  return 1;
}

static int run_lift_modes(const char *state_name,
                          const revac_ta_t *ta,
                          const revac_uav_t *member) {
  int64_t fast[PARAM_D][PARAM_N];
  int64_t fallback[PARAM_D][PARAM_N];
  int64_t checked[PARAM_D][PARAM_N];
  int64_t automatic[PARAM_D][PARAM_N];
  acc_genpp_lift_info_t fast_info;
  acc_genpp_lift_info_t fallback_info;
  acc_genpp_lift_info_t checked_info;
  acc_genpp_lift_info_t automatic_info;
  double started;

  printf("\n[%s]\n", state_name);
  started = monotonic_seconds();
  if (!expect(acc_genpp_memver_exact_lift(
                  fast, &ta->acc_pp, ta->acc, member->handle,
                  &member->witness, ACC_GENPP_LIFT_FORCE_INT64, &fast_info),
              "int64 快速路径计算失败")) {
    return 0;
  }
  printf("int64 快速路径：%.6f 秒\n", monotonic_seconds() - started);

  started = monotonic_seconds();
  if (!expect(acc_genpp_memver_exact_lift(
                  fallback, &ta->acc_pp, ta->acc, member->handle,
                  &member->witness, ACC_GENPP_LIFT_FORCE_FLINT,
                  &fallback_info),
              "强制 FLINT 回退计算失败")) {
    return 0;
  }
  printf("强制 FLINT 回退：%.6f 秒\n", monotonic_seconds() - started);

  started = monotonic_seconds();
  if (!expect(acc_genpp_memver_exact_lift(
                  checked, &ta->acc_pp, ta->acc, member->handle,
                  &member->witness, ACC_GENPP_LIFT_CROSSCHECK, &checked_info),
              "双路径一致性检查失败")) {
    return 0;
  }
  printf("双路径一致性检查：%.6f 秒\n", monotonic_seconds() - started);

  started = monotonic_seconds();
  if (!expect(acc_genpp_memver_exact_lift(
                  automatic, &ta->acc_pp, ta->acc, member->handle,
                  &member->witness, ACC_GENPP_LIFT_AUTO, &automatic_info),
              "自动后端计算失败")) {
    return 0;
  }
  printf("自动选择路径：%.6f 秒\n", monotonic_seconds() - started);

  return expect(fast_info.int64_bound_ok,
                "当前参数的公开上界应允许 int64") &&
         expect(fast_info.backend == ACC_GENPP_LIFT_BACKEND_INT64,
                "强制 int64 未报告 int64 后端") &&
         expect(fallback_info.backend == ACC_GENPP_LIFT_BACKEND_FLINT &&
                    fallback_info.fallback_used,
                "强制回退未报告 FLINT 后端/回退标记") &&
         expect(checked_info.crosscheck_performed &&
                    checked_info.crosscheck_passed,
                "双路径检查标记不完整") &&
         expect(automatic_info.backend == ACC_GENPP_LIFT_BACKEND_INT64 &&
                    !automatic_info.fallback_used,
                "自动模式没有选择安全的 int64 快速路径") &&
         expect(memcmp(fast, fallback, sizeof(fast)) == 0,
                "int64 与强制 FLINT 输出不一致") &&
         expect(memcmp(fast, checked, sizeof(fast)) == 0,
                "int64 与交叉检查输出不一致") &&
         expect(memcmp(fast, automatic, sizeof(fast)) == 0,
                "int64 与自动模式输出不一致");
}

int main(void) {
  revac_ta_t ta;
  revac_public_state_t public_state;
  revac_uav_t revoked;
  revac_uav_t member;
  acc_genpp_update_t update;
  uint8_t revoked_message[PARAM_M * PARAM_N / 8];
  uint8_t member_message[PARAM_M * PARAM_N / 8];
  int ok = 0;

  arith_setup();
  random_init();
  revac_ta_init(&ta);
  revac_uav_init(&revoked);
  revac_uav_init(&member);

  revac_ta_setup_q19(&ta);
  revac_uav_keygen(&revoked, UINT64_C(0x8500000000000001), ta.ta_pk.seed);
  revac_uav_keygen(&member, UINT64_C(0x8500000000000002), ta.ta_pk.seed);
  randombytes(revoked_message, sizeof(revoked_message));
  randombytes(member_message, sizeof(member_message));
  if (!expect(revac_ta_register_uav(&ta, &revoked), "被撤销 UAV 注册失败") ||
      !expect(revac_ta_register_uav(&ta, &member), "存活 UAV 注册失败") ||
      !expect(revac_ta_issue_uav(&ta, &revoked, revoked_message),
              "被撤销 UAV 凭证签发失败") ||
      !expect(revac_ta_issue_uav(&ta, &member, member_message),
              "存活 UAV 凭证签发失败") ||
      !expect(acc_genpp_memver(&ta.acc_pp, ta.acc, member.handle,
                              &member.witness),
              "撤销前 MemVer 失败") ||
      !run_lift_modes("撤销前见证", &ta, &member)) {
    goto cleanup;
  }

  if (!expect(revac_ta_revoke_uav(&ta, revoked.handle, &update),
              "撤销操作失败")) {
    goto cleanup;
  }
  revac_ta_public_state(&public_state, &ta);
  revac_uav_apply_update(&member, &public_state, &update);
  if (!expect(acc_genpp_memver(&ta.acc_pp, ta.acc, member.handle,
                              &member.witness),
              "撤销后 MemVer 失败") ||
      !run_lift_modes("撤销后更新见证", &ta, &member)) {
    goto cleanup;
  }

  printf("\nexact-lift 快速路径、强制回退和双路径一致性测试全部通过。\n");
  ok = 1;

cleanup:
  revac_uav_clear(&member);
  revac_uav_clear(&revoked);
  revac_ta_clear(&ta);
  arith_teardown();
  return ok ? 0 : 1;
}
