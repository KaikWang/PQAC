#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "acc_genpp.h"
#include "arith.h"
#include "random.h"

static int expect(int condition, const char *message) {
  if (!condition) {
    fprintf(stderr, "[失败] %s\n", message);
    return 0;
  }
  return 1;
}

static int witness_equal(const acc_genpp_wit_t *lhs,
                         const acc_genpp_wit_t *rhs) {
  if (!poly_q_vec_d_equal(lhs->left, rhs->left) ||
      !poly_q_vec_d_equal(lhs->v12, rhs->v12) ||
      !poly_q_equal(lhs->tag, rhs->tag)) {
    return 0;
  }
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    if (!poly_q_vec_d_equal(lhs->v2[i], rhs->v2[i])) {
      return 0;
    }
  }
  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
    for (size_t digit = 0; digit < ACC_GENPP_GADGET_DIM; digit++) {
      if (!poly_q_vec_d_equal(lhs->right[bit][digit],
                              rhs->right[bit][digit])) {
        return 0;
      }
    }
  }
  return 1;
}

static int test_scalar_decomposition(uint64_t *tested) {
  const coeff_q bound = (coeff_q)(PARAM_Q / 2);
  coeff_q digits[ACC_GENPP_GADGET_DIM];

  *tested = 0;
  for (coeff_q value = -bound; value <= bound; value++) {
    int64_t reconstructed = 0;

    if (!expect(acc_genpp_balanced_ternary_decompose(digits, value),
                "中心代表元超出 12 位平衡三进制覆盖范围")) {
      return 0;
    }
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      if (!expect(digits[j] >= -1 && digits[j] <= 1,
                  "平衡三进制 digit 不属于 {-1,0,1}")) {
        return 0;
      }
      reconstructed +=
          (int64_t)digits[j] * (int64_t)acc_genpp_gadget_power(j);
    }
    if (!expect(reconstructed == value, "平衡三进制重构不等于输入")) {
      return 0;
    }
    (*tested)++;
  }

  if (!expect(!acc_genpp_balanced_ternary_decompose(
                  digits, (coeff_q)ACC_GENPP_GADGET_MAX_CENTERED + 1),
              "覆盖上界之外的正整数应保留非零最高 carry") ||
      !expect(!acc_genpp_balanced_ternary_decompose(
                  digits, -(coeff_q)ACC_GENPP_GADGET_MAX_CENTERED - 1),
              "覆盖上界之外的负整数应保留非零最高 carry")) {
    return 0;
  }
  return 1;
}

static int test_indicator_and_delta(void) {
  const uint32_t revoked = UINT32_C(0xa5c39e71);
  acc_genpp_pp_t pp;
  acc_genpp_td_t td;
  acc_genpp_state_t state;
  acc_genpp_wit_t revoked_witness, local_delta, cached_delta, untouched;
  acc_genpp_update_t add_update, delete_update;
  acc_genpp_gmat_t captured[ACC_GENPP_DELTA4_LAYERS];
  acc_genpp_gmat_t *capture[ACC_GENPP_HANDLE_BITS] = {0};
  poly_q_vec_d delta4[ACC_GENPP_DELTA4_LAYERS][ACC_GENPP_GADGET_DIM];
  poly_q_vec_d accumulator, after_delete, accumulator_delta;
  int ok = 0;

  acc_genpp_pp_init(&pp);
  acc_genpp_td_init(&td);
  acc_genpp_state_init(&state);
  acc_genpp_wit_init(&revoked_witness);
  acc_genpp_wit_init(&local_delta);
  acc_genpp_wit_init(&cached_delta);
  acc_genpp_wit_init(&untouched);
  poly_q_vec_d_init(accumulator);
  poly_q_vec_d_init(after_delete);
  poly_q_vec_d_init(accumulator_delta);
  for (size_t i = 0; i < ACC_GENPP_DELTA4_LAYERS; i++) {
    acc_genpp_gmat_init(&captured[i]);
    capture[i + 1] = &captured[i];
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_vec_d_init(delta4[i][j]);
    }
  }

  acc_genpp_gen(&pp, &td, &state, accumulator);
  if (!expect(acc_genpp_add(after_delete, &revoked_witness, &add_update,
                            &state, &pp, &td, accumulator, revoked),
              "为待撤销句柄创建成员记录") ||
      !expect(acc_genpp_memver(&pp, after_delete, revoked,
                               &revoked_witness),
              "删除前成员见证通过 MemVer")) {
    goto cleanup;
  }
  poly_q_vec_d_set(accumulator, after_delete);
  if (!expect(acc_genpp_delete_trusted_capture_signed(
                  after_delete, &delete_update, &state, &pp, accumulator,
                  revoked, -1, capture),
              "带负随机符号的 Delete 与 T1..T4 捕获") ||
      !expect(delete_update.op == ACC_GENPP_OP_DEL &&
                  delete_update.handle == revoked &&
                  delete_update.sign == -1,
              "Delete 更新消息绑定撤销句柄与随机符号")) {
    goto cleanup;
  }
  poly_q_vec_d_sub(accumulator_delta, after_delete, accumulator);

  acc_genpp_delta4_from_tails(delta4, captured, revoked, pp.ua);

  for (unsigned bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
    const uint32_t survivor = revoked ^ (UINT32_C(1) << bit);

    acc_genpp_wit_clear(&local_delta);
    acc_genpp_wit_init(&local_delta);
    if (!expect(acc_genpp_memwitup_tails(
                    &local_delta, &pp, survivor, &delete_update,
                    NULL, NULL, 0, NULL),
                "LOCAL EvalFX 应用认证随机符号")) {
      goto cleanup;
    }
    if (!expect(acc_genpp_memver(&pp, accumulator_delta, survivor,
                                 &local_delta),
                "专门化 EvalFX 满足删除后的 MemVer 关系")) {
      goto cleanup;
    }

    if (bit < ACC_GENPP_DELTA4_LAYERS) {
      size_t applied_layers = 0;

      acc_genpp_wit_clear(&cached_delta);
      acc_genpp_wit_init(&cached_delta);
      if (!expect(acc_genpp_memwitup_delta4(
                      &cached_delta, survivor, &delete_update, delta4,
                      &applied_layers),
                  "DELTA4 覆盖首个不同位") ||
          !expect(applied_layers == (size_t)bit + 1,
                  "DELTA4 应用层数等于首个不同位加一") ||
          !expect(witness_equal(&cached_delta, &local_delta),
                  "DELTA4 与 LOCAL EvalFX 逐系数相等")) {
        goto cleanup;
      }
    } else {
      size_t applied_layers = 99;

      acc_genpp_wit_clear(&cached_delta);
      acc_genpp_wit_init(&cached_delta);
      acc_genpp_wit_clear(&untouched);
      acc_genpp_wit_init(&untouched);
      if (!expect(!acc_genpp_memwitup_delta4(
                      &cached_delta, survivor, &delete_update, delta4,
                      &applied_layers),
                  "DELTA4 未覆盖时必须拒绝并回退") ||
          !expect(applied_layers == 0,
                  "DELTA4 未覆盖时不得报告已应用层") ||
          !expect(witness_equal(&cached_delta, &untouched),
                  "DELTA4 未覆盖时不得修改见证")) {
        goto cleanup;
      }
    }
  }

  acc_genpp_wit_clear(&local_delta);
  acc_genpp_wit_init(&local_delta);
  if (!expect(acc_genpp_memwitup_tails(
                  &local_delta, &pp, revoked, &delete_update,
                  NULL, NULL, 0, NULL),
              "被撤销句柄的带符号 LOCAL EvalFX 可执行")) {
    goto cleanup;
  }
  if (!expect(!acc_genpp_memver(&pp, accumulator_delta, revoked,
                                &local_delta),
              "被撤销句柄不能用诚实 EvalFX 增量满足新累加器关系")) {
    goto cleanup;
  }

  ok = 1;

cleanup:
  for (size_t i = 0; i < ACC_GENPP_DELTA4_LAYERS; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_vec_d_clear(delta4[i][j]);
    }
    acc_genpp_gmat_clear(&captured[i]);
  }
  poly_q_vec_d_clear(accumulator_delta);
  poly_q_vec_d_clear(after_delete);
  poly_q_vec_d_clear(accumulator);
  acc_genpp_wit_clear(&untouched);
  acc_genpp_wit_clear(&cached_delta);
  acc_genpp_wit_clear(&local_delta);
  acc_genpp_wit_clear(&revoked_witness);
  acc_genpp_state_clear(&state);
  acc_genpp_td_clear(&td);
  acc_genpp_pp_clear(&pp);
  return ok;
}

int main(void) {
  uint64_t scalar_cases = 0;
  int ok;

  arith_setup();
  random_init();
  ok = test_scalar_decomposition(&scalar_cases) &&
       test_indicator_and_delta();
  arith_teardown();

  if (!ok) {
    return 1;
  }
  printf("balanced base-3 正确性回归通过。\n");
  printf("  穷举中心代表元：%" PRIu64 " 个（-%ld 到 %ld）\n",
         scalar_cases, (long)(PARAM_Q / 2), (long)(PARAM_Q / 2));
  printf("  每个输入均满足：12 digits、digit 属于 {-1,0,1}、整数精确重构。\n");
  printf("  EvalFX/MemVer：首个不同位 j=0..31 全部通过。\n");
  printf("  DELTA4/LOCAL：j=0..3 逐系数相等；j=4..31 安全拒绝且不修改。\n");
  printf("  撤销目标：同句柄诚实更新后 MemVer 拒绝。\n");
  return 0;
}
