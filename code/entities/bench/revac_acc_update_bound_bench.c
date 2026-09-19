#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "acc_genpp.h"
#include "arith.h"
#include "poly_q_sampling.h"
#include "random.h"
#include "revac_context.h"
#include "revac_labrador_bridge.h"
#include "revac_uav.h"

#define DEFAULT_DELETIONS 32u
#define MAX_DELETIONS 2000u
#define RESEARCH_RIGHT_PER_DELETION UINT64_C(8388608)
#define RESEARCH_RIGHT_MAX UINT64_C(1073741824)

static uint64_t splitmix64(uint64_t *state) {
  uint64_t z = (*state += UINT64_C(0x9e3779b97f4a7c15));
  z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
  return z ^ (z >> 31);
}

static size_t env_count(const char *name, size_t fallback) {
  const char *value = getenv(name);
  char *end = NULL;
  unsigned long parsed;

  if (value == NULL || value[0] == '\0') return fallback;
  parsed = strtoul(value, &end, 10);
  if (end == value || *end != '\0' || parsed == 0 || parsed > MAX_DELETIONS) {
    fprintf(stderr, "%s 必须是 [1, %u] 内的整数。\n", name,
            (unsigned)MAX_DELETIONS);
    exit(2);
  }
  return (size_t)parsed;
}

static size_t common_low_prefix(uint32_t lhs, uint32_t rhs) {
  uint32_t difference = lhs ^ rhs;
  size_t prefix = 0;

  if (difference == 0) return ACC_GENPP_HANDLE_BITS;
  while ((difference & 1u) == 0u) {
    prefix++;
    difference >>= 1;
  }
  return prefix;
}

static uint64_t right_norm2(const acc_genpp_wit_t *wit) {
  uint64_t value = 0;

  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
    for (size_t digit = 0; digit < ACC_GENPP_GADGET_DIM; digit++) {
      value += poly_q_vec_d_norm2(wit->right[bit][digit]);
    }
  }
  return value;
}

/* Exact grouping used by the production Labrador importer: concatenate the
 * 32 handle-indexed polynomials for each (gadget digit, module row). */
static uint64_t right_group_max_norm2(const acc_genpp_wit_t *wit,
                                      size_t group_width,
                                      size_t hot_prefix,
                                      int balanced) {
  uint64_t maximum = 0;
  size_t balanced_map[ACC_GENPP_HANDLE_BITS] = {0};
  const size_t balanced_groups = ACC_GENPP_HANDLE_BITS / group_width;

  if (balanced) {
    uint64_t weights[ACC_GENPP_HANDLE_BITS] = {0};
    size_t sizes[ACC_GENPP_HANDLE_BITS] = {0};
    for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
      size_t target = 0;
      while (sizes[target] == group_width) target++;
      for (size_t group = target + 1; group < balanced_groups; group++) {
        if (sizes[group] < group_width && weights[group] < weights[target]) {
          target = group;
        }
      }
      balanced_map[bit] = target;
      sizes[target]++;
      weights[target] += UINT64_C(1) << (ACC_GENPP_HANDLE_BITS - 1 - bit);
    }
  }

  for (size_t digit = 0; digit < ACC_GENPP_GADGET_DIM; digit++) {
    for (size_t row = 0; row < PARAM_D; row++) {
      if (balanced) {
        uint64_t values[ACC_GENPP_HANDLE_BITS] = {0};
        for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
          values[balanced_map[bit]] +=
              poly_q_sq_norm2(wit->right[bit][digit]->entries[row]);
        }
        for (size_t group = 0; group < balanced_groups; group++) {
          if (values[group] > maximum) maximum = values[group];
        }
        continue;
      }
      size_t first = 0;
      while (first < ACC_GENPP_HANDLE_BITS) {
        uint64_t value = 0;
        size_t last = (hot_prefix && first < hot_prefix)
                          ? first + 1
                          : (hot_prefix ? ACC_GENPP_HANDLE_BITS
                                        : first + group_width);
        if (last > ACC_GENPP_HANDLE_BITS) last = ACC_GENPP_HANDLE_BITS;
        for (size_t bit = first; bit < last; bit++) {
          value += poly_q_sq_norm2(wit->right[bit][digit]->entries[row]);
        }
        if (value > maximum) maximum = value;
        first = last;
      }
    }
  }
  return maximum;
}

static int experimental_balanced_groups(void) {
  const char *value = getenv("REVAC_RIGHT_GROUP_BALANCED");
  return value != NULL && value[0] != '\0' && strcmp(value, "0") != 0;
}

static size_t experimental_hot_prefix(void) {
  const char *value = getenv("REVAC_RIGHT_GROUP_HOT_PREFIX");
  char *end = NULL;
  unsigned long parsed;

  if (value == NULL || value[0] == '\0') return 0;
  parsed = strtoul(value, &end, 10);
  if (end == value || *end != '\0' || parsed >= ACC_GENPP_HANDLE_BITS) {
    fprintf(stderr,
            "REVAC_RIGHT_GROUP_HOT_PREFIX 必须在 [0, 句柄位数) 内。\n");
    exit(2);
  }
  return (size_t)parsed;
}

static size_t experimental_group_width(void) {
  const char *value = getenv("REVAC_RIGHT_GROUP_WIDTH");
  char *end = NULL;
  unsigned long parsed;

  if (value == NULL || value[0] == '\0') return ACC_GENPP_HANDLE_BITS;
  parsed = strtoul(value, &end, 10);
  if (end == value || *end != '\0' || parsed == 0 ||
      parsed > ACC_GENPP_HANDLE_BITS ||
      ACC_GENPP_HANDLE_BITS % parsed != 0) {
    fprintf(stderr, "REVAC_RIGHT_GROUP_WIDTH 必须是句柄位数的正因数。\n");
    exit(2);
  }
  return (size_t)parsed;
}

static int fresh_update_u_enabled(void) {
  const char *mode = getenv("REVAC_UPDATE_U_MODE");

  if (mode == NULL || mode[0] == '\0' || strcmp(mode, "fixed") == 0) {
    return 0;
  }
  if (strcmp(mode, "fresh") == 0) return 1;
  fprintf(stderr, "REVAC_UPDATE_U_MODE 必须为 fixed 或 fresh。\n");
  exit(2);
}

static void derive_update_u(poly_q_vec_d out, uint64_t *prng,
                            size_t deletion_index) {
  uint8_t seed[SEED_BYTES];

  for (size_t offset = 0; offset < sizeof(seed); offset += sizeof(uint64_t)) {
    const uint64_t word = splitmix64(prng);
    memcpy(seed + offset, &word, sizeof(word));
  }
  for (size_t i = 0; i < sizeof(deletion_index) && i < sizeof(seed); i++) {
    seed[i] ^= (uint8_t)(deletion_index >> (8 * i));
  }
  poly_q_vec_d_uniform(out, seed, DOMAIN_SEPARATOR_U);
}

static uint64_t refresh_norm2(const acc_genpp_wit_t *wit) {
  uint64_t value = 0;

  for (size_t digit = 0; digit < ACC_GENPP_GADGET_DIM; digit++) {
    value += poly_q_vec_d_norm2(wit->refresh_z[digit]);
  }
  return value;
}

static uint64_t v2_norm2(const acc_genpp_wit_t *wit) {
  uint64_t value = 0;

  for (size_t index = 0; index < ACC_GENPP_TRUNC_K_L; index++) {
    value += poly_q_vec_d_norm2(wit->v2[index]);
  }
  return value;
}

static int64_t centered(coeff_q value) {
  int64_t result = (int64_t)value;
  if (result > PARAM_Q / 2) result -= PARAM_Q;
  return result;
}

static uint64_t right_max_abs(const acc_genpp_wit_t *wit) {
  uint64_t maximum = 0;

  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
    for (size_t digit = 0; digit < ACC_GENPP_GADGET_DIM; digit++) {
      for (size_t row = 0; row < PARAM_D; row++) {
        for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
          int64_t value = centered(poly_q_get_coeff(
              wit->right[bit][digit]->entries[row], coeff));
          uint64_t magnitude = (uint64_t)(value < 0 ? -value : value);
          if (magnitude > maximum) maximum = magnitude;
        }
      }
    }
  }
  return maximum;
}

static void add_right(acc_genpp_wit_t *accumulator,
                      const acc_genpp_wit_t *increment) {
  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
    for (size_t digit = 0; digit < ACC_GENPP_GADGET_DIM; digit++) {
      poly_q_vec_d_add(accumulator->right[bit][digit],
                       accumulator->right[bit][digit],
                       increment->right[bit][digit]);
    }
  }
}

static void negate_right(acc_genpp_wit_t *increment) {
  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
    for (size_t digit = 0; digit < ACC_GENPP_GADGET_DIM; digit++) {
      poly_q_vec_d_neg(increment->right[bit][digit],
                       increment->right[bit][digit]);
    }
  }
}

static int random_update_sign_enabled(void) {
  const char *mode = getenv("REVAC_UPDATE_SIGN_MODE");

  if (mode == NULL || mode[0] == '\0' || strcmp(mode, "positive") == 0) {
    return 0;
  }
  if (strcmp(mode, "rademacher") == 0) return 1;
  fprintf(stderr,
          "REVAC_UPDATE_SIGN_MODE 必须为 positive 或 rademacher。\n");
  exit(2);
}

/* Advance the public accumulator with exactly the image of an EvalFX witness
 * increment under the surviving member's MemVer matrix.  This avoids a
 * second EvalF chain while keeping the benchmark's final relation genuine. */
static void apply_increment_to_accumulator(
    poly_q_vec_d accumulator, const acc_genpp_pp_t *pp,
    uint32_t member_handle, const acc_genpp_wit_t *increment) {
  poly_q_mat_d_d block, gadget;

  poly_q_mat_d_d_init(block);
  poly_q_mat_d_d_init(gadget);
  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
    const int member_bit = (int)((member_handle >> bit) & 1u);
    for (size_t digit = 0; digit < ACC_GENPP_GADGET_DIM; digit++) {
      poly_q_mat_d_d_set(block, pp->Ba.bit[bit].blocks[digit]);
      if (member_bit) {
        poly_q_mat_d_d_zero(gadget);
        for (size_t row = 0; row < PARAM_D; row++) {
          poly_q_set_coeff(gadget->rows[row]->entries[row], 0,
                           (coeff_q)acc_genpp_gadget_power(digit));
        }
        poly_q_mat_d_d_sub(block, block, gadget);
      }
      poly_q_mat_d_d_muladd_vec_d(
          accumulator, block, increment->right[bit][digit]);
    }
  }
  poly_q_mat_d_d_clear(gadget);
  poly_q_mat_d_d_clear(block);
}

static void witness_set(acc_genpp_wit_t *dst, const acc_genpp_wit_t *src) {
  poly_q_vec_d_set(dst->left, src->left);
  poly_q_vec_d_set(dst->v12, src->v12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_vec_d_set(dst->v2[i], src->v2[i]);
  }
  poly_q_set(dst->tag, src->tag);
  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
    for (size_t digit = 0; digit < ACC_GENPP_GADGET_DIM; digit++) {
      poly_q_vec_d_set(dst->right[bit][digit], src->right[bit][digit]);
    }
  }
  memcpy(dst->refresh_r, src->refresh_r, sizeof(dst->refresh_r));
  for (size_t digit = 0; digit < ACC_GENPP_GADGET_DIM; digit++) {
    poly_q_vec_d_set(dst->refresh_z[digit], src->refresh_z[digit]);
  }
}

static double monotonic_seconds(void) {
  struct timespec value;
  clock_gettime(CLOCK_MONOTONIC, &value);
  return (double)value.tv_sec + (double)value.tv_nsec / 1000000000.0;
}

static int online_proof_requested(void) {
  const char *value = getenv("REVAC_UPDATE_BOUND_LABRADOR");
  return value != NULL && value[0] != '\0' && strcmp(value, "0") != 0;
}

static int prove_final_witness(const acc_genpp_pp_t *pp,
                               const poly_q_vec_d accumulator,
                               const acc_genpp_wit_t *witness,
                               uint32_t member_handle,
                               uint32_t deletion_count) {
  static const uint8_t user_seed[SEED_BYTES] = {0x42};
  revac_public_state_t public_state = {0};
  revac_uav_t uav;
  revac_show_proof_t proof;
  revac_show_context_t context;
  poly_q_vec_d link_randomness[2];
  double started, prove_seconds, verify_seconds;
  int ok = 0;

  revac_uav_init(&uav);
  revac_uav_keygen(&uav, UINT64_C(0x500), user_seed);
  uav.handle = member_handle;
  uav.last_seq = deletion_count;
  uav.registered = 1;
  uav.issued = 1;
  witness_set(&uav.witness, witness);
  public_state.acc_pp = pp;
  public_state.acc = accumulator;
  public_state.acc_seq = deletion_count;
  revac_context_compute_pp_hash(public_state.context_pp_hash, pp);

  revac_show_proof_init(&proof);
  revac_show_context_init(&context);
  poly_q_vec_d_init(link_randomness[0]);
  poly_q_vec_d_init(link_randomness[1]);
  revac_same_handle_commitment_generate(
      proof.handle_commitment, link_randomness, &public_state, member_handle);
  revac_show_context_from_public(&context, &public_state, uav.pk.seed);

  started = monotonic_seconds();
  if (!revac_labrador_memver_prove(
          &proof, &public_state, &uav, link_randomness)) {
    fprintf(stderr, "500-window Labrador prove 失败。\n");
    goto cleanup;
  }
  prove_seconds = monotonic_seconds() - started;
  started = monotonic_seconds();
  if (!revac_labrador_memver_verify(&public_state, &context, &proof)) {
    fprintf(stderr, "500-window Labrador verify 失败。\n");
    goto cleanup;
  }
  verify_seconds = monotonic_seconds() - started;
  printf("500-window Labrador：prove=%.6f 秒，verify=%.6f 秒，"
         "payload=%zu 字节（%.2f KiB）。\n",
         prove_seconds, verify_seconds, proof.acc_zk_proof_len,
         (double)proof.acc_zk_proof_len / 1024.0);
  ok = 1;

cleanup:
  poly_q_vec_d_clear(link_randomness[1]);
  poly_q_vec_d_clear(link_randomness[0]);
  revac_show_context_clear(&context);
  revac_show_proof_clear(&proof);
  revac_uav_clear(&uav);
  return ok;
}

static int is_checkpoint(size_t completed, size_t target) {
  static const size_t checkpoints[] = {
      1, 2, 4, 8, 16, 32, 64, 128, 256, 500, 512, 1000, 2000};

  if (completed == target) return 1;
  for (size_t i = 0; i < sizeof(checkpoints) / sizeof(checkpoints[0]); i++) {
    if (completed == checkpoints[i]) return 1;
  }
  return 0;
}

int main(void) {
  const size_t deletions = env_count("REVAC_UPDATE_BOUND_DELETIONS",
                                     DEFAULT_DELETIONS);
  const char *proof_deletions_env =
      getenv("REVAC_UPDATE_BOUND_PROOF_DELETIONS");
  const size_t proof_deletions =
      (proof_deletions_env != NULL && proof_deletions_env[0] != '\0')
          ? env_count("REVAC_UPDATE_BOUND_PROOF_DELETIONS", deletions)
          : deletions;
  const uint32_t member_handle = UINT32_C(0x6a09e667);
  uint64_t prng = UINT64_C(0x243f6a8885a308d3);
  uint64_t update_prng = UINT64_C(0x13198a2e03707344);
  uint64_t sign_prng = UINT64_C(0xa4093822299f31d0);
  const int fresh_update_u = fresh_update_u_enabled();
  const int random_update_sign = random_update_sign_enabled();
  const size_t group_width = experimental_group_width();
  const size_t hot_prefix = experimental_hot_prefix();
  const int balanced_groups = experimental_balanced_groups();
  acc_genpp_pp_t pp;
  acc_genpp_td_t td;
  acc_genpp_state_t state;
  acc_genpp_wit_t cumulative, increment, plain_fresh;
  acc_genpp_update_t add_update;
  poly_q_vec_d accumulator, update_u;
  uint64_t fresh_norm2;
  uint64_t plain_fresh_norm2;
  uint8_t plain_seed[ACC_GENPP_SEED_BYTES];
  long double sum_increment_norm2 = 0.0L;
  long double sum_increment_norm = 0.0L;
  uint64_t max_increment_norm2 = 0;
  uint64_t max_increment_group_norm2 = 0;
  uint64_t prefix_histogram[ACC_GENPP_HANDLE_BITS + 1] = {0};

  arith_setup();
  random_init();
  acc_genpp_pp_init(&pp);
  if (hot_prefix && balanced_groups) {
    fprintf(stderr, "热前缀分组与权重平衡分组不能同时启用。\n");
    return 2;
  }
  acc_genpp_td_init(&td);
  acc_genpp_state_init(&state);
  acc_genpp_wit_init(&cumulative);
  acc_genpp_wit_init(&increment);
  acc_genpp_wit_init(&plain_fresh);
  poly_q_vec_d_init(accumulator);
  poly_q_vec_d_init(update_u);
  acc_genpp_gen(&pp, &td, &state, accumulator);
  if (!acc_genpp_add(accumulator, &cumulative, &add_update, &state, &pp,
                     &td, accumulator, member_handle)) {
    fprintf(stderr, "无法生成基准成员的新鲜见证。\n");
    return 1;
  }
  fresh_norm2 = acc_genpp_wit_norm2(&cumulative);
  acc_genpp_Fk(plain_seed, state.ka, member_handle, accumulator);
  acc_genpp_sampleleft(&plain_fresh, &pp, &td, accumulator, member_handle,
                       plain_seed);
  plain_fresh_norm2 = acc_genpp_wit_norm2(&plain_fresh);

  printf("累加器更新界真实轨迹测试（不修改生产参数）\n");
  printf("q=%ld，gadget=balanced-base-%d/%d 层，成员句柄=0x%08" PRIx32
         "，删除次数=%zu。\n",
         (long)PARAM_Q, ACC_GENPP_GADGET_BASE, ACC_GENPP_GADGET_DIM,
         member_handle, deletions);
  if (proof_deletions != deletions) {
    printf("仅证明参数使用 deletion_count=%zu（用于快速容量筛选）。\n",
           proof_deletions);
  }
  if (balanced_groups) {
    printf("Labrador right 实验分组=按 EvalFX 活跃权重平衡的 %zu x %zu，"
           "研究上界²=%" PRIu64 "。\n",
           ACC_GENPP_HANDLE_BITS / group_width, group_width,
           (uint64_t)RESEARCH_RIGHT_MAX);
  } else if (hot_prefix) {
    printf("Labrador right 实验分组=前 %zu 个热位单列 + %zu 位冷后缀，"
           "研究上界²=%" PRIu64 "。\n",
           hot_prefix, ACC_GENPP_HANDLE_BITS - hot_prefix,
           (uint64_t)RESEARCH_RIGHT_MAX);
  } else {
    printf("Labrador right 实验分组宽度=%zu，研究上界²=%" PRIu64 "。\n",
           group_width, (uint64_t)RESEARCH_RIGHT_MAX);
  }
  printf("新鲜完整见证 norm²=%" PRIu64 "，norm=%.3f（%.4f q）。\n",
         fresh_norm2, sqrt((double)fresh_norm2),
         sqrt((double)fresh_norm2) / (double)PARAM_Q);
  printf("  当前随机化见证分解：left=%" PRIu64 "，v12=%" PRIu64
         "，v2=%" PRIu64 "，right=%" PRIu64 "，refresh-z=%" PRIu64
         "。\n",
         poly_q_vec_d_norm2(cumulative.left),
         poly_q_vec_d_norm2(cumulative.v12), v2_norm2(&cumulative),
         right_norm2(&cumulative), refresh_norm2(&cumulative));
  printf("普通 SampleLeft 对照 norm²=%" PRIu64 "，norm=%.3f（%.4f q）；"
         "该对照只用于定位随机化开销。\n",
         plain_fresh_norm2, sqrt((double)plain_fresh_norm2),
         sqrt((double)plain_fresh_norm2) / (double)PARAM_Q);
  printf("说明：轨迹使用固定种子的唯一伪随机删除句柄；输出是实现实测，"
         "不是对抗性安全上界。\n\n");
  printf("更新向量模式：%s。%s\n\n",
         fresh_update_u ? "每次删除独立种子展开 u_t" : "沿用固定公共 u",
         fresh_update_u
             ? "该模式用于验证中心化增量路线，尚未接入生产协议。"
             : "该模式与当前生产协议一致。");
  printf("更新符号模式：%s。%s\n\n",
         random_update_sign ? "独立 Rademacher 符号" : "固定正号",
         random_update_sign
             ? "仅用于评估可证明平方根增长路线，尚未接入生产协议。"
             : "与当前生产协议一致。");

  for (size_t index = 0; index < deletions; index++) {
    uint32_t revoked;
    size_t prefix;
    uint64_t increment_norm2;
    uint64_t increment_group_norm2;

    do {
      revoked = (uint32_t)splitmix64(&prng);
    } while (revoked == member_handle);
    prefix = common_low_prefix(revoked, member_handle);
    prefix_histogram[prefix]++;

    if (fresh_update_u) {
      derive_update_u(update_u, &update_prng, index + 1);
    } else {
      poly_q_vec_d_set(update_u, pp.ua);
    }
    acc_genpp_evalfx_apply(&increment, &pp.Ba, revoked, member_handle,
                           update_u);
    if (random_update_sign && (splitmix64(&sign_prng) & 1u)) {
      negate_right(&increment);
    }
    increment_norm2 = right_norm2(&increment);
    increment_group_norm2 =
        right_group_max_norm2(&increment, group_width, hot_prefix,
                              balanced_groups);
    if (increment_norm2 > max_increment_norm2) {
      max_increment_norm2 = increment_norm2;
    }
    if (increment_group_norm2 > max_increment_group_norm2) {
      max_increment_group_norm2 = increment_group_norm2;
    }
    sum_increment_norm2 += (long double)increment_norm2;
    sum_increment_norm += sqrt((double)increment_norm2);
    apply_increment_to_accumulator(accumulator, &pp, member_handle,
                                   &increment);
    add_right(&cumulative, &increment);

    if (is_checkpoint(index + 1, deletions)) {
      const uint64_t cumulative_norm2 = right_norm2(&cumulative);
      const uint64_t full_norm2 = acc_genpp_wit_norm2(&cumulative);
      const uint64_t group_norm2 =
          right_group_max_norm2(&cumulative, group_width, hot_prefix,
                                balanced_groups);
      const double rss_projection = sqrt(
          (double)fresh_norm2 + (double)sum_increment_norm2);
      const double triangle_bound = (double)sum_increment_norm;
      uint64_t research_bound =
          (uint64_t)(index + 1) * RESEARCH_RIGHT_PER_DELETION;
      if (research_bound > RESEARCH_RIGHT_MAX) {
        research_bound = RESEARCH_RIGHT_MAX;
      }
      if (!acc_genpp_memver(&pp, accumulator, member_handle, &cumulative)) {
        fprintf(stderr, "删除 %zu 次后的 MemVer 一致性检查失败。\n",
                index + 1);
        return 1;
      }

      printf("  删除 %4zu 次：累计 right norm²=%" PRIu64
             "，完整 norm²=%" PRIu64 "，完整 norm=%.3f（%.4f q），"
             "最大 right 分组 norm²=%" PRIu64 "，right 最大系数=%" PRIu64
             "；RSS 投影=%.3f，"
             "逐次三角界=%.3f。\n",
             index + 1, cumulative_norm2, full_norm2,
             sqrt((double)full_norm2),
             sqrt((double)full_norm2) / (double)PARAM_Q,
             group_norm2, right_max_abs(&cumulative), rss_projection,
             triangle_bound);
      printf("                 500-window 研究界²=%" PRIu64
             "，该轨迹%s（余量 %.3fx）；此项是实验筛选，不是安全证明。\n",
             research_bound, group_norm2 <= research_bound ? "通过" : "失败",
             group_norm2 == 0 ? 0.0 :
                 (double)research_bound / (double)group_norm2);
    }
  }

  printf("\n单次更新统计：均方范数²=%.3Lf，RMS 范数=%.3f，"
         "最大单次范数=%.3f，最大单次 right 分组 norm²=%" PRIu64 "。\n",
         sum_increment_norm2 / (long double)deletions,
         sqrt((double)(sum_increment_norm2 / (long double)deletions)),
         sqrt((double)max_increment_norm2), max_increment_group_norm2);
  printf("低位共同前缀分布：");
  for (size_t prefix = 0; prefix <= ACC_GENPP_HANDLE_BITS; prefix++) {
    if (prefix_histogram[prefix] != 0) {
      printf(" p%zu=%" PRIu64, prefix, prefix_histogram[prefix]);
    }
  }
  printf("\n");
  if (sqrt((double)acc_genpp_wit_norm2(&cumulative)) < (double)PARAM_Q) {
    printf("判读：该固定轨迹的完整见证范数低于 q；这只说明实现轨迹有余量，"
           "将 RSS 作为安全界仍需证明中心化、亚高斯性和自适应条件。\n");
  } else {
    printf("判读：该固定轨迹的完整见证范数已经不满足 beta<q；"
           "即使更新项呈 RSS 增长，也不能据此认证当前 M-SIS 参数。\n");
  }
  if (online_proof_requested() &&
      !prove_final_witness(&pp, accumulator, &cumulative, member_handle,
                           (uint32_t)proof_deletions)) {
    return 1;
  }

  poly_q_vec_d_clear(update_u);
  poly_q_vec_d_clear(accumulator);
  acc_genpp_wit_clear(&plain_fresh);
  acc_genpp_wit_clear(&increment);
  acc_genpp_wit_clear(&cumulative);
  acc_genpp_state_clear(&state);
  acc_genpp_td_clear(&td);
  acc_genpp_pp_clear(&pp);
  arith_teardown();
  return 0;
}
