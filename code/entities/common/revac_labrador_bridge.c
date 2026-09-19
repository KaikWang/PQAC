#include "revac_labrador_bridge.h"

#include "arith.h"
#include "fips202.h"
#include "poly_q_sampling.h"
#include "revac_proof_codec.h"
#include "acc96_exact_lift.h"
#include "acc96_labrador_wire.h"
#include "revac_q96_link_wire.h"

#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/syscall.h>
#include <sys/un.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#ifndef REVAC_LABRADOR_DEFAULT_PYTHON
#define REVAC_LABRADOR_DEFAULT_PYTHON "/usr/bin/python3"
#endif
#ifndef REVAC_LABRADOR_DEFAULT_SDE
#define REVAC_LABRADOR_DEFAULT_SDE "sde64"
#endif
#ifndef REVAC_LABRADOR_DEFAULT_SCRIPT
#define REVAC_LABRADOR_DEFAULT_SCRIPT \
  "../zk_lazer/revac_acc_memver_import_labrador.py"
#endif
#ifndef REVAC_LABRADOR_DEFAULT_WORKER_SCRIPT
#define REVAC_LABRADOR_DEFAULT_WORKER_SCRIPT \
  "../zk_lazer/revac_acc_memver_worker.py"
#endif
#ifndef REVAC_LABRADOR_DEFAULT_Q96_SCRIPT
#define REVAC_LABRADOR_DEFAULT_Q96_SCRIPT \
  "../zk_lazer/revac_acc96_memver_labrador.py"
#endif
#ifndef REVAC_LABRADOR_DEFAULT_Q96_WORKER_SCRIPT
#define REVAC_LABRADOR_DEFAULT_Q96_WORKER_SCRIPT \
  "../zk_lazer/revac_acc96_memver_worker.py"
#endif

#define REVAC_WORKER_MAGIC "RVMVWB03"
#define REVAC_WORKER_MAGIC_BYTES 8
#define REVAC_WORKER_PP_ID_BYTES 32
#define REVAC_WORKER_MAX_FRAME (64U * 1024U * 1024U)

static double revac_labrador_monotonic_seconds(void) {
  struct timespec now;

  (void)clock_gettime(CLOCK_MONOTONIC, &now);
  return (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
}

static void revac_labrador_print_timing(const char *label, double started) {
  printf("[Labrador 桥接计时] %-35s %10.6f 秒\n",
         label, revac_labrador_monotonic_seconds() - started);
}

static const char *revac_env_or_default(const char *name, const char *fallback) {
  const char *value = getenv(name);
  return (value != NULL && value[0] != '\0') ? value : fallback;
}

int revac_labrador_online_enabled(void) {
  const char *enabled = getenv("REVAC_LABRADOR_ONLINE");
  return enabled != NULL &&
         (strcmp(enabled, "1") == 0 ||
          strcmp(enabled, "true") == 0 ||
          strcmp(enabled, "TRUE") == 0 ||
          strcmp(enabled, "yes") == 0);
}

static int revac_labrador_use_sde(void) {
  const char *use_sde = getenv("REVAC_LABRADOR_USE_SDE");
  return use_sde == NULL || strcmp(use_sde, "0") != 0;
}

static const char *revac_lift_backend_name(acc_genpp_lift_backend_t backend) {
  switch (backend) {
    case ACC_GENPP_LIFT_BACKEND_INT64:
      return "int64 快速路径";
    case ACC_GENPP_LIFT_BACKEND_FLINT:
      return "FLINT 精确回退";
    case ACC_GENPP_LIFT_BACKEND_NONE:
    default:
      return "未选择";
  }
}

static int revac_lift_mode_from_env(acc_genpp_lift_mode_t *mode,
                                    int *use_python_reference) {
  const char *value = revac_env_or_default("REVAC_LIFT_BACKEND", "auto");

  *use_python_reference = 0;
  if (strcmp(value, "auto") == 0) {
    *mode = ACC_GENPP_LIFT_AUTO;
  } else if (strcmp(value, "int64") == 0) {
    *mode = ACC_GENPP_LIFT_FORCE_INT64;
  } else if (strcmp(value, "flint") == 0 || strcmp(value, "fallback") == 0) {
    *mode = ACC_GENPP_LIFT_FORCE_FLINT;
  } else if (strcmp(value, "check") == 0 ||
             strcmp(value, "crosscheck") == 0) {
    *mode = ACC_GENPP_LIFT_CROSSCHECK;
  } else if (strcmp(value, "python") == 0) {
    *mode = ACC_GENPP_LIFT_AUTO;
    *use_python_reference = 1;
  } else {
    fprintf(stderr,
            "无效的 REVAC_LIFT_BACKEND=%s；可选 auto、int64、flint、check、python\n",
            value);
    return 0;
  }
  return 1;
}

static void json_lift(FILE *fp, int64_t lift[PARAM_D][PARAM_N]) {
  fprintf(fp, "[");
  for (size_t row = 0; row < PARAM_D; row++) {
    if (row != 0) {
      fprintf(fp, ",");
    }
    fprintf(fp, "[");
    for (size_t k = 0; k < PARAM_N; k++) {
      if (k != 0) {
        fprintf(fp, ",");
      }
      fprintf(fp, "%" PRId64, lift[row][k]);
    }
    fprintf(fp, "]");
  }
  fprintf(fp, "]");
}

static void json_poly(FILE *fp, const poly_q p) {
  fprintf(fp, "[");
  for (size_t i = 0; i < PARAM_N; i++) {
    if (i) {
      fprintf(fp, ",");
    }
    fprintf(fp, "%ld", (long)poly_q_get_coeff_centered(p, i));
  }
  fprintf(fp, "]");
}

static void json_vec_d(FILE *fp, const poly_q_vec_d v) {
  fprintf(fp, "[");
  for (size_t i = 0; i < PARAM_D; i++) {
    if (i) {
      fprintf(fp, ",");
    }
    json_poly(fp, v->entries[i]);
  }
  fprintf(fp, "]");
}

static void json_mat_d_d(FILE *fp, const poly_q_mat_d_d m) {
  fprintf(fp, "[");
  for (size_t r = 0; r < PARAM_D; r++) {
    if (r) {
      fprintf(fp, ",");
    }
    fprintf(fp, "[");
    for (size_t c = 0; c < PARAM_D; c++) {
      if (c) {
        fprintf(fp, ",");
      }
      json_poly(fp, m->rows[r]->entries[c]);
    }
    fprintf(fp, "]");
  }
  fprintf(fp, "]");
}

static int write_memver_json(FILE *fp,
                             const revac_public_state_t *public_state,
                             const poly_q_vec_d acc,
                             const revac_uav_t *uav,
                             const revac_show_proof_t *show_proof,
                             const poly_q_vec_d link_randomness[2]) {
  poly_q_mat_d_d Atrap, A_link, D_link;
  int include_witness = (uav != NULL);
  int use_python_reference = 0;
  int have_c_lift = 0;
  int64_t lift[PARAM_D][PARAM_N];
  int64_t link_lift[PARAM_D][PARAM_N];
  acc_genpp_lift_mode_t lift_mode = ACC_GENPP_LIFT_AUTO;
  acc_genpp_lift_info_t lift_info;

  if (show_proof == NULL || (include_witness && link_randomness == NULL)) {
    return 0;
  }
  if (include_witness &&
      !revac_lift_mode_from_env(&lift_mode, &use_python_reference)) {
    return 0;
  }
  if (include_witness && !use_python_reference) {
    const double lift_started = revac_labrador_monotonic_seconds();
    if (!acc_genpp_memver_exact_lift(lift, public_state->acc_pp, acc,
                                     uav->handle, &uav->witness,
                                     lift_mode, &lift_info)) {
      fprintf(stderr, "C exact-lift 计算或双路径一致性检查失败\n");
      return 0;
    }
    have_c_lift = 1;
    if (!revac_same_handle_exact_lift(
            link_lift, public_state, uav->handle, link_randomness,
            show_proof->handle_commitment)) {
      fprintf(stderr, "same-handle commitment exact-lift 计算失败\n");
      return 0;
    }
    printf("[Lift 计时] C exact-lift（%s）：%.6f 秒\n",
           revac_lift_backend_name(lift_info.backend),
           revac_labrador_monotonic_seconds() - lift_started);
    printf("[Lift 检查] int64 公共上界=%" PRIu64 "，安全=%s，FLINT 回退=%s\n",
           lift_info.public_abs_bound,
           lift_info.int64_bound_ok ? "是" : "否",
           lift_info.fallback_used ? "是" : "否");
    if (lift_info.crosscheck_performed) {
      printf("[Lift 检查] int64 与 FLINT 逐系数一致性：%s\n",
             lift_info.crosscheck_passed ? "通过" : "失败");
    }
  } else if (include_witness) {
    printf("[Lift 后端] 使用 Python 参考实现（仅用于回归/性能对照）\n");
  }

  poly_q_mat_d_d_init(Atrap);
  poly_q_mat_d_d_init(A_link);
  poly_q_mat_d_d_init(D_link);
  poly_q_mat_d_d_uniform(Atrap, public_state->acc_pp->A.seed,
                         DOMAIN_SEPARATOR_A, 0);
  revac_same_handle_public_matrices(A_link, D_link, public_state);

  fprintf(fp, "{\n");
  fprintf(fp,
          "\"params\":{\"n\":%d,\"d\":%d,\"q\":%ld,\"handle_bits\":%d,"
          "\"gadget_base\":%d,\"gadget_dim\":%d,\"gadget_signed\":%d,"
          "\"trunc_k_l\":%d,\"deletion_count\":%u},\n",
          PARAM_N, PARAM_D, (long)PARAM_Q, ACC_GENPP_HANDLE_BITS,
          ACC_GENPP_GADGET_BASE, ACC_GENPP_GADGET_DIM,
          ACC_GENPP_GADGET_SIGNED, ACC_GENPP_TRUNC_K_L,
          (unsigned)public_state->acc_seq);

  if (include_witness) {
    fprintf(fp, "\"x\":%u,\n", uav->handle);
    fprintf(fp, "\"x_bits\":[");
    for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
      if (i) {
        fprintf(fp, ",");
      }
      fprintf(fp, "%u", (unsigned)((uav->handle >> i) & 1u));
    }
    fprintf(fp, "],\n");
  }

  fprintf(fp, "\"Atrap\":");
  json_mat_d_d(fp, Atrap);
  fprintf(fp, ",\n\"link_A\":");
  json_mat_d_d(fp, A_link);
  fprintf(fp, ",\n\"link_D\":");
  json_mat_d_d(fp, D_link);
  fprintf(fp, ",\n\"link_commitment\":");
  json_vec_d(fp, show_proof->handle_commitment);

  fprintf(fp, ",\n\"trunc_B\":[");
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    if (i) {
      fprintf(fp, ",");
    }
    json_mat_d_d(fp, public_state->acc_pp->trunc_B[i]);
  }

  fprintf(fp, "],\n\"Ba\":[");
  for (size_t h = 0; h < ACC_GENPP_HANDLE_BITS; h++) {
    if (h) {
      fprintf(fp, ",");
    }
    fprintf(fp, "[");
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      if (j) {
        fprintf(fp, ",");
      }
      json_mat_d_d(fp, public_state->acc_pp->Ba.bit[h].blocks[j]);
    }
    fprintf(fp, "]");
  }

  fprintf(fp, "],\n\"acc\":");
  json_vec_d(fp, acc);

  if (have_c_lift) {
    fprintf(fp, ",\n\"lift\":");
    json_lift(fp, lift);
    fprintf(fp, ",\n\"link_lift\":");
    json_lift(fp, link_lift);
  }

  if (include_witness) {
    fprintf(fp, ",\n\"witness\":{");
    fprintf(fp, "\"left\":");
    json_vec_d(fp, uav->witness.left);
    fprintf(fp, ",\"v12\":");
    json_vec_d(fp, uav->witness.v12);
    fprintf(fp, ",\"v2\":[");
    for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
      if (i) {
        fprintf(fp, ",");
      }
      json_vec_d(fp, uav->witness.v2[i]);
    }
    fprintf(fp, "],\"tag\":");
    json_poly(fp, uav->witness.tag);
    fprintf(fp, ",\"link_r\":[");
    json_vec_d(fp, link_randomness[0]);
    fprintf(fp, ",");
    json_vec_d(fp, link_randomness[1]);
    fprintf(fp, "]");
    fprintf(fp, ",\"right\":[");
    for (size_t h = 0; h < ACC_GENPP_HANDLE_BITS; h++) {
      if (h) {
        fprintf(fp, ",");
      }
      fprintf(fp, "[");
      for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
        if (j) {
          fprintf(fp, ",");
        }
        json_vec_d(fp, uav->witness.right[h][j]);
      }
      fprintf(fp, "]");
    }
    fprintf(fp, "]}");
  }

  fprintf(fp, "\n}\n");
  poly_q_mat_d_d_clear(D_link);
  poly_q_mat_d_d_clear(A_link);
  poly_q_mat_d_d_clear(Atrap);
  return ferror(fp) == 0;
}

static int make_temp_path(char path[64]) {
  int fd;

  snprintf(path, 64, "/tmp/revac_labrador_%ld_XXXXXX", (long)getpid());
  fd = mkstemp(path);
  if (fd < 0) {
    return -1;
  }
  return fd;
}

static int write_json_temp(char path[64],
                           const revac_public_state_t *public_state,
                           const poly_q_vec_d acc,
                           const revac_uav_t *uav,
                           const revac_show_proof_t *show_proof,
                           const poly_q_vec_d link_randomness[2]) {
  int fd = make_temp_path(path);
  FILE *fp;
  int ok;

  if (fd < 0) {
    return 0;
  }
  fp = fdopen(fd, "w");
  if (fp == NULL) {
    close(fd);
    unlink(path);
    return 0;
  }
  ok = write_memver_json(fp, public_state, acc, uav, show_proof,
                         link_randomness);
  if (fclose(fp) != 0) {
    ok = 0;
  }
  if (!ok) {
    unlink(path);
  }
  return ok;
}

static int write_bytes_temp(char path[64], const uint8_t *bytes, size_t len) {
  int fd = make_temp_path(path);
  size_t off = 0;

  if (fd < 0) {
    return 0;
  }
  while (off < len) {
    ssize_t got = write(fd, bytes + off, len - off);
    if (got < 0) {
      if (errno == EINTR) {
        continue;
      }
      close(fd);
      unlink(path);
      return 0;
    }
    off += (size_t)got;
  }
  if (close(fd) != 0) {
    unlink(path);
    return 0;
  }
  return 1;
}

static int read_file(const char *path, uint8_t **out, size_t *out_len) {
  FILE *fp = fopen(path, "rb");
  long len;
  uint8_t *buf;
  int ok = 0;

  *out = NULL;
  *out_len = 0;
  if (fp == NULL) {
    return 0;
  }
  if (fseek(fp, 0, SEEK_END) != 0) {
    goto cleanup;
  }
  len = ftell(fp);
  if (len < 0 || fseek(fp, 0, SEEK_SET) != 0) {
    goto cleanup;
  }
  buf = malloc((size_t)len);
  if (buf == NULL && len != 0) {
    goto cleanup;
  }
  if (len != 0 && fread(buf, 1, (size_t)len, fp) != (size_t)len) {
    free(buf);
    goto cleanup;
  }
  *out = buf;
  *out_len = (size_t)len;
  ok = 1;

cleanup:
  fclose(fp);
  return ok;
}

static int spawn_wait(char *const argv[]) {
  pid_t pid = fork();
  int status;

  if (pid < 0) {
    return 0;
  }
  if (pid == 0) {
    execv(argv[0], argv);
    _exit(127);
  }
  do {
    if (waitpid(pid, &status, 0) < 0) {
      if (errno == EINTR) {
        continue;
      }
      return 0;
    }
    break;
  } while (1);
  return WIFEXITED(status) && WEXITSTATUS(status) == 0;
}

static int run_labrador(const char *mode,
                        const char *json_path,
                        const char *proof_flag,
                        const char *proof_path) {
  const char *python = revac_env_or_default("REVAC_LABRADOR_PYTHON",
                                            REVAC_LABRADOR_DEFAULT_PYTHON);
  const char *script = revac_env_or_default("REVAC_LABRADOR_SCRIPT",
                                            REVAC_LABRADOR_DEFAULT_SCRIPT);
  const char *sde = revac_env_or_default("REVAC_LABRADOR_SDE",
                                         REVAC_LABRADOR_DEFAULT_SDE);

  if (revac_labrador_use_sde()) {
    char *const argv[] = {
      (char *)sde, "-icx", "--", (char *)python, (char *)script,
      (char *)json_path, "--mode", (char *)mode,
      (char *)proof_flag, (char *)proof_path, NULL
    };
    return spawn_wait(argv);
  }

  {
    char *const argv[] = {
      (char *)python, (char *)script, (char *)json_path, "--mode",
      (char *)mode, (char *)proof_flag, (char *)proof_path, NULL
    };
    return spawn_wait(argv);
  }
}

static int run_labrador_q96_prove(const char *public_path,
                                  const char *prover_path,
                                  const char *link_path,
                                  const char *proof_path) {
  const char *python = revac_env_or_default("REVAC_LABRADOR_PYTHON",
                                            REVAC_LABRADOR_DEFAULT_PYTHON);
  const char *script = revac_env_or_default("REVAC_LABRADOR_Q96_SCRIPT",
                                            REVAC_LABRADOR_DEFAULT_Q96_SCRIPT);
  char *const argv[] = {
      (char *)python, (char *)script, "--public", (char *)public_path,
      "--prover", (char *)prover_path, "--link", (char *)link_path,
      "--prove", "--skip-verify", "--proof-out", (char *)proof_path, NULL};
  return spawn_wait(argv);
}

static int run_labrador_q96_verify(const char *public_path,
                                   const char *instance_path,
                                   const char *proof_path) {
  const char *python = revac_env_or_default("REVAC_LABRADOR_PYTHON",
                                            REVAC_LABRADOR_DEFAULT_PYTHON);
  const char *script = revac_env_or_default("REVAC_LABRADOR_Q96_SCRIPT",
                                            REVAC_LABRADOR_DEFAULT_Q96_SCRIPT);
  char *const argv[] = {
      (char *)python, (char *)script, "--public", (char *)public_path,
      "--instance", (char *)instance_path, "--proof-in",
      (char *)proof_path, NULL};
  return spawn_wait(argv);
}

/*
 * Binary, same-device worker bridge.  The worker is deliberately keyed by a
 * public-parameter identifier.  It never receives a witness except in a local
 * prove request, and never stores that witness after the response is sent.
 */
typedef struct {
  uint8_t *data;
  size_t len;
  size_t cap;
} revac_worker_buf_t;

static void revac_worker_buf_clear(revac_worker_buf_t *buf) {
  free(buf->data);
  memset(buf, 0, sizeof(*buf));
}

static int revac_worker_buf_reserve(revac_worker_buf_t *buf, size_t len) {
  size_t cap;
  uint8_t *grown;

  if (len > REVAC_WORKER_MAX_FRAME || buf->len > REVAC_WORKER_MAX_FRAME - len) {
    return 0;
  }
  if (buf->len + len <= buf->cap) {
    return 1;
  }
  cap = buf->cap == 0 ? 4096 : buf->cap;
  while (cap < buf->len + len) {
    if (cap > REVAC_WORKER_MAX_FRAME / 2) {
      cap = REVAC_WORKER_MAX_FRAME;
      break;
    }
    cap *= 2;
  }
  grown = realloc(buf->data, cap);
  if (grown == NULL) {
    return 0;
  }
  buf->data = grown;
  buf->cap = cap;
  return 1;
}

static int revac_worker_buf_append(revac_worker_buf_t *buf,
                                   const void *data, size_t len) {
  if (!revac_worker_buf_reserve(buf, len)) {
    return 0;
  }
  memcpy(buf->data + buf->len, data, len);
  buf->len += len;
  return 1;
}

static int revac_worker_buf_u32(revac_worker_buf_t *buf, uint32_t value) {
  uint8_t out[4] = {
      (uint8_t)value, (uint8_t)(value >> 8), (uint8_t)(value >> 16),
      (uint8_t)(value >> 24)};
  return revac_worker_buf_append(buf, out, sizeof(out));
}

static int revac_worker_buf_i64(revac_worker_buf_t *buf, int64_t value) {
  uint64_t raw = (uint64_t)value;
  uint8_t out[8];
  for (size_t i = 0; i < sizeof(out); i++) {
    out[i] = (uint8_t)(raw >> (8 * i));
  }
  return revac_worker_buf_append(buf, out, sizeof(out));
}

static int revac_worker_buf_poly(revac_worker_buf_t *buf, const poly_q poly) {
  const size_t encoded_bytes = PARAM_N * sizeof(uint32_t);

  if (!revac_worker_buf_reserve(buf, encoded_bytes)) {
    return 0;
  }
  for (size_t i = 0; i < PARAM_N; i++) {
    const uint32_t value = (uint32_t)(int32_t)poly_q_get_coeff_centered(poly, i);
    uint8_t *dst = buf->data + buf->len + i * sizeof(uint32_t);
    dst[0] = (uint8_t)value;
    dst[1] = (uint8_t)(value >> 8);
    dst[2] = (uint8_t)(value >> 16);
    dst[3] = (uint8_t)(value >> 24);
  }
  buf->len += encoded_bytes;
  return 1;
}

static int revac_worker_buf_vec_d(revac_worker_buf_t *buf,
                                  const poly_q_vec_d vec) {
  for (size_t row = 0; row < PARAM_D; row++) {
    if (!revac_worker_buf_poly(buf, vec->entries[row])) {
      return 0;
    }
  }
  return 1;
}

static int revac_worker_buf_mat_d_d(revac_worker_buf_t *buf,
                                    const poly_q_mat_d_d mat) {
  for (size_t row = 0; row < PARAM_D; row++) {
    for (size_t col = 0; col < PARAM_D; col++) {
      if (!revac_worker_buf_poly(buf, mat->rows[row]->entries[col])) {
        return 0;
      }
    }
  }
  return 1;
}

static void revac_worker_pp_id(
    uint8_t out[REVAC_WORKER_PP_ID_BYTES],
    const revac_public_state_t *public_state) {
  static const uint8_t domain[] = "REVAC-LABRADOR-MEMVER-PP-V3-RANDOM-REFRESH";
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, domain, sizeof(domain) - 1);
  shake256_absorb(&state, public_state->context_pp_hash,
                  sizeof(public_state->context_pp_hash));
  shake256_absorb(&state, public_state->acc_pp->refresh_seed,
                  sizeof(public_state->acc_pp->refresh_seed));
  shake256_absorb(&state,
                  (const uint8_t *)&public_state->acc_pp->refresh_bits,
                  sizeof(public_state->acc_pp->refresh_bits));
  /* context_pp_hash binds A.seed and Ba.  Bind trunc_B here too: the current
     revocation-context hash predates the MemVer worker and omits it. */
  for (size_t block = 0; block < ACC_GENPP_TRUNC_K_L; block++) {
    for (size_t row = 0; row < PARAM_D; row++) {
      for (size_t col = 0; col < PARAM_D; col++) {
        for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
          const uint32_t value = (uint32_t)poly_q_get_coeff(
              public_state->acc_pp->trunc_B[block]->rows[row]->entries[col],
              coeff);
          uint8_t encoded[4] = {(uint8_t)value, (uint8_t)(value >> 8),
                                (uint8_t)(value >> 16),
                                (uint8_t)(value >> 24)};
          shake256_absorb(&state, encoded, sizeof(encoded));
        }
      }
    }
  }
  shake256_finalize(&state);
  shake256_squeeze(out, REVAC_WORKER_PP_ID_BYTES, &state);
}

static int revac_worker_build_public(revac_worker_buf_t *out,
                                     const revac_public_state_t *public_state) {
  poly_q_mat_d_d Atrap, A_link, D_link;
  const acc_genpp_pp_t *pp = public_state->acc_pp;
  int ok = 1;
  const char *refresh_blob_mode = getenv("REVAC_REFRESH_PUBLIC_BLOB");
  const int seed_only =
      refresh_blob_mode == NULL || refresh_blob_mode[0] == '\0' ||
      strcmp(refresh_blob_mode, "explicit") != 0;
  const size_t matrix_count =
      1 + ACC_GENPP_TRUNC_K_L +
      ACC_GENPP_HANDLE_BITS * ACC_GENPP_GADGET_DIM + 2 +
      (seed_only ? 0 :
       ((size_t)pp->refresh_bits + 1) * ACC_GENPP_GADGET_DIM);
  const size_t expected = REVAC_WORKER_MAGIC_BYTES + sizeof(pp->A.seed) +
      sizeof(uint32_t) + sizeof(pp->refresh_seed) + 1 +
      matrix_count * PARAM_D * PARAM_D * PARAM_N * sizeof(uint32_t);

  if (!revac_worker_buf_reserve(out, expected)) {
    return 0;
  }

  poly_q_mat_d_d_init(Atrap);
  poly_q_mat_d_d_init(A_link);
  poly_q_mat_d_d_init(D_link);
  poly_q_mat_d_d_uniform(Atrap, pp->A.seed, DOMAIN_SEPARATOR_A, 0);
  revac_same_handle_public_matrices(A_link, D_link, public_state);
  ok &= revac_worker_buf_append(out, REVAC_WORKER_MAGIC,
                                REVAC_WORKER_MAGIC_BYTES);
  ok &= revac_worker_buf_append(out, pp->A.seed, sizeof(pp->A.seed));
  ok &= revac_worker_buf_u32(out, pp->refresh_bits);
  ok &= revac_worker_buf_append(out, pp->refresh_seed,
                                sizeof(pp->refresh_seed));
  {
    const uint8_t mode = seed_only ? 1u : 0u;
    ok &= revac_worker_buf_append(out, &mode, sizeof(mode));
  }
  ok &= revac_worker_buf_mat_d_d(out, Atrap);
  for (size_t i = 0; ok && i < ACC_GENPP_TRUNC_K_L; i++) {
    ok &= revac_worker_buf_mat_d_d(out, pp->trunc_B[i]);
  }
  for (size_t bit = 0; ok && bit < ACC_GENPP_HANDLE_BITS; bit++) {
    for (size_t gadget = 0; ok && gadget < ACC_GENPP_GADGET_DIM; gadget++) {
      ok &= revac_worker_buf_mat_d_d(out, pp->Ba.bit[bit].blocks[gadget]);
    }
  }
  ok &= revac_worker_buf_mat_d_d(out, A_link);
  ok &= revac_worker_buf_mat_d_d(out, D_link);
  if (!seed_only) {
    acc_genpp_gmat_t refresh;
    acc_genpp_gmat_init(&refresh);
    for (size_t i = 0; ok && i <= pp->refresh_bits; i++) {
      acc_genpp_refresh_matrix_expand(&refresh, pp->refresh_seed, i);
      for (size_t j = 0; ok && j < ACC_GENPP_GADGET_DIM; j++) {
        ok &= revac_worker_buf_mat_d_d(out, refresh.blocks[j]);
      }
    }
    acc_genpp_gmat_clear(&refresh);
  }
  poly_q_mat_d_d_clear(D_link);
  poly_q_mat_d_d_clear(A_link);
  poly_q_mat_d_d_clear(Atrap);
  return ok;
}

static int revac_worker_build_dynamic_verify(
    revac_worker_buf_t *out, const revac_show_context_t *ctx,
    const revac_show_proof_t *proof) {
  const uint8_t *proof_bytes = proof->acc_zk_proof;
  size_t proof_len = proof->acc_zk_proof_len;
  uint8_t *raw_proof = NULL;
  int ok = 0;

  if (proof->acc_zk_kind == REVAC_ACC_ZK_LABRADOR_MEMVER_DEFLATE) {
    if (!revac_deflate_unpack(&raw_proof, &proof_len, proof->acc_zk_proof,
                              proof->acc_zk_proof_len)) {
      return 0;
    }
    proof_bytes = raw_proof;
  } else if (proof->acc_zk_kind != REVAC_ACC_ZK_LABRADOR_MEMVER) {
    return 0;
  }
  if (proof_len == 0 || proof_len > UINT32_MAX) {
    goto cleanup;
  }
  if (getenv("REVAC_LABRADOR_PROFILE") != NULL) {
    fprintf(stderr,
            "[Labrador bridge verify] kind=%u，wire=%zu，raw=%zu，magic=%02x%02x%02x%02x%02x%02x%02x%02x\n",
            (unsigned int)proof->acc_zk_kind, proof->acc_zk_proof_len,
            proof_len,
            proof_len > 0 ? proof_bytes[0] : 0,
            proof_len > 1 ? proof_bytes[1] : 0,
            proof_len > 2 ? proof_bytes[2] : 0,
            proof_len > 3 ? proof_bytes[3] : 0,
            proof_len > 4 ? proof_bytes[4] : 0,
            proof_len > 5 ? proof_bytes[5] : 0,
            proof_len > 6 ? proof_bytes[6] : 0,
            proof_len > 7 ? proof_bytes[7] : 0);
  }
  ok = revac_worker_buf_u32(out, ctx->acc_seq) &&
       revac_worker_buf_vec_d(out, ctx->acc) &&
       revac_worker_buf_vec_d(out, proof->handle_commitment) &&
       revac_worker_buf_u32(out, (uint32_t)proof_len) &&
       revac_worker_buf_append(out, proof_bytes, proof_len);

cleanup:
  free(raw_proof);
  return ok;
}

static int revac_worker_build_dynamic_prove(
    revac_worker_buf_t *out, const revac_public_state_t *public_state,
    revac_uav_t *uav, const revac_show_proof_t *show_proof,
    const poly_q_vec_d link_randomness[2]) {
  int64_t lift[PARAM_D][PARAM_N];
  int64_t link_lift[PARAM_D][PARAM_N];
  acc_genpp_lift_mode_t lift_mode = ACC_GENPP_LIFT_AUTO;
  acc_genpp_lift_info_t lift_info;
  int use_python_reference = 0;
  int lift_cache_hit = 0;
  const char *lift_cache_env = getenv("REVAC_MEMVER_LIFT_CACHE");
  const int lift_cache_enabled =
      lift_cache_env == NULL || lift_cache_env[0] == '\0' ||
      lift_cache_env[0] != '0';
  const double lift_started = revac_labrador_monotonic_seconds();

  if (!revac_lift_mode_from_env(&lift_mode, &use_python_reference) ||
      use_python_reference) {
    return 0;
  }
  if (lift_cache_enabled && lift_mode == ACC_GENPP_LIFT_AUTO &&
      uav->memver_cache_valid &&
      uav->memver_lift_cache_valid &&
      memcmp(uav->memver_lift_cache_key, uav->memver_cache_key,
             REVAC_CONTEXT_HASH_BYTES) == 0) {
    memcpy(lift, uav->memver_lift_cache, sizeof(lift));
    memset(&lift_info, 0, sizeof(lift_info));
    lift_info.backend = ACC_GENPP_LIFT_BACKEND_INT64;
    lift_cache_hit = 1;
  } else {
    if (!acc_genpp_memver_exact_lift(lift, public_state->acc_pp,
                                     public_state->acc, uav->handle,
                                     &uav->witness, lift_mode, &lift_info)) {
      return 0;
    }
    if (lift_cache_enabled && lift_mode == ACC_GENPP_LIFT_AUTO &&
        uav->memver_cache_valid) {
      memcpy(uav->memver_lift_cache, lift, sizeof(lift));
      memcpy(uav->memver_lift_cache_key, uav->memver_cache_key,
             REVAC_CONTEXT_HASH_BYTES);
      uav->memver_lift_cache_valid = 1;
    }
  }
  if (!revac_same_handle_exact_lift(
          link_lift, public_state, uav->handle, link_randomness,
          show_proof->handle_commitment)) {
    return 0;
  }
  {
    const size_t i32_polys =
        2 * PARAM_D + 2 * PARAM_D + 2 * PARAM_D +
        ACC_GENPP_TRUNC_K_L * PARAM_D + 1 +
        ACC_GENPP_HANDLE_BITS * ACC_GENPP_GADGET_DIM * PARAM_D +
        ACC_GENPP_GADGET_DIM * PARAM_D;
    const size_t i64_polys = 2 * PARAM_D;
    const size_t expected = sizeof(uint32_t) * 3 +
        ACC_GENPP_REFRESH_BYTES_MAX +
        i32_polys * PARAM_N * sizeof(uint32_t) +
        i64_polys * PARAM_N * sizeof(uint64_t);
    if (!revac_worker_buf_reserve(out, expected)) {
      return 0;
    }
  }
  printf("[Lift 计时] C exact-lift（%s，MemVer lift 缓存%s）：%.6f 秒\n",
         revac_lift_backend_name(lift_info.backend),
         lift_cache_hit ? "命中" : "未命中",
         revac_labrador_monotonic_seconds() - lift_started);
  if (!lift_cache_hit) {
    printf("[Lift 检查] int64 公共上界=%" PRIu64 "，安全=%s，FLINT 回退=%s\n",
           lift_info.public_abs_bound,
           lift_info.int64_bound_ok ? "是" : "否",
           lift_info.fallback_used ? "是" : "否");
    if (lift_info.crosscheck_performed) {
      printf("[Lift 检查] int64 与 FLINT 逐系数一致性：%s\n",
             lift_info.crosscheck_passed ? "通过" : "失败");
    }
  }
  if (!revac_worker_buf_u32(out, public_state->acc_seq) ||
      !revac_worker_buf_vec_d(out, public_state->acc) ||
      !revac_worker_buf_vec_d(out, show_proof->handle_commitment) ||
      !revac_worker_buf_u32(out, uav->handle) ||
      !revac_worker_buf_vec_d(out, link_randomness[0]) ||
      !revac_worker_buf_vec_d(out, link_randomness[1]) ||
      !revac_worker_buf_vec_d(out, uav->witness.left) ||
      !revac_worker_buf_vec_d(out, uav->witness.v12)) {
    return 0;
  }
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    if (!revac_worker_buf_vec_d(out, uav->witness.v2[i])) {
      return 0;
    }
  }
  if (!revac_worker_buf_poly(out, uav->witness.tag)) {
    return 0;
  }
  for (size_t bit = 0; bit < ACC_GENPP_HANDLE_BITS; bit++) {
    for (size_t gadget = 0; gadget < ACC_GENPP_GADGET_DIM; gadget++) {
      if (!revac_worker_buf_vec_d(out, uav->witness.right[bit][gadget])) {
        return 0;
      }
    }
  }
  if (!revac_worker_buf_u32(out, public_state->acc_pp->refresh_bits) ||
      !revac_worker_buf_append(out, uav->witness.refresh_r,
                               ACC_GENPP_REFRESH_BYTES_MAX)) {
    return 0;
  }
  for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
    if (!revac_worker_buf_vec_d(out, uav->witness.refresh_z[j])) {
      return 0;
    }
  }
  for (size_t row = 0; row < PARAM_D; row++) {
    for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
      if (!revac_worker_buf_i64(out, lift[row][coeff])) {
        return 0;
      }
    }
  }
  for (size_t row = 0; row < PARAM_D; row++) {
    for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
      if (!revac_worker_buf_i64(out, link_lift[row][coeff])) {
        return 0;
      }
    }
  }
  return 1;
}

static int revac_worker_bridge_enabled(void) {
  const char *value = revac_env_or_default("REVAC_LABRADOR_BRIDGE", "worker");
  return strcmp(value, "json") != 0 && strcmp(value, "legacy") != 0;
}

static int revac_worker_socket_path(char out[sizeof(((struct sockaddr_un *)0)->sun_path)]) {
  const char *value = getenv("REVAC_LABRADOR_SOCKET");
  if (value != NULL && value[0] != '\0') {
    if (strlen(value) >= sizeof(((struct sockaddr_un *)0)->sun_path)) {
      return 0;
    }
    strcpy(out, value);
    return 1;
  }
  return snprintf(out, sizeof(((struct sockaddr_un *)0)->sun_path),
                  "/tmp/revac_labrador_memver_%ld.sock", (long)getuid()) > 0;
}

static int revac_worker_connect(const char *socket_path) {
  struct sockaddr_un address;
  int fd = socket(AF_UNIX, SOCK_STREAM, 0);
  if (fd < 0) {
    return -1;
  }
  memset(&address, 0, sizeof(address));
  address.sun_family = AF_UNIX;
  strncpy(address.sun_path, socket_path, sizeof(address.sun_path) - 1);
  if (connect(fd, (const struct sockaddr *)&address, sizeof(address)) != 0) {
    close(fd);
    return -1;
  }
  return fd;
}

static int revac_worker_write_full(int fd, const uint8_t *data, size_t len) {
  size_t offset = 0;
  while (offset < len) {
    ssize_t sent = write(fd, data + offset, len - offset);
    if (sent < 0) {
      if (errno == EINTR) continue;
      return 0;
    }
    offset += (size_t)sent;
  }
  return 1;
}

static int revac_worker_read_full(int fd, uint8_t *data, size_t len) {
  size_t offset = 0;
  while (offset < len) {
    ssize_t got = read(fd, data + offset, len - offset);
    if (got <= 0) {
      if (got < 0 && errno == EINTR) continue;
      return 0;
    }
    offset += (size_t)got;
  }
  return 1;
}

static int revac_worker_send_frame(int fd, const revac_worker_buf_t *frame) {
  uint8_t header[4];
  if (frame->len == 0 || frame->len > REVAC_WORKER_MAX_FRAME) return 0;
  header[0] = (uint8_t)frame->len;
  header[1] = (uint8_t)(frame->len >> 8);
  header[2] = (uint8_t)(frame->len >> 16);
  header[3] = (uint8_t)(frame->len >> 24);
  return revac_worker_write_full(fd, header, sizeof(header)) &&
         revac_worker_write_full(fd, frame->data, frame->len);
}

static int revac_worker_recv_frame(int fd, uint8_t **out, size_t *out_len) {
  uint8_t header[4];
  size_t len;
  uint8_t *data;
  *out = NULL;
  *out_len = 0;
  if (!revac_worker_read_full(fd, header, sizeof(header))) return 0;
  len = (size_t)header[0] | ((size_t)header[1] << 8) |
        ((size_t)header[2] << 16) | ((size_t)header[3] << 24);
  if (len == 0 || len > REVAC_WORKER_MAX_FRAME) return 0;
  data = malloc(len);
  if (data == NULL) return 0;
  if (!revac_worker_read_full(fd, data, len)) {
    free(data);
    return 0;
  }
  *out = data;
  *out_len = len;
  return 1;
}

static int revac_worker_spawn(const char *socket_path) {
  const char *python = revac_env_or_default("REVAC_LABRADOR_PYTHON",
                                            REVAC_LABRADOR_DEFAULT_PYTHON);
  const char *script = getenv("REVAC_LABRADOR_WORKER_SCRIPT");
  const char *sde = revac_env_or_default("REVAC_LABRADOR_SDE",
                                         REVAC_LABRADOR_DEFAULT_SDE);
  char derived_script[1024];

  if (script == NULL || script[0] == '\0') {
    const char *import_script = revac_env_or_default(
        "REVAC_LABRADOR_SCRIPT", REVAC_LABRADOR_DEFAULT_SCRIPT);
    const char *slash = strrchr(import_script, '/');
    if (slash != NULL && (size_t)(slash - import_script) +
                             sizeof("/revac_acc_memver_worker.py") <
                             sizeof(derived_script)) {
      const size_t dir_len = (size_t)(slash - import_script);
      memcpy(derived_script, import_script, dir_len);
      memcpy(derived_script + dir_len, "/revac_acc_memver_worker.py",
             sizeof("/revac_acc_memver_worker.py"));
      script = derived_script;
    } else {
      script = REVAC_LABRADOR_DEFAULT_WORKER_SCRIPT;
    }
  }
  pid_t pid = fork();
  if (pid < 0) return 0;
  if (pid == 0) {
    /*
     * The bridge is intentionally persistent.  Do not let its inherited
     * stdout/stderr keep a caller's ``2>&1 | tee ...`` pipe alive after the C
     * demo exits: otherwise tee never receives EOF and the interactive shell
     * appears to hang even though all tests have completed.  A caller that
     * needs native profiler output can opt in to a separate log file.
     */
    const char *worker_log = getenv("REVAC_LABRADOR_WORKER_LOG");
    int null_fd;
    int log_fd;

    (void)setsid();
    null_fd = open("/dev/null", O_RDONLY);
    if (null_fd >= 0) {
      (void)dup2(null_fd, STDIN_FILENO);
      if (null_fd > STDERR_FILENO) (void)close(null_fd);
    }
    log_fd = open((worker_log != NULL && worker_log[0] != '\0') ?
                      worker_log : "/dev/null",
                  O_WRONLY | O_CREAT | O_APPEND, 0600);
    if (log_fd >= 0) {
      (void)dup2(log_fd, STDOUT_FILENO);
      (void)dup2(log_fd, STDERR_FILENO);
      if (log_fd > STDERR_FILENO) (void)close(log_fd);
    }
    /* A non-interactive shell may keep an extra copy of the pipeline's write
       end above fd 2.  If the persistent worker inherits it, ``tee`` never
       observes EOF after the C program exits.  Close every non-stdio fd before
       exec; the worker creates its Unix socket only after exec. */
#ifdef SYS_close_range
    if (syscall(SYS_close_range, 3U, ~0U, 0U) != 0)
#endif
    {
      long open_max = sysconf(_SC_OPEN_MAX);
      if (open_max < 0) open_max = 1024;
      for (long inherited_fd = 3; inherited_fd < open_max; inherited_fd++) {
        (void)close((int)inherited_fd);
      }
    }
    if (revac_labrador_use_sde()) {
      char *const argv[] = {(char *)sde, "-icx", "--", (char *)python,
                            (char *)script, "--socket", (char *)socket_path,
                            NULL};
      execv(argv[0], argv);
    } else {
      char *const argv[] = {(char *)python, (char *)script, "--socket",
                            (char *)socket_path, NULL};
      execv(argv[0], argv);
    }
    _exit(127);
  }
  return 1;
}

static int revac_worker_open(const char *socket_path) {
  int fd = revac_worker_connect(socket_path);
  if (fd >= 0) return fd;
  if (!revac_worker_spawn(socket_path)) return -1;
  for (size_t retry = 0; retry < 50; retry++) {
    struct timespec delay = {0, 20000000L};
    nanosleep(&delay, NULL);
    fd = revac_worker_connect(socket_path);
    if (fd >= 0) return fd;
  }
  return -1;
}

static int revac_worker_make_request(revac_worker_buf_t *request,
                                     char command,
                                     const uint8_t pp_id[REVAC_WORKER_PP_ID_BYTES],
                                     const revac_worker_buf_t *dynamic,
                                     const revac_worker_buf_t *public_blob) {
  if (!revac_worker_buf_append(request, &command, 1) ||
      !revac_worker_buf_append(request, pp_id, REVAC_WORKER_PP_ID_BYTES) ||
      dynamic->len > UINT32_MAX ||
      !revac_worker_buf_u32(request, (uint32_t)dynamic->len) ||
      !revac_worker_buf_append(request, dynamic->data, dynamic->len)) {
    return 0;
  }
  if (public_blob != NULL) {
    if (public_blob->len > UINT32_MAX ||
        !revac_worker_buf_u32(request, (uint32_t)public_blob->len) ||
        !revac_worker_buf_append(request, public_blob->data, public_blob->len)) {
      return 0;
    }
  }
  return 1;
}

static int revac_worker_request_ex(char mode,
                                   const revac_public_state_t *public_state,
                                   const revac_worker_buf_t *dynamic,
                                   int allow_public_install,
                                   uint8_t **response,
                                   size_t *response_len) {
  char socket_path[sizeof(((struct sockaddr_un *)0)->sun_path)];
  uint8_t pp_id[REVAC_WORKER_PP_ID_BYTES];
  revac_worker_buf_t request = {0}, public_blob = {0};
  uint8_t *reply = NULL;
  size_t reply_len = 0;
  int fd = -1;
  int ok = 0;

  *response = NULL;
  *response_len = 0;
  if (!revac_worker_socket_path(socket_path)) goto cleanup;
  revac_worker_pp_id(pp_id, public_state);
  fd = revac_worker_open(socket_path);
  if (fd < 0 || !revac_worker_make_request(&request, mode, pp_id, dynamic, NULL) ||
      !revac_worker_send_frame(fd, &request) ||
      !revac_worker_recv_frame(fd, &reply, &reply_len)) goto cleanup;
  close(fd);
  fd = -1;
  if (reply_len == 1 && reply[0] == 'N') {
    if (!allow_public_install) goto cleanup;
    free(reply);
    reply = NULL;
    revac_worker_buf_clear(&request);
    if (!revac_worker_build_public(&public_blob, public_state)) goto cleanup;
    fd = revac_worker_open(socket_path);
    if (fd < 0 ||
        !revac_worker_make_request(&request, mode == 'P' ? 'I' : 'J', pp_id,
                                   dynamic, &public_blob) ||
        !revac_worker_send_frame(fd, &request) ||
        !revac_worker_recv_frame(fd, &reply, &reply_len)) goto cleanup;
    close(fd);
    fd = -1;
  }
  if (reply_len < 2 || reply[0] != 'S') {
    if (reply_len > 1 && reply[0] == 'F') {
      fprintf(stderr, "Labrador worker 拒绝请求：%.*s\n",
              (int)(reply_len - 1), (const char *)(reply + 1));
    }
    goto cleanup;
  }
  if (mode == 'P') {
    printf("[Labrador worker] 常驻原生 Pack 参数缓存：%s\n",
           reply[1] != 0 ? "命中（无反序列化）" : "首次生成");
  }
  memmove(reply, reply + 2, reply_len - 2);
  *response = reply;
  *response_len = reply_len - 2;
  reply = NULL;
  ok = 1;

cleanup:
  if (fd >= 0) close(fd);
  free(reply);
  revac_worker_buf_clear(&request);
  revac_worker_buf_clear(&public_blob);
  return ok;
}

static int revac_worker_request(char mode,
                                const revac_public_state_t *public_state,
                                const revac_worker_buf_t *dynamic,
                                uint8_t **response, size_t *response_len) {
  return revac_worker_request_ex(mode, public_state, dynamic, 1,
                                 response, response_len);
}

/* q96 uses an independent resident worker because its ring, wire ABI and
 * 52-constraint statement are deliberately separated from the q19 profile. */
static int revac_q96_worker_socket_path(
    char out[sizeof(((struct sockaddr_un *)0)->sun_path)]) {
  const char *value = getenv("REVAC_LABRADOR_Q96_SOCKET");
  if (value != NULL && value[0] != '\0') {
    if (strlen(value) >= sizeof(((struct sockaddr_un *)0)->sun_path)) return 0;
    strcpy(out, value);
    return 1;
  }
  return snprintf(out, sizeof(((struct sockaddr_un *)0)->sun_path),
                  "/tmp/revac_labrador_q96_%ld.sock", (long)getuid()) > 0;
}

static int revac_q96_worker_spawn(const char *socket_path) {
  const char *python = revac_env_or_default("REVAC_LABRADOR_PYTHON",
                                            REVAC_LABRADOR_DEFAULT_PYTHON);
  const char *script = revac_env_or_default(
      "REVAC_LABRADOR_Q96_WORKER_SCRIPT",
      REVAC_LABRADOR_DEFAULT_Q96_WORKER_SCRIPT);
  pid_t pid = fork();
  if (pid < 0) return 0;
  if (pid == 0) {
    const char *worker_log = getenv("REVAC_LABRADOR_Q96_WORKER_LOG");
    int null_fd, log_fd;
    (void)setsid();
    null_fd = open("/dev/null", O_RDONLY);
    if (null_fd >= 0) {
      (void)dup2(null_fd, STDIN_FILENO);
      if (null_fd > STDERR_FILENO) (void)close(null_fd);
    }
    log_fd = open((worker_log != NULL && worker_log[0] != '\0') ?
                      worker_log : "/dev/null",
                  O_WRONLY | O_CREAT | O_APPEND, 0600);
    if (log_fd >= 0) {
      (void)dup2(log_fd, STDOUT_FILENO);
      (void)dup2(log_fd, STDERR_FILENO);
      if (log_fd > STDERR_FILENO) (void)close(log_fd);
    }
#ifdef SYS_close_range
    if (syscall(SYS_close_range, 3U, ~0U, 0U) != 0)
#endif
    {
      long open_max = sysconf(_SC_OPEN_MAX);
      if (open_max < 0) open_max = 1024;
      for (long inherited_fd = 3; inherited_fd < open_max; inherited_fd++)
        (void)close((int)inherited_fd);
    }
    {
      char *const argv[] = {(char *)python, (char *)script, "--socket",
                            (char *)socket_path, NULL};
      execv(argv[0], argv);
    }
    _exit(127);
  }
  return 1;
}

static int revac_q96_worker_open(const char *socket_path) {
  int fd = revac_worker_connect(socket_path);
  if (fd >= 0) return fd;
  if (!revac_q96_worker_spawn(socket_path)) return -1;
  for (size_t retry = 0; retry < 100; retry++) {
    struct timespec delay = {0, 20000000L};
    nanosleep(&delay, NULL);
    fd = revac_worker_connect(socket_path);
    if (fd >= 0) return fd;
  }
  return -1;
}

static void revac_q96_worker_pp_id(uint8_t out[REVAC_WORKER_PP_ID_BYTES],
                                   const uint8_t *public_wire,
                                   size_t public_len) {
  keccak_state state;
  shake256_init(&state);
  shake256_absorb(&state, public_wire, public_len);
  shake256_finalize(&state);
  shake256_squeeze(out, REVAC_WORKER_PP_ID_BYTES, &state);
}

static int revac_q96_worker_request(
    char mode, const uint8_t *public_wire, size_t public_len,
    const revac_worker_buf_t *dynamic,
    uint8_t **response, size_t *response_len, int *cache_hit) {
  char socket_path[sizeof(((struct sockaddr_un *)0)->sun_path)];
  uint8_t pp_id[REVAC_WORKER_PP_ID_BYTES];
  revac_worker_buf_t request = {0}, public_blob = {0};
  uint8_t *reply = NULL;
  size_t reply_len = 0;
  int fd = -1, ok = 0;

  *response = NULL;
  *response_len = 0;
  if (cache_hit != NULL) *cache_hit = 0;
  if (!revac_q96_worker_socket_path(socket_path)) goto cleanup;
  revac_q96_worker_pp_id(pp_id, public_wire, public_len);
  fd = revac_q96_worker_open(socket_path);
  if (fd < 0 || !revac_worker_make_request(
                    &request, mode, pp_id, dynamic, NULL) ||
      !revac_worker_send_frame(fd, &request) ||
      !revac_worker_recv_frame(fd, &reply, &reply_len)) goto cleanup;
  close(fd);
  fd = -1;
  if (reply_len == 1 && reply[0] == 'N') {
    char install_mode;

    /* A compact prover-witness reuse miss is expected after a worker restart,
       accumulator transition, or UAV switch.  Its caller retries with the
       complete private payload; never misinterpret it as a verifier install. */
    if (mode == 'R') goto cleanup;
    if (mode == 'P') {
      install_mode = 'I';
    } else if (mode == 'V') {
      install_mode = 'J';
    } else if (mode == 'U') {
      install_mode = 'K';
    } else if (mode == 'G') {
      install_mode = 'L';
    } else {
      goto cleanup;
    }
    free(reply);
    reply = NULL;
    revac_worker_buf_clear(&request);
    if (!revac_worker_buf_append(&public_blob, public_wire, public_len))
      goto cleanup;
    fd = revac_q96_worker_open(socket_path);
    if (fd < 0 || !revac_worker_make_request(
                      &request, install_mode, pp_id,
                      dynamic, &public_blob) ||
        !revac_worker_send_frame(fd, &request) ||
        !revac_worker_recv_frame(fd, &reply, &reply_len)) goto cleanup;
    close(fd);
    fd = -1;
  }
  if (reply_len < 2 || reply[0] != 'S') {
    if (reply_len > 1 && reply[0] == 'F')
      fprintf(stderr, "q96 Labrador worker 拒绝请求：%.*s\n",
              (int)(reply_len - 1), (const char *)(reply + 1));
    goto cleanup;
  }
  if (cache_hit != NULL) *cache_hit = reply[1] != 0;
  memmove(reply, reply + 2, reply_len - 2);
  *response = reply;
  *response_len = reply_len - 2;
  reply = NULL;
  ok = 1;

cleanup:
  if (fd >= 0) close(fd);
  free(reply);
  revac_worker_buf_clear(&request);
  revac_worker_buf_clear(&public_blob);
  return ok;
}

int revac_labrador_q96_prover_prewarm(
    const revac_public_state_t *public_state, const revac_uav_t *uav) {
  acc96_exact_lift_t *lift = NULL;
  uint8_t *public_wire = NULL, *prover_wire = NULL, *link_wire = NULL;
  uint8_t *reply = NULL;
  size_t reply_len = 0;
  unsigned lift_bits = 0;
  revac_worker_buf_t dynamic = {0};
  poly_q_vec_d commitment, link_randomness[2];
  int commitment_initialized = 0;
  int cache_hit = 0;
  int ok = 0;

  if (!revac_labrador_online_enabled() || public_state == NULL || uav == NULL ||
      public_state->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      public_state->acc96_ctx == NULL || public_state->acc96_pp == NULL ||
      public_state->acc96_acc == NULL || uav->acc96_witness == NULL ||
      uav->last_seq != public_state->acc96_seq)
    return 0;
  poly_q_vec_d_init(commitment);
  for (size_t i = 0; i < 2; i++) poly_q_vec_d_init(link_randomness[i]);
  commitment_initialized = 1;
  revac_same_handle_commitment_generate(
      commitment, link_randomness, public_state, uav->handle);
  lift = malloc(sizeof(*lift));
  public_wire = malloc(ACC96_LABRADOR_PUBLIC_BYTES);
  prover_wire = malloc(ACC96_LABRADOR_PROVER_BYTES);
  link_wire = malloc(REVAC_Q96_LINK_WIRE_BYTES);
  if (lift == NULL || public_wire == NULL || prover_wire == NULL ||
      link_wire == NULL ||
      !acc96_memver_exact_lift(lift, &lift_bits, public_state->acc96_pp,
                               public_state->acc96_acc, uav->handle,
                               uav->acc96_witness, public_state->acc96_ctx) ||
      !acc96_labrador_public_encode(
          public_wire, ACC96_LABRADOR_PUBLIC_BYTES,
          public_state->acc96_public_seed, public_state->acc96_pp) ||
      !acc96_labrador_prover_encode(
          prover_wire, ACC96_LABRADOR_PROVER_BYTES,
          public_state->acc96_seq, uav->handle, public_state->acc96_acc,
          uav->acc96_witness, lift, public_state->acc96_ctx) ||
      !revac_q96_link_wire_encode(
          link_wire, REVAC_Q96_LINK_WIRE_BYTES, public_state, uav->handle,
          link_randomness, commitment) ||
      !revac_worker_buf_append(
          &dynamic, prover_wire, ACC96_LABRADOR_PROVER_BYTES) ||
      !revac_worker_buf_append(
          &dynamic, link_wire, REVAC_Q96_LINK_WIRE_BYTES) ||
      !revac_q96_worker_request(
          'U', public_wire, ACC96_LABRADOR_PUBLIC_BYTES, &dynamic,
          &reply, &reply_len, &cache_hit) || reply_len != 0)
    goto cleanup;
  printf("[q96 Labrador worker] UAV prover 预热：%s\n",
         cache_hit ? "缓存已存在并刷新 witness" : "statement/参数首次构造");
  ok = 1;

cleanup:
  free(reply);
  revac_worker_buf_clear(&dynamic);
  free(link_wire);
  free(prover_wire);
  free(public_wire);
  free(lift);
  if (commitment_initialized) {
    for (size_t i = 0; i < 2; i++) poly_q_vec_d_clear(link_randomness[i]);
    poly_q_vec_d_clear(commitment);
  }
  return ok;
}

int revac_labrador_q96_verifier_prewarm(
    const revac_public_state_t *public_state) {
  uint8_t *public_wire = NULL, *instance_wire = NULL, *reply = NULL;
  size_t reply_len = 0;
  revac_worker_buf_t dynamic = {0};
  poly_q_vec_d commitment;
  int commitment_initialized = 0;
  int cache_hit = 0;
  int ok = 0;

  if (!revac_labrador_online_enabled() || public_state == NULL ||
      public_state->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      public_state->acc96_ctx == NULL || public_state->acc96_pp == NULL ||
      public_state->acc96_acc == NULL)
    return 0;
  poly_q_vec_d_init(commitment);
  poly_q_vec_d_zero(commitment);
  commitment_initialized = 1;
  public_wire = malloc(ACC96_LABRADOR_PUBLIC_BYTES);
  instance_wire = malloc(REVAC_Q96_INSTANCE_WIRE_BYTES);
  if (public_wire == NULL || instance_wire == NULL ||
      !acc96_labrador_public_encode(
          public_wire, ACC96_LABRADOR_PUBLIC_BYTES,
          public_state->acc96_public_seed, public_state->acc96_pp) ||
      !revac_q96_instance_wire_encode(
          instance_wire, REVAC_Q96_INSTANCE_WIRE_BYTES,
          public_state, commitment) ||
      !revac_worker_buf_append(
          &dynamic, instance_wire, REVAC_Q96_INSTANCE_WIRE_BYTES) ||
      !revac_q96_worker_request(
          'G', public_wire, ACC96_LABRADOR_PUBLIC_BYTES, &dynamic,
          &reply, &reply_len, &cache_hit) || reply_len != 0)
    goto cleanup;
  printf("[q96 Labrador worker] GCS verifier 预热：%s\n",
         cache_hit ? "缓存已存在" : "statement/参数首次构造");
  ok = 1;

cleanup:
  free(reply);
  revac_worker_buf_clear(&dynamic);
  free(instance_wire);
  free(public_wire);
  if (commitment_initialized) poly_q_vec_d_clear(commitment);
  return ok;
}

struct revac_labrador_prove_async {
  pthread_t thread;
  const revac_public_state_t *public_state;
  revac_uav_t *uav;
  const revac_show_proof_t *show_proof;
  const poly_q_vec_d *link_randomness;
  revac_worker_buf_t dynamic;
  uint8_t *worker_proof;
  size_t worker_proof_len;
  double prepare_started;
  double worker_started;
  int prepare_ok;
  int prepare_in_thread;
  int request_ok;
  int thread_started;
  int q96_mode;
  revac_show_proof_t q96_proof;
  int q96_proof_initialized;
};

static int revac_labrador_async_prepare_enabled(void) {
  const char *value = getenv("REVAC_SHOW_PARALLEL_PREPARE");

  /* This only overlaps independent work on disjoint q and qshow objects; it
     does not alter the relation, transcript, proof bytes, or worker format. */
  return value == NULL || value[0] == '\0' || value[0] != '0';
}

static void *revac_labrador_prove_async_main(void *opaque) {
  revac_labrador_prove_async_t *job = opaque;

  if (job->q96_mode) {
    job->request_ok = revac_labrador_q96_memver_prove(
        &job->q96_proof, job->public_state, job->uav,
        job->link_randomness);
    return NULL;
  }

  if (job->prepare_in_thread) {
    job->prepare_ok = revac_worker_build_dynamic_prove(
        &job->dynamic, job->public_state, job->uav,
        job->show_proof, job->link_randomness);
    revac_labrador_print_timing(
        "异步线程 C exact-lift + 动态见证二进制编码",
        job->prepare_started);
  }
  if (!job->prepare_ok) {
    return NULL;
  }

  /* Never expand/install public matrices on this thread: that path touches
     the in-process R_q arithmetic and must not overlap qshow arithmetic.
     A cold worker answers N; finish then installs public state serially while
     reusing this job's already encoded exact lift and witness. */
  job->request_ok = revac_worker_request_ex(
      'P', job->public_state, &job->dynamic, 0,
      &job->worker_proof, &job->worker_proof_len);
  return NULL;
}

int revac_labrador_memver_prove_async_start(
    revac_labrador_prove_async_t **job_out,
    const revac_show_proof_t *proof,
    const revac_public_state_t *public_state,
    revac_uav_t *uav,
    const poly_q_vec_d link_randomness[2]) {
  revac_labrador_prove_async_t *job;
  double started;

  if (job_out == NULL || proof == NULL || public_state == NULL ||
      uav == NULL ||
      (public_state->acc_profile == REVAC_ACC_PROFILE_Q19 &&
       !revac_worker_bridge_enabled()) ||
      (public_state->acc_profile != REVAC_ACC_PROFILE_Q19 &&
       public_state->acc_profile != REVAC_ACC_PROFILE_Q96)) {
    return 0;
  }
  *job_out = NULL;
  job = calloc(1, sizeof(*job));
  if (job == NULL) return 0;
  started = revac_labrador_monotonic_seconds();
  job->public_state = public_state;
  job->uav = uav;
  job->show_proof = proof;
  job->link_randomness = link_randomness;
  job->prepare_started = started;
  job->q96_mode = public_state->acc_profile == REVAC_ACC_PROFILE_Q96;
  if (job->q96_mode) {
    /* Only the handle commitment is shared with the accumulator relation.
       Keep a private proof shell so the worker thread and signature prover
       never write the same C object. */
    memset(&job->q96_proof, 0, sizeof(job->q96_proof));
    poly_q_vec_d_init(job->q96_proof.handle_commitment);
    poly_q_vec_d_set(job->q96_proof.handle_commitment,
                     proof->handle_commitment);
    job->q96_proof.acc_zk_kind = REVAC_ACC_ZK_NONE;
    job->q96_proof_initialized = 1;
    job->prepare_ok = 1;
  }
  job->prepare_in_thread = revac_labrador_async_prepare_enabled();
  if (!job->q96_mode && !job->prepare_in_thread) {
    job->prepare_ok = revac_worker_build_dynamic_prove(
        &job->dynamic, public_state, uav, proof, link_randomness);
    if (!job->prepare_ok) {
      revac_worker_buf_clear(&job->dynamic);
      free(job);
      return 0;
    }
    revac_labrador_print_timing(
        "串行 C exact-lift + 动态见证二进制编码", started);
  }
  job->worker_started = revac_labrador_monotonic_seconds();
  if (pthread_create(&job->thread, NULL, revac_labrador_prove_async_main,
                     job) != 0) {
    if (job->q96_proof_initialized) {
      free(job->q96_proof.acc_zk_proof);
      poly_q_vec_d_clear(job->q96_proof.handle_commitment);
    }
    revac_worker_buf_clear(&job->dynamic);
    free(job);
    return 0;
  }
  job->thread_started = 1;
  *job_out = job;
  return 1;
}

int revac_labrador_memver_prove_async_finish(
    revac_labrador_prove_async_t *job,
    revac_show_proof_t *proof) {
  int ok = 0;

  if (job == NULL || proof == NULL) return 0;
  if (job->thread_started) {
    (void)pthread_join(job->thread, NULL);
  }
  if (job->q96_mode) {
    if (job->request_ok && job->q96_proof.acc_zk_proof != NULL &&
        job->q96_proof.acc_zk_proof_len != 0) {
      ok = revac_show_proof_set_acc_zk(
          proof, job->q96_proof.acc_zk_kind,
          job->q96_proof.acc_zk_proof,
          job->q96_proof.acc_zk_proof_len);
    }
    free(job->q96_proof.acc_zk_proof);
    job->q96_proof.acc_zk_proof = NULL;
    if (job->q96_proof_initialized) {
      poly_q_vec_d_clear(job->q96_proof.handle_commitment);
    }
    revac_worker_buf_clear(&job->dynamic);
    free(job);
    return ok;
  }
  if (job->prepare_ok && !job->request_ok) {
    /* The usual reason is a cold worker replying N.  Signature generation is
       complete now, so it is safe to install/expand public matrices on this
       thread while reusing the already encoded exact lift and witness. */
    free(job->worker_proof);
    job->worker_proof = NULL;
    job->worker_proof_len = 0;
    job->worker_started = revac_labrador_monotonic_seconds();
    job->request_ok = revac_worker_request_ex(
        'P', job->public_state, &job->dynamic, 1,
        &job->worker_proof, &job->worker_proof_len);
  }
  if (job->request_ok && job->worker_proof_len != 0) {
    ok = revac_show_proof_set_acc_zk(
        proof, REVAC_ACC_ZK_LABRADOR_MEMVER,
        job->worker_proof, job->worker_proof_len);
    revac_labrador_print_timing(
        "异步 Unix socket + worker 动态绑定/prove",
        job->worker_started);
  }
  free(job->worker_proof);
  revac_worker_buf_clear(&job->dynamic);
  free(job);
  return ok;
}

int revac_labrador_memver_prove(revac_show_proof_t *proof,
                                const revac_public_state_t *public_state,
                                revac_uav_t *uav,
                                const poly_q_vec_d link_randomness[2]) {
  char json_path[64] = {0};
  char proof_path[64] = {0};
  uint8_t *proof_bytes = NULL;
  size_t proof_len = 0;
  int proof_fd;
  double step_started;
  double total_started = revac_labrador_monotonic_seconds();
  int ok = 0;

  if (revac_worker_bridge_enabled()) {
    revac_worker_buf_t dynamic = {0};
    uint8_t *worker_proof = NULL;
    size_t worker_proof_len = 0;
    step_started = revac_labrador_monotonic_seconds();
    if (revac_worker_build_dynamic_prove(&dynamic, public_state, uav, proof,
                                         link_randomness)) {
      revac_labrador_print_timing("C exact-lift + 动态见证二进制编码",
                                  step_started);
      step_started = revac_labrador_monotonic_seconds();
    }
    if (dynamic.len != 0 &&
        revac_worker_request('P', public_state, &dynamic, &worker_proof,
                             &worker_proof_len) &&
        worker_proof_len != 0) {
      ok = revac_show_proof_set_acc_zk(proof, REVAC_ACC_ZK_LABRADOR_MEMVER,
                                       worker_proof, worker_proof_len);
      revac_labrador_print_timing("Unix socket + worker 动态绑定/prove",
                                  step_started);
      free(worker_proof);
      revac_worker_buf_clear(&dynamic);
      revac_labrador_print_timing("Labrador MemVer 证明桥接总计", total_started);
      return ok;
    }
    free(worker_proof);
    revac_worker_buf_clear(&dynamic);
    fprintf(stderr, "Labrador 二进制 worker 不可用，回退到 JSON 子进程桥接。\n");
  }

  step_started = revac_labrador_monotonic_seconds();
  if (!write_json_temp(json_path, public_state, public_state->acc, uav,
                       proof, link_randomness)) {
    goto cleanup;
  }
  revac_labrador_print_timing("C 导出 MemVer JSON", step_started);
  proof_fd = make_temp_path(proof_path);
  if (proof_fd < 0) {
    goto cleanup;
  }
  close(proof_fd);

  step_started = revac_labrador_monotonic_seconds();
  if (!run_labrador("prove", json_path, "--proof-out", proof_path)) {
    goto cleanup;
  }
  revac_labrador_print_timing("Python/Labrador 子进程总计", step_started);
  step_started = revac_labrador_monotonic_seconds();
  if (!read_file(proof_path, &proof_bytes, &proof_len)) {
    goto cleanup;
  }
  revac_labrador_print_timing("C 读取序列化 proof", step_started);
  step_started = revac_labrador_monotonic_seconds();
  ok = revac_show_proof_set_acc_zk(proof, REVAC_ACC_ZK_LABRADOR_MEMVER,
                                   proof_bytes, proof_len);
  revac_labrador_print_timing("proof 写入 show proof", step_started);

cleanup:
  free(proof_bytes);
  if (json_path[0]) {
    unlink(json_path);
  }
  if (proof_path[0]) {
    unlink(proof_path);
  }
  revac_labrador_print_timing("Labrador MemVer 证明桥接总计", total_started);
  return ok;
}

int revac_labrador_memver_verify(
                                 const revac_public_state_t *public_state,
                                 const revac_show_context_t *ctx,
                                 const revac_show_proof_t *proof) {
  char json_path[64] = {0};
  char proof_path[64] = {0};
  uint8_t *raw_proof = NULL;
  const uint8_t *proof_bytes = proof->acc_zk_proof;
  size_t proof_len = proof->acc_zk_proof_len;
  int ok = 0;

  if (proof->acc_zk_proof == NULL || proof->acc_zk_proof_len == 0) {
    return 0;
  }
  if (proof->acc_zk_kind != REVAC_ACC_ZK_LABRADOR_MEMVER_DEFLATE &&
      proof->acc_zk_kind != REVAC_ACC_ZK_LABRADOR_MEMVER) {
    return 0;
  }
  if (revac_worker_bridge_enabled()) {
    revac_worker_buf_t dynamic = {0};
    uint8_t *reply = NULL;
    size_t reply_len = 0;
    const double started = revac_labrador_monotonic_seconds();
    double worker_started;
    /* The binary builder handles compressed and uncompressed envelopes itself.
       Do not inflate here and then immediately inflate the same proof again. */
    if (revac_worker_build_dynamic_verify(&dynamic, ctx, proof)) {
      revac_labrador_print_timing("GCS verify 请求二进制编码", started);
      worker_started = revac_labrador_monotonic_seconds();
    } else {
      worker_started = started;
    }
    if (dynamic.len != 0 &&
        revac_worker_request('V', public_state, &dynamic, &reply, &reply_len) &&
        reply_len == 1) {
      ok = reply[0] == 1;
      revac_labrador_print_timing("Unix socket + worker 动态绑定/verify",
                                  worker_started);
      free(reply);
      revac_worker_buf_clear(&dynamic);
      return ok;
    }
    free(reply);
    revac_worker_buf_clear(&dynamic);
    fprintf(stderr, "Labrador 二进制 worker 不可用，回退到 JSON 子进程桥接。\n");
  }
  /* Only the legacy JSON/file bridge needs a persistent raw proof view. */
  if (proof->acc_zk_kind == REVAC_ACC_ZK_LABRADOR_MEMVER_DEFLATE) {
    if (!revac_deflate_unpack(&raw_proof, &proof_len, proof->acc_zk_proof,
                              proof->acc_zk_proof_len)) {
      return 0;
    }
    proof_bytes = raw_proof;
  } else {
    proof_bytes = proof->acc_zk_proof;
    proof_len = proof->acc_zk_proof_len;
  }
  if (!write_json_temp(json_path, public_state, ctx->acc, NULL, proof,
                       NULL)) {
    goto cleanup;
  }
  if (!write_bytes_temp(proof_path, proof_bytes, proof_len)) {
    goto cleanup;
  }
  ok = run_labrador("verify", json_path, "--proof-in", proof_path);

cleanup:
  free(raw_proof);
  if (json_path[0]) {
    unlink(json_path);
  }
  if (proof_path[0]) {
    unlink(proof_path);
  }
  return ok;
}

int revac_labrador_q96_memver_prove(
    revac_show_proof_t *proof, const revac_public_state_t *public_state,
    const revac_uav_t *uav, const poly_q_vec_d link_randomness[2]) {
  acc96_exact_lift_t *lift = NULL;
  uint8_t *public_wire = NULL, *prover_wire = NULL, *link_wire = NULL;
  uint8_t *proof_bytes = NULL;
  size_t proof_len = 0;
  revac_worker_buf_t dynamic = {0};
  unsigned lift_bits = 0;
  char public_path[64] = {0}, prover_path[64] = {0};
  char link_path[64] = {0}, proof_path[64] = {0};
  int proof_fd = -1;
  int worker_cache_hit = 0;
  int ok = 0;
  const char *profile_value = getenv("REVAC_LABRADOR_PROFILE");
  const int profile = profile_value != NULL && profile_value[0] != '\0' &&
                      profile_value[0] != '0';
  const char *disable_reuse_value =
      getenv("REVAC_LABRADOR_Q96_DISABLE_WITNESS_REUSE");
  const int disable_witness_reuse =
      disable_reuse_value != NULL && disable_reuse_value[0] != '\0' &&
      disable_reuse_value[0] != '0';
  double stage_started = revac_labrador_monotonic_seconds();
  static const uint8_t reuse_magic[8] = {
      'R', 'V', '8', '8', 'R', 'U', '0', '1'};

  if (proof == NULL || public_state == NULL || uav == NULL ||
      link_randomness == NULL ||
      public_state->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      public_state->acc96_ctx == NULL || public_state->acc96_pp == NULL ||
      public_state->acc96_acc == NULL || uav->acc96_witness == NULL ||
      uav->last_seq != public_state->acc96_seq)
    return 0;
  public_wire = malloc(ACC96_LABRADOR_PUBLIC_BYTES);
  link_wire = malloc(REVAC_Q96_LINK_WIRE_BYTES);
  if (public_wire == NULL || link_wire == NULL ||
      !acc96_labrador_public_encode(public_wire,
                                    ACC96_LABRADOR_PUBLIC_BYTES,
                                    public_state->acc96_public_seed,
                                    public_state->acc96_pp) ||
      !revac_q96_link_wire_encode(
          link_wire, REVAC_Q96_LINK_WIRE_BYTES, public_state, uav->handle,
          link_randomness, proof->handle_commitment))
    goto cleanup;
  if (profile)
    revac_labrador_print_timing("q96 公参标识 + 新 link opening 编码",
                                stage_started);
  if (strcmp(revac_env_or_default("REVAC_LABRADOR_Q96_BRIDGE", "worker"),
             "legacy") != 0) {
    uint8_t packed_accumulator[ACC96_POLY_PACKED_BYTES];

    /* Fast steady-state path: the mode-0600 worker already owns the last
       private witness.  Bind reuse to the exact public state before sending
       only the fresh q19 same-handle opening. */
    stage_started = revac_labrador_monotonic_seconds();
    if (!disable_witness_reuse &&
        acc96_poly_pack(packed_accumulator, public_state->acc96_acc) ==
            ACC96_POLY_PACKED_BYTES &&
        revac_worker_buf_append(&dynamic, reuse_magic,
                                sizeof(reuse_magic)) &&
        revac_worker_buf_u32(&dynamic, public_state->acc96_seq) &&
        revac_worker_buf_u32(&dynamic, uav->handle) &&
        revac_worker_buf_append(&dynamic, packed_accumulator,
                                sizeof(packed_accumulator)) &&
        revac_worker_buf_append(&dynamic, link_wire,
                                REVAC_Q96_LINK_WIRE_BYTES) &&
        revac_q96_worker_request(
            'R', public_wire, ACC96_LABRADOR_PUBLIC_BYTES, &dynamic,
            &proof_bytes, &proof_len, &worker_cache_hit) && proof_len != 0) {
      printf("[q96 Labrador worker] 私密见证缓存：命中\n");
      ok = revac_show_proof_set_acc_zk(
          proof, REVAC_ACC_ZK_Q96_LABRADOR_MEMVER, proof_bytes, proof_len);
      if (profile)
        revac_labrador_print_timing("q96 缓存见证请求 + 新 proof",
                                    stage_started);
      goto cleanup;
    }
    free(proof_bytes);
    proof_bytes = NULL;
    proof_len = 0;
    revac_worker_buf_clear(&dynamic);

    /* Cold/restart/update fallback: reconstruct the exact lift and install a
       complete private witness.  A successful proof refreshes worker reuse. */
    stage_started = revac_labrador_monotonic_seconds();
    lift = malloc(sizeof(*lift));
    prover_wire = malloc(ACC96_LABRADOR_PROVER_BYTES);
    if (lift == NULL || prover_wire == NULL ||
        !acc96_memver_exact_lift(lift, &lift_bits, public_state->acc96_pp,
                                 public_state->acc96_acc, uav->handle,
                                 uav->acc96_witness,
                                 public_state->acc96_ctx) ||
        !acc96_labrador_prover_encode(
            prover_wire, ACC96_LABRADOR_PROVER_BYTES,
            public_state->acc96_seq, uav->handle, public_state->acc96_acc,
            uav->acc96_witness, lift, public_state->acc96_ctx))
      goto cleanup;
    if (profile)
      revac_labrador_print_timing("q96 exact-lift + 完整见证编码",
                                  stage_started);
    stage_started = revac_labrador_monotonic_seconds();
    if (!revac_worker_buf_append(&dynamic, prover_wire,
                                 ACC96_LABRADOR_PROVER_BYTES) ||
        !revac_worker_buf_append(&dynamic, link_wire,
                                 REVAC_Q96_LINK_WIRE_BYTES) ||
        !revac_q96_worker_request(
            'P', public_wire, ACC96_LABRADOR_PUBLIC_BYTES, &dynamic,
            &proof_bytes, &proof_len, &worker_cache_hit) || proof_len == 0)
      goto cleanup;
    printf("[q96 Labrador worker] statement/Pack 参数缓存：%s\n",
           worker_cache_hit ? "命中" : "首次构造");
    ok = revac_show_proof_set_acc_zk(
        proof, REVAC_ACC_ZK_Q96_LABRADOR_MEMVER, proof_bytes, proof_len);
    if (profile)
      revac_labrador_print_timing("q96 完整见证请求 + 新 proof",
                                  stage_started);
    goto cleanup;
  }
  lift = malloc(sizeof(*lift));
  prover_wire = malloc(ACC96_LABRADOR_PROVER_BYTES);
  if (lift == NULL || prover_wire == NULL ||
      !acc96_memver_exact_lift(lift, &lift_bits, public_state->acc96_pp,
                               public_state->acc96_acc, uav->handle,
                               uav->acc96_witness,
                               public_state->acc96_ctx) ||
      !acc96_labrador_prover_encode(
          prover_wire, ACC96_LABRADOR_PROVER_BYTES,
          public_state->acc96_seq, uav->handle, public_state->acc96_acc,
          uav->acc96_witness, lift, public_state->acc96_ctx))
    goto cleanup;
  if (!write_bytes_temp(public_path, public_wire,
                        ACC96_LABRADOR_PUBLIC_BYTES) ||
      !write_bytes_temp(prover_path, prover_wire,
                        ACC96_LABRADOR_PROVER_BYTES) ||
      !write_bytes_temp(link_path, link_wire, REVAC_Q96_LINK_WIRE_BYTES))
    goto cleanup;
  proof_fd = make_temp_path(proof_path);
  if (proof_fd < 0) goto cleanup;
  close(proof_fd);
  proof_fd = -1;
  if (!run_labrador_q96_prove(public_path, prover_path, link_path,
                              proof_path) ||
      !read_file(proof_path, &proof_bytes, &proof_len) || proof_len == 0)
    goto cleanup;
  ok = revac_show_proof_set_acc_zk(
      proof, REVAC_ACC_ZK_Q96_LABRADOR_MEMVER, proof_bytes, proof_len);

cleanup:
  if (proof_fd >= 0) close(proof_fd);
  free(proof_bytes);
  free(link_wire);
  free(prover_wire);
  free(public_wire);
  free(lift);
  revac_worker_buf_clear(&dynamic);
  if (public_path[0]) unlink(public_path);
  if (prover_path[0]) unlink(prover_path);
  if (link_path[0]) unlink(link_path);
  if (proof_path[0]) unlink(proof_path);
  return ok;
}

int revac_labrador_q96_memver_verify(
    const revac_public_state_t *public_state,
    const revac_show_context_t *ctx,
    const revac_show_proof_t *proof) {
  uint8_t *public_wire = NULL, *instance_wire = NULL, *raw_proof = NULL;
  const uint8_t *proof_bytes;
  size_t proof_len;
  revac_worker_buf_t dynamic = {0};
  uint8_t *worker_reply = NULL;
  size_t worker_reply_len = 0;
  int worker_cache_hit = 0;
  char public_path[64] = {0}, instance_path[64] = {0}, proof_path[64] = {0};
  int ok = 0;

  if (public_state == NULL || ctx == NULL || proof == NULL ||
      public_state->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      public_state->acc96_pp == NULL || public_state->acc96_acc == NULL ||
      ctx->acc_profile != REVAC_ACC_PROFILE_Q96 ||
      ctx->acc_seq != public_state->acc96_seq ||
      !revac_show_context_matches_public(ctx, public_state) ||
      proof->acc_zk_proof == NULL || proof->acc_zk_proof_len == 0 ||
      (proof->acc_zk_kind != REVAC_ACC_ZK_Q96_LABRADOR_MEMVER &&
       proof->acc_zk_kind != REVAC_ACC_ZK_Q96_LABRADOR_MEMVER_DEFLATE))
    return 0;
  proof_bytes = proof->acc_zk_proof;
  proof_len = proof->acc_zk_proof_len;
  if (proof->acc_zk_kind == REVAC_ACC_ZK_Q96_LABRADOR_MEMVER_DEFLATE) {
    if (!revac_deflate_unpack(&raw_proof, &proof_len, proof_bytes, proof_len))
      return 0;
    proof_bytes = raw_proof;
  }
  public_wire = malloc(ACC96_LABRADOR_PUBLIC_BYTES);
  instance_wire = malloc(REVAC_Q96_INSTANCE_WIRE_BYTES);
  if (public_wire == NULL || instance_wire == NULL ||
      !acc96_labrador_public_encode(public_wire,
                                    ACC96_LABRADOR_PUBLIC_BYTES,
                                    public_state->acc96_public_seed,
                                    public_state->acc96_pp) ||
      !revac_q96_instance_wire_encode(
          instance_wire, REVAC_Q96_INSTANCE_WIRE_BYTES, public_state,
          proof->handle_commitment))
    goto cleanup;
  if (strcmp(revac_env_or_default("REVAC_LABRADOR_Q96_BRIDGE", "worker"),
             "legacy") != 0) {
    if (proof_len > UINT32_MAX ||
        !revac_worker_buf_append(&dynamic, instance_wire,
                                 REVAC_Q96_INSTANCE_WIRE_BYTES) ||
        !revac_worker_buf_u32(&dynamic, (uint32_t)proof_len) ||
        !revac_worker_buf_append(&dynamic, proof_bytes, proof_len) ||
        !revac_q96_worker_request(
            'V', public_wire, ACC96_LABRADOR_PUBLIC_BYTES, &dynamic,
            &worker_reply, &worker_reply_len, &worker_cache_hit) ||
        worker_reply_len != 1)
      goto cleanup;
    printf("[q96 Labrador worker] GCS statement 缓存：%s\n",
           worker_cache_hit ? "命中" : "首次构造");
    ok = worker_reply[0] == 1;
    goto cleanup;
  }
  if (!write_bytes_temp(public_path, public_wire,
                        ACC96_LABRADOR_PUBLIC_BYTES) ||
      !write_bytes_temp(instance_path, instance_wire,
                        REVAC_Q96_INSTANCE_WIRE_BYTES) ||
      !write_bytes_temp(proof_path, proof_bytes, proof_len))
    goto cleanup;
  ok = run_labrador_q96_verify(public_path, instance_path, proof_path);

cleanup:
  free(raw_proof);
  free(worker_reply);
  revac_worker_buf_clear(&dynamic);
  free(instance_wire);
  free(public_wire);
  if (public_path[0]) unlink(public_path);
  if (instance_path[0]) unlink(instance_path);
  if (proof_path[0]) unlink(proof_path);
  return ok;
}
