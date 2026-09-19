#include "acc_genpp.h"
#include "poly_q_sampling.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

static uint64_t acc_lift_public_bound(int *fits_int64) {
  const __uint128_t coeff_bound = (__uint128_t)(PARAM_Q / 2);
  const __uint128_t conv_bound =
      (__uint128_t)PARAM_N * coeff_bound * coeff_bound;
  const __uint128_t matrix_products =
      (__uint128_t)PARAM_D +
      (__uint128_t)ACC_GENPP_TRUNC_K_L * PARAM_D +
      (__uint128_t)ACC_GENPP_HANDLE_BITS * ACC_GENPP_GADGET_DIM * PARAM_D +
      (__uint128_t)ACC_GENPP_GADGET_DIM * PARAM_D;
  __uint128_t trunc_scale_sum = 0;
  __uint128_t gadget_scale_sum = 0;
  __uint128_t scale = ACC_GENPP_TRUNC_B_POW_L;
  __uint128_t power = 1;
  __uint128_t bound;

  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    trunc_scale_sum += scale;
    scale *= ACC_GENPP_TRUNC_B;
  }
  for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
    gadget_scale_sum += power;
    power *= ACC_GENPP_GADGET_BASE;
  }

  /* Includes left, the public accumulator coefficient, every matrix-vector
     convolution, every tag*v2 term, and every x_h*G_j*right term. */
  bound = 2 * coeff_bound + matrix_products * conv_bound +
          trunc_scale_sum * conv_bound +
          (__uint128_t)ACC_GENPP_HANDLE_BITS * gadget_scale_sum * coeff_bound;
  *fits_int64 = bound <= (__uint128_t)INT64_MAX;
  return bound > (__uint128_t)UINT64_MAX ? UINT64_MAX : (uint64_t)bound;
}

static void acc_lift_i64_poly_muladd(
    int64_t out[PARAM_N],
    const poly_q lhs,
    const poly_q rhs,
    int64_t scale) {
  int64_t lhs_coeffs[PARAM_N];
  int64_t rhs_coeffs[PARAM_N];

  for (size_t i = 0; i < PARAM_N; i++) {
    lhs_coeffs[i] = poly_q_get_coeff_centered(lhs, i);
    rhs_coeffs[i] = poly_q_get_coeff_centered(rhs, i);
  }
  for (size_t i = 0; i < PARAM_N; i++) {
    const int64_t scaled_a = lhs_coeffs[i] * scale;
    if (scaled_a == 0) {
      continue;
    }
    for (size_t j = 0; j < PARAM_N - i; j++) {
      out[i + j] += scaled_a * rhs_coeffs[j];
    }
    for (size_t j = PARAM_N - i; j < PARAM_N; j++) {
      out[i + j - PARAM_N] -= scaled_a * rhs_coeffs[j];
    }
  }
}

static void acc_lift_i64_mat_vec_muladd(
    int64_t out[PARAM_N],
    const poly_q_mat_d_d matrix,
    size_t row,
    const poly_q_vec_d vector,
    int64_t scale) {
  for (size_t column = 0; column < PARAM_D; column++) {
    acc_lift_i64_poly_muladd(out, matrix->rows[row]->entries[column],
                             vector->entries[column], scale);
  }
}

static int acc_lift_int64(
    int64_t lift[PARAM_D][PARAM_N],
    const acc_genpp_pp_t *pp,
    const poly_q_vec_d acc,
    uint32_t handle,
    const acc_genpp_wit_t *witness) {
  int64_t expression[PARAM_D][PARAM_N];
  poly_q_mat_d_d Atrap;
  acc_genpp_gmat_t refresh_matrix;
  int row_ok[PARAM_D];
  const char *parallel_env = getenv("REVAC_LIFT_PARALLEL");
  const int parallel_rows =
      parallel_env == NULL || parallel_env[0] == '\0' ||
      parallel_env[0] != '0';

  poly_q_mat_d_d_init(Atrap);
  acc_genpp_gmat_init(&refresh_matrix);
  poly_q_mat_d_d_uniform(Atrap, pp->A.seed, DOMAIN_SEPARATOR_A, 0);
  if (pp->refresh_bits != 0) {
    acc_genpp_refresh_matrix_compose(&refresh_matrix, pp,
                                     witness->refresh_r);
  }

  /* Every output row is independent.  Keeping the complete row pipeline in
     one worker avoids repeatedly creating OpenMP teams in the handle/gadget
     loops and retains the row accumulator in the same core's caches. */
#ifdef _OPENMP
#pragma omp parallel for num_threads(PARAM_D) if(parallel_rows) schedule(static)
#endif
  for (size_t row = 0; row < PARAM_D; row++) {
    int64_t trunc_scale = ACC_GENPP_TRUNC_B_POW_L;
    row_ok[row] = 1;
    for (size_t k = 0; k < PARAM_N; k++) {
      expression[row][k] =
          poly_q_get_coeff_centered(witness->left->entries[row], k);
    }
    acc_lift_i64_mat_vec_muladd(expression[row], Atrap, row,
                                witness->v12, 1);
    for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
      acc_lift_i64_mat_vec_muladd(expression[row], pp->trunc_B[i], row,
                                  witness->v2[i], -1);
      acc_lift_i64_poly_muladd(expression[row], witness->tag,
                               witness->v2[i]->entries[row], trunc_scale);
      trunc_scale *= ACC_GENPP_TRUNC_B;
    }

    for (size_t h = 0; h < ACC_GENPP_HANDLE_BITS; h++) {
      const int handle_bit = (int)((handle >> h) & 1u);
      for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
        const int64_t gadget_power = (int64_t)acc_genpp_gadget_power(j);
        acc_lift_i64_mat_vec_muladd(expression[row],
                                    pp->Ba.bit[h].blocks[j], row,
                                    witness->right[h][j], 1);
        if (handle_bit) {
          for (size_t k = 0; k < PARAM_N; k++) {
            expression[row][k] -=
                gadget_power * poly_q_get_coeff_centered(
                                   witness->right[h][j]->entries[row], k);
          }
        }
      }
    }

    if (pp->refresh_bits != 0) {
      for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
        acc_lift_i64_mat_vec_muladd(expression[row],
                                    refresh_matrix.blocks[j], row,
                                    witness->refresh_z[j], 1);
      }
    }
    for (size_t k = 0; k < PARAM_N; k++) {
      const int64_t difference =
          poly_q_get_coeff_centered(acc->entries[row], k) -
          expression[row][k];
      if (difference % PARAM_Q != 0) {
        row_ok[row] = 0;
        break;
      }
      lift[row][k] = difference / PARAM_Q;
    }
  }
  acc_genpp_gmat_clear(&refresh_matrix);
  poly_q_mat_d_d_clear(Atrap);
  for (size_t row = 0; row < PARAM_D; row++) {
    if (!row_ok[row]) {
      return 0;
    }
  }
  return 1;
}

static void acc_lift_fmpz_poly_from_poly_q(fmpz_poly_t out,
                                           const poly_q input) {
  fmpz_poly_zero(out);
  for (size_t i = 0; i < PARAM_N; i++) {
    fmpz_poly_set_coeff_si(out, (slong)i,
                           (slong)poly_q_get_coeff_centered(input, i));
  }
}

static void acc_lift_fmpz_poly_muladd(
    fmpz_t out[PARAM_N],
    fmpz_poly_t lhs,
    fmpz_poly_t rhs,
    fmpz_poly_t product,
    fmpz_t low,
    fmpz_t high,
    fmpz_t term,
    const poly_q lhs_q,
    const poly_q rhs_q,
    slong scale) {
  acc_lift_fmpz_poly_from_poly_q(lhs, lhs_q);
  acc_lift_fmpz_poly_from_poly_q(rhs, rhs_q);
  fmpz_poly_mul(product, lhs, rhs);
  for (size_t k = 0; k < PARAM_N; k++) {
    fmpz_poly_get_coeff_fmpz(low, product, (slong)k);
    fmpz_poly_get_coeff_fmpz(high, product, (slong)(k + PARAM_N));
    fmpz_sub(term, low, high);
    fmpz_addmul_si(out[k], term, scale);
  }
}

static void acc_lift_fmpz_mat_vec_muladd(
    fmpz_t out[PARAM_N],
    fmpz_poly_t lhs,
    fmpz_poly_t rhs,
    fmpz_poly_t product,
    fmpz_t low,
    fmpz_t high,
    fmpz_t term,
    const poly_q_mat_d_d matrix,
    size_t row,
    const poly_q_vec_d vector,
    slong scale) {
  for (size_t column = 0; column < PARAM_D; column++) {
    acc_lift_fmpz_poly_muladd(out, lhs, rhs, product, low, high, term,
                              matrix->rows[row]->entries[column],
                              vector->entries[column], scale);
  }
}

static int acc_lift_flint(
    int64_t lift[PARAM_D][PARAM_N],
    const acc_genpp_pp_t *pp,
    const poly_q_vec_d acc,
    uint32_t handle,
    const acc_genpp_wit_t *witness) {
  fmpz_t expression[PARAM_D][PARAM_N];
  fmpz_poly_t lhs, rhs, product;
  fmpz_t low, high, term, difference;
  poly_q_mat_d_d Atrap;
  slong trunc_scale = ACC_GENPP_TRUNC_B_POW_L;
  int ok = 0;
  acc_genpp_gmat_t refresh_matrix;

  fmpz_poly_init(lhs);
  fmpz_poly_init(rhs);
  fmpz_poly_init(product);
  fmpz_init(low);
  fmpz_init(high);
  fmpz_init(term);
  fmpz_init(difference);
  poly_q_mat_d_d_init(Atrap);
  acc_genpp_gmat_init(&refresh_matrix);
  poly_q_mat_d_d_uniform(Atrap, pp->A.seed, DOMAIN_SEPARATOR_A, 0);
  for (size_t row = 0; row < PARAM_D; row++) {
    for (size_t k = 0; k < PARAM_N; k++) {
      fmpz_init(expression[row][k]);
      fmpz_set_si(expression[row][k],
                  (slong)poly_q_get_coeff_centered(
                      witness->left->entries[row], k));
    }
    acc_lift_fmpz_mat_vec_muladd(expression[row], lhs, rhs, product,
                                 low, high, term, Atrap, row,
                                 witness->v12, 1);
  }

  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    for (size_t row = 0; row < PARAM_D; row++) {
      acc_lift_fmpz_mat_vec_muladd(expression[row], lhs, rhs, product,
                                   low, high, term, pp->trunc_B[i], row,
                                   witness->v2[i], -1);
      acc_lift_fmpz_poly_muladd(expression[row], lhs, rhs, product,
                                low, high, term, witness->tag,
                                witness->v2[i]->entries[row], trunc_scale);
    }
    trunc_scale *= ACC_GENPP_TRUNC_B;
  }

  for (size_t h = 0; h < ACC_GENPP_HANDLE_BITS; h++) {
    const int handle_bit = (int)((handle >> h) & 1u);
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      const slong gadget_power = (slong)acc_genpp_gadget_power(j);
      for (size_t row = 0; row < PARAM_D; row++) {
        acc_lift_fmpz_mat_vec_muladd(expression[row], lhs, rhs, product,
                                     low, high, term,
                                     pp->Ba.bit[h].blocks[j], row,
                                     witness->right[h][j], 1);
        if (handle_bit) {
          for (size_t k = 0; k < PARAM_N; k++) {
            fmpz_set_si(term, (slong)poly_q_get_coeff_centered(
                                  witness->right[h][j]->entries[row], k));
            fmpz_submul_si(expression[row][k], term, gadget_power);
          }
        }
      }
    }
  }


  if (pp->refresh_bits != 0) {
    acc_genpp_refresh_matrix_compose(&refresh_matrix, pp,
                                     witness->refresh_r);
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      for (size_t row = 0; row < PARAM_D; row++) {
        acc_lift_fmpz_mat_vec_muladd(expression[row], lhs, rhs, product,
                                     low, high, term,
                                     refresh_matrix.blocks[j], row,
                                     witness->refresh_z[j], 1);
      }
    }
  }

  for (size_t row = 0; row < PARAM_D; row++) {
    for (size_t k = 0; k < PARAM_N; k++) {
      fmpz_set_si(difference,
                  (slong)poly_q_get_coeff_centered(acc->entries[row], k));
      fmpz_sub(difference, difference, expression[row][k]);
      if (!fmpz_divisible_si(difference, PARAM_Q)) {
        goto cleanup;
      }
      fmpz_divexact_si(difference, difference, PARAM_Q);
      if (!fmpz_fits_si(difference)) {
        goto cleanup;
      }
      lift[row][k] = (int64_t)fmpz_get_si(difference);
    }
  }
  ok = 1;

cleanup:
  for (size_t row = 0; row < PARAM_D; row++) {
    for (size_t k = 0; k < PARAM_N; k++) {
      fmpz_clear(expression[row][k]);
    }
  }
  poly_q_mat_d_d_clear(Atrap);
  acc_genpp_gmat_clear(&refresh_matrix);
  fmpz_clear(difference);
  fmpz_clear(term);
  fmpz_clear(high);
  fmpz_clear(low);
  fmpz_poly_clear(product);
  fmpz_poly_clear(rhs);
  fmpz_poly_clear(lhs);
  return ok;
}

int acc_genpp_memver_exact_lift(
    int64_t lift[PARAM_D][PARAM_N],
    const acc_genpp_pp_t *pp,
    const poly_q_vec_d acc,
    uint32_t handle,
    const acc_genpp_wit_t *witness,
    acc_genpp_lift_mode_t mode,
    acc_genpp_lift_info_t *info) {
  acc_genpp_lift_info_t local_info;
  int64_t fast_lift[PARAM_D][PARAM_N];
  int64_t flint_lift[PARAM_D][PARAM_N];

  memset(&local_info, 0, sizeof(local_info));
  local_info.public_abs_bound =
      acc_lift_public_bound(&local_info.int64_bound_ok);
  if (lift == NULL || pp == NULL || acc == NULL || witness == NULL) {
    return 0;
  }

  switch (mode) {
    case ACC_GENPP_LIFT_AUTO:
      if (local_info.int64_bound_ok) {
        if (!acc_lift_int64(lift, pp, acc, handle, witness)) {
          return 0;
        }
        local_info.backend = ACC_GENPP_LIFT_BACKEND_INT64;
      } else {
        if (!acc_lift_flint(lift, pp, acc, handle, witness)) {
          return 0;
        }
        local_info.backend = ACC_GENPP_LIFT_BACKEND_FLINT;
        local_info.fallback_used = 1;
      }
      break;
    case ACC_GENPP_LIFT_FORCE_INT64:
      if (!local_info.int64_bound_ok ||
          !acc_lift_int64(lift, pp, acc, handle, witness)) {
        return 0;
      }
      local_info.backend = ACC_GENPP_LIFT_BACKEND_INT64;
      break;
    case ACC_GENPP_LIFT_FORCE_FLINT:
      if (!acc_lift_flint(lift, pp, acc, handle, witness)) {
        return 0;
      }
      local_info.backend = ACC_GENPP_LIFT_BACKEND_FLINT;
      local_info.fallback_used = 1;
      break;
    case ACC_GENPP_LIFT_CROSSCHECK:
      if (!local_info.int64_bound_ok ||
          !acc_lift_int64(fast_lift, pp, acc, handle, witness) ||
          !acc_lift_flint(flint_lift, pp, acc, handle, witness)) {
        return 0;
      }
      local_info.crosscheck_performed = 1;
      local_info.crosscheck_passed =
          memcmp(fast_lift, flint_lift, sizeof(fast_lift)) == 0;
      if (!local_info.crosscheck_passed) {
        if (info != NULL) {
          *info = local_info;
        }
        return 0;
      }
      memcpy(lift, fast_lift, sizeof(fast_lift));
      local_info.backend = ACC_GENPP_LIFT_BACKEND_INT64;
      break;
    default:
      return 0;
  }

  if (info != NULL) {
    *info = local_info;
  }
  return 1;
}
