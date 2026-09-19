#include "acc_genpp.h"
#include "acc_truncated_params_generated.h"

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "covariance.h"
#include "fips202.h"
#include "poly_q_sampling.h"
#include "randombytes.h"
#include "random.h"

_Static_assert(ACC_TRUNC_PARAM_Q == PARAM_Q,
               "accumulator truncated-sampler modulus mismatch");
_Static_assert(ACC_TRUNC_PARAM_N == PARAM_N && ACC_TRUNC_PARAM_D == PARAM_D,
               "accumulator truncated-sampler ring mismatch");
_Static_assert(ACC_TRUNC_PARAM_B == ACC_GENPP_TRUNC_B &&
                   ACC_TRUNC_PARAM_K == ACC_GENPP_TRUNC_K &&
                   ACC_TRUNC_PARAM_L == ACC_GENPP_TRUNC_L,
               "accumulator truncated-sampler gadget mismatch");

typedef struct {
  poly_q_mat_d_d blocks[ACC_GENPP_GADGET_DIM][ACC_GENPP_GADGET_DIM];
} acc_genpp_ginv_mat_t;

uint64_t acc_genpp_gadget_power(size_t j) {
  uint64_t power = 1;

  for (size_t i = 0; i < j; i++) {
    power *= ACC_GENPP_GADGET_BASE;
  }
  return power;
}

/*
 * Extract one digit in {-1, 0, 1} and divide the remaining centered integer
 * by three.  Subtracting the selected digit makes the division exact even for
 * negative C operands.
 */
static coeff_q acc_genpp_balanced_ternary_digit(coeff_q *value) {
  coeff_q digit = *value % ACC_GENPP_GADGET_BASE;

  if (digit > 1) {
    digit -= ACC_GENPP_GADGET_BASE;
  } else if (digit < -1) {
    digit += ACC_GENPP_GADGET_BASE;
  }
  *value = (*value - digit) / ACC_GENPP_GADGET_BASE;
  return digit;
}

int acc_genpp_balanced_ternary_decompose(
    coeff_q digits[ACC_GENPP_GADGET_DIM], coeff_q value) {
  if (digits == NULL) {
    return 0;
  }
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    digits[i] = acc_genpp_balanced_ternary_digit(&value);
  }
  return value == 0;
}

static int acc_genpp_handle_bit(uint32_t x, size_t pos) {
  return (int)((x >> pos) & 1u);
}

/*
 * EvalFX visits handle bits from least to most significant and stops after the
 * first mismatch.  A mismatch at bit p consumes T_1,...,T_{p+1}; T_32 is
 * handled as the gadget identity and is never stored.  Equal selectors use
 * every stored tail.  The function is local to the UAV and does not expose the
 * resulting prefix length on the wire.
 */
static size_t acc_genpp_evalfx_last_tail(uint32_t selector, uint32_t point) {
  uint32_t difference = selector ^ point;
  size_t first_difference = 0;

  if (difference == 0) {
    return ACC_GENPP_HANDLE_BITS - 1;
  }
  while ((difference & 1u) == 0u) {
    first_difference++;
    difference >>= 1;
  }
  if (first_difference + 1 >= ACC_GENPP_HANDLE_BITS) {
    return ACC_GENPP_HANDLE_BITS - 1;
  }
  return first_difference + 1;
}

static void acc_genpp_gadget_block(poly_q_mat_d_d out, size_t j) {
  poly_q_mat_d_d_zero(out);
  for (size_t i = 0; i < PARAM_D; i++) {
    poly_q_set_coeff(out->rows[i]->entries[i], 0, (coeff_q)acc_genpp_gadget_power(j));
  }
}

static void acc_genpp_Di(acc_genpp_gmat_t *out, const acc_genpp_Ba_t *B, size_t i, int bit) {
  poly_q_mat_d_d tmp;
  poly_q_mat_d_d_init(tmp);

  for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
    if (bit) {
      poly_q_mat_d_d_set(out->blocks[j], B->bit[i].blocks[j]);
    } else {
      acc_genpp_gadget_block(tmp, j);
      poly_q_mat_d_d_sub(out->blocks[j], tmp, B->bit[i].blocks[j]);
    }
  }

  poly_q_mat_d_d_clear(tmp);
}

static void acc_genpp_ginv_mat_init(acc_genpp_ginv_mat_t *m) {
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_mat_d_d_init(m->blocks[i][j]);
    }
  }
}

static void acc_genpp_ginv_mat_clear(acc_genpp_ginv_mat_t *m) {
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_mat_d_d_clear(m->blocks[i][j]);
    }
  }
}

static void acc_genpp_ginv_vec(poly_q_vec_d out[ACC_GENPP_GADGET_DIM], const poly_q_vec_d in) {
  for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
    poly_q_vec_d_zero(out[j]);
  }

  for (size_t row = 0; row < PARAM_D; row++) {
    for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
      coeff_q c = poly_q_get_coeff_centered(in->entries[row], coeff);
      for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
        coeff_q digit = acc_genpp_balanced_ternary_digit(&c);
        if (digit != 0) {
          poly_q_set_coeff(out[j]->entries[row], coeff, digit);
        }
      }
      assert(c == 0);
    }
  }
}

static void acc_genpp_ginv_mat(acc_genpp_ginv_mat_t *out, const acc_genpp_gmat_t *in) {
  for (size_t r = 0; r < ACC_GENPP_GADGET_DIM; r++) {
    for (size_t c = 0; c < ACC_GENPP_GADGET_DIM; c++) {
      poly_q_mat_d_d_zero(out->blocks[r][c]);
    }
  }

  for (size_t cblock = 0; cblock < ACC_GENPP_GADGET_DIM; cblock++) {
    for (size_t r = 0; r < PARAM_D; r++) {
      for (size_t c = 0; c < PARAM_D; c++) {
        for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
          coeff_q v = poly_q_get_coeff_centered(in->blocks[cblock]->rows[r]->entries[c], coeff);
          for (size_t digit_index = 0; digit_index < ACC_GENPP_GADGET_DIM; digit_index++) {
            coeff_q digit = acc_genpp_balanced_ternary_digit(&v);
            if (digit != 0) {
              poly_q_set_coeff(out->blocks[digit_index][cblock]->rows[r]->entries[c],
                               coeff, digit);
            }
          }
          assert(v == 0);
        }
      }
    }
  }
}

static void acc_genpp_gmat_mul_ginv_serial(
    acc_genpp_gmat_t *out,
    const acc_genpp_gmat_t *lhs,
    const acc_genpp_ginv_mat_t *rhs) {
  poly_q_mat_d_d acc, tmp;
  poly_q_mat_d_d_init(acc);
  poly_q_mat_d_d_init(tmp);

  for (size_t c = 0; c < ACC_GENPP_GADGET_DIM; c++) {
    poly_q_mat_d_d_zero(acc);
    for (size_t k = 0; k < ACC_GENPP_GADGET_DIM; k++) {
      poly_q_mat_d_d_mul_mat_d_d(tmp, lhs->blocks[k], rhs->blocks[k][c]);
      poly_q_mat_d_d_add(acc, acc, tmp);
    }
    poly_q_mat_d_d_set(out->blocks[c], acc);
  }

  poly_q_mat_d_d_clear(tmp);
  poly_q_mat_d_d_clear(acc);
}

size_t acc_genpp_parallel_threads(void) {
  size_t workers = 1;

#ifdef _OPENMP
  const char *configured = getenv("REVAC_ACC_THREADS");
  int runtime_max = omp_get_max_threads();

  if (runtime_max > 1) {
    workers = (size_t)runtime_max;
  }
  if (configured != NULL && configured[0] != '\0') {
    char *end = NULL;
    unsigned long requested = strtoul(configured, &end, 10);

    if (end != configured && *end == '\0' && requested > 0) {
      workers = (size_t)requested;
      if (runtime_max > 0 && workers > (size_t)runtime_max) {
        workers = (size_t)runtime_max;
      }
    }
  }
#endif
  if (workers > ACC_GENPP_GADGET_DIM) {
    workers = ACC_GENPP_GADGET_DIM;
  }
  return workers;
}

/* Compute one output gadget column without the global TMP/ACC scratch used by
   poly_q_mat_d_d_mul_mat_d_d().  Distinct output columns therefore have no
   writable state in common and are safe to evaluate concurrently. */
static void acc_genpp_gmat_mul_ginv_column(
    acc_genpp_gmat_t *out,
    const acc_genpp_gmat_t *lhs,
    const acc_genpp_ginv_mat_t *rhs,
    size_t output_block) {
  poly_q acc, product;

  poly_q_init(acc);
  poly_q_init(product);
  for (size_t row = 0; row < PARAM_D; row++) {
    for (size_t col = 0; col < PARAM_D; col++) {
      poly_q_zero(acc);
      for (size_t gadget = 0; gadget < ACC_GENPP_GADGET_DIM; gadget++) {
        for (size_t inner = 0; inner < PARAM_D; inner++) {
          poly_q_mul(product,
                     lhs->blocks[gadget]->rows[row]->entries[inner],
                     rhs->blocks[gadget][output_block]
                         ->rows[inner]
                         ->entries[col]);
          poly_q_add(acc, acc, product);
        }
      }
      poly_q_set(out->blocks[output_block]->rows[row]->entries[col], acc);
    }
  }
  poly_q_clear(product);
  poly_q_clear(acc);
}

static void acc_genpp_gmat_mul_ginv(acc_genpp_gmat_t *out,
                                    const acc_genpp_gmat_t *lhs,
                                    const acc_genpp_ginv_mat_t *rhs) {
  const size_t workers = acc_genpp_parallel_threads();

  if (workers <= 1) {
    acc_genpp_gmat_mul_ginv_serial(out, lhs, rhs);
    return;
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(workers)
  for (int output_block = 0; output_block < ACC_GENPP_GADGET_DIM;
       output_block++) {
    acc_genpp_gmat_mul_ginv_column(out, lhs, rhs, (size_t)output_block);
  }
#else
  for (size_t output_block = 0; output_block < ACC_GENPP_GADGET_DIM;
       output_block++) {
    acc_genpp_gmat_mul_ginv_column(out, lhs, rhs, output_block);
  }
#endif
}

static void acc_genpp_gmat_mul_gvec(poly_q_vec_d out,
                                    const acc_genpp_gmat_t *lhs,
                                    const poly_q_vec_d rhs[ACC_GENPP_GADGET_DIM]) {
  poly_q_vec_d_zero(out);
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_mat_d_d_muladd_vec_d(out, lhs->blocks[i], rhs[i]);
  }
}

void acc_genpp_refresh_matrix_expand(
    acc_genpp_gmat_t *out, const uint8_t master[ACC_GENPP_SEED_BYTES],
    size_t matrix_index) {
  uint8_t block_seed[ACC_GENPP_SEED_BYTES];
  keccak_state state;

  for (size_t block = 0; block < ACC_GENPP_GADGET_DIM; block++) {
    shake256_init(&state);
    shake256_absorb(&state, master, ACC_GENPP_SEED_BYTES);
    shake256_absorb(&state, (const uint8_t *)&matrix_index,
                    sizeof(matrix_index));
    shake256_absorb(&state, (const uint8_t *)&block, sizeof(block));
    shake256_finalize(&state);
    shake256_squeeze(block_seed, sizeof(block_seed), &state);
    poly_q_mat_d_d_uniform(out->blocks[block], block_seed,
                           DOMAIN_SEPARATOR_DS, 0);
  }
}

void acc_genpp_refresh_matrix_compose(
    acc_genpp_gmat_t *out, const acc_genpp_pp_t *pp,
    const uint8_t r[ACC_GENPP_REFRESH_BYTES_MAX]) {
  acc_genpp_gmat_t expanded;

  acc_genpp_gmat_init(&expanded);
  acc_genpp_refresh_matrix_expand(out, pp->refresh_seed, 0);
  for (size_t i = 0; i < pp->refresh_bits; i++) {
    if (((r[i >> 3] >> (i & 7u)) & 1u) == 0) continue;
    acc_genpp_refresh_matrix_expand(&expanded, pp->refresh_seed, i + 1);
    for (size_t block = 0; block < ACC_GENPP_GADGET_DIM; block++) {
      poly_q_mat_d_d_add(out->blocks[block], out->blocks[block],
                         expanded.blocks[block]);
    }
  }
  acc_genpp_gmat_clear(&expanded);
}

static void acc_genpp_refresh_sample(acc_genpp_wit_t *wit,
                                     const acc_genpp_pp_t *pp,
                                     const uint8_t seed[ACC_GENPP_SEED_BYTES]) {
  keccak_state state;
  uint8_t buf[SHAKE256_RATE];
  size_t pos = sizeof(buf);

  memset(wit->refresh_r, 0, sizeof(wit->refresh_r));
  shake256(wit->refresh_r, (pp->refresh_bits + 7) / 8, seed,
           ACC_GENPP_SEED_BYTES);
  if ((pp->refresh_bits & 7u) != 0) {
    wit->refresh_r[pp->refresh_bits >> 3] &=
        (uint8_t)((1u << (pp->refresh_bits & 7u)) - 1u);
  }
  shake256_init(&state);
  shake256_absorb(&state, seed, ACC_GENPP_SEED_BYTES);
  shake256_absorb(&state, (const uint8_t *)"refresh-z", 9);
  shake256_finalize(&state);
  for (size_t block = 0; block < ACC_GENPP_GADGET_DIM; block++) {
    poly_q_vec_d_zero(wit->refresh_z[block]);
    for (size_t row = 0; row < PARAM_D; row++) {
      for (size_t coeff = 0; coeff < PARAM_N; coeff++) {
        uint8_t sample;
        do {
          if (pos == sizeof(buf)) {
            shake256_squeeze(buf, sizeof(buf), &state);
            pos = 0;
          }
          sample = buf[pos++] & 3u;
        } while (sample == 3u);
        poly_q_set_coeff(wit->refresh_z[block]->entries[row], coeff,
                         (coeff_q)sample - 1);
      }
    }
  }
}

static void acc_genpp_sampleleft_randomized(
    acc_genpp_wit_t *e, const acc_genpp_pp_t *pp, const acc_genpp_td_t *td,
    const poly_q_vec_d u, uint32_t x,
    const uint8_t ra[ACC_GENPP_SEED_BYTES]) {
  acc_genpp_gmat_t Ar;
  poly_q_vec_d contribution, target;

  if (!ACC_GENPP_RANDOMIZED_WITNESS_REFRESH) {
    memset(e->refresh_r, 0, sizeof(e->refresh_r));
    for (size_t block = 0; block < ACC_GENPP_GADGET_DIM; block++) {
      poly_q_vec_d_zero(e->refresh_z[block]);
    }
    acc_genpp_sampleleft(e, pp, td, u, x, ra);
    return;
  }

  acc_genpp_gmat_init(&Ar);
  poly_q_vec_d_init(contribution);
  poly_q_vec_d_init(target);
  acc_genpp_refresh_sample(e, pp, ra);
  acc_genpp_refresh_matrix_compose(&Ar, pp, e->refresh_r);
  acc_genpp_gmat_mul_gvec(contribution, &Ar, e->refresh_z);
  poly_q_vec_d_sub(target, u, contribution);
  acc_genpp_sampleleft(e, pp, td, target, x, ra);
  poly_q_vec_d_clear(target);
  poly_q_vec_d_clear(contribution);
  acc_genpp_gmat_clear(&Ar);
}

static void acc_genpp_sample_small_vec_d(poly_q_vec_d out,
                                         const uint8_t seed[ACC_GENPP_SEED_BYTES],
                                         uint32_t x,
                                         size_t bit,
                                         size_t gadget) {
  keccak_state state;
  uint8_t buf[SHAKE256_RATE];
  size_t pos = sizeof(buf);

  poly_q_vec_d_zero(out);
  shake256_init(&state);
  shake256_absorb(&state, seed, ACC_GENPP_SEED_BYTES);
  shake256_absorb(&state, (const uint8_t *)&x, sizeof(x));
  shake256_absorb(&state, (const uint8_t *)&bit, sizeof(bit));
  shake256_absorb(&state, (const uint8_t *)&gadget, sizeof(gadget));
  shake256_finalize(&state);

  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_N; j++) {
      uint8_t sample;
      do {
        if (pos >= sizeof(buf)) {
          shake256_squeeze(buf, sizeof(buf), &state);
          pos = 0;
        }
        sample = buf[pos++] & 3u;
      } while (sample == 3u);
      poly_q_set_coeff(out->entries[i], j, (coeff_q)sample - 1);
    }
  }
}

static void acc_genpp_mul_Mx_right(poly_q_vec_d out,
                                   const acc_genpp_pp_t *pp,
                                   uint32_t x,
                                   const poly_q_vec_d right[ACC_GENPP_HANDLE_BITS][ACC_GENPP_GADGET_DIM]) {
  poly_q_mat_d_d Gblock, Mblock;

  poly_q_vec_d_zero(out);
  poly_q_mat_d_d_init(Gblock);
  poly_q_mat_d_d_init(Mblock);

  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    int xbit = acc_genpp_handle_bit(x, i);
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_mat_d_d_set(Mblock, pp->Ba.bit[i].blocks[j]);
      if (xbit) {
        acc_genpp_gadget_block(Gblock, j);
        poly_q_mat_d_d_sub(Mblock, Mblock, Gblock);
      }
      poly_q_mat_d_d_muladd_vec_d(out, Mblock, right[i][j]);
    }
  }

  poly_q_mat_d_d_clear(Mblock);
  poly_q_mat_d_d_clear(Gblock);
}

static void acc_trunc_vec_gaussian_sqrt_s4sq_sgsq(poly_q_vec_d out) {
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_N; j++) {
      poly_q_set_coeff(out->entries[i], j, SampleZ(0, ACC_TRUNC_SQRT_S4SQ_SGSQ));
    }
  }
}

static void acc_trunc_compute_covariance(poly_real_mat_2d_2d S,
                                         const poly_q_mat_d_d RRstar[2][2]) {
  poly_real finv, tmp;

  poly_real_init(finv);
  poly_real_init(tmp);

  for (size_t i = 0; i < 2 * PARAM_D; i++) {
    for (size_t j = 0; j < 2 * PARAM_D; j++) {
      poly_real_from_poly_q(S->rows[i]->entries[j],
          RRstar[i / PARAM_D][j / PARAM_D]->rows[i % PARAM_D]->entries[j % PARAM_D]);
      poly_real_mul_scalar(S->rows[i]->entries[j], S->rows[i]->entries[j], ACC_TRUNC_SGINVSQ_S4INVSQ);
    }
  }
  for (size_t i = 0; i < PARAM_D; i++) {
    poly_real_add_constant(S->rows[i]->entries[i], ACC_TRUNC_S1SQ);
  }
  for (size_t i = PARAM_D; i < 2 * PARAM_D; i++) {
    poly_real_add_constant(S->rows[i]->entries[i], ACC_TRUNC_S3SQ);
  }

  for (size_t i = 2 * PARAM_D - 1; i > PARAM_D - 1; i--) {
    poly_real_invert(finv, S->rows[i]->entries[i]);
    for (size_t j = 0; j < i; j++) {
      poly_real_mul(S->rows[j]->entries[i], S->rows[j]->entries[i], finv);
    }
    for (size_t j = 0; j < i; j++) {
      for (size_t k = 0; k < i; k++) {
        poly_real_mul(tmp, S->rows[j]->entries[i], S->rows[i]->entries[k]);
        poly_real_sub(S->rows[j]->entries[k], S->rows[j]->entries[k], tmp);
      }
    }
  }

  poly_real_clear(tmp);
  poly_real_clear(finv);
}

static void acc_trunc_sample_perturb(
    poly_q_vec_d pL[ACC_GENPP_TRUNC_L],
    poly_q_vec_d p12,
    poly_q_vec_d p2[ACC_GENPP_TRUNC_K_L],
    const poly_real m_sG_t_tstar,
    const poly_q_mat_d_d R[2][ACC_GENPP_TRUNC_K_L],
    const poly_real_mat_2d_2d S) {
  poly_real tmp_sub, tmp_mul, tmp_sample;
  poly_real_vec_2d c;
  poly_real S_d[PARAM_D][PARAM_D];

  poly_real_init(tmp_sub);
  poly_real_init(tmp_mul);
  poly_real_init(tmp_sample);
  poly_real_vec_2d_init(c);
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_real_init(S_d[i][j]);
    }
  }

  arb_poly_set(tmp_sub, m_sG_t_tstar);
  poly_real_add_constant(tmp_sub, (coeff_real)ACC_TRUNC_S2SQ);
  arb_poly_zero(tmp_mul);
  for (size_t i = 1; i < ACC_GENPP_TRUNC_L; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_q_samplefz(pL[i]->entries[j], tmp_sub, tmp_mul);
    }
  }

  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    acc_trunc_vec_gaussian_sqrt_s4sq_sgsq(p2[i]);
  }

  poly_q_mat_d_d_mul_vec_d(pL[0], R[0][0], p2[0]);
  poly_q_mat_d_d_mul_vec_d(p12, R[1][0], p2[0]);
  for (size_t i = 1; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_mat_d_d_muladd_vec_d(pL[0], R[0][i], p2[i]);
    poly_q_mat_d_d_muladd_vec_d(p12, R[1][i], p2[i]);
  }

  for (size_t i = 0; i < PARAM_D; i++) {
    poly_real_from_poly_q(c->entries[i], pL[0]->entries[i]);
    poly_real_from_poly_q(c->entries[i + PARAM_D], p12->entries[i]);
    poly_real_mul_scalar(c->entries[i], c->entries[i], ACC_TRUNC_NEGSGSQ_DIV_S4SQ_SGSQ);
    poly_real_mul_scalar(c->entries[i + PARAM_D], c->entries[i + PARAM_D], ACC_TRUNC_NEGSGSQ_DIV_S4SQ_SGSQ);
  }

  for (size_t i = 2 * PARAM_D; i > PARAM_D; i--) {
    poly_q_samplefz(p12->entries[i - 1 - PARAM_D], S->rows[i - 1]->entries[i - 1], c->entries[i - 1]);
    poly_real_from_poly_q(tmp_sample, p12->entries[i - 1 - PARAM_D]);
    poly_real_sub(tmp_sub, tmp_sample, c->entries[i - 1]);
    for (size_t j = 0; j < i - 1; j++) {
      poly_real_mul(tmp_mul, tmp_sub, S->rows[j]->entries[i - 1]);
      poly_real_add(c->entries[j], c->entries[j], tmp_mul);
    }
  }

  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      arb_poly_set(S_d[i][j], S->rows[i]->entries[j]);
    }
    poly_real_add(S_d[i][i], S_d[i][i], m_sG_t_tstar);
  }

  for (size_t i = PARAM_D; i > 0; i--) {
    poly_q_samplefz(pL[0]->entries[i - 1], S_d[i - 1][i - 1], c->entries[i - 1]);
    poly_real_invert(tmp_sub, S_d[i - 1][i - 1]);
    for (size_t j = 0; j < i - 1; j++) {
      poly_real_mul(S_d[j][i - 1], S_d[j][i - 1], tmp_sub);
    }
    for (size_t j = 0; j < i - 1; j++) {
      for (size_t k = 0; k < i - 1; k++) {
        poly_real_mul(tmp_mul, S_d[j][i - 1], S_d[i - 1][k]);
        poly_real_sub(S_d[j][k], S_d[j][k], tmp_mul);
      }
    }
    poly_real_from_poly_q(tmp_sample, pL[0]->entries[i - 1]);
    poly_real_sub(tmp_sub, tmp_sample, c->entries[i - 1]);
    for (size_t j = 0; j < i - 1; j++) {
      poly_real_mul(tmp_mul, tmp_sub, S_d[j][i - 1]);
      poly_real_add(c->entries[j], c->entries[j], tmp_mul);
    }
  }

  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_real_clear(S_d[i][j]);
    }
  }
  poly_real_vec_2d_clear(c);
  poly_real_clear(tmp_sample);
  poly_real_clear(tmp_mul);
  poly_real_clear(tmp_sub);
}

static void acc_trunc_sample_klein(poly_q_vec_d v[ACC_GENPP_TRUNC_K], const poly_q_vec_d w) {
  for (int64_t i = PARAM_N * PARAM_D * ACC_GENPP_TRUNC_K; i > 0; i--) {
    int64_t i_1 = (i - 1) / (PARAM_N * PARAM_D);
    int64_t i_2 = (i - 1) % (PARAM_N * PARAM_D);
    int64_t i_2_1 = i_2 / PARAM_N;
    int64_t i_2_2 = i_2 % PARAM_N;
    coeff_q w_coeff = poly_q_get_coeff(w->entries[i_2_1], (size_t)i_2_2);
    double di;

    if (i == PARAM_N * PARAM_D * ACC_GENPP_TRUNC_K) {
      di = acc_trunc_neg_scaled_gadget_gso[0][i_1] * ((double)w_coeff);
    } else {
      coeff_q v_coeff = poly_q_get_coeff_centered(v[0]->entries[i_2_1], (size_t)i_2_2);
      di = acc_trunc_neg_scaled_gadget_gso[0][i_1] * ((double)(w_coeff + v_coeff));
      for (size_t j = 1; j < ACC_GENPP_TRUNC_K; j++) {
        v_coeff = poly_q_get_coeff_centered(v[j]->entries[i_2_1], (size_t)i_2_2);
        di += acc_trunc_neg_scaled_gadget_gso[j][i_1] * ((double)v_coeff);
      }
    }

    int64_t zi = SampleZ(di, acc_trunc_klein_widths[i_1]);
    for (size_t j = 0; j < ACC_GENPP_TRUNC_K; j++) {
      coeff_q v_coeff = poly_q_get_coeff(v[j]->entries[i_2_1], (size_t)i_2_2);
      v_coeff += zi * acc_trunc_gadget_basis[j][i_1];
      poly_q_set_coeff(v[j]->entries[i_2_1], (size_t)i_2_2, v_coeff);
    }
  }
}

static void acc_trunc_sampler(
    poly_q_vec_d v11,
    poly_q_vec_d v12,
    poly_q_vec_d v2[ACC_GENPP_TRUNC_K_L],
    const poly_q_mat_d_d R[2][ACC_GENPP_TRUNC_K_L],
    const poly_q_mat_d_d A,
    const poly_q_mat_d_d B[ACC_GENPP_TRUNC_K_L],
    const poly_q_vec_d u,
    const poly_q tag,
    const poly_real_mat_2d_2d S) {
  poly_q_vec_d pL[ACC_GENPP_TRUNC_L], p12, p2[ACC_GENPP_TRUNC_K_L];
  poly_q_vec_d w, tmp, y[ACC_GENPP_TRUNC_K];
  poly_q taginv;
  poly_real m_sG_t_tstar;
  int64_t bexpi;

  poly_q_init(taginv);
  poly_real_init(m_sG_t_tstar);
  poly_q_vec_d_init(w);
  poly_q_vec_d_init(tmp);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_L; i++) {
    poly_q_vec_d_init(pL[i]);
    poly_q_vec_d_init(y[i]);
  }
  poly_q_vec_d_init(p12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_vec_d_init(p2[i]);
    poly_q_vec_d_init(y[ACC_GENPP_TRUNC_L + i]);
  }

  poly_q_conjugate(taginv, tag);
  poly_q_mul(taginv, taginv, tag);
  poly_real_from_poly_q(m_sG_t_tstar, taginv);
  poly_real_mul_scalar(m_sG_t_tstar, m_sG_t_tstar, (coeff_real)ACC_TRUNC_NEG_SGSQ);

  acc_trunc_sample_perturb(pL, p12, p2, m_sG_t_tstar, R, S);

  bexpi = ACC_GENPP_TRUNC_B;
  for (size_t i = 1; i < ACC_GENPP_TRUNC_L; i++) {
    poly_q_vec_d_mul_scalar(tmp, pL[i], bexpi);
    poly_q_vec_d_add(pL[0], pL[0], tmp);
    bexpi *= ACC_GENPP_TRUNC_B;
  }

  poly_q_invert(taginv, tag);
  poly_q_vec_d_sub(w, u, pL[0]);
  poly_q_mat_d_d_mulsub_vec_d(w, A, p12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_mat_d_d_muladd_vec_d(w, B[i], p2[i]);
  }
  poly_q_vec_d_mul_poly(w, w, taginv);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_vec_d_mul_scalar(tmp, p2[i], bexpi);
    poly_q_vec_d_sub(w, w, tmp);
    bexpi *= ACC_GENPP_TRUNC_B;
  }

  acc_trunc_sample_klein(y, w);
  poly_q_vec_d_add(y[0], y[0], w);

  poly_q_vec_d_mul_poly(y[0], y[0], tag);
  poly_q_vec_d_add(pL[0], pL[0], y[0]);
  bexpi = ACC_GENPP_TRUNC_B;
  for (size_t i = 1; i < ACC_GENPP_TRUNC_L; i++) {
    poly_q_vec_d_mul_poly(y[i], y[i], tag);
    poly_q_vec_d_mul_scalar(y[i], y[i], bexpi);
    poly_q_vec_d_add(pL[0], pL[0], y[i]);
    bexpi *= ACC_GENPP_TRUNC_B;
  }

  poly_q_mat_d_d_mul_vec_d(v11, R[0][0], y[ACC_GENPP_TRUNC_L]);
  for (size_t i = 1; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_mat_d_d_muladd_vec_d(v11, R[0][i], y[ACC_GENPP_TRUNC_L + i]);
  }
  poly_q_vec_d_add(v11, v11, pL[0]);

  poly_q_mat_d_d_mul_vec_d(v12, R[1][0], y[ACC_GENPP_TRUNC_L]);
  for (size_t i = 1; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_mat_d_d_muladd_vec_d(v12, R[1][i], y[ACC_GENPP_TRUNC_L + i]);
  }
  poly_q_vec_d_add(v12, v12, p12);

  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_vec_d_add(v2[i], p2[i], y[ACC_GENPP_TRUNC_L + i]);
  }

  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_vec_d_clear(p2[i]);
    poly_q_vec_d_clear(y[ACC_GENPP_TRUNC_L + i]);
  }
  poly_q_vec_d_clear(p12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_L; i++) {
    poly_q_vec_d_clear(pL[i]);
    poly_q_vec_d_clear(y[i]);
  }
  poly_q_vec_d_clear(tmp);
  poly_q_vec_d_clear(w);
  poly_real_clear(m_sG_t_tstar);
  poly_q_clear(taginv);
}

static void acc_genpp_sample_Ba(acc_genpp_Ba_t *B, const uint8_t seed[SEED_BYTES]) {
  uint8_t block_seed[SEED_BYTES];
  keccak_state state;

  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      shake256_init(&state);
      shake256_absorb(&state, seed, SEED_BYTES);
      shake256_absorb(&state, (const uint8_t *)&i, sizeof(i));
      shake256_absorb(&state, (const uint8_t *)&j, sizeof(j));
      shake256_finalize(&state);
      shake256_squeeze(block_seed, SEED_BYTES, &state);
      poly_q_mat_d_d_uniform(B->bit[i].blocks[j], block_seed, DOMAIN_SEPARATOR_DS, 0);
    }
  }
}

static int acc_genpp_dict_find(const acc_genpp_state_t *st, uint32_t x) {
  for (size_t i = 0; i < ACC_GENPP_MAX_DICT_ENTRIES; i++) {
    if (st->dict[i].used && st->dict[i].handle == x) {
      return (int)i;
    }
  }
  return -1;
}

int acc_genpp_update_sign(const acc_genpp_update_t *update) {
  if (update == NULL) {
    return 0;
  }
  if (update->sign == 0 || update->sign == 1) {
    return 1;
  }
  return update->sign == -1 ? -1 : 0;
}

int acc_genpp_update_sign_is_explicit(const acc_genpp_update_t *update) {
  return update != NULL && (update->sign == 1 || update->sign == -1);
}

static int acc_genpp_dict_put(acc_genpp_state_t *st, uint32_t x, acc_genpp_op_t op, uint32_t index) {
  int pos = acc_genpp_dict_find(st, x);
  if (pos < 0) {
    for (size_t i = 0; i < ACC_GENPP_MAX_DICT_ENTRIES; i++) {
      if (!st->dict[i].used) {
        pos = (int)i;
        st->dict[i].used = 1;
        st->dict[i].handle = x;
        break;
      }
    }
  }
  if (pos < 0) {
    return 0;
  }
  st->dict[pos].op = op;
  st->dict[pos].index = index;
  if (op != ACC_GENPP_OP_DEL) {
    st->dict[pos].update_sign = 0;
  }
  return 1;
}

static void acc_genpp_bml_push(acc_genpp_state_t *st, const acc_genpp_update_t *upmsg) {
  size_t pos;

  if (st->bml_len < ACC_GENPP_BML_CAPACITY) {
    pos = (st->bml_start + st->bml_len) % ACC_GENPP_BML_CAPACITY;
    st->bml_len++;
  } else {
    pos = st->bml_start;
    st->bml_start = (st->bml_start + 1) % ACC_GENPP_BML_CAPACITY;
  }
  st->bml[pos] = *upmsg;
}

void acc_genpp_gmat_init(acc_genpp_gmat_t *m) {
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_mat_d_d_init(m->blocks[i]);
  }
}

void acc_genpp_gmat_clear(acc_genpp_gmat_t *m) {
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_mat_d_d_clear(m->blocks[i]);
  }
}

void acc_genpp_Ba_init(acc_genpp_Ba_t *B) {
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    acc_genpp_gmat_init(&B->bit[i]);
  }
}

void acc_genpp_Ba_clear(acc_genpp_Ba_t *B) {
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    acc_genpp_gmat_clear(&B->bit[i]);
  }
}

void acc_genpp_state_init(acc_genpp_state_t *st) {
  memset(st->ka, 0, sizeof(st->ka));
  st->ctr = 0;
  poly_q_vec_d_init(st->acc0);
  memset(st->dict, 0, sizeof(st->dict));
  memset(st->bml, 0, sizeof(st->bml));
  st->bml_start = 0;
  st->bml_len = 0;
}

void acc_genpp_state_clear(acc_genpp_state_t *st) {
  poly_q_vec_d_clear(st->acc0);
}

void acc_genpp_pp_init(acc_genpp_pp_t *pp) {
  for (size_t i = 0; i < PARAM_KH; i++) {
    poly_q_mat_d_d_init(pp->A.B[i]);
  }
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_mat_d_d_init(pp->trunc_B[i]);
  }
  memset(pp->A.seed, 0, sizeof(pp->A.seed));
  acc_genpp_Ba_init(&pp->Ba);
  poly_q_vec_d_init(pp->ua);
  memset(pp->refresh_seed, 0, sizeof(pp->refresh_seed));
  pp->refresh_bits = ACC_GENPP_RANDOMIZED_WITNESS_REFRESH
                         ? ACC_GENPP_REFRESH_BITS_MIN
                         : 0;
}

void acc_genpp_pp_clear(acc_genpp_pp_t *pp) {
  for (size_t i = 0; i < PARAM_KH; i++) {
    poly_q_mat_d_d_clear(pp->A.B[i]);
  }
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_mat_d_d_clear(pp->trunc_B[i]);
  }
  acc_genpp_Ba_clear(&pp->Ba);
  poly_q_vec_d_clear(pp->ua);
}

void acc_genpp_td_init(acc_genpp_td_t *td) {
  for (size_t i = 0; i < PARAM_KH; i++) {
    for (size_t j = 0; j < 2; j++) {
      poly_q_mat_d_d_init(td->RA.R[j][i]);
    }
  }
  poly_real_mat_2d_2d_init(td->RA.S);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    for (size_t j = 0; j < 2; j++) {
      poly_q_mat_d_d_init(td->trunc_R[j][i]);
    }
  }
  poly_real_mat_2d_2d_init(td->trunc_S);
}

void acc_genpp_td_clear(acc_genpp_td_t *td) {
  for (size_t i = 0; i < PARAM_KH; i++) {
    for (size_t j = 0; j < 2; j++) {
      poly_q_mat_d_d_clear(td->RA.R[j][i]);
    }
  }
  poly_real_mat_2d_2d_clear(td->RA.S);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    for (size_t j = 0; j < 2; j++) {
      poly_q_mat_d_d_clear(td->trunc_R[j][i]);
    }
  }
  poly_real_mat_2d_2d_clear(td->trunc_S);
}

void acc_genpp_wit_init(acc_genpp_wit_t *wit) {
  poly_q_vec_d_init(wit->left);
  poly_q_vec_d_init(wit->v12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_vec_d_init(wit->v2[i]);
  }
  poly_q_init(wit->tag);
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_vec_d_init(wit->right[i][j]);
    }
  }
  memset(wit->refresh_r, 0, sizeof(wit->refresh_r));
  for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
    poly_q_vec_d_init(wit->refresh_z[j]);
  }
}

void acc_genpp_wit_clear(acc_genpp_wit_t *wit) {
  poly_q_vec_d_clear(wit->left);
  poly_q_vec_d_clear(wit->v12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_vec_d_clear(wit->v2[i]);
  }
  poly_q_clear(wit->tag);
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_vec_d_clear(wit->right[i][j]);
    }
  }
  for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
    poly_q_vec_d_clear(wit->refresh_z[j]);
  }
}

size_t acc_genpp_acc_size_bytes(void) {
  return (PARAM_D * PARAM_N * PARAM_Q_BITLEN + 7) / 8;
}

size_t acc_genpp_wit_compact_size_bytes(void) {
  const size_t vec_bits = PARAM_D * PARAM_N * PARAM_Q_BITLEN;
  const size_t tag_bits = PARAM_N;
  const size_t vec_count = 2 + ACC_GENPP_TRUNC_K_L;
  return (vec_count * vec_bits + tag_bits + 7) / 8;
}

size_t acc_genpp_wit_full_size_bytes(void) {
  const size_t vec_bits = PARAM_D * PARAM_N * PARAM_Q_BITLEN;
  const size_t tag_bits = PARAM_N;
  const size_t vec_count = 2 + ACC_GENPP_TRUNC_K_L +
      ACC_GENPP_HANDLE_BITS * ACC_GENPP_GADGET_DIM +
      (ACC_GENPP_RANDOMIZED_WITNESS_REFRESH ? ACC_GENPP_GADGET_DIM : 0);
  const size_t refresh_bits = ACC_GENPP_RANDOMIZED_WITNESS_REFRESH
      ? ACC_GENPP_REFRESH_BITS_MAX : 0;
  return (vec_count * vec_bits + tag_bits + refresh_bits + 7) / 8;
}

uint64_t acc_genpp_wit_norm2(const acc_genpp_wit_t *wit) {
  uint64_t norm = 0;

  norm += poly_q_vec_d_norm2(wit->left);
  norm += poly_q_vec_d_norm2(wit->v12);
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    norm += poly_q_vec_d_norm2(wit->v2[i]);
  }
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      norm += poly_q_vec_d_norm2(wit->right[i][j]);
    }
  }
  if (ACC_GENPP_RANDOMIZED_WITNESS_REFRESH) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      norm += poly_q_vec_d_norm2(wit->refresh_z[j]);
    }
    for (size_t i = 0; i < ACC_GENPP_REFRESH_BYTES_MAX; i++) {
      norm += (uint64_t)__builtin_popcount((unsigned)wit->refresh_r[i]);
    }
  }
  return norm;
}

void acc_genpp_Fk(uint8_t r[ACC_GENPP_SEED_BYTES],
                  const uint8_t ka[ACC_GENPP_SEED_BYTES],
                  uint32_t x,
                  const poly_q_vec_d acc) {
  keccak_state state;
  coeff_q coeff;

  shake256_init(&state);
  shake256_absorb(&state, ka, ACC_GENPP_SEED_BYTES);
  shake256_absorb(&state, (const uint8_t *)&x, sizeof(x));
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_N; j++) {
      coeff = poly_q_get_coeff(acc->entries[i], j);
      shake256_absorb(&state, (const uint8_t *)&coeff, sizeof(coeff));
    }
  }
  shake256_finalize(&state);
  shake256_squeeze(r, ACC_GENPP_SEED_BYTES, &state);
}

static void acc_genpp_evalf_to(acc_genpp_gmat_t *out,
                               const acc_genpp_Ba_t *B,
                               uint32_t selector,
                               size_t stop_index,
                               acc_genpp_gmat_t *capture[ACC_GENPP_HANDLE_BITS]) {
  acc_genpp_gmat_t cur, D;
  acc_genpp_ginv_mat_t inv;

  acc_genpp_gmat_init(&cur);
  acc_genpp_gmat_init(&D);
  acc_genpp_ginv_mat_init(&inv);

  acc_genpp_Di(&cur, B, ACC_GENPP_HANDLE_BITS - 1,
               acc_genpp_handle_bit(selector, ACC_GENPP_HANDLE_BITS - 1));
  if (capture != NULL && capture[ACC_GENPP_HANDLE_BITS - 1] != NULL) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_mat_d_d_set(capture[ACC_GENPP_HANDLE_BITS - 1]->blocks[j],
                         cur.blocks[j]);
    }
  }

  for (size_t off = 1; off < ACC_GENPP_HANDLE_BITS - stop_index; off++) {
    size_t i = ACC_GENPP_HANDLE_BITS - 1 - off;
    acc_genpp_ginv_mat(&inv, &cur);
    acc_genpp_Di(&D, B, i, acc_genpp_handle_bit(selector, i));
    acc_genpp_gmat_mul_ginv(&cur, &D, &inv);
    if (capture != NULL && capture[i] != NULL) {
      for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
        poly_q_mat_d_d_set(capture[i]->blocks[j], cur.blocks[j]);
      }
    }
  }

  for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
    poly_q_mat_d_d_set(out->blocks[j], cur.blocks[j]);
  }

  acc_genpp_ginv_mat_clear(&inv);
  acc_genpp_gmat_clear(&D);
  acc_genpp_gmat_clear(&cur);
}

void acc_genpp_evalf(acc_genpp_gmat_t *out,
                     const acc_genpp_Ba_t *B,
                     uint32_t selector) {
  acc_genpp_evalf_to(out, B, selector, 0, NULL);
}

void acc_genpp_evalf_tail1(acc_genpp_gmat_t *out,
                           const acc_genpp_Ba_t *B,
                           uint32_t selector) {
  acc_genpp_evalf_to(out, B, selector, 1, NULL);
}

void acc_genpp_evalf_capture(
    acc_genpp_gmat_t *out,
    const acc_genpp_Ba_t *B,
    uint32_t selector,
    acc_genpp_gmat_t *capture[ACC_GENPP_HANDLE_BITS]) {
  acc_genpp_evalf_to(out, B, selector, 0, capture);
}

int acc_genpp_evalfx_apply_tails(
    acc_genpp_wit_t *delta,
    const acc_genpp_Ba_t *B,
    uint32_t selector,
    uint32_t point,
    const poly_q_vec_d u,
    const uint8_t *tail_indices,
    const acc_genpp_gmat_t *tail_values,
    size_t tail_count,
    size_t *local_composites) {
  const size_t last_tail = acc_genpp_evalfx_last_tail(selector, point);
  acc_genpp_gmat_t tail[ACC_GENPP_HANDLE_BITS], cur, D;
  poly_q_vec_d ginv_u[ACC_GENPP_GADGET_DIM];
  poly_q_vec_d tmp_in[ACC_GENPP_GADGET_DIM], tmp_out[ACC_GENPP_GADGET_DIM];
  acc_genpp_ginv_mat_t inv;
  int prefix_equal = 1;
  size_t lower = 0;
  size_t checkpoint = 0;
  size_t composites = 0;

  if ((tail_count != 0 && (tail_indices == NULL || tail_values == NULL)) ||
      tail_count >= ACC_GENPP_HANDLE_BITS) {
    return 0;
  }
  for (size_t i = 0; i < tail_count; i++) {
    if (tail_indices[i] == 0 ||
        tail_indices[i] >= ACC_GENPP_HANDLE_BITS - 1 ||
        (i != 0 && tail_indices[i] <= tail_indices[i - 1])) {
      return 0;
    }
  }
  if (local_composites != NULL) {
    *local_composites = 0;
  }

  /* Retain only the tails consumed by this UAV's handle prefix. */
  for (size_t i = 1; i <= last_tail; i++) {
    acc_genpp_gmat_init(&tail[i]);
  }
  acc_genpp_gmat_init(&cur);
  acc_genpp_gmat_init(&D);
  acc_genpp_ginv_mat_init(&inv);
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_vec_d_init(ginv_u[i]);
    poly_q_vec_d_init(tmp_in[i]);
    poly_q_vec_d_init(tmp_out[i]);
  }

  /*
   * Fill T_1,...,T_last_tail interval by interval.  A published T_high is
   * an anchor for the missing values immediately below it.  T_31=D_31 is the
   * final implicit anchor and costs no matrix composite.
   */
  while (lower < last_tail) {
    size_t high;

    if (checkpoint < tail_count) {
      high = tail_indices[checkpoint];
      for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
        poly_q_mat_d_d_set(cur.blocks[j], tail_values[checkpoint].blocks[j]);
      }
      checkpoint++;
    } else {
      high = ACC_GENPP_HANDLE_BITS - 1;
      acc_genpp_Di(&cur, B, high, acc_genpp_handle_bit(selector, high));
    }

    if (high <= lower) {
      goto fail;
    }
    if (high <= last_tail) {
      for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
        poly_q_mat_d_d_set(tail[high].blocks[j], cur.blocks[j]);
      }
    }
    for (size_t i = high - 1; i > lower; i--) {
      acc_genpp_ginv_mat(&inv, &cur);
      acc_genpp_Di(&D, B, i, acc_genpp_handle_bit(selector, i));
      acc_genpp_gmat_mul_ginv(&cur, &D, &inv);
      composites++;
      if (i <= last_tail) {
        for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
          poly_q_mat_d_d_set(tail[i].blocks[j], cur.blocks[j]);
        }
      }
    }
    lower = high;
  }

  acc_genpp_ginv_vec(ginv_u, u);
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_vec_d_zero(delta->right[i][j]);
    }
  }
  poly_q_vec_d_zero(delta->left);

  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    if (!prefix_equal) {
      break;
    }

    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_vec_d_set(tmp_in[j], ginv_u[j]);
    }

    if (i + 1 < ACC_GENPP_HANDLE_BITS) {
      acc_genpp_ginv_mat(&inv, &tail[i + 1]);
      for (size_t r = 0; r < ACC_GENPP_GADGET_DIM; r++) {
        poly_q_vec_d_zero(tmp_out[r]);
        for (size_t c = 0; c < ACC_GENPP_GADGET_DIM; c++) {
          poly_q_mat_d_d_muladd_vec_d(tmp_out[r], inv.blocks[r][c], tmp_in[c]);
        }
      }
    } else {
      for (size_t c = 0; c < ACC_GENPP_GADGET_DIM; c++) {
        poly_q_vec_d_set(tmp_out[c], tmp_in[c]);
      }
    }

    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      if (acc_genpp_handle_bit(selector, i)) {
        poly_q_vec_d_set(delta->right[i][j], tmp_out[j]);
      } else {
        poly_q_vec_d_neg(delta->right[i][j], tmp_out[j]);
      }
    }

    if (acc_genpp_handle_bit(selector, i) != acc_genpp_handle_bit(point, i)) {
      prefix_equal = 0;
    }
  }

  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_vec_d_clear(tmp_out[i]);
    poly_q_vec_d_clear(tmp_in[i]);
    poly_q_vec_d_clear(ginv_u[i]);
  }
  acc_genpp_ginv_mat_clear(&inv);
  acc_genpp_gmat_clear(&D);
  acc_genpp_gmat_clear(&cur);
  for (size_t i = 1; i <= last_tail; i++) {
    acc_genpp_gmat_clear(&tail[i]);
  }
  if (local_composites != NULL) {
    *local_composites = composites;
  }
  return 1;

fail:
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_vec_d_clear(tmp_out[i]);
    poly_q_vec_d_clear(tmp_in[i]);
    poly_q_vec_d_clear(ginv_u[i]);
  }
  acc_genpp_ginv_mat_clear(&inv);
  acc_genpp_gmat_clear(&D);
  acc_genpp_gmat_clear(&cur);
  for (size_t i = 1; i <= last_tail; i++) {
    acc_genpp_gmat_clear(&tail[i]);
  }
  return 0;
}

void acc_genpp_delta_from_tails(
    poly_q_vec_d delta[][ACC_GENPP_GADGET_DIM], size_t layers,
    const acc_genpp_gmat_t tails[], uint32_t selector,
    const poly_q_vec_d u) {
  poly_q_vec_d ginv_u[ACC_GENPP_GADGET_DIM];
  acc_genpp_ginv_mat_t inv;

  acc_genpp_ginv_mat_init(&inv);
  for (size_t gadget = 0; gadget < ACC_GENPP_GADGET_DIM; gadget++) {
    poly_q_vec_d_init(ginv_u[gadget]);
  }
  acc_genpp_ginv_vec(ginv_u, u);

  for (size_t layer = 0; layer < layers; layer++) {
    acc_genpp_ginv_mat(&inv, &tails[layer]);
    for (size_t row_gadget = 0; row_gadget < ACC_GENPP_GADGET_DIM;
         row_gadget++) {
      poly_q_vec_d_zero(delta[layer][row_gadget]);
      for (size_t col_gadget = 0; col_gadget < ACC_GENPP_GADGET_DIM;
           col_gadget++) {
        poly_q_mat_d_d_muladd_vec_d(delta[layer][row_gadget],
                                    inv.blocks[row_gadget][col_gadget],
                                    ginv_u[col_gadget]);
      }
      if (!acc_genpp_handle_bit(selector, layer)) {
        poly_q_vec_d_neg(delta[layer][row_gadget],
                         delta[layer][row_gadget]);
      }
    }
  }

  for (size_t gadget = 0; gadget < ACC_GENPP_GADGET_DIM; gadget++) {
    poly_q_vec_d_clear(ginv_u[gadget]);
  }
  acc_genpp_ginv_mat_clear(&inv);
}

void acc_genpp_delta4_from_tails(
    poly_q_vec_d delta[ACC_GENPP_DELTA4_LAYERS][ACC_GENPP_GADGET_DIM],
    const acc_genpp_gmat_t tails[ACC_GENPP_DELTA4_LAYERS],
    uint32_t selector,
    const poly_q_vec_d u) {
  acc_genpp_delta_from_tails(delta, ACC_GENPP_DELTA4_LAYERS, tails,
                            selector, u);
}

int acc_genpp_memwitup_delta(
    acc_genpp_wit_t *sx, uint32_t member_handle,
    const acc_genpp_update_t *update,
    const poly_q_vec_d delta[][ACC_GENPP_GADGET_DIM], size_t layers,
    size_t *applied_layers) {
  uint32_t difference;
  size_t first_difference = 0;
  const int sign = acc_genpp_update_sign(update);

  if (applied_layers != NULL) {
    *applied_layers = 0;
  }
  if (sx == NULL || update == NULL || delta == NULL || layers == 0 ||
      layers > ACC_GENPP_HANDLE_BITS ||
      update->op != ACC_GENPP_OP_DEL || sign == 0) {
    return 0;
  }
  difference = update->handle ^ member_handle;
  if (difference == 0) {
    return 0;
  }
  while ((difference & 1U) == 0U &&
         first_difference < ACC_GENPP_HANDLE_BITS) {
    first_difference++;
    difference >>= 1;
  }
  if (first_difference >= layers) {
    return 0;
  }

  for (size_t layer = 0; layer <= first_difference; layer++) {
    for (size_t gadget = 0; gadget < ACC_GENPP_GADGET_DIM; gadget++) {
      if (sign > 0) {
        poly_q_vec_d_add(sx->right[layer][gadget],
                         sx->right[layer][gadget], delta[layer][gadget]);
      } else {
        poly_q_vec_d_sub(sx->right[layer][gadget],
                         sx->right[layer][gadget], delta[layer][gadget]);
      }
    }
  }
  if (applied_layers != NULL) {
    *applied_layers = first_difference + 1;
  }
  return 1;
}

int acc_genpp_memwitup_delta4(
    acc_genpp_wit_t *sx,
    uint32_t member_handle,
    const acc_genpp_update_t *update,
    const poly_q_vec_d delta[ACC_GENPP_DELTA4_LAYERS]
                            [ACC_GENPP_GADGET_DIM],
    size_t *applied_layers) {
  return acc_genpp_memwitup_delta(sx, member_handle, update, delta,
                                  ACC_GENPP_DELTA4_LAYERS, applied_layers);
}

void acc_genpp_evalfx_apply(acc_genpp_wit_t *delta,
                            const acc_genpp_Ba_t *B,
                            uint32_t selector,
                            uint32_t point,
                            const poly_q_vec_d u) {
  (void)acc_genpp_evalfx_apply_tails(delta, B, selector, point, u,
                                     NULL, NULL, 0, NULL);
}

void acc_genpp_sampleleft(acc_genpp_wit_t *e,
                          const acc_genpp_pp_t *pp,
                          const acc_genpp_td_t *td,
                          const poly_q_vec_d u,
                          uint32_t x,
                          const uint8_t ra[ACC_GENPP_SEED_BYTES]) {
  poly_q_mat_d_d A;
  poly_q_vec_d right_contrib, target;

  poly_q_mat_d_d_init(A);
  poly_q_vec_d_init(right_contrib);
  poly_q_vec_d_init(target);

  poly_q_zero(e->tag);
  poly_q_set_coeff(e->tag, 0, 1);
  poly_q_mat_d_d_uniform(A, pp->A.seed, DOMAIN_SEPARATOR_A, 0);

  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      acc_genpp_sample_small_vec_d(e->right[i][j], ra, x, i, j);
    }
  }

  acc_genpp_mul_Mx_right(right_contrib, pp, x, e->right);
  poly_q_vec_d_sub(target, u, right_contrib);
  acc_trunc_sampler(e->left, e->v12, e->v2, td->trunc_R, A, pp->trunc_B, target, e->tag, td->trunc_S);

  poly_q_vec_d_clear(target);
  poly_q_vec_d_clear(right_contrib);
  poly_q_mat_d_d_clear(A);
}

void acc_genpp_gen(acc_genpp_pp_t *pp,
                   acc_genpp_td_t *td,
                   acc_genpp_state_t *st,
                   poly_q_vec_d A0) {
  uint8_t seed[SEED_BYTES], rseed[SEED_BYTES];
  poly_q_mat_d_d A, RRstar[2][2];

  sep_keygen(&pp->A, &td->RA);
  randombytes(seed, sizeof(seed));
  acc_genpp_sample_Ba(&pp->Ba, seed);
  memset(pp->refresh_seed, 0, sizeof(pp->refresh_seed));
  if (ACC_GENPP_RANDOMIZED_WITNESS_REFRESH) {
    randombytes(pp->refresh_seed, sizeof(pp->refresh_seed));
    const char *bits_env = getenv("REVAC_RANDOMIZED_REFRESH_BITS");
    unsigned long bits = bits_env != NULL ? strtoul(bits_env, NULL, 10) :
                                           ACC_GENPP_REFRESH_BITS_MIN;
    pp->refresh_bits = bits == ACC_GENPP_REFRESH_BITS_MAX ?
                       ACC_GENPP_REFRESH_BITS_MAX : ACC_GENPP_REFRESH_BITS_MIN;
  } else {
    pp->refresh_bits = 0;
  }

  poly_q_mat_d_d_init(A);
  poly_q_mat_d_d_init(RRstar[0][0]);
  poly_q_mat_d_d_init(RRstar[0][1]);
  poly_q_mat_d_d_init(RRstar[1][0]);
  poly_q_mat_d_d_init(RRstar[1][1]);

  poly_q_mat_d_d_uniform(A, pp->A.seed, DOMAIN_SEPARATOR_A, 0);
  do {
    uint32_t kappa = 0;
    randombytes(rseed, sizeof(rseed));
    for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
      for (size_t j = 0; j < 2; j++) {
        poly_q_mat_d_d_binomial(td->trunc_R[j][i], rseed, kappa++,
                                DOMAIN_SEPARATOR_R);
      }
    }
  } while (sk_sq_spectral_norm_columns(
      RRstar, ACC_GENPP_TRUNC_K_L, td->trunc_R,
      ACC_TRUNC_R_MAX_SQ_SPECTRAL_NORM));

  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_mat_d_d_mul_mat_d_d(pp->trunc_B[i], A, td->trunc_R[1][i]);
    poly_q_mat_d_d_add(pp->trunc_B[i], pp->trunc_B[i], td->trunc_R[0][i]);
  }
  acc_trunc_compute_covariance(td->trunc_S, RRstar);

  randombytes(st->ka, sizeof(st->ka));
  poly_q_vec_d_uniform(pp->ua, seed, DOMAIN_SEPARATOR_U);
  poly_q_vec_d_uniform(A0, seed, DOMAIN_SEPARATOR_DS);
  poly_q_vec_d_set(st->acc0, A0);
  st->ctr = 0;
  memset(st->dict, 0, sizeof(st->dict));

  poly_q_mat_d_d_clear(RRstar[1][1]);
  poly_q_mat_d_d_clear(RRstar[1][0]);
  poly_q_mat_d_d_clear(RRstar[0][1]);
  poly_q_mat_d_d_clear(RRstar[0][0]);
  poly_q_mat_d_d_clear(A);
}

int acc_genpp_add(poly_q_vec_d A_out,
                  acc_genpp_wit_t *sx,
                  acc_genpp_update_t *upmsg,
                  acc_genpp_state_t *st,
                  const acc_genpp_pp_t *pp,
                  const acc_genpp_td_t *td,
                  const poly_q_vec_d A,
                  uint32_t x) {
  uint8_t ra[ACC_GENPP_SEED_BYTES];
  int pos = acc_genpp_dict_find(st, x);

  if (pos >= 0 && st->dict[pos].op != ACC_GENPP_OP_EMPTY) {
    return 0;
  }

  acc_genpp_Fk(ra, st->ka, x, A);
  acc_genpp_sampleleft_randomized(sx, pp, td, A, x, ra);
  if (!acc_genpp_dict_put(st, x, ACC_GENPP_OP_ADD, st->ctr)) {
    return 0;
  }

  poly_q_vec_d_set(A_out, A);
  upmsg->seq = st->ctr;
  upmsg->op = ACC_GENPP_OP_ADD;
  upmsg->handle = 0;
  upmsg->sign = 0;
  return 1;
}

static int acc_genpp_delete_capture_impl(
    poly_q_vec_d A_out,
    acc_genpp_update_t *upmsg,
    acc_genpp_state_t *st,
    const acc_genpp_pp_t *pp,
    const poly_q_vec_d A,
    uint32_t x,
    const acc_genpp_wit_t *sx,
    int require_member_witness,
    int sign,
    acc_genpp_gmat_t *capture[ACC_GENPP_HANDLE_BITS]) {
  acc_genpp_gmat_t B1x;
  poly_q_vec_d ginv_ua[ACC_GENPP_GADGET_DIM], delta;
  int pos;

  if (st->ctr >= ACC_GENPP_MAX_DELETIONS || (sign != 1 && sign != -1)) {
    return 0;
  }
  pos = acc_genpp_dict_find(st, x);
  if (pos < 0 || st->dict[pos].op != ACC_GENPP_OP_ADD) {
    return 0;
  }
  if (require_member_witness &&
      (sx == NULL || !acc_genpp_memver(pp, A, x, sx))) {
    return 0;
  }

  acc_genpp_gmat_init(&B1x);
  poly_q_vec_d_init(delta);
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_vec_d_init(ginv_ua[i]);
  }

  acc_genpp_evalf_capture(&B1x, &pp->Ba, x, capture);
  acc_genpp_ginv_vec(ginv_ua, pp->ua);
  acc_genpp_gmat_mul_gvec(delta, &B1x, ginv_ua);
  if (sign > 0) {
    poly_q_vec_d_add(A_out, A, delta);
  } else {
    poly_q_vec_d_sub(A_out, A, delta);
  }

  st->dict[pos].op = ACC_GENPP_OP_DEL;
  st->dict[pos].index = st->ctr + 1;
  st->dict[pos].update_sign = (int8_t)sign;
  st->ctr++;
  upmsg->seq = st->ctr;
  upmsg->op = ACC_GENPP_OP_DEL;
  upmsg->handle = x;
  upmsg->sign = (int8_t)sign;
  acc_genpp_bml_push(st, upmsg);

  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_vec_d_clear(ginv_ua[i]);
  }
  poly_q_vec_d_clear(delta);
  acc_genpp_gmat_clear(&B1x);
  return 1;
}

int acc_genpp_delete_capture(
    poly_q_vec_d A_out,
    acc_genpp_update_t *upmsg,
    acc_genpp_state_t *st,
    const acc_genpp_pp_t *pp,
    const acc_genpp_td_t *td,
    const poly_q_vec_d A,
    uint32_t x,
    const acc_genpp_wit_t *sx,
    acc_genpp_gmat_t *capture[ACC_GENPP_HANDLE_BITS]) {
  (void)td;
  return acc_genpp_delete_capture_impl(A_out, upmsg, st, pp, A, x, sx, 1, 1,
                                       capture);
}

int acc_genpp_delete_trusted_capture(
    poly_q_vec_d A_out,
    acc_genpp_update_t *upmsg,
    acc_genpp_state_t *st,
    const acc_genpp_pp_t *pp,
    const poly_q_vec_d A,
    uint32_t x,
    acc_genpp_gmat_t *capture[ACC_GENPP_HANDLE_BITS]) {
  return acc_genpp_delete_capture_impl(A_out, upmsg, st, pp, A, x, NULL, 0, 1,
                                       capture);
}

int acc_genpp_delete_trusted_capture_signed(
    poly_q_vec_d A_out, acc_genpp_update_t *upmsg, acc_genpp_state_t *st,
    const acc_genpp_pp_t *pp, const poly_q_vec_d A, uint32_t x, int sign,
    acc_genpp_gmat_t *capture[ACC_GENPP_HANDLE_BITS]) {
  return acc_genpp_delete_capture_impl(A_out, upmsg, st, pp, A, x, NULL, 0,
                                       sign, capture);
}

int acc_genpp_delete(poly_q_vec_d A_out,
                     acc_genpp_update_t *upmsg,
                     acc_genpp_state_t *st,
                     const acc_genpp_pp_t *pp,
                     const acc_genpp_td_t *td,
                     const poly_q_vec_d A,
                     uint32_t x,
                     const acc_genpp_wit_t *sx) {
  return acc_genpp_delete_capture(A_out, upmsg, st, pp, td, A, x, sx, NULL);
}

int acc_genpp_memwitup_tails(
    acc_genpp_wit_t *sx,
    const acc_genpp_pp_t *pp,
    uint32_t x,
    const acc_genpp_update_t *upmsg,
    const uint8_t *tail_indices,
    const acc_genpp_gmat_t *tail_values,
    size_t tail_count,
    size_t *local_composites) {
  acc_genpp_wit_t delta;
  int ok;
  const int sign = acc_genpp_update_sign(upmsg);

  if (upmsg->op != ACC_GENPP_OP_DEL) {
    if (local_composites != NULL) {
      *local_composites = 0;
    }
    return 1;
  }
  if (sign == 0) {
    return 0;
  }

  acc_genpp_wit_init(&delta);
  ok = acc_genpp_evalfx_apply_tails(&delta, &pp->Ba, upmsg->handle, x,
                                     pp->ua, tail_indices, tail_values,
                                     tail_count, local_composites);
  if (!ok) {
    acc_genpp_wit_clear(&delta);
    return 0;
  }
  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      if (sign > 0) {
        poly_q_vec_d_add(sx->right[i][j], sx->right[i][j],
                         delta.right[i][j]);
      } else {
        poly_q_vec_d_sub(sx->right[i][j], sx->right[i][j],
                         delta.right[i][j]);
      }
    }
  }
  acc_genpp_wit_clear(&delta);
  return 1;
}

void acc_genpp_memwitup(acc_genpp_wit_t *sx,
                        const acc_genpp_pp_t *pp,
                        uint32_t x,
                        const acc_genpp_update_t *upmsg) {
  (void)acc_genpp_memwitup_tails(sx, pp, x, upmsg, NULL, NULL, 0, NULL);
}

int acc_genpp_bml_get_missing(const acc_genpp_state_t *st,
                              uint32_t last_seq,
                              acc_genpp_update_t *out,
                              size_t out_cap,
                              size_t *out_len) {
  size_t count = 0;
  uint32_t oldest_seq;

  if (out_len == NULL || last_seq > st->ctr) {
    return 0;
  }
  *out_len = 0;
  if (last_seq == st->ctr) {
    return 1;
  }
  if (st->bml_len == 0) {
    return 0;
  }

  oldest_seq = st->ctr - (uint32_t)st->bml_len + 1;
  if (last_seq + 1 < oldest_seq) {
    return 0;
  }

  for (size_t i = 0; i < st->bml_len; i++) {
    size_t pos = (st->bml_start + i) % ACC_GENPP_BML_CAPACITY;
    const acc_genpp_update_t *msg = &st->bml[pos];

    if (msg->seq <= last_seq) {
      continue;
    }
    if (count >= out_cap) {
      return 0;
    }
    if (out != NULL) {
      out[count] = *msg;
    }
    count++;
  }
  *out_len = count;
  return 1;
}

int acc_genpp_memver(const acc_genpp_pp_t *pp,
                     const poly_q_vec_d A,
                     uint32_t x,
                     const acc_genpp_wit_t *sx) {
  poly_q_vec_d lhs, refresh_contrib;
  poly_q_mat_d_d Atrap, Gblock, Mblock;
  poly_q tag_times_bexpi;
  int ok;

  poly_q_vec_d_init(lhs);
  poly_q_vec_d_init(refresh_contrib);
  poly_q_mat_d_d_init(Atrap);
  poly_q_mat_d_d_init(Gblock);
  poly_q_mat_d_d_init(Mblock);
  poly_q_init(tag_times_bexpi);

  poly_q_vec_d_set(lhs, sx->left);

  poly_q_mat_d_d_uniform(Atrap, pp->A.seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_mat_d_d_muladd_vec_d(lhs, Atrap, sx->v12);

  coeff_q bexpi = ACC_GENPP_TRUNC_B_POW_L;
  for (size_t i = 0; i < ACC_GENPP_TRUNC_K_L; i++) {
    poly_q_mat_d_d_set(Mblock, pp->trunc_B[i]);
    poly_q_mul_scalar(tag_times_bexpi, sx->tag, bexpi);
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_q_sub(Mblock->rows[j]->entries[j], Mblock->rows[j]->entries[j], tag_times_bexpi);
    }
    poly_q_mat_d_d_mulsub_vec_d(lhs, Mblock, sx->v2[i]);
    bexpi *= ACC_GENPP_TRUNC_B;
  }

  for (size_t i = 0; i < ACC_GENPP_HANDLE_BITS; i++) {
    int xbit = acc_genpp_handle_bit(x, i);
    for (size_t j = 0; j < ACC_GENPP_GADGET_DIM; j++) {
      poly_q_mat_d_d_set(Mblock, pp->Ba.bit[i].blocks[j]);
      if (xbit) {
        acc_genpp_gadget_block(Gblock, j);
        poly_q_mat_d_d_sub(Mblock, Mblock, Gblock);
      }
      poly_q_mat_d_d_muladd_vec_d(lhs, Mblock, sx->right[i][j]);
    }
  }

  if (pp->refresh_bits != 0) {
    acc_genpp_gmat_t Ar;
    acc_genpp_gmat_init(&Ar);
    acc_genpp_refresh_matrix_compose(&Ar, pp, sx->refresh_r);
    acc_genpp_gmat_mul_gvec(refresh_contrib, &Ar, sx->refresh_z);
    poly_q_vec_d_add(lhs, lhs, refresh_contrib);
    acc_genpp_gmat_clear(&Ar);
  }

  ok = poly_q_vec_d_equal(lhs, A) && (acc_genpp_wit_norm2(sx) <= ACC_GENPP_WIT_BETA_SQ);
  poly_q_clear(tag_times_bexpi);
  poly_q_mat_d_d_clear(Mblock);
  poly_q_mat_d_d_clear(Gblock);
  poly_q_mat_d_d_clear(Atrap);
  poly_q_vec_d_clear(lhs);
  poly_q_vec_d_clear(refresh_contrib);
  return ok;
}

int acc_genpp_memwitsync_apply(acc_genpp_wit_t *sx,
                               const acc_genpp_pp_t *pp,
                               const acc_genpp_td_t *td,
                               const acc_genpp_state_t *st,
                               const poly_q_vec_d A,
                               uint32_t x,
                               acc_genpp_update_apply_fn apply_update,
                               void *apply_opaque) {
  poly_q_vec_d acc_i;
  acc_genpp_gmat_t B1y;
  poly_q_vec_d ginv_ua[ACC_GENPP_GADGET_DIM], delta;
  uint8_t ra[ACC_GENPP_SEED_BYTES];
  int pos = acc_genpp_dict_find(st, x);
  int ok = 1;
  uint32_t add_index;

  if (pos < 0 || st->dict[pos].op != ACC_GENPP_OP_ADD) {
    return 0;
  }
  add_index = st->dict[pos].index;

  poly_q_vec_d_init(acc_i);
  poly_q_vec_d_init(delta);
  acc_genpp_gmat_init(&B1y);
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_vec_d_init(ginv_ua[i]);
  }

  poly_q_vec_d_set(acc_i, st->acc0);
  acc_genpp_ginv_vec(ginv_ua, pp->ua);
  for (size_t i = 0; i < ACC_GENPP_MAX_DICT_ENTRIES; i++) {
    if (st->dict[i].used && st->dict[i].op == ACC_GENPP_OP_DEL &&
        st->dict[i].index > 0 && st->dict[i].index <= add_index) {
      acc_genpp_evalf(&B1y, &pp->Ba, st->dict[i].handle);
      acc_genpp_gmat_mul_gvec(delta, &B1y, ginv_ua);
      if (st->dict[i].update_sign >= 0) {
        poly_q_vec_d_add(acc_i, acc_i, delta);
      } else {
        poly_q_vec_d_sub(acc_i, acc_i, delta);
      }
    }
  }

  acc_genpp_Fk(ra, st->ka, x, acc_i);
  acc_genpp_sampleleft_randomized(sx, pp, td, acc_i, x, ra);
  for (size_t i = 0; i < ACC_GENPP_MAX_DICT_ENTRIES; i++) {
    if (st->dict[i].used && st->dict[i].op == ACC_GENPP_OP_DEL &&
        st->dict[i].index > add_index && st->dict[i].index <= st->ctr) {
      acc_genpp_update_t msg = {
          st->dict[i].index, ACC_GENPP_OP_DEL, st->dict[i].handle,
          st->dict[i].update_sign};

      if (apply_update != NULL) {
        ok = apply_update(apply_opaque, sx, pp, x, &msg);
      } else {
        ok = acc_genpp_memwitup_tails(sx, pp, x, &msg, NULL, NULL, 0,
                                      NULL);
      }
      if (!ok) {
        break;
      }
    }
  }

  (void)A;
  for (size_t i = 0; i < ACC_GENPP_GADGET_DIM; i++) {
    poly_q_vec_d_clear(ginv_ua[i]);
  }
  acc_genpp_gmat_clear(&B1y);
  poly_q_vec_d_clear(delta);
  poly_q_vec_d_clear(acc_i);
  return ok;
}

int acc_genpp_memwitsync(acc_genpp_wit_t *sx,
                         const acc_genpp_pp_t *pp,
                         const acc_genpp_td_t *td,
                         const acc_genpp_state_t *st,
                         const poly_q_vec_d A,
                         uint32_t x) {
  return acc_genpp_memwitsync_apply(sx, pp, td, st, A, x, NULL, NULL);
}
