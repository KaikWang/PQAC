#include "acc96_trapdoor.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include "acc96_truncated_params_generated.h"
#include "fips202.h"

_Static_assert(ACC_TRUNC_PARAM_N == ACC96_N,
               "q96 trapdoor degree mismatch");
_Static_assert(ACC_TRUNC_PARAM_K - ACC_TRUNC_PARAM_L ==
                   ACC96_TRUNC_RIGHT,
               "q96 trapdoor dimension mismatch");

static void absorb_u32_le(keccak_state *state, uint32_t value) {
  const uint8_t encoded[4] = {
      (uint8_t)value, (uint8_t)(value >> 8),
      (uint8_t)(value >> 16), (uint8_t)(value >> 24)};
  shake256_absorb(state, encoded, sizeof(encoded));
}

/* Centered binomial B_1: bit_0-bit_1, hence probabilities 1/4,1/2,1/4.
 * This is the trapdoor distribution used by the reference implementation;
 * uniform ternary sampling has larger variance and invalidates its spectral
 * rejection calibration. */
static void sample_binomial_b1(acc96_poly_t *out,
                               const uint8_t seed[ACC96_SEED_BYTES],
                               uint32_t counter, uint32_t row,
                               uint32_t column) {
  keccak_state state;
  uint8_t byte;
  acc96_poly_zero(out);
  shake256_init(&state);
  shake256_absorb(&state, seed, ACC96_SEED_BYTES);
  absorb_u32_le(&state, UINT32_C(0x52393654));
  absorb_u32_le(&state, counter);
  absorb_u32_le(&state, row);
  absorb_u32_le(&state, column);
  shake256_finalize(&state);
  for (size_t i = 0; i < ACC96_N; i++) {
    shake256_squeeze(&byte, 1, &state);
    const int value = (int)(byte & 1U) - (int)((byte >> 1) & 1U);
    if (value < 0)
      acc96_coeff_from_u128(&out->coeff[i], acc96_modulus() - 1);
    else
      acc96_coeff_from_u128(&out->coeff[i], (acc96_u128)value);
  }
}

static void fft(double complex values[ACC96_N]) {
  for (size_t i = 1, j = 0; i < ACC96_N; i++) {
    size_t bit = ACC96_N >> 1;
    for (; (j & bit) != 0; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      const double complex swap = values[i];
      values[i] = values[j];
      values[j] = swap;
    }
  }
  for (size_t length = 2; length <= ACC96_N; length <<= 1) {
    const double angle = -2.0 * M_PI / (double)length;
    const double complex step = cos(angle) + I * sin(angle);
    for (size_t start = 0; start < ACC96_N; start += length) {
      double complex root = 1.0;
      for (size_t offset = 0; offset < length / 2; offset++) {
        const double complex even = values[start + offset];
        const double complex odd = values[start + offset + length / 2] * root;
        values[start + offset] = even + odd;
        values[start + offset + length / 2] = even - odd;
        root *= step;
      }
    }
  }
}

double acc96_trapdoor_squared_spectral_norm(
    const acc96_trapdoor_t *trapdoor, const acc96_ctx_t *ctx) {
  double *energy = calloc(ACC96_N, sizeof(*energy));
  double complex *values = malloc(ACC96_N * sizeof(*values));
  double maximum = 0;
  if (energy == NULL || values == NULL) abort();
  for (size_t row = 0; row < 2; row++) {
    for (size_t i = 0; i < ACC96_N; i++) energy[i] = 0;
    for (size_t column = 0; column < ACC96_TRUNC_RIGHT; column++) {
      for (size_t i = 0; i < ACC96_N; i++) {
        const double coefficient = (double)acc96_coeff_to_centered(
            &trapdoor->r[row][column].coeff[i], ctx);
        const double angle = M_PI * (double)i / (double)ACC96_N;
        values[i] = coefficient * (cos(angle) + I * sin(angle));
      }
      fft(values);
      for (size_t i = 0; i < ACC96_N; i++) {
        const double real = creal(values[i]);
        const double imag = cimag(values[i]);
        energy[i] += real * real + imag * imag;
      }
    }
    for (size_t i = 0; i < ACC96_N; i++)
      if (energy[i] > maximum) maximum = energy[i];
  }
  free(values);
  free(energy);
  return maximum;
}

int acc96_trapdoor_public_relation_holds(
    const acc96_trapdoor_t *trapdoor, const acc96_relation_pp_t *pp,
    const acc96_ctx_t *ctx) {
  acc96_poly_t *product = malloc(sizeof(*product));
  acc96_poly_t *expected = malloc(sizeof(*expected));
  int ok = product != NULL && expected != NULL;
  if (!ok) {
    free(expected);
    free(product);
    return 0;
  }
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT && ok; i++) {
    acc96_poly_mul(product, &pp->atrap, &trapdoor->r[1][i], ctx);
    acc96_poly_add(expected, product, &trapdoor->r[0][i]);
    ok = acc96_poly_equal(expected, &pp->trunc_b[i]);
  }
  free(expected);
  free(product);
  return ok;
}

int acc96_trapdoor_generate(acc96_trapdoor_t *trapdoor,
                            acc96_relation_pp_t *pp,
                            const uint8_t seed[ACC96_SEED_BYTES],
                            const acc96_ctx_t *ctx) {
  acc96_poly_t *product;
  if (trapdoor == NULL || pp == NULL || seed == NULL || ctx == NULL ||
      !ctx->initialized)
    return 0;
  for (uint32_t counter = 0; counter < UINT32_MAX; counter++) {
    for (size_t row = 0; row < 2; row++)
      for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++)
        sample_binomial_b1(&trapdoor->r[row][i], seed, counter,
                           (uint32_t)row, (uint32_t)i);
    trapdoor->squared_spectral_norm =
        acc96_trapdoor_squared_spectral_norm(trapdoor, ctx);
    if (isfinite(trapdoor->squared_spectral_norm) &&
        trapdoor->squared_spectral_norm <=
            ACC96_TRAPDOOR_NUMERIC_SAFETY *
                ACC_TRUNC_R_MAX_SQ_SPECTRAL_NORM) {
      trapdoor->rejection_counter = counter;
      break;
    }
    if (counter == UINT32_MAX - 1) return 0;
  }
  product = malloc(sizeof(*product));
  if (product == NULL) return 0;
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++) {
    acc96_poly_mul(product, &pp->atrap, &trapdoor->r[1][i], ctx);
    acc96_poly_add(&pp->trunc_b[i], product, &trapdoor->r[0][i]);
  }
  free(product);
  return acc96_trapdoor_public_relation_holds(trapdoor, pp, ctx);
}
