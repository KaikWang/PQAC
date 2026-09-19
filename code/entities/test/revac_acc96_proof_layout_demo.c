#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "acc96_proof_layout.h"

static uint64_t splitmix64(uint64_t *state) {
  uint64_t value = (*state += UINT64_C(0x9e3779b97f4a7c15));
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

static int check_value(acc96_u128 canonical, const acc96_ctx_t *ctx) {
  acc96_coeff_t coefficient;
  int32_t limbs[ACC96_PROOF_LIMBS];
  acc96_coeff_from_u128(&coefficient, canonical);
  if (!acc96_proof_limb_decompose(limbs, &coefficient, ctx)) return 0;
  for (size_t i = 0; i < ACC96_PROOF_LIMBS; i++)
    if (limbs[i] < -ACC96_PROOF_LIMB_BOUND ||
        limbs[i] > ACC96_PROOF_LIMB_BOUND)
      return 0;
  return acc96_proof_limb_reconstruct(limbs) ==
         acc96_coeff_to_centered(&coefficient, ctx);
}

static int check_window_value(int64_t signed_value, const acc96_ctx_t *ctx) {
  acc96_coeff_t coefficient;
  int32_t limbs[ACC96_WINDOW_PROOF_LIMBS];
  acc96_i128 canonical = signed_value;
  if (canonical < 0) canonical += (acc96_i128)acc96_modulus();
  acc96_coeff_from_u128(&coefficient, (acc96_u128)canonical);
  return acc96_window_proof_limb_decompose(limbs, &coefficient, ctx) &&
         acc96_window_proof_limb_reconstruct(limbs) == signed_value;
}

static int check_polyphase_product(const acc96_ctx_t *ctx) {
  static const uint64_t moduli[ACC96_LIMBS] = {
      ACC96_P0, ACC96_P1, ACC96_P2};
  acc96_poly_t *left = calloc(1, sizeof(*left));
  acc96_poly_t *right = calloc(1, sizeof(*right));
  acc96_poly_t *product = calloc(1, sizeof(*product));
  int64_t (*left_phase)[ACC96_LABRADOR_BASE_DEGREE] =
      calloc(ACC96_LABRADOR_PHASES, sizeof(*left_phase));
  int64_t (*right_phase)[ACC96_LABRADOR_BASE_DEGREE] =
      calloc(ACC96_LABRADOR_PHASES, sizeof(*right_phase));
  int64_t (*expected)[ACC96_LABRADOR_BASE_DEGREE] =
      calloc(ACC96_LABRADOR_PHASES, sizeof(*expected));
  int ok = left != NULL && right != NULL && product != NULL &&
           left_phase != NULL && right_phase != NULL && expected != NULL;
  if (!ok) goto cleanup;

  for (size_t phase = 0; phase < ACC96_LABRADOR_PHASES; phase++) {
    for (size_t coefficient = 0;
         coefficient < ACC96_LABRADOR_BASE_DEGREE; coefficient++) {
      const size_t source =
          acc96_proof_polyphase_source_index(phase, coefficient);
      const int64_t left_value = (int64_t)((source * 5U + 1U) % 5U) - 2;
      const int64_t right_value = (int64_t)((source * 7U + 3U) % 5U) - 2;
      acc96_i128 left_canonical = left_value;
      acc96_i128 right_canonical = right_value;
      if (left_canonical < 0) left_canonical += (acc96_i128)acc96_modulus();
      if (right_canonical < 0) right_canonical += (acc96_i128)acc96_modulus();
      acc96_coeff_from_u128(&left->coeff[source],
                            (acc96_u128)left_canonical);
      acc96_coeff_from_u128(&right->coeff[source],
                            (acc96_u128)right_canonical);
      left_phase[phase][coefficient] = left_value;
      right_phase[phase][coefficient] = right_value;
    }
  }

  /* In R[Y]/(Y^256+1), X^(r+s) contributes an additional Y whenever
   * r+s >= 16.  The following is the exact 16-phase product formula. */
  for (size_t left_r = 0; left_r < ACC96_LABRADOR_PHASES; left_r++) {
    for (size_t right_r = 0; right_r < ACC96_LABRADOR_PHASES; right_r++) {
      const size_t output_phase =
          (left_r + right_r) % ACC96_LABRADOR_PHASES;
      const size_t phase_carry =
          (left_r + right_r) / ACC96_LABRADOR_PHASES;
      for (size_t left_i = 0; left_i < ACC96_LABRADOR_BASE_DEGREE;
           left_i++) {
        for (size_t right_i = 0; right_i < ACC96_LABRADOR_BASE_DEGREE;
             right_i++) {
          size_t exponent = left_i + right_i + phase_carry;
          int64_t term = left_phase[left_r][left_i] *
                         right_phase[right_r][right_i];
          if (exponent >= ACC96_LABRADOR_BASE_DEGREE) {
            exponent -= ACC96_LABRADOR_BASE_DEGREE;
            term = -term;
          }
          expected[output_phase][exponent] += term;
        }
      }
    }
  }

  acc96_poly_mul(product, left, right, ctx);
  for (size_t phase = 0; phase < ACC96_LABRADOR_PHASES && ok; phase++) {
    for (size_t coefficient = 0;
         coefficient < ACC96_LABRADOR_BASE_DEGREE && ok; coefficient++) {
      const size_t source =
          acc96_proof_polyphase_source_index(phase, coefficient);
      for (size_t lane = 0; lane < ACC96_LIMBS; lane++) {
        int64_t reduced = expected[phase][coefficient] %
                          (int64_t)moduli[lane];
        if (reduced < 0) reduced += (int64_t)moduli[lane];
        if (product->coeff[source].limb[lane] != (uint64_t)reduced) {
          ok = 0;
          break;
        }
      }
    }
  }

cleanup:
  free(expected);
  free(right_phase);
  free(left_phase);
  free(product);
  free(right);
  free(left);
  return ok;
}

int main(void) {
  acc96_ctx_t ctx;
  acc96_relation_witness_t *witness = calloc(1, sizeof(*witness));
  acc96_poly_t *polynomial = calloc(1, sizeof(*polynomial));
  int32_t (*poly_limbs)[ACC96_N] =
      calloc(ACC96_PROOF_LIMBS, sizeof(*poly_limbs));
  int32_t (*window_limbs)[ACC96_N] =
      calloc(ACC96_WINDOW_PROOF_LIMBS, sizeof(*window_limbs));
  uint64_t state = UINT64_C(0x243f6a8885a308d3);
  if (witness == NULL || polynomial == NULL || poly_limbs == NULL ||
      window_limbs == NULL ||
      !acc96_ctx_init(&ctx))
    return 1;
  if (ACC96_LABRADOR_PHASES != 16U || !check_polyphase_product(&ctx))
    return 1;
  if (ACC96_PROOF_WITNESS_POLYS != 169U ||
      ACC96_PROOF_NON_TAG_POLYS != 168U ||
      ACC96_PROOF_SLOT_TAG != 24U)
    return 1;
  for (size_t slot = 0; slot < ACC96_PROOF_WITNESS_POLYS; slot++) {
    acc96_poly_t *poly = acc96_proof_witness_poly_at_mut(witness, slot);
    if (poly == NULL ||
        ((slot == ACC96_PROOF_SLOT_TAG) != acc96_proof_slot_is_tag(slot)))
      return 1;
    poly->coeff[0].limb[0] = (uint32_t)(slot + 1U);
  }
  for (size_t slot = 0; slot < ACC96_PROOF_WITNESS_POLYS; slot++) {
    const acc96_poly_t *poly = acc96_proof_witness_poly_at(witness, slot);
    if (poly == NULL || poly->coeff[0].limb[0] != slot + 1U) return 1;
  }
  if (acc96_proof_right_slot(23, 5) != 168U ||
      acc96_proof_right_slot(24, 0) != ACC96_PROOF_WITNESS_POLYS)
    return 1;
  {
    int32_t lift_limbs[ACC96_LIFT_PROOF_LIMBS];
    const acc96_i128 lift_edges[] = {
        0, 1, -1, (acc96_i128)ACC96_LIFT_COEFF_BOUND,
        -(acc96_i128)ACC96_LIFT_COEFF_BOUND};
    for (size_t i = 0; i < sizeof(lift_edges) / sizeof(lift_edges[0]); i++)
      if (!acc96_lift_proof_limb_decompose(lift_limbs, lift_edges[i]) ||
          acc96_lift_proof_limb_reconstruct(lift_limbs) != lift_edges[i])
        return 1;
  }
  const acc96_u128 q = acc96_modulus();
  const acc96_u128 edges[] = {0, 1, q - 1, q / 2, q / 2 + 1};
  for (size_t i = 0; i < sizeof(edges) / sizeof(edges[0]); i++)
    if (!check_value(edges[i], &ctx)) return 1;
  for (size_t i = 0; i < 10000; i++) {
    const acc96_u128 raw =
        (((acc96_u128)splitmix64(&state) << 64) | splitmix64(&state)) % q;
    if (!check_value(raw, &ctx)) return 1;
    if (i < ACC96_N)
      acc96_coeff_from_u128(&polynomial->coeff[i], raw);
  }
  if (!acc96_proof_poly_decompose(poly_limbs, polynomial, &ctx) ||
      acc96_proof_poly_max_chunk_norm2(
          (const int32_t (*)[ACC96_N])poly_limbs) >
          ACC96_PROOF_LIMB_CHUNK_NORM2_BOUND)
    return 1;
  if (ACC96_PROOF_LIMB_CHUNK_NORM2_BOUND > UINT64_MAX) return 1;
  {
    const int64_t window_edges[] = {
        0, 1, -1, ACC96_WINDOW_COEFF_BOUND,
        -ACC96_WINDOW_COEFF_BOUND, ACC96_WINDOW_UPDATE_COEFF_BOUND,
        -ACC96_WINDOW_UPDATE_COEFF_BOUND};
    for (size_t i = 0; i < sizeof(window_edges) / sizeof(window_edges[0]); i++)
      if (!check_window_value(window_edges[i], &ctx)) return 1;
    for (size_t i = 0; i < ACC96_N; i++) {
      const uint64_t span =
          UINT64_C(2) * (uint64_t)ACC96_WINDOW_COEFF_BOUND + 1U;
      int64_t value = (int64_t)(splitmix64(&state) % span) -
                      ACC96_WINDOW_COEFF_BOUND;
      acc96_i128 canonical = value;
      if (canonical < 0) canonical += (acc96_i128)q;
      acc96_coeff_from_u128(&polynomial->coeff[i], (acc96_u128)canonical);
    }
    if (!acc96_window_proof_poly_decompose(window_limbs, polynomial, &ctx) ||
        acc96_window_proof_poly_max_chunk_norm2(
            (const int32_t (*)[ACC96_N])window_limbs) >
            ACC96_WINDOW_PROOF_CHUNK_NORM2_BOUND)
      return 1;
    {
      acc96_coeff_t rejected;
      acc96_coeff_from_u128(&rejected,
                            (acc96_u128)ACC96_WINDOW_COEFF_BOUND + 1);
      if (acc96_window_proof_limb_decompose(
              window_limbs[0], &rejected, &ctx))
        return 1;
    }
  }
  printf("q88 proof layout passed: generic 4 x 28-bit and windowed "
         "5 x radix-8191 layouts; per-update bound=%lld; "
         "10000-update gate=%lld; degree-256 compact chunk norm^2 <= 2^33; "
         "16-phase negacyclic product preserved.\n",
         (long long)ACC96_WINDOW_UPDATE_COEFF_BOUND,
         (long long)ACC96_WINDOW_COEFF_BOUND);
  free(window_limbs);
  free(poly_limbs);
  free(polynomial);
  free(witness);
  return 0;
}
