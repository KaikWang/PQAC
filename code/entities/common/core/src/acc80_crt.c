#include "acc80_crt.h"

#include <stdlib.h>
#include <string.h>

static const uint64_t ACC80_MODS[ACC80_LIMBS] = {
    ACC80_P0, ACC80_P1, ACC80_P2};
static const uint64_t ACC80_GARNER_INVERSES[ACC80_LIMBS] = {
    UINT64_C(0), ACC80_GARNER_INV1, ACC80_GARNER_INV2};
static const uint64_t ACC80_TWO64_MODS[ACC80_LIMBS] = {
    UINT64_C(196588), UINT64_C(15138488), UINT64_C(138064159)};
static const uint64_t ACC80_MUL_RECIPROCALS[ACC80_LIMBS] = {
    UINT64_C(17180655636), UINT64_C(34363408712), UINT64_C(34371276513)};

static uint64_t add_mod(uint64_t a, uint64_t b, uint64_t mod) {
  const uint64_t sum = a + b;
  return sum >= mod ? sum - mod : sum;
}

static uint64_t sub_mod(uint64_t a, uint64_t b, uint64_t mod) {
  return a >= b ? a - b : mod - (b - a);
}

static uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t mod) {
  /* Every CRT prime is below 2^31 and all lane values are canonical, hence
   * (mod - 1)^2 < 2^64.  A 128-bit dividend only routes this hot NTT path
   * through libgcc's software __umodti3 without adding any range. */
  return (a * b) % mod;
}

static inline uint64_t mul_mod_precomputed(uint64_t a, uint64_t b,
                                           uint64_t mod,
                                           uint64_t reciprocal) {
  const uint64_t product = a * b;
  const uint64_t quotient = (uint64_t)(
      ((acc80_u128)product * reciprocal) >> 64);
  uint64_t remainder = product - quotient * mod;
  /* reciprocal=floor(2^64/mod) and product<mod^2<2^64 imply that the
   * approximate quotient is at most one below floor(product/mod). */
  if (remainder >= mod) remainder -= mod;
  return remainder;
}

static uint64_t u128_mod_lane(acc80_u128 value, size_t lane) {
  const uint64_t mod = ACC80_MODS[lane];
  const uint64_t low = (uint64_t)value;
  const uint64_t high = (uint64_t)(value >> 64);
  return add_mod(mul_mod(high % mod, ACC80_TWO64_MODS[lane], mod),
                 low % mod, mod);
}

static uint64_t pow_mod(uint64_t base, uint64_t exponent, uint64_t mod) {
  uint64_t result = 1;
  while (exponent != 0) {
    if ((exponent & 1U) != 0) result = mul_mod(result, base, mod);
    base = mul_mod(base, base, mod);
    exponent >>= 1;
  }
  return result;
}

static uint64_t inv_mod(uint64_t value, uint64_t mod) {
  /* Both CRT factors are prime. */
  return pow_mod(value, mod - 2, mod);
}

static uint64_t root_2n(uint64_t mod) {
  const uint64_t exponent = (mod - 1) / (2U * ACC80_N);
  for (uint64_t candidate = 2; candidate < 100000; candidate++) {
    const uint64_t root = pow_mod(candidate, exponent, mod);
    if (pow_mod(root, ACC80_N, mod) == mod - 1 &&
        pow_mod(root, 2U * ACC80_N, mod) == 1) {
      return root;
    }
  }
  return 0;
}

acc80_u128 acc80_modulus(void) {
  return (acc80_u128)ACC80_P0 * ACC80_P1 * ACC80_P2;
}

int acc80_ctx_init(acc80_ctx_t *ctx) {
  if (ctx == NULL) return 0;
  memset(ctx, 0, sizeof(*ctx));
  for (size_t lane = 0; lane < ACC80_LIMBS; lane++) {
    const uint64_t mod = ACC80_MODS[lane];
    ctx->psi[lane] = root_2n(mod);
    if (ctx->psi[lane] == 0) return 0;
    ctx->psi_inv[lane] = inv_mod(ctx->psi[lane], mod);
    ctx->omega[lane] = mul_mod(ctx->psi[lane], ctx->psi[lane], mod);
    ctx->omega_inv[lane] = inv_mod(ctx->omega[lane], mod);
    ctx->n_inv[lane] = inv_mod(ACC80_N, mod);
    ctx->mul_reciprocal[lane] = ACC80_MUL_RECIPROCALS[lane];
    for (size_t stage = 1; stage <= ACC80_LOG_N; stage++) {
      const size_t len = (size_t)1U << stage;
      ctx->stage_root[lane][stage] =
          pow_mod(ctx->omega[lane], ACC80_N / len, mod);
      ctx->stage_root_inv[lane][stage] =
          pow_mod(ctx->omega_inv[lane], ACC80_N / len, mod);
    }
  }
  acc80_u128 prefix = ACC80_P0;
  ctx->garner_inv[0] = 0;
  for (size_t lane = 1; lane < ACC80_LIMBS; lane++) {
    ctx->garner_inv[lane] = inv_mod((uint64_t)(prefix % ACC80_MODS[lane]),
                                    ACC80_MODS[lane]);
    prefix *= ACC80_MODS[lane];
  }
  ctx->initialized = 1;
  return 1;
}

static void ntt(uint64_t values[ACC80_N], const uint64_t stage_root[],
                uint64_t n_inv, uint64_t mod, uint64_t reciprocal,
                int inverse) {
  for (size_t i = 1, j = 0; i < ACC80_N; i++) {
    size_t bit = ACC80_N >> 1;
    for (; (j & bit) != 0; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      const uint64_t swap = values[i];
      values[i] = values[j];
      values[j] = swap;
    }
  }
  for (size_t len = 2, stage = 1; len <= ACC80_N; len <<= 1, stage++) {
    const uint64_t step = stage_root[stage];
    for (size_t start = 0; start < ACC80_N; start += len) {
      uint64_t twiddle = 1;
      for (size_t offset = 0; offset < len / 2; offset++) {
        const uint64_t even = values[start + offset];
        const uint64_t odd = mul_mod_precomputed(
            values[start + offset + len / 2], twiddle, mod, reciprocal);
        values[start + offset] = add_mod(even, odd, mod);
        values[start + offset + len / 2] = sub_mod(even, odd, mod);
        twiddle = mul_mod_precomputed(twiddle, step, mod, reciprocal);
      }
    }
  }
  if (inverse) {
    for (size_t i = 0; i < ACC80_N; i++) {
      values[i] = mul_mod_precomputed(values[i], n_inv, mod, reciprocal);
    }
  }
}

void acc80_poly_zero(acc80_poly_t *out) {
  memset(out, 0, sizeof(*out));
}

void acc80_poly_set(acc80_poly_t *out, const acc80_poly_t *in) {
  memcpy(out, in, sizeof(*out));
}

void acc80_poly_add(acc80_poly_t *out, const acc80_poly_t *lhs,
                    const acc80_poly_t *rhs) {
  for (size_t i = 0; i < ACC80_N; i++) {
    for (size_t lane = 0; lane < ACC80_LIMBS; lane++) {
      out->coeff[i].limb[lane] = add_mod(lhs->coeff[i].limb[lane],
                                         rhs->coeff[i].limb[lane],
                                         ACC80_MODS[lane]);
    }
  }
}

void acc80_poly_sub(acc80_poly_t *out, const acc80_poly_t *lhs,
                    const acc80_poly_t *rhs) {
  for (size_t i = 0; i < ACC80_N; i++) {
    for (size_t lane = 0; lane < ACC80_LIMBS; lane++) {
      out->coeff[i].limb[lane] = sub_mod(lhs->coeff[i].limb[lane],
                                         rhs->coeff[i].limb[lane],
                                         ACC80_MODS[lane]);
    }
  }
}

void acc80_poly_neg(acc80_poly_t *out, const acc80_poly_t *in) {
  for (size_t i = 0; i < ACC80_N; i++) {
    for (size_t lane = 0; lane < ACC80_LIMBS; lane++) {
      out->coeff[i].limb[lane] = in->coeff[i].limb[lane] == 0
          ? 0 : ACC80_MODS[lane] - in->coeff[i].limb[lane];
    }
  }
}

int acc80_poly_equal(const acc80_poly_t *lhs, const acc80_poly_t *rhs) {
  return memcmp(lhs, rhs, sizeof(*lhs)) == 0;
}

void acc80_poly_mul_workspace(acc80_poly_t *out, const acc80_poly_t *lhs,
                              const acc80_poly_t *rhs,
                              const acc80_ctx_t *ctx, uint64_t *left,
                              uint64_t *right) {
  if (out == NULL || lhs == NULL || rhs == NULL || left == NULL ||
      right == NULL || ctx == NULL || !ctx->initialized)
    abort();
  for (size_t lane = 0; lane < ACC80_LIMBS; lane++) {
    const uint64_t mod = ACC80_MODS[lane];
    const uint64_t reciprocal = ctx->mul_reciprocal[lane];
    uint64_t twist = 1;
    for (size_t i = 0; i < ACC80_N; i++) {
      left[i] = mul_mod_precomputed(lhs->coeff[i].limb[lane], twist,
                                    mod, reciprocal);
      right[i] = mul_mod_precomputed(rhs->coeff[i].limb[lane], twist,
                                     mod, reciprocal);
      twist = mul_mod_precomputed(twist, ctx->psi[lane], mod, reciprocal);
    }
    ntt(left, ctx->stage_root[lane], ctx->n_inv[lane], mod, reciprocal, 0);
    ntt(right, ctx->stage_root[lane], ctx->n_inv[lane], mod, reciprocal, 0);
    for (size_t i = 0; i < ACC80_N; i++) {
      left[i] = mul_mod_precomputed(left[i], right[i], mod, reciprocal);
    }
    ntt(left, ctx->stage_root_inv[lane], ctx->n_inv[lane], mod,
        reciprocal, 1);
    twist = 1;
    for (size_t i = 0; i < ACC80_N; i++) {
      out->coeff[i].limb[lane] =
          mul_mod_precomputed(left[i], twist, mod, reciprocal);
      twist = mul_mod_precomputed(twist, ctx->psi_inv[lane], mod,
                                  reciprocal);
    }
  }
}

void acc80_poly_mul(acc80_poly_t *out, const acc80_poly_t *lhs,
                    const acc80_poly_t *rhs, const acc80_ctx_t *ctx) {
  uint64_t *left = malloc(ACC80_N * sizeof(*left));
  uint64_t *right = malloc(ACC80_N * sizeof(*right));
  if (left == NULL || right == NULL || ctx == NULL || !ctx->initialized) {
    free(right);
    free(left);
    abort();
  }
  acc80_poly_mul_workspace(out, lhs, rhs, ctx, left, right);
  free(right);
  free(left);
}

void acc80_coeff_from_u128(acc80_coeff_t *out, acc80_u128 value) {
  for (size_t lane = 0; lane < ACC80_LIMBS; lane++) {
    out->limb[lane] = u128_mod_lane(value, lane);
  }
}

acc80_u128 acc80_coeff_to_u128(const acc80_coeff_t *value,
                               const acc80_ctx_t *ctx) {
  acc80_u128 reconstructed = value->limb[0];
  acc80_u128 prefix = ACC80_P0;
  for (size_t lane = 1; lane < ACC80_LIMBS; lane++) {
    const uint64_t mod = ACC80_MODS[lane];
    /* Before lane 1 the value is below 2^32; before lane 2 it is below
     * P0*P1 < 2^64.  The cast is therefore exact and avoids a software
     * 128-bit remainder in every CRT-to-integer reconstruction. */
    const uint64_t reduced = (uint64_t)reconstructed % mod;
    const uint64_t difference = sub_mod(
        value->limb[lane], reduced, mod);
    const uint64_t inverse = ctx != NULL ? ctx->garner_inv[lane]
                                         : ACC80_GARNER_INVERSES[lane];
    const uint64_t multiplier = mul_mod(difference, inverse, mod);
    reconstructed += prefix * multiplier;
    prefix *= mod;
  }
  return reconstructed;
}

acc80_i128 acc80_coeff_to_centered(const acc80_coeff_t *value,
                                   const acc80_ctx_t *ctx) {
  const acc80_u128 canonical = acc80_coeff_to_u128(value, ctx);
  const acc80_u128 modulus = acc80_modulus();
  return canonical > modulus / 2 ? (acc80_i128)(canonical - modulus)
                                 : (acc80_i128)canonical;
}

static int balanced_decompose_integer(
    int32_t out[ACC80_GADGET_DIGITS], acc80_i128 integer) {
  acc80_i128 remaining = integer;
  for (size_t i = 0; i < ACC80_GADGET_DIGITS; i++) {
    acc80_i128 digit = remaining % (acc80_i128)ACC80_GADGET_BASE;
    if (digit < 0) digit += ACC80_GADGET_BASE;
    if (digit > ACC80_GADGET_BOUND) digit -= ACC80_GADGET_BASE;
    out[i] = (int32_t)digit;
    remaining = (remaining - digit) / (acc80_i128)ACC80_GADGET_BASE;
  }
  return remaining == 0;
}

int acc80_balanced_decompose_canonical(
    int32_t out[ACC80_GADGET_DIGITS], const acc80_coeff_t *value,
    const acc80_ctx_t *ctx) {
  return balanced_decompose_integer(out,
                                    acc80_coeff_to_centered(value, ctx));
}

static acc80_i128 floor_div(acc80_i128 numerator, acc80_i128 denominator) {
  acc80_i128 quotient = numerator / denominator;
  const acc80_i128 remainder = numerator % denominator;
  if (remainder < 0) quotient--;
  return quotient;
}

static acc80_i128 ceil_div(acc80_i128 numerator, acc80_i128 denominator) {
  return -floor_div(-numerator, denominator);
}

static acc80_u128 abs_i128(acc80_i128 value) {
  return value < 0 ? (acc80_u128)(-value) : (acc80_u128)value;
}

static int digits_lex_less(const int32_t lhs[ACC80_GADGET_DIGITS],
                           const int32_t rhs[ACC80_GADGET_DIGITS]) {
  for (size_t offset = 0; offset < ACC80_GADGET_DIGITS; offset++) {
    const size_t i = ACC80_GADGET_DIGITS - 1 - offset;
    if (lhs[i] != rhs[i]) return lhs[i] < rhs[i];
  }
  return 0;
}

int acc80_balanced_decompose_min_norm(
    int32_t out[ACC80_GADGET_DIGITS], const acc80_coeff_t *value,
    const acc80_ctx_t *ctx) {
  const acc80_i128 centered = acc80_coeff_to_centered(value, ctx);
  const acc80_i128 modulus = (acc80_i128)acc80_modulus();
  acc80_i128 radius = 0;
  acc80_i128 power = 1;
  for (size_t i = 0; i < ACC80_GADGET_DIGITS; i++) {
    radius += (acc80_i128)ACC80_GADGET_BOUND * power;
    power *= (acc80_i128)ACC80_GADGET_BASE;
  }

  const acc80_i128 first = ceil_div(-radius - centered, modulus);
  const acc80_i128 last = floor_div(radius - centered, modulus);
  if (first == last) {
    return balanced_decompose_integer(out, centered + first * modulus);
  }
  acc80_u128 best_norm2 = 0;
  acc80_u128 best_l1 = 0;
  acc80_u128 best_representative = 0;
  unsigned best_weight = 0;
  int found = 0;

  for (acc80_i128 multiple = first; multiple <= last; multiple++) {
    int32_t candidate[ACC80_GADGET_DIGITS];
    const acc80_i128 representative = centered + multiple * modulus;
    acc80_u128 norm2 = 0;
    acc80_u128 l1 = 0;
    unsigned weight = 0;
    if (!balanced_decompose_integer(candidate, representative)) continue;
    for (size_t i = 0; i < ACC80_GADGET_DIGITS; i++) {
      const acc80_i128 digit = candidate[i];
      const acc80_u128 magnitude = abs_i128(digit);
      norm2 += magnitude * magnitude;
      l1 += magnitude;
      weight += digit != 0;
    }
    const acc80_u128 representative_magnitude = abs_i128(representative);
    const int better =
        !found || norm2 < best_norm2 ||
        (norm2 == best_norm2 && l1 < best_l1) ||
        (norm2 == best_norm2 && l1 == best_l1 && weight < best_weight) ||
        (norm2 == best_norm2 && l1 == best_l1 && weight == best_weight &&
         representative_magnitude < best_representative) ||
        (norm2 == best_norm2 && l1 == best_l1 && weight == best_weight &&
         representative_magnitude == best_representative &&
         digits_lex_less(candidate, out));
    if (better) {
      memcpy(out, candidate, sizeof(candidate));
      best_norm2 = norm2;
      best_l1 = l1;
      best_weight = weight;
      best_representative = representative_magnitude;
      found = 1;
    }
  }
  return found;
}

int acc80_balanced_decompose(int32_t out[ACC80_GADGET_DIGITS],
                             const acc80_coeff_t *value,
                             const acc80_ctx_t *ctx) {
  return acc80_balanced_decompose_min_norm(out, value, ctx);
}

size_t acc80_poly_pack(uint8_t out[ACC80_POLY_PACKED_BYTES],
                       const acc80_poly_t *poly) {
  size_t offset = 0;
  for (size_t coeff = 0; coeff < ACC80_N; coeff++) {
    acc80_u128 value = acc80_coeff_to_u128(&poly->coeff[coeff], NULL);
    for (size_t byte = 0; byte < ACC80_COEFF_PACKED_BYTES; byte++) {
      out[offset++] = (uint8_t)(value >> (8 * byte));
    }
  }
  return offset;
}

int acc80_poly_unpack(acc80_poly_t *poly, const uint8_t *in, size_t in_len) {
  size_t offset = 0;
  if (poly == NULL || in == NULL || in_len != ACC80_POLY_PACKED_BYTES) return 0;
  for (size_t coeff = 0; coeff < ACC80_N; coeff++) {
    acc80_u128 value = 0;
    for (size_t byte = 0; byte < ACC80_COEFF_PACKED_BYTES; byte++) {
      value |= (acc80_u128)in[offset++] << (8 * byte);
    }
    if (value >= acc80_modulus()) {
      acc80_poly_zero(poly);
      return 0;
    }
    acc80_coeff_from_u128(&poly->coeff[coeff], value);
  }
  return 1;
}
