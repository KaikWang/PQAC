#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "acc80_crt.h"

static uint64_t splitmix64(uint64_t *state) {
  uint64_t value = (*state += UINT64_C(0x9e3779b97f4a7c15));
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

static int check_sparse_product(const acc80_poly_t *lhs,
                                const acc80_poly_t *rhs,
                                const acc80_poly_t *product,
                                const acc80_ctx_t *ctx) {
  acc80_i128 expected[ACC80_N] = {0};
  const size_t positions[] = {
      0, 1, 2, 3, 4, 5, 6, 7,
      ACC80_N - 8, ACC80_N - 7, ACC80_N - 6, ACC80_N - 5,
      ACC80_N - 4, ACC80_N - 3, ACC80_N - 2, ACC80_N - 1};
  for (size_t i = 0; i < sizeof(positions) / sizeof(positions[0]); i++) {
    const size_t lhs_degree = positions[i];
    const acc80_i128 a =
        acc80_coeff_to_centered(&lhs->coeff[lhs_degree], ctx);
    for (size_t j = 0; j < sizeof(positions) / sizeof(positions[0]); j++) {
      const size_t rhs_degree = positions[j];
      const acc80_i128 b =
          acc80_coeff_to_centered(&rhs->coeff[rhs_degree], ctx);
      const size_t degree = lhs_degree + rhs_degree;
      if (degree < ACC80_N) expected[degree] += a * b;
      else expected[degree - ACC80_N] -= a * b;
    }
  }
  const acc80_i128 modulus = (acc80_i128)acc80_modulus();
  for (size_t i = 0; i < ACC80_N; i++) {
    acc80_i128 reduced = expected[i] % modulus;
    if (reduced < 0) reduced += modulus;
    if (acc80_coeff_to_u128(&product->coeff[i], ctx) !=
        (acc80_u128)reduced) return 0;
  }
  return 1;
}

int main(void) {
  acc80_ctx_t ctx;
  acc80_poly_t *lhs = calloc(1, sizeof(*lhs));
  acc80_poly_t *rhs = calloc(1, sizeof(*rhs));
  acc80_poly_t *product = calloc(1, sizeof(*product));
  acc80_poly_t *decoded = calloc(1, sizeof(*decoded));
  uint8_t *wire = malloc(ACC80_POLY_PACKED_BYTES);
  uint64_t random_state = UINT64_C(0x243f6a8885a308d3);
  int ok = lhs != NULL && rhs != NULL && product != NULL && decoded != NULL &&
           wire != NULL &&
           acc80_ctx_init(&ctx);

  const size_t positions[] = {
      0, 1, 2, 3, 4, 5, 6, 7,
      ACC80_N - 8, ACC80_N - 7, ACC80_N - 6, ACC80_N - 5,
      ACC80_N - 4, ACC80_N - 3, ACC80_N - 2, ACC80_N - 1};
  for (size_t i = 0;
       ok && i < sizeof(positions) / sizeof(positions[0]); i++) {
    const int64_t a = (int64_t)(splitmix64(&random_state) % 2001) - 1000;
    const int64_t b = (int64_t)(splitmix64(&random_state) % 2001) - 1000;
    const acc80_i128 modulus = (acc80_i128)acc80_modulus();
    acc80_coeff_from_u128(&lhs->coeff[positions[i]],
        (acc80_u128)(a < 0 ? modulus + a : a));
    acc80_coeff_from_u128(&rhs->coeff[positions[i]],
        (acc80_u128)(b < 0 ? modulus + b : b));
  }
  if (ok) {
    acc80_poly_mul(product, lhs, rhs, &ctx);
    ok = check_sparse_product(lhs, rhs, product, &ctx);
  }
  for (size_t trial = 0; ok && trial < 10000; trial++) {
    acc80_coeff_t value;
    int32_t digits[ACC80_GADGET_DIGITS];
    const acc80_u128 raw =
        (((acc80_u128)splitmix64(&random_state) << 64) |
         splitmix64(&random_state)) % acc80_modulus();
    acc80_coeff_from_u128(&value, raw);
    if (acc80_coeff_to_u128(&value, &ctx) != raw ||
        !acc80_balanced_decompose(digits, &value, &ctx)) ok = 0;
  }
  if (ok) {
    ok = acc80_poly_pack(wire, product) == ACC80_POLY_PACKED_BYTES &&
         acc80_poly_unpack(decoded, wire, ACC80_POLY_PACKED_BYTES) &&
         acc80_poly_equal(product, decoded);
  }
  if (ok) {
    for (size_t byte = 0; byte < ACC80_COEFF_PACKED_BYTES; byte++) {
      wire[byte] = UINT8_C(0xff);
    }
    ok = !acc80_poly_unpack(decoded, wire, ACC80_POLY_PACKED_BYTES);
  }

  free(wire);
  free(decoded);
  free(product);
  free(rhs);
  free(lhs);
  if (!ok) {
    fprintf(stderr, "q96 CRT/NTT differential test failed.\n");
    return 1;
  }
  printf("q96 CRT/NTT core passed: N=%u, three 32-bit limbs, base=%lu, "
         "digits=%u, canonical accumulator wire=%.2f KiB.\n", ACC80_N,
         (unsigned long)ACC80_GADGET_BASE, ACC80_GADGET_DIGITS,
         ACC80_POLY_PACKED_BYTES / 1024.0);
  return 0;
}
