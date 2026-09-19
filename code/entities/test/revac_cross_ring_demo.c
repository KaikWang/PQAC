#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "revac_cross_ring.h"

static int64_t centered_mod_q(int64_t value) {
  value %= PARAM_Q;
  if (value < 0) value += PARAM_Q;
  if (value > PARAM_Q / 2) value -= PARAM_Q;
  return value;
}

int main(void) {
  poly_q a;
  poly_q b;
  poly_q product;
  int64_t embedded_a[ACC96_N];
  int64_t embedded_b[ACC96_N];
  int64_t embedded_product[ACC96_N];
  int64_t integer_product[PARAM_N] = {0};
  int ok = 1;

  arith_setup();
  random_init();
  poly_q_init(a);
  poly_q_init(b);
  poly_q_init(product);
  /* Use full-width deterministic centered coefficients.  The scalar sampler in
   * poly_q_sampling.c is intentionally private to the matrix/vector APIs. */
  for (size_t coefficient = 0; coefficient < PARAM_N; coefficient++) {
    const int64_t av =
        (int64_t)((UINT64_C(1103515245) * (coefficient + 17U) + 12345U) %
                  PARAM_Q) -
        PARAM_Q / 2;
    const int64_t bv =
        (int64_t)((UINT64_C(2654435761) * (coefficient + 29U) + 67891U) %
                  PARAM_Q) -
        PARAM_Q / 2;
    poly_q_set_coeff(a, coefficient, av);
    poly_q_set_coeff(b, coefficient, bv);
  }
  poly_q_mul(product, a, b);
  revac_q19_poly_embed_q96_i64(embedded_a, a);
  revac_q19_poly_embed_q96_i64(embedded_b, b);
  revac_q19_poly_embed_q96_i64(embedded_product, product);

  for (size_t i = 0; i < PARAM_N; i++) {
    for (size_t j = 0; j < PARAM_N; j++) {
      const size_t degree = i + j;
      const int64_t term = embedded_a[16U * i] * embedded_b[16U * j];
      if (degree < PARAM_N)
        integer_product[degree] += term;
      else
        integer_product[degree - PARAM_N] -= term;
    }
  }
  for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
    if ((coefficient % 16U) != 0U &&
        (embedded_a[coefficient] != 0 || embedded_b[coefficient] != 0 ||
         embedded_product[coefficient] != 0)) {
      ok = 0;
      break;
    }
  }
  for (size_t coefficient = 0; coefficient < PARAM_N && ok; coefficient++) {
    if (centered_mod_q(integer_product[coefficient]) !=
        embedded_product[16U * coefficient])
      ok = 0;
  }

  poly_q_clear(product);
  poly_q_clear(b);
  poly_q_clear(a);
  arith_teardown();
  if (!ok) {
    fprintf(stderr, "q19 -> degree-4096 乘法嵌入检查失败\n");
    return 1;
  }
  printf("跨环嵌入通过：X^i -> X^(16i) 保持 q19 negacyclic 乘法，"
         "可让 q96 MemVer 与匿名签名承诺共用 32 个句柄比特。\n");
  return 0;
}
