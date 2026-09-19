#include "revac_cross_ring.h"

#include <string.h>

void revac_q19_coeffs_embed_q96_i64(int64_t out[ACC96_N],
                                    const int64_t in[PARAM_N]) {
  memset(out, 0, ACC96_N * sizeof(*out));
  for (size_t coefficient = 0; coefficient < PARAM_N; coefficient++)
    out[16U * coefficient] = in[coefficient];
}

void revac_q19_poly_embed_q96_i64(int64_t out[ACC96_N],
                                  const poly_q in) {
  memset(out, 0, ACC96_N * sizeof(*out));
  for (size_t coefficient = 0; coefficient < PARAM_N; coefficient++)
    out[16U * coefficient] =
        poly_q_get_coeff_centered(in, coefficient);
}
