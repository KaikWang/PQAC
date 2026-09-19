#ifndef REVAC_CROSS_RING_H
#define REVAC_CROSS_RING_H

#include <stdint.h>

#include "acc96_crt.h"
#include "arith.h"

#if ACC96_N != (16U * PARAM_N)
#error "q19-to-q96 embedding requires ACC96_N = 16 * PARAM_N"
#endif

/* Inject R_q19 = Z_q19[X]/(X^256+1) into the degree-4096 proof ring by
 * X^i -> X^(16i).  This preserves addition and negacyclic multiplication;
 * q19 modular wrap is represented separately by the existing exact lift. */
void revac_q19_poly_embed_q96_i64(int64_t out[ACC96_N],
                                  const poly_q in);

void revac_q19_coeffs_embed_q96_i64(int64_t out[ACC96_N],
                                    const int64_t in[PARAM_N]);

#endif
