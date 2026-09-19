#ifndef ACC96_TRAPDOOR_H
#define ACC96_TRAPDOOR_H

#include <stdint.h>

#include "acc96_relation.h"

typedef struct {
  acc96_poly_t r[2][ACC96_TRUNC_RIGHT];
  double squared_spectral_norm;
  uint32_t rejection_counter;
} acc96_trapdoor_t;

/* The production theorem bound remains ACC_TRUNC_R_MAX_...; the double FFT
 * pre-filter accepts only below 99% of it.  This large guard dominates normal
 * roundoff, while an Arb/ACB interval recomputation remains a release check. */
#define ACC96_TRAPDOOR_NUMERIC_SAFETY 0.99

/* Generate the q96 truncated-gadget trapdoor, reject until its canonical
 * embedding satisfies the generated spectral bound, and replace pp->trunc_b
 * by A*R_1+R_0.  pp must first be expanded by acc96_relation_pp_expand(). */
int acc96_trapdoor_generate(acc96_trapdoor_t *trapdoor,
                            acc96_relation_pp_t *pp,
                            const uint8_t seed[ACC96_SEED_BYTES],
                            const acc96_ctx_t *ctx);

double acc96_trapdoor_squared_spectral_norm(
    const acc96_trapdoor_t *trapdoor,
    const acc96_ctx_t *ctx);

int acc96_trapdoor_public_relation_holds(
    const acc96_trapdoor_t *trapdoor,
    const acc96_relation_pp_t *pp,
    const acc96_ctx_t *ctx);

#endif
