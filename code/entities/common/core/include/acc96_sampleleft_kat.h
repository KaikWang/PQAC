#ifndef ACC96_SAMPLELEFT_KAT_H
#define ACC96_SAMPLELEFT_KAT_H

#include "acc96_accumulator.h"

/* Algebraic scaffold for the q96 truncated SampleLeft implementation.
 * It uses an exact balanced-base-13 gadget preimage but omits Gaussian
 * perturbation and randomized Klein sampling.  It is suitable for relation
 * KATs only and MUST NOT be used by production Issue. */
int acc96_sampleleft_algebraic_kat(
    acc96_relation_witness_t *out,
    const acc96_relation_pp_t *pp,
    const acc96_trapdoor_t *trapdoor,
    const acc96_poly_t *accumulator,
    uint32_t handle,
    void *opaque,
    const acc96_ctx_t *ctx);

#endif
