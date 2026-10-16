#ifndef REVAC_ACC96_SAMPLELEFT_H
#define REVAC_ACC96_SAMPLELEFT_H

#include "acc96_accumulator.h"

typedef struct acc96_sampleleft_ctx acc96_sampleleft_ctx_t;

/* Build the q96 truncated-gadget sampler's immutable covariance/Klein state.
 * Construction validates all interval inversions and is intended once per TA
 * Setup.  thread_count==0 selects the available OpenMP worker count, capped
 * by the number of independent coefficient chains. */
acc96_sampleleft_ctx_t *acc96_sampleleft_ctx_new(
    const acc96_trapdoor_t *trapdoor, const acc96_ctx_t *ctx,
    long precision_bits, size_t thread_count);
void acc96_sampleleft_ctx_free(acc96_sampleleft_ctx_t *sampler);

/* Production q96 SampleLeft callback.  It samples the right-hand hiding
 * vector, subtracts its public contribution, applies the certified
 * perturbation distribution and finishes with exact-rational Klein. Success
 * also requires acc96_relation_fresh_witness_check; on any failure out is
 * cleared. There is no implicit retry or coefficient clipping. */
int acc96_sampleleft_secure(acc96_relation_witness_t *out,
                            const acc96_relation_pp_t *pp,
                            const acc96_trapdoor_t *trapdoor,
                            const acc96_poly_t *accumulator,
                            uint32_t handle, void *opaque,
                            const acc96_ctx_t *ctx);

size_t acc96_sampleleft_thread_count(const acc96_sampleleft_ctx_t *sampler);
double acc96_sampleleft_maximum_arb_radius(
    const acc96_sampleleft_ctx_t *sampler);

#endif
