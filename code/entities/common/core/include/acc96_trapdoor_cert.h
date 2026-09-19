#ifndef REVAC_ACC96_TRAPDOOR_CERT_H
#define REVAC_ACC96_TRAPDOOR_CERT_H

#include "acc96_trapdoor.h"

/* Recompute the canonical-embedding spectral norm with ACB ball arithmetic.
 * Returns one only when every interval is provably below the published bound.
 * This is a setup/release certification path, not the fast rejection filter. */
int acc96_trapdoor_spectral_certify_acb(
    const acc96_trapdoor_t *trapdoor, const acc96_ctx_t *ctx,
    long precision_bits, double *certified_upper_bound);

#endif
