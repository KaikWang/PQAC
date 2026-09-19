#ifndef REVAC_ACC96_PERTURBATION_CERT_H
#define REVAC_ACC96_PERTURBATION_CERT_H

#include "acc96_trapdoor.h"

typedef struct {
  long precision_bits;
  double minimum_second_pivot_lower_bound;
  double minimum_adjusted_schur_lower_bound;
  double independent_variance_lower_bound;
  double maximum_row_energy_upper_bound;
} acc96_perturbation_covariance_cert_t;

/* Certify, with outward-rounded ACB/Arb arithmetic, the covariance objects
 * used by the rank-one q96 truncated-gadget perturbation sampler.  In every
 * canonical embedding this checks
 *
 *   C = diag(s1^2,s3^2) + c R R*,  c < 0,
 *   C_11 > 0,
 *   C_00 - sG^2 - |C_01|^2/C_11 > 0,
 *   s2^2 - sG^2 > 0.
 *
 * The last two expressions are precisely the covariances consumed after the
 * Schur elimination in SamplePerturb.  This is a setup/release check, not an
 * Issue/show hot-path operation. */
int acc96_perturbation_covariance_certify_acb(
    const acc96_trapdoor_t *trapdoor, const acc96_ctx_t *ctx,
    long precision_bits, acc96_perturbation_covariance_cert_t *certificate);

#endif
