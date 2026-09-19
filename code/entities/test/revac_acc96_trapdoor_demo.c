#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "acc96_trapdoor.h"
#include "acc96_trapdoor_cert.h"
#include "acc96_perturbation_cert.h"
#include "acc96_truncated_params_generated.h"

static double seconds(void) {
  struct timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);
  return now.tv_sec + now.tv_nsec / 1e9;
}

int main(void) {
  const uint8_t pp_seed[ACC96_SEED_BYTES] = {0x96, 0x50, 0x50, 0x01};
  const uint8_t trapdoor_seed[ACC96_SEED_BYTES] = {0x96, 0x54, 0x44, 0x01};
  acc96_ctx_t ctx;
  acc96_relation_pp_t *pp = calloc(1, sizeof(*pp));
  acc96_trapdoor_t *trapdoor = calloc(1, sizeof(*trapdoor));
  if (pp == NULL || trapdoor == NULL || !acc96_ctx_init(&ctx) ||
      !acc96_relation_pp_expand(pp, pp_seed, &ctx))
    return 1;
  const double started = seconds();
  if (!acc96_trapdoor_generate(trapdoor, pp, trapdoor_seed, &ctx) ||
      trapdoor->squared_spectral_norm >
          ACC96_TRAPDOOR_NUMERIC_SAFETY *
              ACC_TRUNC_R_MAX_SQ_SPECTRAL_NORM ||
      !acc96_trapdoor_public_relation_holds(trapdoor, pp, &ctx))
    return 1;
  const double acb_started = seconds();
  double acb_upper = 0;
  if (!acc96_trapdoor_spectral_certify_acb(trapdoor, &ctx, 160,
                                           &acb_upper))
    return 1;
  const double acb_time = seconds() - acb_started;
  const double covariance_started = seconds();
  acc96_perturbation_covariance_cert_t covariance_certificate;
  if (!acc96_perturbation_covariance_certify_acb(
          trapdoor, &ctx, 192, &covariance_certificate))
    return 1;
  const double covariance_time = seconds() - covariance_started;
  const acc96_coeff_t saved = trapdoor->r[0][0].coeff[0];
  acc96_coeff_from_u128(&trapdoor->r[0][0].coeff[0], 2);
  double rejected_upper = 0;
  if (acc96_trapdoor_spectral_certify_acb(trapdoor, &ctx, 160,
                                          &rejected_upper))
    return 1;
  trapdoor->r[0][0].coeff[0] = saved;
  acc96_trapdoor_t *adversarial = calloc(1, sizeof(*adversarial));
  if (adversarial == NULL) return 1;
  for (size_t column = 0; column < ACC96_TRUNC_RIGHT; column++)
    for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++)
      acc96_coeff_from_u128(
          &adversarial->r[0][column].coeff[coefficient], 1);
  acc96_perturbation_covariance_cert_t rejected_covariance;
  if (acc96_perturbation_covariance_certify_acb(
          adversarial, &ctx, 192, &rejected_covariance))
    return 1;
  free(adversarial);
  printf("q96 truncated trapdoor passed: spectral norm^2 %.3f <= "
         "guarded %.3f (theorem %.3f), "
         "rejections %u, setup %.3f s, ACB-160 upper %.3f in %.3f s, "
         "covariance ACB-192 pivots %.3e / %.3e / %.3e in %.3f s; "
         "B=A*R1+R0 checked in all 24 columns; malformed and "
         "non-positive covariance trapdoors rejected.\n",
         trapdoor->squared_spectral_norm,
         ACC96_TRAPDOOR_NUMERIC_SAFETY *
             ACC_TRUNC_R_MAX_SQ_SPECTRAL_NORM,
         ACC_TRUNC_R_MAX_SQ_SPECTRAL_NORM,
         trapdoor->rejection_counter, acb_started - started,
         acb_upper, acb_time,
         covariance_certificate.minimum_second_pivot_lower_bound,
         covariance_certificate.minimum_adjusted_schur_lower_bound,
         covariance_certificate.independent_variance_lower_bound,
         covariance_time);
  free(trapdoor);
  free(pp);
  return 0;
}
