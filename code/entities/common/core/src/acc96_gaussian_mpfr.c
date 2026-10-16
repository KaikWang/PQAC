#include "acc96_gaussian_mpfr.h"

#include <mpfr.h>
#include <limits.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "randombytes.h"

static uint64_t screen_hits, screen_fallbacks, screen_audits;
__attribute__((destructor)) static void screen_report(void) {
  const char *mode = getenv("REVAC_GAUSSIAN_CERTIFIED_SCREEN");
  if ((mode && mode[0] == '2') || getenv("REVAC_GAUSSIAN_SCREEN_REPORT"))
    fprintf(stderr, "[certified Gaussian screen] hits=%llu fallbacks=%llu reference_audits=%llu\n",
      (unsigned long long)screen_hits, (unsigned long long)screen_fallbacks,
      (unsigned long long)screen_audits);
}

typedef struct acc96_gaussian_mpfr_workspace {
  struct acc96_gaussian_mpfr_workspace *screen;
  mpfr_t center_lo, center_hi;
  mpfr_t sigma_lo, sigma_hi;
  mpfr_t difference_lo, difference_hi;
  mpfr_t square_lo, square_hi;
  mpfr_t probability_lo, probability_hi;
  mpfr_t temporary_lo, temporary_hi;
  mpfr_t pi_lo, pi_hi;
  mpfr_t alpha_lo, alpha_hi;
  mpfr_t cached_variance;
  int cached_variance_valid;
  mpz_t threshold_lo, threshold_hi;
  mpz_t uniform;
  char sigma_decimal[160];
  char proposal_sigma[ACC96_GAUSSIAN_PROPOSAL_CACHE][160];
  uint32_t proposal_weight[ACC96_GAUSSIAN_PROPOSAL_CACHE]
                          [2 * ACC96_GAUSSIAN_SUPPORT + 1];
  uint32_t proposal_total[ACC96_GAUSSIAN_PROPOSAL_CACHE];
  uint32_t proposal_count;
  uint32_t active_proposal;
  mpfr_prec_t precision;
  uint8_t random_buffer[4096];
  size_t random_offset;
} acc96_gaussian_mpfr_workspace_t;

static void workspace_set_precision(acc96_gaussian_mpfr_workspace_t *workspace,
                                    mpfr_prec_t precision) {
#define SET_PRECISION(name) mpfr_set_prec(workspace->name, precision)
  SET_PRECISION(center_lo);
  SET_PRECISION(center_hi);
  SET_PRECISION(sigma_lo);
  SET_PRECISION(sigma_hi);
  SET_PRECISION(difference_lo);
  SET_PRECISION(difference_hi);
  SET_PRECISION(square_lo);
  SET_PRECISION(square_hi);
  SET_PRECISION(probability_lo);
  SET_PRECISION(probability_hi);
  SET_PRECISION(temporary_lo);
  SET_PRECISION(temporary_hi);
  SET_PRECISION(pi_lo);
  SET_PRECISION(pi_hi);
  SET_PRECISION(alpha_lo);
  SET_PRECISION(alpha_hi);
#undef SET_PRECISION
  workspace->precision = precision;
  workspace->sigma_decimal[0] = '\0';
  workspace->cached_variance_valid = 0;
  mpfr_const_pi(workspace->pi_lo, MPFR_RNDD);
  mpfr_const_pi(workspace->pi_hi, MPFR_RNDU);
}

static acc96_gaussian_mpfr_workspace_t *workspace_new(void) {
  acc96_gaussian_mpfr_workspace_t *workspace = calloc(1, sizeof(*workspace));
  if (workspace == NULL) return NULL;
  mpfr_inits2(ACC96_GAUSSIAN_MPFR_PRECISION_BITS,
              workspace->center_lo, workspace->center_hi,
              workspace->sigma_lo, workspace->sigma_hi,
              workspace->difference_lo, workspace->difference_hi,
              workspace->square_lo, workspace->square_hi,
              workspace->probability_lo, workspace->probability_hi,
              workspace->temporary_lo, workspace->temporary_hi,
              workspace->pi_lo, workspace->pi_hi,
              workspace->alpha_lo, workspace->alpha_hi, workspace->cached_variance,
              (mpfr_ptr)0);
  mpz_inits(workspace->threshold_lo, workspace->threshold_hi,
            workspace->uniform, NULL);
  workspace->precision = ACC96_GAUSSIAN_MPFR_PRECISION_BITS;
  mpfr_const_pi(workspace->pi_lo, MPFR_RNDD);
  mpfr_const_pi(workspace->pi_hi, MPFR_RNDU);
  workspace->random_offset = sizeof(workspace->random_buffer);
  return workspace;
}

static void workspace_free(acc96_gaussian_mpfr_workspace_t *workspace) {
  if (workspace == NULL) return;
  workspace_free(workspace->screen);
  mpz_clears(workspace->uniform, workspace->threshold_hi,
             workspace->threshold_lo, NULL);
  mpfr_clears(workspace->cached_variance, workspace->alpha_hi, workspace->alpha_lo,
              workspace->pi_hi, workspace->pi_lo,
              workspace->temporary_hi, workspace->temporary_lo,
              workspace->probability_hi, workspace->probability_lo,
              workspace->square_hi, workspace->square_lo,
              workspace->difference_hi, workspace->difference_lo,
              workspace->sigma_hi, workspace->sigma_lo,
              workspace->center_hi, workspace->center_lo,
              (mpfr_ptr)0);
  free(workspace);
}

void acc96_gaussian_mpfr_stats_clear(acc96_gaussian_mpfr_stats_t *stats) {
  if (stats == NULL) return;
  workspace_free((acc96_gaussian_mpfr_workspace_t *)stats->workspace);
  stats->workspace = NULL;
}

static void workspace_random(acc96_gaussian_mpfr_workspace_t *workspace,
                             uint8_t *out, size_t length) {
  while (length > 0) {
    size_t available;
    size_t take;
    if (workspace->random_offset == sizeof(workspace->random_buffer)) {
      randombytes(workspace->random_buffer, sizeof(workspace->random_buffer));
      workspace->random_offset = 0;
    }
    available = sizeof(workspace->random_buffer) - workspace->random_offset;
    take = length < available ? length : available;
    memcpy(out, workspace->random_buffer + workspace->random_offset, take);
    workspace->random_offset += take;
    out += take;
    length -= take;
  }
}

static uint32_t uniform_below(acc96_gaussian_mpfr_workspace_t *workspace,
                              uint32_t bound) {
  const uint64_t range = (uint64_t)UINT32_MAX + 1;
  const uint64_t limit = range - range % bound;
  uint32_t value;
  do {
    uint8_t bytes[sizeof(value)];
    workspace_random(workspace, bytes, sizeof(bytes));
    memcpy(&value, bytes, sizeof(value));
  } while ((uint64_t)value >= limit);
  return value % bound;
}

static int weighted_support_sample(acc96_gaussian_mpfr_workspace_t *workspace,
                                   int64_t *value, uint32_t *weight) {
  const uint32_t table = workspace->active_proposal;
  uint32_t choice;
  if (table >= workspace->proposal_count ||
      workspace->proposal_total[table] == 0)
    return 0;
  choice = uniform_below(workspace, workspace->proposal_total[table]);
  for (uint32_t index = 0; index < 2 * ACC96_GAUSSIAN_SUPPORT + 1;
       index++) {
    const uint32_t candidate_weight =
        workspace->proposal_weight[table][index];
    if (choice < candidate_weight) {
      *value = (int64_t)index - ACC96_GAUSSIAN_SUPPORT;
      *weight = candidate_weight;
      return 1;
    }
    choice -= candidate_weight;
  }
  return 0;
}

/* q_z is proportional to ceil(S rho_max(z)), where rho_max takes the maximum
 * over every fractional center in [0,1).  Hence weight[z]/S >= rho(z-c) for
 * every valid center c and acceptance S*rho/weight never exceeds one. */
static int prepare_proposal_table(
    acc96_gaussian_mpfr_workspace_t *workspace, const char *sigma_decimal) {
  uint32_t table;
  uint64_t total = 0;
  for (table = 0; table < workspace->proposal_count; table++) {
    if (strcmp(workspace->proposal_sigma[table], sigma_decimal) == 0) {
      workspace->active_proposal = table;
      return 1;
    }
  }
  if (workspace->proposal_count >= ACC96_GAUSSIAN_PROPOSAL_CACHE)
    return 0;
  table = workspace->proposal_count;
  for (int64_t z = -ACC96_GAUSSIAN_SUPPORT;
       z <= ACC96_GAUSSIAN_SUPPORT; z++) {
    const uint64_t distance = z < 0 ? (uint64_t)-z
                              : z > 1 ? (uint64_t)(z - 1)
                                      : 0;
    uint64_t value;
    mpfr_set_ui(workspace->temporary_lo, distance, MPFR_RNDD);
    mpfr_sqr(workspace->temporary_lo, workspace->temporary_lo, MPFR_RNDD);
    mpfr_mul(workspace->temporary_lo, workspace->temporary_lo,
             workspace->alpha_lo, MPFR_RNDD);
    mpfr_neg(workspace->temporary_lo, workspace->temporary_lo, MPFR_RNDU);
    mpfr_exp(workspace->temporary_lo, workspace->temporary_lo, MPFR_RNDU);
    mpfr_mul_ui(workspace->temporary_lo, workspace->temporary_lo,
                ACC96_GAUSSIAN_PROPOSAL_SCALE, MPFR_RNDU);
    value = mpfr_get_ui(workspace->temporary_lo, MPFR_RNDU);
    if (value == 0) value = 1;
    if (value > UINT32_MAX || total + value > UINT32_MAX) return 0;
    workspace->proposal_weight[table]
                              [z + ACC96_GAUSSIAN_SUPPORT] = (uint32_t)value;
    total += value;
  }
  if (strlen(sigma_decimal) >= sizeof(workspace->proposal_sigma[table]))
    return 0;
  strcpy(workspace->proposal_sigma[table], sigma_decimal);
  workspace->proposal_total[table] = (uint32_t)total;
  workspace->active_proposal = table;
  workspace->proposal_count++;
  return 1;
}

static int prepare_distribution_interval(
    acc96_gaussian_mpfr_workspace_t *workspace,
    const mpz_t center_numerator, const mpz_t center_denominator,
    const char *sigma_decimal) {
  workspace->cached_variance_valid = 0;
  mpfr_set_z(workspace->temporary_lo, center_numerator, MPFR_RNDD);
  mpfr_set_z(workspace->temporary_hi, center_numerator, MPFR_RNDU);
  mpfr_set_z(workspace->difference_lo, center_denominator, MPFR_RNDD);
  mpfr_set_z(workspace->difference_hi, center_denominator, MPFR_RNDU);
  mpfr_div(workspace->center_lo, workspace->temporary_lo,
           workspace->difference_hi, MPFR_RNDD);
  mpfr_div(workspace->center_hi, workspace->temporary_hi,
           workspace->difference_lo, MPFR_RNDU);

  if (strcmp(workspace->sigma_decimal, sigma_decimal) != 0) {
    if (strlen(sigma_decimal) >= sizeof(workspace->sigma_decimal) ||
        mpfr_set_str(workspace->sigma_lo, sigma_decimal, 10, MPFR_RNDD) != 0 ||
        mpfr_set_str(workspace->sigma_hi, sigma_decimal, 10, MPFR_RNDU) != 0 ||
        mpfr_sgn(workspace->sigma_lo) <= 0) {
      return 0;
    }
    strcpy(workspace->sigma_decimal, sigma_decimal);
    mpfr_sqr(workspace->temporary_lo, workspace->sigma_lo, MPFR_RNDD);
    mpfr_sqr(workspace->temporary_hi, workspace->sigma_hi, MPFR_RNDU);
    mpfr_div(workspace->alpha_lo, workspace->pi_lo,
             workspace->temporary_hi, MPFR_RNDD);
    mpfr_div(workspace->alpha_hi, workspace->pi_hi,
             workspace->temporary_lo, MPFR_RNDU);
  }
  return prepare_proposal_table(workspace, sigma_decimal);
}

/* Enclose p = exp(-pi (candidate-center)^2/sigma^2) * W/(M*w).
 * The exact center and decimal sigma parameter were enclosed once for this
 * SampleZ invocation by prepare_distribution_interval(). */
static int acceptance_interval(
    acc96_gaussian_mpfr_workspace_t *workspace, int64_t candidate,
    uint32_t proposal_weight) {

  mpfr_set_si(workspace->difference_lo, candidate, MPFR_RNDD);
  mpfr_set_si(workspace->difference_hi, candidate, MPFR_RNDU);
  mpfr_sub(workspace->difference_lo, workspace->difference_lo,
           workspace->center_hi, MPFR_RNDD);
  mpfr_sub(workspace->difference_hi, workspace->difference_hi,
           workspace->center_lo, MPFR_RNDU);

  if (mpfr_sgn(workspace->difference_lo) >= 0) {
    mpfr_sqr(workspace->square_lo, workspace->difference_lo, MPFR_RNDD);
    mpfr_sqr(workspace->square_hi, workspace->difference_hi, MPFR_RNDU);
  } else if (mpfr_sgn(workspace->difference_hi) <= 0) {
    mpfr_sqr(workspace->square_lo, workspace->difference_hi, MPFR_RNDD);
    mpfr_sqr(workspace->square_hi, workspace->difference_lo, MPFR_RNDU);
  } else {
    mpfr_set_zero(workspace->square_lo, 0);
    mpfr_sqr(workspace->temporary_lo, workspace->difference_lo, MPFR_RNDU);
    mpfr_sqr(workspace->temporary_hi, workspace->difference_hi, MPFR_RNDU);
    mpfr_max(workspace->square_hi, workspace->temporary_lo,
             workspace->temporary_hi, MPFR_RNDU);
  }

  mpfr_mul(workspace->temporary_lo, workspace->square_lo,
           workspace->alpha_lo, MPFR_RNDD);
  mpfr_mul(workspace->temporary_hi, workspace->square_hi,
           workspace->alpha_hi, MPFR_RNDU);
  mpfr_neg(workspace->temporary_lo, workspace->temporary_lo, MPFR_RNDU);
  mpfr_neg(workspace->temporary_hi, workspace->temporary_hi, MPFR_RNDD);
  /* Negation reverses the interval: -upper is the lower exponent. */
  mpfr_swap(workspace->temporary_lo, workspace->temporary_hi);
  mpfr_exp(workspace->probability_lo, workspace->temporary_lo, MPFR_RNDD);
  mpfr_exp(workspace->probability_hi, workspace->temporary_hi, MPFR_RNDU);
  mpfr_mul_ui(workspace->probability_lo, workspace->probability_lo,
              ACC96_GAUSSIAN_PROPOSAL_SCALE, MPFR_RNDD);
  mpfr_mul_ui(workspace->probability_hi, workspace->probability_hi,
              ACC96_GAUSSIAN_PROPOSAL_SCALE, MPFR_RNDU);
  mpfr_div_ui(workspace->probability_lo, workspace->probability_lo,
              proposal_weight, MPFR_RNDD);
  mpfr_div_ui(workspace->probability_hi, workspace->probability_hi,
              proposal_weight, MPFR_RNDU);
  if (mpfr_sgn(workspace->probability_lo) < 0 ||
      mpfr_cmp_ui(workspace->probability_hi, 1) > 0 ||
      mpfr_cmp(workspace->probability_lo, workspace->probability_hi) > 0)
    return 0;
  mpfr_mul_2ui(workspace->probability_lo, workspace->probability_lo,
               ACC96_GAUSSIAN_UNIFORM_BITS, MPFR_RNDD);
  mpfr_mul_2ui(workspace->probability_hi, workspace->probability_hi,
               ACC96_GAUSSIAN_UNIFORM_BITS, MPFR_RNDU);
  mpfr_get_z(workspace->threshold_lo, workspace->probability_lo, MPFR_RNDU);
  mpfr_get_z(workspace->threshold_hi, workspace->probability_hi, MPFR_RNDU);
  return 1;
}

static int prepare_point_distribution(
    acc96_gaussian_mpfr_workspace_t *workspace,
    mpfr_srcptr center, mpfr_srcptr variance) {
  if (mpfr_sgn(variance) <= 0 || !mpfr_number_p(center) ||
      !mpfr_number_p(variance))
    return 0;
  mpfr_set(workspace->center_lo, center, MPFR_RNDD);
  mpfr_set(workspace->center_hi, center, MPFR_RNDU);
  if (workspace->cached_variance_valid && mpfr_equal_p(workspace->cached_variance, variance))
    return 1;
  workspace->sigma_decimal[0] = '\0'; /* Point alpha cannot serve a rational sigma cache. */
  mpfr_sqrt(workspace->sigma_lo, variance, MPFR_RNDD);
  mpfr_sqrt(workspace->sigma_hi, variance, MPFR_RNDU);
  if (mpfr_sgn(workspace->sigma_lo) <= 0) return 0;
  mpfr_sqr(workspace->temporary_lo, workspace->sigma_lo, MPFR_RNDD);
  mpfr_sqr(workspace->temporary_hi, workspace->sigma_hi, MPFR_RNDU);
  mpfr_div(workspace->alpha_lo, workspace->pi_lo,
           workspace->temporary_hi, MPFR_RNDD);
  mpfr_div(workspace->alpha_hi, workspace->pi_hi,
           workspace->temporary_lo, MPFR_RNDU);
  mpfr_set_prec(workspace->cached_variance, mpfr_get_prec(variance));
  mpfr_set(workspace->cached_variance, variance, MPFR_RNDN); /* Exact same-precision copy. */
  workspace->cached_variance_valid = 1;
  return 1;
}

int acc96_gaussian_mpfr_sample_center_variance(
    int64_t *sample, mpfr_srcptr center, mpfr_srcptr variance,
    void *opaque) {
  acc96_gaussian_mpfr_stats_t *stats = opaque;
  acc96_gaussian_mpfr_workspace_t *workspace = NULL;
  int transient_workspace = 0;
  uint8_t random[ACC96_GAUSSIAN_UNIFORM_BITS / 8];
  int64_t base;
  uint32_t radius;
  uint32_t proposal_bound;
  int ok = 0;

  if (sample == NULL || center == NULL || variance == NULL ||
      !mpfr_fits_slong_p(center, MPFR_RNDD))
    return 0;
  if (stats != NULL) {
    workspace = (acc96_gaussian_mpfr_workspace_t *)stats->workspace;
    if (workspace == NULL) {
      workspace = workspace_new();
      stats->workspace = workspace;
    }
  } else {
    workspace = workspace_new();
    transient_workspace = 1;
  }
  if (workspace == NULL) return 0;
  if (stats != NULL) stats->calls++;
  if (!prepare_point_distribution(workspace, center, variance)) goto cleanup;

  base = (int64_t)mpfr_get_si(center, MPFR_RNDD);
  mpfr_mul_ui(workspace->temporary_hi, workspace->sigma_hi, 12,
              MPFR_RNDU);
  if (!mpfr_fits_uint_p(workspace->temporary_hi, MPFR_RNDU)) goto cleanup;
  radius = (uint32_t)mpfr_get_ui(workspace->temporary_hi, MPFR_RNDU);
  if (radius > (UINT32_MAX - 1U) / 2U - 2U) goto cleanup;
  radius += 2U;
  proposal_bound = 2U * radius + 1U;

  const char *screen_mode = getenv("REVAC_GAUSSIAN_CERTIFIED_SCREEN");
  const int use_screen = screen_mode == NULL || screen_mode[0] != '0';
  const int audit_screen = screen_mode != NULL && screen_mode[0] == '2';
  const int report_screen = audit_screen || getenv("REVAC_GAUSSIAN_SCREEN_REPORT") != NULL;
  if (use_screen) {
    if (workspace->screen == NULL) {
      workspace->screen = workspace_new();
      if (workspace->screen == NULL) goto cleanup;
      workspace_set_precision(workspace->screen, 64);
    }
    if (!prepare_point_distribution(workspace->screen, center, variance))
      goto cleanup;
  }

  for (;;) {
    const uint32_t offset = uniform_below(workspace, proposal_bound);
    int64_t candidate;
    if (offset <= radius) {
      if (base < INT64_MIN + (int64_t)(radius - offset)) goto cleanup;
      candidate = base - (int64_t)(radius - offset);
    } else {
      if (base > INT64_MAX - (int64_t)(offset - radius)) goto cleanup;
      candidate = base + (int64_t)(offset - radius);
    }
    if (stats != NULL) stats->proposals++;
    workspace_random(workspace, random, sizeof(random));
    mpz_import(workspace->uniform, sizeof(random), 1, 1, 1, 0, random);
    int screened = 0;
    if (use_screen) {
      /* RNDD/RNDU encloses the SAME exact p. A separated uniform value
       * proves the decision; an overlapping interval falls back unchanged.
       * Proposal support, source MPFR values and all RNG draws are original. */
      if (acceptance_interval(workspace->screen, candidate,
                              ACC96_GAUSSIAN_PROPOSAL_SCALE)) {
        if (mpz_cmp(workspace->uniform, workspace->screen->threshold_lo) < 0)
          screened = 1;
        else if (mpz_cmp(workspace->uniform, workspace->screen->threshold_hi) >= 0)
          screened = -1;
      } /* A failed coarse enclosure also falls back to the original path. */
      if (report_screen) __atomic_fetch_add(screened ? &screen_hits : &screen_fallbacks, 1, __ATOMIC_RELAXED);
      if (screened != 0 && !audit_screen) {
        if (stats != NULL && workspace->precision > stats->maximum_precision_bits)
          stats->maximum_precision_bits = (uint32_t)workspace->precision;
        if (screened == 1) {
          if (stats != NULL) stats->certified_accepts++;
          *sample = candidate; ok = 1; goto cleanup;
        }
        if (stats != NULL) stats->certified_rejects++;
        continue;
      }
    }
    for (;;) {
      /* A uniform proposal has constant mass, so accepting with rho is
       * obtained by setting the existing envelope ratio M/weight to one. */
      if (!acceptance_interval(workspace, candidate,
                               ACC96_GAUSSIAN_PROPOSAL_SCALE))
        goto cleanup;
      if (stats != NULL && workspace->precision > stats->maximum_precision_bits)
        stats->maximum_precision_bits = (uint32_t)workspace->precision;
      if (mpz_cmp(workspace->uniform, workspace->threshold_lo) < 0) {
        if (audit_screen && screened == -1) abort();
        if (audit_screen && screened) __atomic_fetch_add(&screen_audits, 1, __ATOMIC_RELAXED);
        if (stats != NULL) stats->certified_accepts++;
        *sample = candidate;
        ok = 1;
        goto cleanup;
      }
      if (mpz_cmp(workspace->uniform, workspace->threshold_hi) >= 0) {
        if (audit_screen && screened == 1) abort();
        if (audit_screen && screened) __atomic_fetch_add(&screen_audits, 1, __ATOMIC_RELAXED);
        if (stats != NULL) stats->certified_rejects++;
        break;
      }
      if (workspace->precision >= 6144) goto cleanup;
      workspace_set_precision(workspace, workspace->precision * 2);
      if (stats != NULL) stats->interval_refinements++;
      if (!prepare_point_distribution(workspace, center, variance))
        goto cleanup;
    }
  }

cleanup:
  if (transient_workspace) workspace_free(workspace);
  return ok;
}

int acc96_gaussian_mpfr_sample_rational(
    int64_t *sample, const mpz_t center_numerator,
    const mpz_t center_denominator, const char *sigma_decimal,
    void *opaque) {
  acc96_gaussian_mpfr_stats_t *stats = opaque;
  acc96_gaussian_mpfr_workspace_t *workspace = NULL;
  int transient_workspace = 0;
  uint8_t random[ACC96_GAUSSIAN_UNIFORM_BITS / 8];
  int ok = 0;

  if (sample == NULL || center_numerator == NULL ||
      center_denominator == NULL || mpz_sgn(center_denominator) <= 0 ||
      mpz_sgn(center_numerator) < 0 ||
      mpz_cmp(center_numerator, center_denominator) >= 0 ||
      sigma_decimal == NULL) {
    return 0;
  }
  if (stats != NULL) {
    workspace = (acc96_gaussian_mpfr_workspace_t *)stats->workspace;
    if (workspace == NULL) {
      workspace = workspace_new();
      stats->workspace = workspace;
    }
  } else {
    workspace = workspace_new();
    transient_workspace = 1;
  }
  if (workspace == NULL) return 0;
  if (stats != NULL) stats->calls++;
  if (!prepare_distribution_interval(workspace, center_numerator,
                                     center_denominator, sigma_decimal))
    goto cleanup;

  for (;;) {
    int64_t candidate;
    uint32_t proposal_weight;
    if (!weighted_support_sample(workspace, &candidate, &proposal_weight))
      goto cleanup;
    if (stats != NULL) stats->proposals++;
    workspace_random(workspace, random, sizeof(random));
    mpz_import(workspace->uniform, sizeof(random), 1, 1, 1, 0, random);
    for (;;) {
      if (!acceptance_interval(workspace, candidate, proposal_weight))
        goto cleanup;
      if (stats != NULL && workspace->precision > stats->maximum_precision_bits)
        stats->maximum_precision_bits = (uint32_t)workspace->precision;
      if (mpz_cmp(workspace->uniform, workspace->threshold_lo) < 0) {
        if (stats != NULL) stats->certified_accepts++;
        *sample = candidate;
        ok = 1;
        goto cleanup;
      }
      if (mpz_cmp(workspace->uniform, workspace->threshold_hi) >= 0) {
        if (stats != NULL) stats->certified_rejects++;
        break;
      }
      if (workspace->precision >= 6144) goto cleanup;
      workspace_set_precision(workspace, workspace->precision * 2);
      if (stats != NULL) stats->interval_refinements++;
      if (!prepare_distribution_interval(workspace, center_numerator,
                                         center_denominator, sigma_decimal))
        goto cleanup;
    }
  }

cleanup:
  if (transient_workspace) workspace_free(workspace);
  return ok;
}
