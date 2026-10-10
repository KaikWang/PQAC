#ifndef REVAC_LABRADOR_PROFILE_H
#define REVAC_LABRADOR_PROFILE_H

#include <stddef.h>
#include <stdlib.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Runtime-only helpers used by the revocable-credential integration.  The
 * profiler and the optional OpenMP team size do not alter transcripts,
 * parameters, or rejection decisions. */
static inline int revac_labrador_profile_enabled(void) {
  const char *value = getenv("REVAC_LABRADOR_PROFILE");
  return value != NULL && value[0] != '\0' && value[0] != '0';
}

static inline int revac_labrador_profile_detail_enabled(void) {
  const char *value = getenv("REVAC_LABRADOR_PROFILE_DETAIL");
  return revac_labrador_profile_enabled() &&
         value != NULL && value[0] != '\0' && value[0] != '0';
}

static inline double revac_labrador_profile_seconds(void) {
  struct timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);
  return (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
}

static inline int revac_labrador_threads(void) {
#ifdef _OPENMP
  const char *value = getenv("REVAC_LABRADOR_THREADS");
  char *end = NULL;
  long requested = 8;

  if (value != NULL && value[0] != '\0') {
    requested = strtol(value, &end, 10);
    if (end == value || *end != '\0' || requested < 1) {
      requested = 8;
    }
  }
  if (requested > omp_get_max_threads()) {
    requested = omp_get_max_threads();
  }
  return (int)requested;
#else
  return 1;
#endif
}

/* Keep the coarse-grained fold parallelism from recursively creating a new
 * OpenMP team in lower-level vector helpers.  Nested teams are particularly
 * expensive for refresh(), which operates on many 32-polynomial chunks. */
static inline int revac_labrador_parallel_work(size_t work,
                                               size_t min_work) {
#ifdef _OPENMP
  return revac_labrador_threads() > 1 && work >= min_work &&
         !omp_in_parallel();
#else
  (void)work;
  (void)min_work;
  return 0;
#endif
}

static inline int revac_labrador_team_threads(size_t independent_work) {
  const int configured = revac_labrador_threads();
  if (independent_work == 0) return 1;
  return independent_work < (size_t)configured ? (int)independent_work
                                                : configured;
}

/* LNP rejection checks only read immutable mask/shift polynomials and reduce
 * integer products.  They are independent of the transcript-sensitive
 * reducers guarded by revac_labrador_threads(), so they get a separate,
 * opt-in team size.  Keep the default serial: this adds no OpenMP team
 * creation to production unless the experiment explicitly enables it. */
static inline int revac_labrador_lnp_norm_threads(size_t work) {
#ifdef _OPENMP
  const char *value = getenv("REVAC_LABRADOR_LNP_NORM_THREADS");
  char *end = NULL;
  long requested = 1;

  if (value != NULL && value[0] != '\0') {
    requested = strtol(value, &end, 10);
    if (end == value || *end != '\0' || requested < 1)
      requested = 1;
  }
  if (requested > omp_get_max_threads())
    requested = omp_get_max_threads();
  if (work != 0 && (size_t)requested > work)
    requested = (long)work;
  return (int)requested;
#else
  (void)work;
  return 1;
#endif
}

/* Long CRT conversion helpers used to be serial upstream.  Keep their
 * experimental block parallelism opt-in independently so randomized Pack
 * failures can be attributed to one kernel without disabling unrelated,
 * known-disjoint work.  REVAC_LABRADOR_CONVERSION_PARALLEL enables all five
 * helpers; a kernel-specific variable overrides the group setting. */
static inline int revac_labrador_conversion_parallel_enabled(
    const char *kernel_variable) {
#ifdef _OPENMP
  const char *value = getenv(kernel_variable);
  if (value == NULL || value[0] == '\0')
    value = getenv("REVAC_LABRADOR_CONVERSION_PARALLEL");
  return value != NULL && value[0] != '\0' && value[0] != '0';
#else
  (void)kernel_variable;
  return 0;
#endif
}

/* The garbage-term reducers use approximate CRT/NTT accumulators whose
 * current implementation is not reentrant across independent rows.  Keep
 * them serial by default; the switch exists only for controlled diagnostics. */
static inline int revac_labrador_garbage_parallel_enabled(void) {
  const char *value = getenv("REVAC_LABRADOR_GARBAGE_PARALLEL");
  return value != NULL && value[0] != '\0' && value[0] != '0';
}

/* Rotation aggregation writes each degree-sized output block independently.
 * It is therefore safe to split large blocks across a team, unlike the
 * approximate garbage reducers above.  Keep an opt-out for regression and
 * single-core measurements. */
static inline int revac_labrador_rotation_parallel_enabled(void) {
#ifdef _OPENMP
  const char *value = getenv("REVAC_LABRADOR_ROTATION_PARALLEL");
  return value == NULL || value[0] == '\0' || value[0] != '0';
#else
  return 0;
#endif
}

static inline int revac_labrador_rotation_threads(size_t work) {
#ifdef _OPENMP
  const char *value = getenv("REVAC_LABRADOR_ROTATION_THREADS");
  char *end = NULL;
  long requested = revac_labrador_threads();

  if (value != NULL && value[0] != '\0') {
    requested = strtol(value, &end, 10);
    if (end == value || *end != '\0' || requested < 1)
      requested = 1;
  }
  if (requested > omp_get_max_threads())
    requested = omp_get_max_threads();
  if (work != 0 && (size_t)requested > work)
    requested = (long)work;
  return (int)requested;
#else
  (void)work;
  return 1;
#endif
}

/* Clearing a long CRT vector touches six independent prime projections.
 * Parallelising only large vectors avoids team overhead in small rounds. */
static inline int revac_labrador_zero_parallel_enabled(void) {
#ifdef _OPENMP
  const char *value = getenv("REVAC_LABRADOR_ZERO_PARALLEL");
  return value != NULL && value[0] != '\0' && value[0] != '0';
#else
  return 0;
#endif
}

/* A commitment has K independent CRT projections.  Parallelising those
 * projections is deterministic and avoids the memory-contention regression
 * seen when whole commitments were scheduled concurrently.  Restrict its
 * use to large, non-nested calls; callers with aliased accumulation scratch
 * remain serial. */
static inline int revac_labrador_commit_parallel_enabled(void) {
#ifdef _OPENMP
  const char *value = getenv("REVAC_LABRADOR_COMMIT_PARALLEL");
  return value != NULL && value[0] != '\0' && value[0] != '0';
#else
  return 0;
#endif
}

/* The two LNP section pairs call polynomial kernels which share native
 * implementation state.  Running the pairs concurrently has produced rare
 * proofs that fail a later tail-reduction round, especially after a long
 * rejection path.  Keep this diagnostic switch opt-in until every callee is
 * proven reentrant; correctness must not depend on scheduler interleavings. */
static inline int revac_labrador_lnp_sections_parallel_enabled(void) {
#ifdef _OPENMP
  const char *value = getenv("REVAC_LABRADOR_LNP_SECTIONS_PARALLEL");
  return value != NULL && value[0] != '\0' && value[0] != '0';
#else
  return 0;
#endif
}

/* Coarse Labradoodle tasks operate on disjoint output objects and remain the
 * main safe source of fold parallelism.  Keep an opt-out for deterministic
 * diagnostics and for hosts where team creation outweighs the work. */
static inline int revac_labrador_ldd_parallel_enabled(void) {
#ifdef _OPENMP
  const char *value = getenv("REVAC_LABRADOR_LDD_PARALLEL");
  return value == NULL || value[0] == '\0' || value[0] != '0';
#else
  return 0;
#endif
}

/* A kernel-specific value overrides the coarse LDD switch.  This keeps
 * correctness isolation cheap: one suspect fold can be serialized while the
 * independent projection and aggregation stages remain parallel. */
static inline int revac_labrador_ldd_kernel_parallel_enabled(
    const char *kernel_variable) {
#ifdef _OPENMP
  const char *value = getenv(kernel_variable);
  if (value != NULL && value[0] != '\0')
    return value[0] != '0';
  return revac_labrador_ldd_parallel_enabled();
#else
  (void)kernel_variable;
  return 0;
#endif
}

#endif
