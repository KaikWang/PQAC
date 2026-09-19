#ifndef REVAC_LABRADOR_BRIDGE_H
#define REVAC_LABRADOR_BRIDGE_H

#include "revocable_ac.h"

typedef struct revac_labrador_prove_async revac_labrador_prove_async_t;

int revac_labrador_online_enabled(void);

int revac_labrador_memver_prove(revac_show_proof_t *proof,
                                const revac_public_state_t *public_state,
                                revac_uav_t *uav,
                                const poly_q_vec_d link_randomness[2]);

/* Prepare the exact lift synchronously, then submit only the immutable
   binary request to the external worker.  A failed finish is safe to retry
   through revac_labrador_memver_prove (e.g. worker cold-cache miss). */
int revac_labrador_memver_prove_async_start(
    revac_labrador_prove_async_t **job,
    const revac_show_proof_t *proof,
    const revac_public_state_t *public_state,
    revac_uav_t *uav,
    const poly_q_vec_d link_randomness[2]);

int revac_labrador_memver_prove_async_finish(
    revac_labrador_prove_async_t *job,
    revac_show_proof_t *proof);

int revac_labrador_memver_verify(
                                 const revac_public_state_t *public_state,
                                 const revac_show_context_t *ctx,
                                 const revac_show_proof_t *proof);

int revac_labrador_q96_memver_prove(
    revac_show_proof_t *proof, const revac_public_state_t *public_state,
    const revac_uav_t *uav, const poly_q_vec_d link_randomness[2]);

int revac_labrador_q96_memver_verify(
    const revac_public_state_t *public_state,
    const revac_show_context_t *ctx,
    const revac_show_proof_t *proof);

/* Move immutable q96 statement/parameter construction out of the first
   latency-sensitive show.  These calls never create or accept a proof. */
int revac_labrador_q96_prover_prewarm(
    const revac_public_state_t *public_state, const revac_uav_t *uav);
int revac_labrador_q96_verifier_prewarm(
    const revac_public_state_t *public_state);

#endif
