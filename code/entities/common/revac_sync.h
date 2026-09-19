#ifndef REVAC_SYNC_H
#define REVAC_SYNC_H

#include <stddef.h>
#include <stdint.h>

#include "revac_context.h"

#define REVAC_SYNC_WIRE_VERSION REVAC_PARAM_SYNC_WIRE_VERSION
#define REVAC_SYNC_WITNESS_VEC_COUNT                                  \
  ((size_t)2 + ACC_GENPP_TRUNC_K_L +                                  \
   (size_t)ACC_GENPP_HANDLE_BITS * ACC_GENPP_GADGET_DIM +              \
   (ACC_GENPP_RANDOMIZED_WITNESS_REFRESH ? ACC_GENPP_GADGET_DIM : 0))
#define REVAC_SYNC_REFRESH_BITS                                       \
  (ACC_GENPP_RANDOMIZED_WITNESS_REFRESH ? ACC_GENPP_REFRESH_BITS_MAX : 0)
#define REVAC_SYNC_WITNESS_BYTES                                      \
  ((REVAC_SYNC_WITNESS_VEC_COUNT * PARAM_D * PARAM_N *                \
        PARAM_Q_BITLEN +                                              \
    PARAM_N + REVAC_SYNC_REFRESH_BITS + 7) /                          \
   8)
#define REVAC_SYNC_WIRE_HEADER_BYTES 156

/* q96 recovery is deliberately a rare TA-to-member fallback rather than a
   public broadcast.  A freshly sampled short coefficient is sent once as a
   canonical signed int32 and reconstructed into all three CRT lanes; the
   normal signed q96 BML update remains small. */
#define REVAC_Q96_SYNC_WITNESS_POLYS                                  \
  ((size_t)2 + ACC96_TRUNC_RIGHT + 1 +                               \
   (size_t)ACC96_HANDLE_BITS * ACC96_GADGET_DIGITS)
#define REVAC_Q96_SYNC_WITNESS_BYTES                                  \
  (REVAC_Q96_SYNC_WITNESS_POLYS * ACC96_N * 4U)
#define REVAC_Q96_SYNC_WIRE_HEADER_BYTES 156

int revac_sync_request_create(revac_sync_request_t *request,
                              uint32_t handle,
                              uint32_t last_seq,
                              const revac_public_state_t *public_state);

int revac_sync_response_pack(uint8_t **wire,
                             size_t *wire_len,
                             const revac_sync_request_t *request,
                             uint32_t acc_seq,
                             const poly_q_vec_d acc,
                             const acc_genpp_wit_t *witness,
                             const sep_sk_t *ta_sk,
                             const sep_pk_t *ta_pk);

int revac_sync_response_unpack(acc_genpp_wit_t *witness,
                               uint32_t *response_seq,
                               const uint8_t *wire,
                               size_t wire_len,
                               const revac_sync_request_t *request,
                               const revac_public_state_t *public_state);

int revac_q96_sync_response_pack(
    uint8_t **wire, size_t *wire_len,
    const revac_sync_request_t *request, uint32_t acc_seq,
    const acc96_poly_t *accumulator,
    const acc96_relation_witness_t *witness,
    const acc96_ctx_t *ctx,
    const sep_sk_t *ta_sk, const sep_pk_t *ta_pk);

int revac_q96_sync_response_unpack(
    acc96_relation_witness_t *witness, uint32_t *response_seq,
    const uint8_t *wire, size_t wire_len,
    const revac_sync_request_t *request,
    const revac_public_state_t *public_state);

void revac_sync_witness_set(acc_genpp_wit_t *dst,
                            const acc_genpp_wit_t *src);

#endif
