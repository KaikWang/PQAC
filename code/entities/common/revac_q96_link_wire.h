#ifndef REVAC_Q96_LINK_WIRE_H
#define REVAC_Q96_LINK_WIRE_H

#include <stddef.h>
#include <stdint.h>

#include "revocable_ac.h"

#define REVAC_Q96_LINK_MAGIC "RV88LK02"
#define REVAC_Q96_INSTANCE_MAGIC "RV88IN02"
#define REVAC_Q96_LINK_MAGIC_BYTES 8U
#define REVAC_Q96_LINK_DECOMP_BASE 256
#define REVAC_Q96_LINK_DECOMP_LIMBS 4U
#define REVAC_Q96_LINK_POLY_BYTES (PARAM_N * 4U)
#define REVAC_Q96_LINK_WIRE_BYTES                                      \
  (REVAC_Q96_LINK_MAGIC_BYTES + 4U + 3U * REVAC_CONTEXT_HASH_BYTES +   \
   (PARAM_D + 2U * PARAM_D +                                           \
    PARAM_D * REVAC_Q96_LINK_DECOMP_LIMBS) *                           \
       REVAC_Q96_LINK_POLY_BYTES)
#define REVAC_Q96_INSTANCE_WIRE_BYTES                                  \
  (REVAC_Q96_LINK_MAGIC_BYTES + 4U + 4U + ACC96_POLY_PACKED_BYTES +    \
   3U * REVAC_CONTEXT_HASH_BYTES +                                     \
   PARAM_D * REVAC_Q96_LINK_POLY_BYTES)

/* Local-only input for the q88 LaBRADOR prover.  It carries the public
 * q19 handle commitment followed by its private opening and exact-lift
 * decomposition.  It must never be sent to TA/GCS or placed in a show
 * packet. */
int revac_q96_link_wire_encode(
    uint8_t *out, size_t out_len,
    const revac_public_state_t *public_state, uint32_t handle,
    const poly_q_vec_d randomness[2], const poly_q_vec_d commitment);

/* Public, per-show verifier input.  The commitment is shared with the
 * signature proof; no q88 witness, link opening, or exact lift is included. */
int revac_q96_instance_wire_encode(
    uint8_t *out, size_t out_len, const revac_public_state_t *public_state,
    const poly_q_vec_d commitment);

#endif
