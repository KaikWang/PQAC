#ifndef ACC96_DELTA4_CODEC_H
#define ACC96_DELTA4_CODEC_H

#include <stddef.h>
#include <stdint.h>

#include "acc96_relation.h"

/* Lossless q88 delta-prefix transport format.  Historical acc96 names are retained
 * for API compatibility.  Coefficients are centered, ZigZag mapped, and
 * bit-packed in independent degree-256 blocks. */
#define ACC96_DELTA4_CODEC_MAGIC "RV88D4B1"
#define ACC96_DELTA4_CODEC_MAGIC_BYTES 8U
#define ACC96_DELTA4_CODEC_VERSION 1U
#define ACC96_DELTA4_CODEC_BLOCK_COEFFS 256U
#define ACC96_DELTA4_CODEC_HEADER_BYTES 32U

/* Encode the first layers of DELTA4.  A receiver whose first differing handle
 * bit is j only needs layers=j+1.  The returned allocation belongs to the
 * caller. */
int acc96_delta4_pack_alloc(uint8_t **out, size_t *out_len,
                            const acc96_delta4_t *delta, size_t layers,
                            const acc96_ctx_t *ctx);

/* Decode atomically: on failure, delta and layers_out are unchanged.  Missing
 * suffix layers are set to zero on success. */
int acc96_delta4_unpack(acc96_delta4_t *delta, size_t *layers_out,
                        const uint8_t *encoded, size_t encoded_len,
                        const acc96_ctx_t *ctx);

/* Encode/decode one independently fetchable layer of the six-layer prefix.
 * The codec header byte 28 carries the zero-based layer index.  A layer
 * packet is authenticated by the signed six-layer manifest at the protocol
 * level; these helpers only provide lossless canonical packing. */
int acc96_delta6_pack_layer_alloc(uint8_t **out, size_t *out_len,
                                  const acc96_delta6_t *delta, size_t layer,
                                  const acc96_ctx_t *ctx);
int acc96_delta6_unpack_layer(acc96_delta6_t *delta, size_t *layer_out,
                              const uint8_t *encoded, size_t encoded_len,
                              const acc96_ctx_t *ctx);

#endif
