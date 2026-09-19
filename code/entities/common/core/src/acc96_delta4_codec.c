#include "acc96_delta4_codec.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

#define ACC96_DELTA4_CODEC_BLOCKS_PER_POLY \
  (ACC96_N / ACC96_DELTA4_CODEC_BLOCK_COEFFS)

static void store_u32(uint8_t out[4], uint32_t value) {
  out[0] = (uint8_t)value;
  out[1] = (uint8_t)(value >> 8);
  out[2] = (uint8_t)(value >> 16);
  out[3] = (uint8_t)(value >> 24);
}

static uint32_t load_u32(const uint8_t in[4]) {
  return (uint32_t)in[0] | ((uint32_t)in[1] << 8) |
         ((uint32_t)in[2] << 16) | ((uint32_t)in[3] << 24);
}

static acc96_u128 zigzag(acc96_i128 value) {
  if (value < 0) return (acc96_u128)(-value) * 2U - 1U;
  return (acc96_u128)value * 2U;
}

static unsigned bit_width(acc96_u128 value) {
  unsigned width = 0;
  while (value != 0) {
    width++;
    value >>= 1;
  }
  return width;
}

static void write_bits(uint8_t *out, size_t *bit_position,
                       acc96_u128 value, unsigned width) {
  while (width != 0) {
    const unsigned offset = (unsigned)(*bit_position & 7U);
    const unsigned available = 8U - offset;
    const unsigned take = width < available ? width : available;
    const unsigned mask = (1U << take) - 1U;
    out[*bit_position >> 3] |=
        (uint8_t)(((unsigned)value & mask) << offset);
    value >>= take;
    *bit_position += take;
    width -= take;
  }
}

static acc96_u128 read_bits(const uint8_t *in, size_t *bit_position,
                            unsigned width) {
  acc96_u128 value = 0;
  unsigned written = 0;
  while (width != 0) {
    const unsigned offset = (unsigned)(*bit_position & 7U);
    const unsigned available = 8U - offset;
    const unsigned take = width < available ? width : available;
    const unsigned mask = (1U << take) - 1U;
    const unsigned chunk = (in[*bit_position >> 3] >> offset) & mask;
    value |= (acc96_u128)chunk << written;
    *bit_position += take;
    written += take;
    width -= take;
  }
  return value;
}

static int coefficient_from_zigzag(acc96_coeff_t *out, acc96_u128 encoded) {
  const acc96_u128 modulus = acc96_modulus();
  const int negative = (int)(encoded & 1U);
  const acc96_u128 magnitude = negative ? (encoded + 1U) / 2U : encoded / 2U;
  if (magnitude > modulus / 2U) return 0;
  acc96_coeff_from_u128(out,
      negative && magnitude != 0 ? modulus - magnitude : magnitude);
  return 1;
}

int acc96_delta4_pack_alloc(uint8_t **out, size_t *out_len,
                            const acc96_delta4_t *delta, size_t layers,
                            const acc96_ctx_t *ctx) {
  const size_t blocks_per_layer =
      ACC96_GADGET_DIGITS * ACC96_DELTA4_CODEC_BLOCKS_PER_POLY;
  size_t descriptor_count;
  size_t payload_bits = 0;
  size_t payload_bytes;
  size_t total;
  uint8_t *encoded;
  uint8_t *widths;
  uint8_t *payload;
  size_t descriptor = 0;
  size_t bit_position = 0;

  if (out == NULL || out_len == NULL) return 0;
  *out = NULL;
  *out_len = 0;
  if (delta == NULL || ctx == NULL || !ctx->initialized || layers == 0 ||
      layers > ACC96_DELTA4_LAYERS ||
      ACC96_N % ACC96_DELTA4_CODEC_BLOCK_COEFFS != 0)
    return 0;
  descriptor_count = layers * blocks_per_layer;
  for (size_t layer = 0; layer < layers; layer++) {
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++) {
      for (size_t block = 0; block < ACC96_DELTA4_CODEC_BLOCKS_PER_POLY;
           block++) {
        unsigned width = 0;
        const size_t first = block * ACC96_DELTA4_CODEC_BLOCK_COEFFS;
        for (size_t i = 0; i < ACC96_DELTA4_CODEC_BLOCK_COEFFS; i++) {
          const acc96_i128 centered = acc96_coeff_to_centered(
              &delta->right[layer][digit].coeff[first + i], ctx);
          const unsigned current = bit_width(zigzag(centered));
          if (current > width) width = current;
        }
        if (width > ACC96_MODULUS_BITS) return 0;
        if (payload_bits > SIZE_MAX -
                               (size_t)width *
                                   ACC96_DELTA4_CODEC_BLOCK_COEFFS)
          return 0;
        payload_bits += (size_t)width * ACC96_DELTA4_CODEC_BLOCK_COEFFS;
      }
    }
  }
  payload_bytes = (payload_bits + 7U) / 8U;
  if (descriptor_count > UINT32_MAX || payload_bytes > UINT32_MAX ||
      ACC96_DELTA4_CODEC_HEADER_BYTES > SIZE_MAX - descriptor_count ||
      ACC96_DELTA4_CODEC_HEADER_BYTES + descriptor_count >
          SIZE_MAX - payload_bytes)
    return 0;
  total = ACC96_DELTA4_CODEC_HEADER_BYTES + descriptor_count + payload_bytes;
  if (total > UINT32_MAX) return 0;
  encoded = calloc(total, 1);
  if (encoded == NULL) return 0;
  memcpy(encoded, ACC96_DELTA4_CODEC_MAGIC,
         ACC96_DELTA4_CODEC_MAGIC_BYTES);
  encoded[8] = ACC96_DELTA4_CODEC_VERSION;
  encoded[9] = (uint8_t)layers;
  encoded[10] = (uint8_t)ACC96_GADGET_DIGITS;
  encoded[11] = 8U; /* log2(256 coefficients per block). */
  store_u32(encoded + 12, ACC96_N);
  store_u32(encoded + 16, (uint32_t)descriptor_count);
  store_u32(encoded + 20, (uint32_t)payload_bytes);
  store_u32(encoded + 24, (uint32_t)total);
  widths = encoded + ACC96_DELTA4_CODEC_HEADER_BYTES;
  payload = widths + descriptor_count;
  for (size_t layer = 0; layer < layers; layer++) {
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++) {
      for (size_t block = 0; block < ACC96_DELTA4_CODEC_BLOCKS_PER_POLY;
           block++) {
        unsigned width = 0;
        const size_t first = block * ACC96_DELTA4_CODEC_BLOCK_COEFFS;
        for (size_t i = 0; i < ACC96_DELTA4_CODEC_BLOCK_COEFFS; i++) {
          const acc96_i128 centered = acc96_coeff_to_centered(
              &delta->right[layer][digit].coeff[first + i], ctx);
          const unsigned current = bit_width(zigzag(centered));
          if (current > width) width = current;
        }
        widths[descriptor++] = (uint8_t)width;
        for (size_t i = 0; i < ACC96_DELTA4_CODEC_BLOCK_COEFFS; i++) {
          const acc96_i128 centered = acc96_coeff_to_centered(
              &delta->right[layer][digit].coeff[first + i], ctx);
          write_bits(payload, &bit_position, zigzag(centered), width);
        }
      }
    }
  }
  if (descriptor != descriptor_count || bit_position != payload_bits) {
    free(encoded);
    return 0;
  }
  *out = encoded;
  *out_len = total;
  return 1;
}

int acc96_delta4_unpack(acc96_delta4_t *delta, size_t *layers_out,
                        const uint8_t *encoded, size_t encoded_len,
                        const acc96_ctx_t *ctx) {
  const size_t blocks_per_layer =
      ACC96_GADGET_DIGITS * ACC96_DELTA4_CODEC_BLOCKS_PER_POLY;
  acc96_delta4_t *decoded = NULL;
  size_t layers;
  size_t descriptor_count;
  size_t payload_bytes;
  size_t payload_bits = 0;
  const uint8_t *widths;
  const uint8_t *payload;
  size_t descriptor = 0;
  size_t bit_position = 0;
  int ok = 0;

  if (delta == NULL || encoded == NULL || ctx == NULL || !ctx->initialized ||
      encoded_len < ACC96_DELTA4_CODEC_HEADER_BYTES ||
      memcmp(encoded, ACC96_DELTA4_CODEC_MAGIC,
             ACC96_DELTA4_CODEC_MAGIC_BYTES) != 0 ||
      encoded[8] != ACC96_DELTA4_CODEC_VERSION ||
      encoded[10] != ACC96_GADGET_DIGITS || encoded[11] != 8U ||
      load_u32(encoded + 12) != ACC96_N || load_u32(encoded + 24) != encoded_len ||
      load_u32(encoded + 28) != 0)
    return 0;
  layers = encoded[9];
  descriptor_count = load_u32(encoded + 16);
  payload_bytes = load_u32(encoded + 20);
  if (layers == 0 || layers > ACC96_DELTA4_LAYERS ||
      descriptor_count != layers * blocks_per_layer ||
      ACC96_DELTA4_CODEC_HEADER_BYTES + descriptor_count > encoded_len ||
      payload_bytes != encoded_len - ACC96_DELTA4_CODEC_HEADER_BYTES -
                           descriptor_count)
    return 0;
  widths = encoded + ACC96_DELTA4_CODEC_HEADER_BYTES;
  payload = widths + descriptor_count;
  for (size_t i = 0; i < descriptor_count; i++) {
    if (widths[i] > ACC96_MODULUS_BITS ||
        payload_bits > SIZE_MAX -
                           (size_t)widths[i] *
                               ACC96_DELTA4_CODEC_BLOCK_COEFFS)
      return 0;
    payload_bits +=
        (size_t)widths[i] * ACC96_DELTA4_CODEC_BLOCK_COEFFS;
  }
  if ((payload_bits + 7U) / 8U != payload_bytes) return 0;
  if ((payload_bits & 7U) != 0 && payload_bytes != 0) {
    const unsigned used = (unsigned)(payload_bits & 7U);
    if ((payload[payload_bytes - 1] & (uint8_t)(0xffU << used)) != 0)
      return 0;
  }
  decoded = calloc(1, sizeof(*decoded));
  if (decoded == NULL) return 0;
  for (size_t layer = 0; layer < layers; layer++) {
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++) {
      for (size_t block = 0; block < ACC96_DELTA4_CODEC_BLOCKS_PER_POLY;
           block++) {
        const unsigned width = widths[descriptor++];
        const size_t first = block * ACC96_DELTA4_CODEC_BLOCK_COEFFS;
        for (size_t i = 0; i < ACC96_DELTA4_CODEC_BLOCK_COEFFS; i++) {
          const acc96_u128 value = read_bits(payload, &bit_position, width);
          if (!coefficient_from_zigzag(
                  &decoded->right[layer][digit].coeff[first + i], value))
            goto cleanup;
        }
      }
    }
  }
  if (descriptor != descriptor_count || bit_position != payload_bits)
    goto cleanup;
  memcpy(delta, decoded, sizeof(*delta));
  if (layers_out != NULL) *layers_out = layers;
  ok = 1;
cleanup:
  free(decoded);
  return ok;
}

int acc96_delta6_pack_layer_alloc(uint8_t **out, size_t *out_len,
                                  const acc96_delta6_t *delta, size_t layer,
                                  const acc96_ctx_t *ctx) {
  acc96_delta4_t one;
  if (out == NULL || out_len == NULL || delta == NULL ||
      layer >= ACC96_DELTA6_LAYERS) {
    return 0;
  }
  memset(&one, 0, sizeof(one));
  memcpy(one.right[0], delta->right[layer], sizeof(one.right[0]));
  if (!acc96_delta4_pack_alloc(out, out_len, &one, 1, ctx)) return 0;
  /* Byte 28 was reserved in v1; using it for the layer index preserves the
     old four-layer prefix format while making each new packet self-describing. */
  (*out)[28] = (uint8_t)layer;
  return 1;
}

int acc96_delta6_unpack_layer(acc96_delta6_t *delta, size_t *layer_out,
                              const uint8_t *encoded, size_t encoded_len,
                              const acc96_ctx_t *ctx) {
  uint8_t *copy;
  acc96_delta4_t one;
  size_t decoded_layers = 0;
  size_t layer;
  if (delta == NULL || encoded == NULL || ctx == NULL || encoded_len < 32 ||
      encoded[28] >= ACC96_DELTA6_LAYERS) return 0;
  layer = encoded[28];
  copy = malloc(encoded_len);
  if (copy == NULL) return 0;
  memcpy(copy, encoded, encoded_len);
  copy[28] = 0;
  memset(&one, 0, sizeof(one));
  if (!acc96_delta4_unpack(&one, &decoded_layers, copy, encoded_len, ctx) ||
      decoded_layers != 1) {
    free(copy);
    return 0;
  }
  memcpy(delta->right[layer], one.right[0], sizeof(one.right[0]));
  if (layer_out != NULL) *layer_out = layer;
  free(copy);
  return 1;
}
