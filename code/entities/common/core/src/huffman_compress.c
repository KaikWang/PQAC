/*************************************************
 * Huffman compression / decompression for discrete
 * Gaussian samples (z1, z2, z3) in the Show proof.
 *
 * Encoding per coefficient:
 *   [1 bit: sign (1=negative)]
 *   [Huffman code for bucket k of |v|]
 *   [(max(0, k-1)) bits: offset within bucket]
 *
 * Auto-generated tables in huffman_compress.h
 *************************************************/

#include "huffman_compress.h"
#include <string.h>

/* ================================================
 * Bit-level write buffer
 * ================================================ */

typedef struct {
    uint8_t *buf;
    size_t   byte_pos;
    int      bit_off;   /* 0..7, next write position (0=MSB of current byte) */
} bw_t;

static void bw_init(bw_t *bw, uint8_t *buf) {
    bw->buf  = buf;
    bw->byte_pos = 0;
    bw->bit_off  = 0;
    bw->buf[0]   = 0;
}

static void bw_write_bit(bw_t *bw, int bit) {
    if (bw->bit_off == 8) { bw->bit_off = 0; bw->byte_pos++; bw->buf[bw->byte_pos] = 0; }
    bw->buf[bw->byte_pos] |= (uint8_t)((bit & 1) << (7 - bw->bit_off));
    bw->bit_off++;
}

/* write nbits from value, MSB first */
static void bw_write_bits(bw_t *bw, uint64_t value, int nbits) {
    for (int i = nbits - 1; i >= 0; i--)
        bw_write_bit(bw, (int)((value >> i) & 1));
}

/* returns total bytes written */
static size_t bw_finish(bw_t *bw) {
    return bw->bit_off > 0 ? bw->byte_pos + 1 : bw->byte_pos;
}

/* ================================================
 * Bit-level read buffer
 * ================================================ */

typedef struct {
    const uint8_t *buf;
    size_t   byte_pos;
    size_t   byte_limit;
    int      bit_off;
    int      error;
} br_t;

static void br_init(br_t *br, const uint8_t *buf, size_t len) {
    br->buf   = buf;
    br->byte_pos   = 0;
    br->byte_limit = len;
    br->bit_off    = 0;
    br->error      = 0;
}

static int br_read_bit(br_t *br) {
    if (br->byte_pos >= br->byte_limit) {
        br->error = 1;
        return 0;
    }
    int b = (br->buf[br->byte_pos] >> (7 - br->bit_off)) & 1;
    if (++br->bit_off == 8) { br->bit_off = 0; br->byte_pos++; }
    return b;
}

/* read nbits, MSB first */
static uint64_t br_read_bits(br_t *br, int nbits) {
    uint64_t v = 0;
    for (int i = 0; i < nbits; i++)
        v = (v << 1) | (uint64_t)br_read_bit(br);
    return v;
}

/* ================================================
 * Bucket mapping
 * ================================================ */

static inline int coeff_to_bucket(uint64_t abs_v) {
    if (abs_v == 0) return 0;
#if defined(__GNUC__) || defined(__clang__)
    return 64 - __builtin_clzll(abs_v);
#else
    int k = 1;
    while (abs_v >>= 1) k++;
    return k;
#endif
}

static inline uint64_t bucket_offset(uint64_t abs_v, int k) {
    if (k <= 1) return 0;
    return abs_v - (1ULL << (k - 1));
}

static inline int offset_bitcount(int k) {
    return (k <= 1) ? 0 : (k - 1);
}

/* ================================================
 * Pack (compress)
 * ================================================ */

static size_t huffman_pack(
    uint8_t *out,
    const int64_t *coeffs, size_t count,
    const uint32_t *codes, const uint8_t *lens,
    int max_bucket, int esc_bucket)
{
    bw_t bw;
    bw_init(&bw, out);
    for (size_t i = 0; i < count; i++) {
        int64_t v  = coeffs[i];
        int sign   = (v < 0);
        /* careful: -INT64_MIN overflows; clamp to INT64_MAX in that case */
        uint64_t av = (v == INT64_MIN) ? ((uint64_t)INT64_MAX + 1) :
                      (v < 0)          ? (uint64_t)(-v) : (uint64_t)v;

        int k = coeff_to_bucket(av);
        if (k > max_bucket) k = esc_bucket;
        if (lens[k] == 0) return 0;

        bw_write_bit(&bw, sign);               /* sign */
        bw_write_bits(&bw, codes[k], lens[k]); /* Huffman bucket code */

        if (k <= max_bucket) {
            int ob = offset_bitcount(k);
            if (ob > 0) bw_write_bits(&bw, bucket_offset(av, k), ob);
        } else {
            /* escape: 48 raw bits of |v| (covers all tail possibilities) */
            bw_write_bits(&bw, av >> 32, 16);
            bw_write_bits(&bw, av & 0xFFFFFFFFULL, 32);
        }
    }
    return bw_finish(&bw);
}

/* ================================================
 * Unpack (decompress)
 * ================================================ */

static size_t huffman_unpack(
    int64_t *out,
    const uint8_t *in, size_t packed_bytes, size_t count,
    const uint32_t *codes, const uint8_t *lens,
    int max_bucket, int esc_bucket)
{
    br_t br;
    br_init(&br, in, packed_bytes);

    /* find max code length for decode loop */
    int max_len = 0;
    for (int i = 0; i <= esc_bucket; i++)
        if (lens[i] > max_len) max_len = lens[i];

    for (size_t i = 0; i < count; i++) {
        int sign = br_read_bit(&br);
        if (br.error) return 0;

        /* decode Huffman: read bit-by-bit until we match a code */
        int k = -1;
        uint64_t cand = 0;
        for (int b = 0; b < max_len && k < 0; b++) {
            cand = (cand << 1) | (uint64_t)br_read_bit(&br);
            /* compare with all codes of length (b+1) */
            int tlen = b + 1;
            for (int j = 0; j <= esc_bucket; j++) {
                if (lens[j] == tlen && codes[j] == cand) {
                    k = j;
                    break;
                }
            }
        }
        if (k < 0 || (k > max_bucket && k != esc_bucket) || br.error)
            return 0;

        /* reconstruct absolute value */
        uint64_t av;
        if (k == esc_bucket) {
            /* escape: read 48 raw bits */
            av = br_read_bits(&br, 48);
        } else if (k <= 1) {
            av = (uint64_t)k;
        } else {
            int ob = offset_bitcount(k);
            av = (1ULL << (k - 1)) + br_read_bits(&br, ob);
        }
        if (br.error || av > (uint64_t)INT64_MAX || (sign && av == 0))
            return 0;

        /* reconstruct signed value */
        int64_t sv = (int64_t)av;
        if (sign && av > 0) sv = -sv;
        /* Note: if av == 0 and sign == 1, we still output 0 (av=0 means v=0) */
        out[i] = sv;
    }
    {
        size_t used_bits = br.byte_pos * 8 + (size_t)br.bit_off;
        size_t used_bytes = (used_bits + 7) / 8;
        if (br.error || used_bytes != packed_bytes) return 0;
        if (used_bits % 8 != 0 &&
            (in[packed_bytes - 1] & ((1U << (8 - used_bits % 8)) - 1U)) != 0)
            return 0;
    }
    return count;
}

/* ================================================
 * Per-sigma public API
 * ================================================ */

size_t huffman_pack_z1(uint8_t *out, const int64_t *coeffs, size_t count) {
    return huffman_pack(out, coeffs, count,
        huffman_z1_code, huffman_z1_len,
        HUFFMAN_Z1_MAX_BUCKET, HUFFMAN_Z1_ESCAPE_BUCKET);
}
size_t huffman_pack_z2(uint8_t *out, const int64_t *coeffs, size_t count) {
    return huffman_pack(out, coeffs, count,
        huffman_z2_code, huffman_z2_len,
        HUFFMAN_Z2_MAX_BUCKET, HUFFMAN_Z2_ESCAPE_BUCKET);
}
size_t huffman_pack_z3(uint8_t *out, const int64_t *coeffs, size_t count) {
    return huffman_pack(out, coeffs, count,
        huffman_z3_code, huffman_z3_len,
        HUFFMAN_Z3_MAX_BUCKET, HUFFMAN_Z3_ESCAPE_BUCKET);
}

size_t huffman_unpack_z1(int64_t *out, const uint8_t *in, size_t packed_bytes, size_t count) {
    return huffman_unpack(out, in, packed_bytes, count,
        huffman_z1_code, huffman_z1_len,
        HUFFMAN_Z1_MAX_BUCKET, HUFFMAN_Z1_ESCAPE_BUCKET);
}
size_t huffman_unpack_z2(int64_t *out, const uint8_t *in, size_t packed_bytes, size_t count) {
    return huffman_unpack(out, in, packed_bytes, count,
        huffman_z2_code, huffman_z2_len,
        HUFFMAN_Z2_MAX_BUCKET, HUFFMAN_Z2_ESCAPE_BUCKET);
}
size_t huffman_unpack_z3(int64_t *out, const uint8_t *in, size_t packed_bytes, size_t count) {
    return huffman_unpack(out, in, packed_bytes, count,
        huffman_z3_code, huffman_z3_len,
        HUFFMAN_Z3_MAX_BUCKET, HUFFMAN_Z3_ESCAPE_BUCKET);
}

/* ================================================
 * Compact bit-packing for uniform-in-R_q elements
 *
 * Each coefficient is in [0, q_show-1] or centered
 * in [-q_show/2, q_show/2].  We store the centered
 * value as a signed integer using 58 bits.
 *
 * Encoding: 58 bits/coeff, MSB first, contiguous.
 * ================================================ */

size_t compact_pack_uniform(uint8_t *out, const int64_t *coeffs, size_t count) {
    bw_t bw;
    bw_init(&bw, out);
    for (size_t i = 0; i < count; i++) {
        /* coeffs[i] is from poly_qshow_get_coeff(): non-negative, in [0, q-1].
         * q < 2^58, so 58 bits are lossless. */
        uint64_t v = (uint64_t)coeffs[i] & 0x3FFFFFFFFFFFFFFULL; /* 58-bit mask */
        bw_write_bits(&bw, v, COMPACT_UNIFORM_BITS);
    }
    return bw_finish(&bw);
}

size_t compact_unpack_uniform(int64_t *out, const uint8_t *in, size_t packed_bytes, size_t count) {
    br_t br;
    br_init(&br, in, packed_bytes);
    for (size_t i = 0; i < count; i++) {
        uint64_t v = br_read_bits(&br, COMPACT_UNIFORM_BITS);
        /* Value is non-negative [0, q-1], no sign extension needed */
        out[i] = (int64_t)v;
    }
    return count;
}

/* ================================================
 * Compact bit-packing for challenge polynomial c
 *
 * Coefficients are in [-rho, rho] with rho=8.
 * Encoded as (value + 8) in [0,16] using 5 bits.
 * ================================================ */

size_t compact_pack_challenge(uint8_t *out, const int64_t *coeffs, size_t count) {
    bw_t bw;
    bw_init(&bw, out);
    for (size_t i = 0; i < count; i++) {
        uint64_t v = (uint64_t)(coeffs[i] + COMPACT_CHALLENGE_OFFSET);
        bw_write_bits(&bw, v, COMPACT_CHALLENGE_BITS);
    }
    return bw_finish(&bw);
}

size_t compact_unpack_challenge(int64_t *out, const uint8_t *in, size_t packed_bytes, size_t count) {
    br_t br;
    br_init(&br, in, packed_bytes);
    for (size_t i = 0; i < count; i++) {
        uint64_t v = br_read_bits(&br, COMPACT_CHALLENGE_BITS);
        out[i] = (int64_t)v - COMPACT_CHALLENGE_OFFSET;
    }
    return count;
}
