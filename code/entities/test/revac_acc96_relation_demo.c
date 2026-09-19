#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <zlib.h>

#include "acc96_delta4_codec.h"
#include "acc96_relation.h"

static double seconds(void) {
  struct timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);
  return now.tv_sec + now.tv_nsec / 1e9;
}

static size_t benchmark_repetitions(void) {
  const char *value = getenv("REVAC_Q88_DELTA6_BENCH_REPETITIONS");
  char *end = NULL;
  unsigned long parsed;
  if (value == NULL || *value == '\0') return 0;
  parsed = strtoul(value, &end, 10);
  if (end == value || *end != '\0' || parsed == 0 || parsed > 50) return 0;
  return (size_t)parsed;
}

static int benchmark_delta6_updates(const acc96_relation_pp_t *pp,
                                    const acc96_poly_t *base_accumulator,
                                    const acc96_relation_witness_t *base,
                                    uint32_t member,
                                    const acc96_ctx_t *ctx) {
  const size_t repetitions = benchmark_repetitions();
  if (repetitions == 0) return 1;
  printf("q88 DELTA6 update benchmark (%zu repetitions; raw seconds):\n",
         repetitions);
  for (size_t bit = 0; bit <= ACC96_DELTA6_LAYERS; bit++) {
    acc96_delta6_t *delta = calloc(1, sizeof(*delta));
    acc96_poly_t *post_accumulator = malloc(sizeof(*post_accumulator));
    acc96_relation_update_t canonical_update;
    if (delta == NULL || post_accumulator == NULL) {
      free(delta);
      free(post_accumulator);
      return 0;
    }
    acc96_poly_set(post_accumulator, base_accumulator);
    if (!acc96_relation_delete_delta6(
            post_accumulator, &canonical_update, delta, pp,
            member ^ (UINT32_C(1) << bit), 1, ctx)) {
      free(delta);
      free(post_accumulator);
      return 0;
    }
    printf("  first-difference %zu (%s):", bit,
           bit < ACC96_DELTA6_LAYERS ? "DELTA6" : "LOCAL");
    for (size_t run = 0; run < repetitions; run++) {
      acc96_relation_witness_t *candidate = malloc(sizeof(*candidate));
      acc96_relation_witness_t *expected = malloc(sizeof(*expected));
      size_t applied = 0;
      double started;
      double elapsed;
      int ok;
      if (candidate == NULL || expected == NULL) {
        free(candidate);
        free(expected);
        free(delta);
        free(post_accumulator);
        return 0;
      }
      memcpy(candidate, base, sizeof(*candidate));
      memcpy(expected, base, sizeof(*expected));
      started = seconds();
      if (bit < ACC96_DELTA6_LAYERS) {
        ok = acc96_relation_memwitup_delta6(
            candidate, member, &canonical_update, delta, bit + 1, &applied);
      } else {
        ok = acc96_relation_memwitup(candidate, pp, member,
                                     &canonical_update, ctx);
      }
      elapsed = seconds() - started;
      if (!ok || (bit < ACC96_DELTA6_LAYERS && applied != bit + 1) ||
          !acc96_relation_memwitup(expected, pp, member, &canonical_update,
                                   ctx) ||
          memcmp(candidate, expected, sizeof(*candidate)) != 0 ||
          !acc96_relation_memver(pp, post_accumulator, member, candidate, 0,
                                ctx)) {
        free(candidate);
        free(expected);
        free(delta);
        free(post_accumulator);
        return 0;
      }
      printf(" %.6f", elapsed);
      free(candidate);
      free(expected);
    }
    putchar('\n');
    free(delta);
    free(post_accumulator);
  }
  return 1;
}

static unsigned u128_bits(acc96_u128 value) {
  unsigned bits = 0;
  while (value != 0) {
    bits++;
    value >>= 1;
  }
  return bits == 0 ? 1 : bits;
}

static void print_delta4_layer_sizes(const acc96_delta4_t *delta,
                                     const acc96_ctx_t *ctx) {
  const size_t coefficients = ACC96_GADGET_DIGITS * ACC96_N;
  double covered_weighted_kib = 0;
  double cumulative_kib = 0;
  for (size_t layer = 0; layer < ACC96_DELTA4_LAYERS; layer++) {
    acc96_u128 maximum = 0;
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++) {
      for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++) {
        const acc96_i128 centered = acc96_coeff_to_centered(
            &delta->right[layer][digit].coeff[coefficient], ctx);
        const acc96_u128 magnitude = centered < 0
                                         ? (acc96_u128)(-centered)
                                         : (acc96_u128)centered;
        if (magnitude > maximum) maximum = magnitude;
      }
    }
    /* ZigZag maps -m/+m to at most 2m, so this is an exact fixed-width
       representation for every coefficient in the layer. */
    const unsigned bits = u128_bits(maximum * 2);
    const double layer_kib =
        ((double)coefficients * bits + 7.0) / (8.0 * 1024.0);
    cumulative_kib += layer_kib;
    covered_weighted_kib += cumulative_kib / (double)(UINT64_C(1) <<
                                                       (layer + 1));
    printf("q88 DELTA4 layer %zu: max |coefficient|=%llu, "
           "exact ZigZag width %u bit, fixed payload %.2f KiB, "
           "cumulative prefix %.2f KiB.\n",
           layer,
           (unsigned long long)maximum, bits, layer_kib, cumulative_kib);
  }
  printf("q88 receiver-specific DELTA prefix: expected payload %.2f KiB "
         "over all random survivors (6.25%% use LOCAL); %.2f KiB conditioned "
         "on DELTA4 coverage.\n",
         covered_weighted_kib, covered_weighted_kib / 0.9375);
}

static int test_delta4_codec(const acc96_delta4_t *delta,
                             const acc96_ctx_t *ctx) {
  acc96_delta4_t *decoded = calloc(1, sizeof(*decoded));
  double expected_all_kib = 0.0;
  double covered_probability = 0.0;
  if (decoded == NULL) return 0;
  for (size_t layers = 1; layers <= ACC96_DELTA4_LAYERS; layers++) {
    uint8_t *encoded = NULL;
    size_t encoded_len = 0;
    size_t decoded_layers = 0;
    const double probability = 1.0 / (double)(UINT64_C(1) << layers);
    if (!acc96_delta4_pack_alloc(&encoded, &encoded_len, delta, layers, ctx) ||
        !acc96_delta4_unpack(decoded, &decoded_layers, encoded, encoded_len,
                             ctx) ||
        decoded_layers != layers ||
        memcmp(decoded->right, delta->right,
               layers * sizeof(delta->right[0])) != 0) {
      free(encoded);
      free(decoded);
      return 0;
    }
    for (size_t layer = layers; layer < ACC96_DELTA4_LAYERS; layer++) {
      acc96_delta4_t zero = {0};
      if (memcmp(&decoded->right[layer], &zero.right[layer],
                 sizeof(decoded->right[layer])) != 0) {
        free(encoded);
        free(decoded);
        return 0;
      }
    }
    expected_all_kib += probability * encoded_len / 1024.0;
    covered_probability += probability;
    printf("q88 DELTA4 block codec prefix %zu: %.2f KiB.\n",
           layers, encoded_len / 1024.0);
    if (layers == ACC96_DELTA4_LAYERS) {
      uLongf compressed_len = compressBound((uLong)encoded_len);
      uint8_t *compressed = malloc((size_t)compressed_len);
      if (compressed != NULL &&
          compress2(compressed, &compressed_len, encoded,
                    (uLong)encoded_len, Z_BEST_COMPRESSION) == Z_OK) {
        printf("q88 DELTA4 optional DEFLATE probe: %.2f KiB (%.1f%% of "
               "block codec).\n",
               (double)compressed_len / 1024.0,
               100.0 * (double)compressed_len / (double)encoded_len);
      }
      free(compressed);
    }
    /* A width above the q88 modulus size must fail before mutating output. */
    if (layers == ACC96_DELTA4_LAYERS) {
      acc96_delta4_t saved;
      memcpy(&saved, decoded, sizeof(saved));
      encoded[ACC96_DELTA4_CODEC_HEADER_BYTES] =
          (uint8_t)(ACC96_MODULUS_BITS + 1U);
      if (acc96_delta4_unpack(decoded, NULL, encoded, encoded_len, ctx) ||
          memcmp(decoded, &saved, sizeof(saved)) != 0) {
        free(encoded);
        free(decoded);
        return 0;
      }
    }
    free(encoded);
  }
  printf("q88 DELTA4 block codec receiver-specific expectation: "
         "%.2f KiB over all survivors; %.2f KiB conditioned on 4-layer "
         "coverage.\n",
         expected_all_kib, expected_all_kib / covered_probability);
  free(decoded);
  return 1;
}

static long double poly_chunk_norm2(const acc96_poly_t *poly,
                                    size_t chunk,
                                    const acc96_ctx_t *ctx) {
  long double norm = 0;
  const size_t first = chunk * 256;
  for (size_t i = first; i < first + 256; i++) {
    const long double value =
        (long double)acc96_coeff_to_centered(&poly->coeff[i], ctx);
    norm += value * value;
  }
  return norm;
}

static long double max_backend_chunk_norm2(
    const acc96_relation_witness_t *witness, const acc96_ctx_t *ctx) {
  long double maximum = 0;
#define ACCOUNT_POLY(poly_) do {                                           \
    for (size_t chunk_ = 0; chunk_ < ACC96_N / 256; chunk_++) {            \
      const long double norm_ = poly_chunk_norm2((poly_), chunk_, ctx);    \
      if (norm_ > maximum) maximum = norm_;                                \
    }                                                                      \
  } while (0)
  ACCOUNT_POLY(&witness->left);
  ACCOUNT_POLY(&witness->v12);
  ACCOUNT_POLY(&witness->tag);
  for (size_t i = 0; i < ACC96_TRUNC_RIGHT; i++)
    ACCOUNT_POLY(&witness->v2[i]);
  for (size_t bit = 0; bit < ACC96_HANDLE_BITS; bit++)
    for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++)
      ACCOUNT_POLY(&witness->right.right[bit][digit]);
#undef ACCOUNT_POLY
  return maximum;
}

int main(void) {
  static const uint8_t pp_seed[ACC96_SEED_BYTES] = {
      0x51, 0x39, 0x36, 0x2d, 0x70, 0x70, 0x2d, 0x76,
      0x31, 0x00, 0x73, 0x68, 0x61, 0x6b, 0x65, 0x32,
      0x35, 0x36, 0x00, 0xa5, 0x5a, 0x18, 0x42, 0x99,
      0x01, 0x7f, 0x33, 0xc4, 0x0d, 0xee, 0x20, 0x26};
  static const uint8_t witness_seed[ACC96_SEED_BYTES] = {
      0x77, 0x69, 0x74, 0x6e, 0x65, 0x73, 0x73, 0x2d,
      0x6b, 0x61, 0x74, 0x2d, 0x71, 0x39, 0x36, 0x00,
      0x19, 0x88, 0x73, 0x44, 0xfa, 0x21, 0x38, 0x91,
      0xcc, 0x05, 0x62, 0x7a, 0x10, 0xbe, 0xef, 0x01};
  const uint32_t member = UINT32_C(0xa5c39e71);
  const uint32_t revoked[] = {
      UINT32_C(0xa5c39e70), UINT32_C(0xa5c39e73),
      UINT32_C(0xa5c39e75)};
  acc96_ctx_t ctx;
  acc96_relation_pp_t *pp = calloc(1, sizeof(*pp));
  acc96_relation_witness_t *witness = calloc(1, sizeof(*witness));
  acc96_poly_t *accumulator = calloc(1, sizeof(*accumulator));
  if (pp == NULL || witness == NULL || accumulator == NULL ||
      !acc96_ctx_init(&ctx) || !acc96_relation_pp_expand(pp, pp_seed, &ctx))
    return 1;

  acc96_relation_make_kat_instance(accumulator, witness, pp, member,
                                   witness_seed, &ctx);
  const long double initial_norm = acc96_relation_witness_norm2(witness, &ctx);
  if (!acc96_relation_memver(pp, accumulator, member, witness,
                            initial_norm + 1, &ctx)) {
    fprintf(stderr, "q88 initial MemVer relation failed.\n");
    return 1;
  }

  /* Compare exact single-checkpoint broadcasts.  The deepest covered member
     exercises the maximum receiver work for each T2/T4/T8 choice. */
  {
    static const size_t checkpoint_bits[] = {2, 4, 8};
    const uint32_t selector = UINT32_C(0x6d52a9c0);
    for (size_t c = 0; c < sizeof(checkpoint_bits) /
                                sizeof(checkpoint_bits[0]); c++) {
      const size_t checkpoint_bit = checkpoint_bits[c];
      const uint32_t receiver =
          selector ^ (UINT32_C(1) << (checkpoint_bit - 1));
      acc96_gmat_t *checkpoint = malloc(sizeof(*checkpoint));
      acc96_evalfx_delta_t *local_delta = malloc(sizeof(*local_delta));
      acc96_evalfx_delta_t *checkpoint_delta =
          malloc(sizeof(*checkpoint_delta));
      double build_started, build_elapsed, local_started, local_elapsed;
      double resume_started, resume_elapsed;
      if (checkpoint == NULL || local_delta == NULL ||
          checkpoint_delta == NULL)
        return 1;
      build_started = seconds();
      if (!acc96_evalf_tail_checkpoint(checkpoint, &pp->ba, selector,
                                       checkpoint_bit, &ctx))
        return 1;
      build_elapsed = seconds() - build_started;
      local_started = seconds();
      acc96_evalfx(local_delta, &pp->ba, selector, receiver, &pp->ua, &ctx);
      local_elapsed = seconds() - local_started;
      resume_started = seconds();
      if (!acc96_evalfx_from_tail(checkpoint_delta, &pp->ba, selector,
                                  receiver, &pp->ua, checkpoint_bit,
                                  checkpoint, &ctx))
        return 1;
      resume_elapsed = seconds() - resume_started;
      if (memcmp(local_delta, checkpoint_delta, sizeof(*local_delta)) != 0) {
        fprintf(stderr, "q88 T%zu checkpoint/local mismatch.\n",
                checkpoint_bit);
        return 1;
      }
      printf("q88 T%zu checkpoint equality passed: TA standalone build %.3f s, "
             "UAV LOCAL %.3f s, resume %.3f s, exact payload %.2f KiB, "
             "uniform-handle coverage %.3f%%.\n",
             checkpoint_bit, build_elapsed, local_elapsed, resume_elapsed,
             (ACC96_GADGET_DIGITS * ACC96_N * ACC96_MODULUS_BITS / 8.0) /
                 1024.0,
             100.0 * (1.0 - 1.0 / (double)(UINT64_C(1) << checkpoint_bit)));
      free(checkpoint_delta);
      free(local_delta);
      free(checkpoint);
    }
    {
      acc96_gmat_t *checkpoint = malloc(sizeof(*checkpoint));
      acc96_evalfx_delta_t *local_delta = malloc(sizeof(*local_delta));
      acc96_evalfx_delta_t *checkpoint_delta =
          malloc(sizeof(*checkpoint_delta));
      if (checkpoint == NULL || local_delta == NULL ||
          checkpoint_delta == NULL ||
          !acc96_evalf_tail_checkpoint(checkpoint, &pp->ba, selector,
                                       ACC96_TAIL4_BIT, &ctx))
        return 1;
      for (size_t bit = 0; bit < ACC96_TAIL4_BIT; bit++) {
        const uint32_t receiver = selector ^ (UINT32_C(1) << bit);
        acc96_evalfx(local_delta, &pp->ba, selector, receiver, &pp->ua, &ctx);
        if (!acc96_evalfx_from_tail(checkpoint_delta, &pp->ba, selector,
                                    receiver, &pp->ua, ACC96_TAIL4_BIT,
                                    checkpoint, &ctx) ||
            memcmp(local_delta, checkpoint_delta,
                   sizeof(*local_delta)) != 0) {
          fprintf(stderr, "q88 T4 first-difference bit %zu mismatch.\n", bit);
          return 1;
        }
      }
      if (acc96_evalfx_from_tail(
              checkpoint_delta, &pp->ba, selector,
              selector ^ (UINT32_C(1) << ACC96_TAIL4_BIT), &pp->ua,
              ACC96_TAIL4_BIT, checkpoint, &ctx)) {
          fprintf(stderr, "q88 T4 accepted an uncovered first difference.\n");
        return 1;
      }
      free(checkpoint_delta);
      free(local_delta);
      free(checkpoint);
      printf("q88 T4 exhaustive prefix check passed: bits 0..3 exact, "
             "bit 4 safely rejected.\n");
    }

    /* DELTA6 selection: every first-difference bit 0..5 must agree with
       LOCAL, while bit 6 is deliberately rejected before mutating state. */
    {
      acc96_relation_witness_t *base6 = calloc(1, sizeof(*base6));
      acc96_delta6_t *delta6 = calloc(1, sizeof(*delta6));
      size_t layer_packet_total = 0;
      if (delta6 == NULL || base6 == NULL)
        return 1;
      memcpy(base6, witness, sizeof(*base6));
      for (size_t bit = 0; bit < ACC96_DELTA6_LAYERS; bit++) {
        acc96_relation_witness_t *local6 = calloc(1, sizeof(*local6));
        acc96_relation_witness_t *fast6 = calloc(1, sizeof(*fast6));
        acc96_relation_update_t update6 = {
            .sequence = (uint32_t)(bit + 1),
            .revoked_handle = member ^ (UINT32_C(1) << bit)};
        size_t used6 = 0;
        if (local6 == NULL || fast6 == NULL) return 1;
        if (!acc96_relation_delta6(delta6, pp, update6.revoked_handle, &ctx))
          return 1;
        memcpy(local6, base6, sizeof(*local6));
        memcpy(fast6, base6, sizeof(*fast6));
        if (!acc96_relation_memwitup(local6, pp, member,
                                     &update6, &ctx) ||
            !acc96_relation_memwitup_delta6(fast6, member, &update6,
                                            delta6, bit + 1, &used6) ||
            used6 != bit + 1 || memcmp(local6, fast6, sizeof(*local6)) != 0) {
          free(local6);
          free(fast6);
          fprintf(stderr, "q88 DELTA6 bit %zu/local mismatch.\n", bit);
          return 1;
        }
        free(local6);
        free(fast6);
      }
      for (size_t bit = ACC96_DELTA6_LAYERS;
           bit < ACC96_HANDLE_BITS; bit++) {
        acc96_relation_witness_t *saved = calloc(1, sizeof(*saved));
        acc96_relation_update_t outside = {
            .sequence = (uint32_t)(bit + 1),
            .revoked_handle = member ^ (UINT32_C(1) << bit)};
        size_t used6 = 77;
        if (saved == NULL) return 1;
        memcpy(saved, base6, sizeof(*saved));
        if (acc96_relation_memwitup_delta6(base6, member, &outside, delta6,
                                           ACC96_DELTA6_LAYERS, &used6) ||
            used6 != 0 || memcmp(base6, saved, sizeof(*saved)) != 0) {
          free(saved);
          fprintf(stderr, "q88 DELTA6 uncovered bit %zu mutated witness.\n",
                  bit);
          return 1;
        }
        free(saved);
      }
      for (size_t layer = 0; layer < ACC96_DELTA6_LAYERS; layer++) {
        uint8_t *wire = NULL;
        size_t wire_len = 0, decoded_layer = 99;
        acc96_delta6_t *decoded = calloc(1, sizeof(*decoded));
        if (!acc96_delta6_pack_layer_alloc(&wire, &wire_len, delta6, layer,
                                           &ctx) ||
            decoded == NULL ||
            !acc96_delta6_unpack_layer(decoded, &decoded_layer, wire,
                                       wire_len, &ctx) ||
            decoded_layer != layer ||
            memcmp(decoded->right[layer], delta6->right[layer],
                   sizeof(decoded->right[layer])) != 0) {
          free(wire);
          free(decoded);
          fprintf(stderr, "q88 DELTA6 layer codec mismatch.\n");
          return 1;
        }
        layer_packet_total += wire_len;
        printf("q88 DELTA6 layer %zu packet: %zu bytes (%.2f KiB).\n",
               layer, wire_len, wire_len / 1024.0);
        {
          acc96_delta6_t *saved = malloc(sizeof(*saved));
          if (saved == NULL) return 1;
          memcpy(saved, decoded, sizeof(*saved));
          wire[28] = (uint8_t)ACC96_DELTA6_LAYERS;
          if (acc96_delta6_unpack_layer(decoded, NULL, wire, wire_len, &ctx) ||
              memcmp(saved, decoded, sizeof(*saved)) != 0) {
            free(saved);
            free(wire);
            free(decoded);
            fprintf(stderr, "q88 DELTA6 invalid layer index mutated output.\n");
            return 1;
          }
          free(saved);
        }
        free(wire);
        free(decoded);
      }
      printf("q88 DELTA6 six-layer packets: %zu bytes total (%.2f KiB), "
             "%.2f bytes mean (%.2f KiB).\n",
             layer_packet_total, layer_packet_total / 1024.0,
             (double)layer_packet_total / ACC96_DELTA6_LAYERS,
             (double)layer_packet_total /
                 (ACC96_DELTA6_LAYERS * 1024.0));
      free(base6);
      free(delta6);
      printf("q88 DELTA6 exhaustive selector check passed: bits 0..5 exact, "
             "bits 6..23 safely rejected without mutation.\n");
    }
    if (!benchmark_delta6_updates(pp, accumulator, witness, member, &ctx)) {
      fprintf(stderr, "q88 DELTA6 update benchmark failed.\n");
      return 1;
    }
  }

  double delete_total = 0;
  double delta_build_total = 0;
  double combined_delete_total = 0;
  double local_update_total = 0;
  double delta_update_total = 0;
  for (size_t i = 0; i < sizeof(revoked) / sizeof(revoked[0]); i++) {
    acc96_relation_update_t update;
    acc96_delta4_t *delta4 = malloc(sizeof(*delta4));
    acc96_delta4_t *combined_delta4 = malloc(sizeof(*combined_delta4));
    acc96_poly_t *combined_accumulator = malloc(sizeof(*combined_accumulator));
    acc96_poly_t *tail_accumulator = malloc(sizeof(*tail_accumulator));
    acc96_gmat_t *tail4 = malloc(sizeof(*tail4));
    acc96_relation_witness_t *local = malloc(sizeof(*local));
    acc96_relation_witness_t *fast = malloc(sizeof(*fast));
    size_t applied_layers = 0;
    if (delta4 == NULL || combined_delta4 == NULL ||
        combined_accumulator == NULL || tail_accumulator == NULL ||
        tail4 == NULL || local == NULL || fast == NULL)
      return 1;
    acc96_poly_set(combined_accumulator, accumulator);
    acc96_poly_set(tail_accumulator, accumulator);
    double started = seconds();
    if (!acc96_relation_delta4(delta4, pp, revoked[i], &ctx)) return 1;
    delta_build_total += seconds() - started;
    if (i == 0) {
      print_delta4_layer_sizes(delta4, &ctx);
      {
        size_t zero_count = 0;
        size_t wide37_count = 0;
        size_t total_count = (size_t)ACC96_DELTA4_LAYERS *
                             ACC96_GADGET_DIGITS * ACC96_N;
        for (size_t layer = 0; layer < ACC96_DELTA4_LAYERS; layer++)
          for (size_t digit = 0; digit < ACC96_GADGET_DIGITS; digit++)
            for (size_t coefficient = 0; coefficient < ACC96_N; coefficient++)
              {
                const acc96_i128 centered = acc96_coeff_to_centered(
                    &delta4->right[layer][digit].coeff[coefficient], &ctx);
                const acc96_u128 magnitude = centered < 0
                    ? (acc96_u128)(-centered) : (acc96_u128)centered;
                zero_count += magnitude == 0;
                wide37_count += magnitude > (UINT64_C(1) << 35);
              }
        printf("q88 DELTA4 zero coefficients: %zu/%zu (%.3f%%); "
               "37-bit-or-wider magnitudes: %zu/%zu (%.3f%%).\n",
               zero_count, total_count,
               100.0 * (double)zero_count / (double)total_count,
               wide37_count, total_count,
               100.0 * (double)wide37_count / (double)total_count);
      }
      if (!test_delta4_codec(delta4, &ctx)) {
        fprintf(stderr, "q88 DELTA4 block codec regression failed.\n");
        return 1;
      }
    }
    memcpy(local, witness, sizeof(*local));
    memcpy(fast, witness, sizeof(*fast));
    started = seconds();
    acc96_relation_delete(accumulator, &update, pp, revoked[i],
                          (uint32_t)i + 1, &ctx);
    delete_total += seconds() - started;
    {
      acc96_relation_update_t combined_update;
      started = seconds();
      if (!acc96_relation_delete_delta4(
              combined_accumulator, &combined_update, combined_delta4, pp,
              revoked[i], (uint32_t)i + 1, &ctx))
        return 1;
      combined_delete_total += seconds() - started;
      if (!acc96_poly_equal(accumulator, combined_accumulator) ||
          memcmp(delta4, combined_delta4, sizeof(*delta4)) != 0 ||
          memcmp(&update, &combined_update, sizeof(update)) != 0) {
        fprintf(stderr, "q88 fused Delete+DELTA4 mismatch.\n");
        return 1;
      }
    }
    {
      acc96_relation_update_t tail_update;
      if (!acc96_relation_delete_tail4(
              tail_accumulator, &tail_update, tail4, pp, revoked[i],
              (uint32_t)i + 1, &ctx) ||
          !acc96_poly_equal(accumulator, tail_accumulator) ||
          memcmp(&update, &tail_update, sizeof(update)) != 0) {
        fprintf(stderr, "q88 fused Delete+T4 mismatch.\n");
        return 1;
      }
    }
    started = seconds();
    if (!acc96_relation_memwitup(local, pp, member, &update, &ctx))
      return 1;
    local_update_total += seconds() - started;
    started = seconds();
    if (!acc96_relation_memwitup_delta4(fast, member, &update, delta4,
                                        &applied_layers))
      return 1;
    delta_update_total += seconds() - started;
    if (applied_layers != i + 1 || memcmp(local, fast, sizeof(*local)) != 0)
      return 1;
    memcpy(witness, fast, sizeof(*witness));
    free(fast);
    free(local);
    free(combined_accumulator);
    free(tail4);
    free(tail_accumulator);
    free(combined_delta4);
    free(delta4);
    if (!acc96_relation_memver(pp, accumulator, member, witness, 0, &ctx)) {
      fprintf(stderr, "q88 MemVer preservation failed after deletion %zu.\n",
              i + 1);
      return 1;
    }
  }

  {
    acc96_delta4_t *delta4 = malloc(sizeof(*delta4));
    acc96_relation_witness_t *saved = malloc(sizeof(*saved));
    const acc96_relation_update_t uncovered = {
        .sequence = 4, .revoked_handle = member ^ UINT32_C(0x10)};
    size_t applied_layers = 99;
    if (delta4 == NULL || saved == NULL ||
        !acc96_relation_delta4(delta4, pp, uncovered.revoked_handle, &ctx))
      return 1;
    memcpy(saved, witness, sizeof(*saved));
    if (acc96_relation_memwitup_delta4(witness, member, &uncovered, delta4,
                                       &applied_layers) ||
        applied_layers != 0 || memcmp(saved, witness, sizeof(*saved)) != 0) {
      fprintf(stderr, "q88 DELTA4 uncovered-prefix mutation.\n");
      return 1;
    }
    free(saved);
    free(delta4);
  }

  /* Both relation and norm gates are active in production. */
  const long double final_norm = acc96_relation_witness_norm2(witness, &ctx);
  const long double max_chunk_norm = max_backend_chunk_norm2(witness, &ctx);
  if (!acc96_relation_memver(pp, accumulator, member, witness,
                            final_norm + 1, &ctx) ||
      acc96_relation_memver(pp, accumulator, member, witness,
                            final_norm * 0.5L, &ctx)) {
    fprintf(stderr, "q88 norm gate failed.\n");
    return 1;
  }
  acc96_coeff_from_u128(&witness->left.coeff[0],
      (acc96_coeff_to_u128(&witness->left.coeff[0], &ctx) + 1) %
      acc96_modulus());
  if (acc96_relation_memver(pp, accumulator, member, witness, 0, &ctx)) {
    fprintf(stderr, "q88 tampered witness was accepted.\n");
    return 1;
  }

  printf("q88 MemVer/Delete/MemWitUp and DELTA4/LOCAL equality passed: "
         "%zu deletions, "
         "Delete %.3f + DELTA4 %.3f s/update, fused %.3f s/update, "
         "LOCAL %.3f s/update, DELTA4 add %.6f s/update, "
         "initial norm^2 %.0Lf, final norm^2 %.0Lf, "
         "max degree-256 chunk norm^2 %.0Lf, "
         "accumulator wire %.2f KiB, witness projection %.2f MiB.\n",
         sizeof(revoked) / sizeof(revoked[0]),
         delete_total / (sizeof(revoked) / sizeof(revoked[0])),
         delta_build_total / (sizeof(revoked) / sizeof(revoked[0])),
         combined_delete_total / (sizeof(revoked) / sizeof(revoked[0])),
         local_update_total / (sizeof(revoked) / sizeof(revoked[0])),
         delta_update_total / (sizeof(revoked) / sizeof(revoked[0])),
         initial_norm, final_norm, max_chunk_norm,
         ACC96_POLY_PACKED_BYTES / 1024.0,
         ((2 + ACC96_TRUNC_RIGHT +
           ACC96_HANDLE_BITS * ACC96_GADGET_DIGITS) *
              ACC96_POLY_PACKED_BYTES + ACC96_N / 8.0) / 1048576.0);
  free(accumulator);
  free(witness);
  free(pp);
  return 0;
}
