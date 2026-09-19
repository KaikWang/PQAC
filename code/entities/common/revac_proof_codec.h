#ifndef REVAC_PROOF_CODEC_H
#define REVAC_PROOF_CODEC_H

#include <stddef.h>
#include <stdint.h>

#include "revocable_ac.h"

typedef struct {
  size_t magic;
  size_t tA;
  size_t hints;
  size_t tB;
  size_t z3;
  size_t h;
  size_t t1;
  size_t challenge_and_counter;
  size_t z1;
  size_t z2;
  size_t total;
} revac_show_signature_size_breakdown_t;

size_t revac_show_signature_raw_size_bytes(void);
size_t revac_show_signature_compact_bound_bytes(void);
int revac_show_signature_compact_size_breakdown(
    const show_proof_t *proof,
    revac_show_signature_size_breakdown_t *sizes);
size_t revac_show_proof_wire_size_bytes(const revac_show_proof_t *proof);

int revac_deflate_pack(uint8_t **out, size_t *out_len,
                       const uint8_t *in, size_t in_len);
int revac_deflate_unpack(uint8_t **out, size_t *out_len,
                         const uint8_t *in, size_t in_len);

int revac_show_signature_pack_compact(uint8_t **out, size_t *out_len,
                                      const show_proof_t *proof);
int revac_show_signature_unpack_compact(show_proof_t *proof,
                                        const uint8_t *in, size_t in_len);

int revac_show_proof_pack_compact(uint8_t **out, size_t *out_len,
                                  const revac_show_proof_t *proof,
                                  int deflate_acc_zk);
int revac_show_proof_unpack_compact(revac_show_proof_t *proof,
                                    const uint8_t *in, size_t in_len);

#endif
