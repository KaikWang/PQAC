#ifndef REVOCABLE_AC_H
#define REVOCABLE_AC_H

#include <stdint.h>
#include <stddef.h>

#include "acc_genpp.h"
#include "acc96_accumulator.h"
#include "acc96_exact_lift.h"
#include "osig.h"
#include "show.h"

#define REVAC_MAX_UAVS ACC_GENPP_MAX_DICT_ENTRIES
#define REVAC_GCS_CERT_VALID 1
#define REVAC_SHOW_NONCE_BYTES REVAC_PARAM_SHOW_NONCE_BYTES
#define REVAC_SAME_HANDLE_COMMITMENT_BYTES \
  ((PARAM_D * PARAM_N * PARAM_Q_BITLEN + 7) / 8)
#define REVAC_CONTEXT_HASH_BYTES REVAC_PARAM_HASH_BYTES
#define REVAC_TA_SYNC_DELTA_LAYERS REVAC_PARAM_PRIVATE_DELTA_LAYERS
#define REVAC_SYNC_NONCE_BYTES REVAC_PARAM_SYNC_NONCE_BYTES
#define REVAC_Q96_DELTA_LAYERS ACC96_DELTA6_LAYERS

typedef enum {
  REVAC_ACC_PROFILE_Q19 = 0,
  REVAC_ACC_PROFILE_Q96 = 1
} revac_acc_profile_t;

typedef struct acc96_sampleleft_ctx acc96_sampleleft_ctx_t;

typedef enum {
  REVAC_CONTEXT_PROFILE_LOCAL = 0,
  REVAC_CONTEXT_PROFILE_DELTA4 = 1
} revac_context_profile_t;

/* One q19 public BML broadcast record.  The small revocation update and its
   mandatory, authenticated DELTA4 payload have one shared lifetime. */
typedef struct {
  acc_genpp_update_t update;
  uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t delta_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t *delta_data;
  size_t delta_data_len;
  uint8_t *signature;
  size_t signature_len;
  int used;
} revac_bml_entry_t;

/* TA-private, reusable DELTA6 material for one revocation event.  Unlike the
   public DELTA4 packet this cache is not signed or exported.  It is retained
   for the accumulator epoch until the safe baseline watermark passes seq. */
typedef struct {
  uint32_t seq;
  uint32_t handle;
  int8_t sign;
  uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t delta_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t *delta_data;
  size_t delta_data_len;
  int used;
} revac_sync_delta_cache_entry_t;

/* One TA-owned q96 BML record.  The production wire carries the exact T4
   EvalF checkpoint; delta_data is a generic authenticated payload name. */
typedef struct {
  acc96_relation_update_t update;
  uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t accumulator_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t delta_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t *delta_data;
  size_t delta_data_len;
  uint8_t *signature;
  size_t signature_len;
  int used;
} revac_q96_bml_entry_t;

/* Public q88 DELTA6 manifest plus six independently fetchable layer blobs.
 * The existing q96_bml records remain as a compatibility T4 path; new UAVs
 * use this authenticated, selective-prefix representation. */
typedef struct {
  acc96_relation_update_t update;
  uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t accumulator_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t layer_hash[REVAC_Q96_DELTA_LAYERS][REVAC_CONTEXT_HASH_BYTES];
  uint32_t layer_len[REVAC_Q96_DELTA_LAYERS];
  uint8_t *layer_data[REVAC_Q96_DELTA_LAYERS];
  uint8_t *signature;
  size_t signature_len;
  int used;
} revac_q96_delta6_entry_t;

typedef struct {
  revac_context_profile_t profile;
  uint32_t seq;
  uint32_t handle;
  int8_t sign;
  poly_q_vec_d delta[ACC_GENPP_DELTA4_LAYERS][ACC_GENPP_GADGET_DIM];
} revac_delta_context_t;

typedef enum {
  REVAC_ACC_ZK_NONE = 0,
  REVAC_ACC_ZK_LABRADOR_MEMVER = 1,
  REVAC_ACC_ZK_LABRADOR_MEMVER_DEFLATE = 2,
  REVAC_ACC_ZK_Q96_LABRADOR_MEMVER = 3,
  REVAC_ACC_ZK_Q96_LABRADOR_MEMVER_DEFLATE = 4
} revac_acc_zk_kind_t;

typedef struct {
  uint64_t id;
  uint32_t handle;
  user_pk_t pk;
  int used;
  int issued;
  int revoked;
  /* Accumulator deletion sequence at issuance.  Updates with seq greater
     than this baseline must be replayed by MemWitSync. */
  uint32_t baseline_seq;
} revac_uav_record_t;

typedef struct {
  sep_pk_t ta_pk;
  sep_sk_t ta_sk;
  acc_genpp_pp_t acc_pp;
  acc_genpp_td_t acc_td;
  acc_genpp_state_t acc_st;
  poly_q_vec_d acc;
  uint8_t context_pp_hash[REVAC_CONTEXT_HASH_BYTES];
  /* The revocable-credential BML is a ring of complete, signed DELTA4
     broadcasts.  The core accumulator retains its own small update log. */
  revac_bml_entry_t bml[ACC_GENPP_BML_CAPACITY];
  revac_q96_bml_entry_t q96_bml[ACC_GENPP_BML_CAPACITY];
  revac_q96_delta6_entry_t q96_delta6_bml[ACC_GENPP_BML_CAPACITY];
  /* Private sync acceleration is decoupled from BML and can retain every
     deletion in the current accumulator epoch. */
  revac_sync_delta_cache_entry_t
      sync_delta_cache[ACC_GENPP_MAX_DELETIONS];
  revac_uav_record_t uavs[REVAC_MAX_UAVS];
  /* Default q96 accumulator backend.  Pointers keep the legacy q19 entity
     footprint unchanged until production Setup allocates q96 state. */
  revac_acc_profile_t acc_profile;
  acc96_ctx_t acc96_ctx;
  acc96_accumulator_t *acc96_state;
  /* Public seed retained so a remote prover/verifier expands the fixed q96
     matrices locally instead of transferring the multi-megabyte pp. */
  uint8_t acc96_public_seed[ACC96_SEED_BYTES];
  /* Remark 4.3 refresh-window identifier.  Setup starts at window zero; this
     implementation fails closed at the 50000-deletion boundary until an
     authenticated refresh/reissue operation is implemented. */
  uint32_t acc96_window_id;
  acc96_sampleleft_fn acc96_sampler;
  void *acc96_sampler_opaque;
  /* Non-NULL only for the Setup-owned production q96 sampler cache. */
  acc96_sampleleft_ctx_t *acc96_secure_sampler;
  long double acc96_witness_norm_sq_bound;
} revac_ta_t;

/* Public, read-only state installed on UAV/GCS devices.  It deliberately
   contains no TA signing key, accumulator trapdoor, PRF key, or member DB. */
typedef struct {
  revac_acc_profile_t acc_profile;
  const sep_pk_t *ta_pk;
  const acc_genpp_pp_t *acc_pp;
  const __poly_q_vec_d *acc;
  uint32_t acc_seq;
  uint8_t context_pp_hash[REVAC_CONTEXT_HASH_BYTES];
  uint32_t window_id;
  /* Explicit digest of the accumulator/credential parameter set.  It is
     currently equal to context_pp_hash, but remains a separate protocol field
     so a future refresh can authenticate independently versioned parameters. */
  uint8_t parameter_digest[REVAC_CONTEXT_HASH_BYTES];
  /* TA-owned read-only BML records.  A cross-device decoder owns equivalent
     records after authenticating their DELTA4 packets. */
  const revac_bml_entry_t *bml[ACC_GENPP_BML_CAPACITY];
  size_t bml_len;
  const revac_q96_bml_entry_t *q96_bml[ACC_GENPP_BML_CAPACITY];
  size_t q96_bml_len;
  const acc96_ctx_t *acc96_ctx;
  const acc96_relation_pp_t *acc96_pp;
  const acc96_poly_t *acc96_acc;
  uint8_t acc96_public_seed[ACC96_SEED_BYTES];
  uint32_t acc96_seq;
  long double acc96_witness_norm_sq_bound;
} revac_public_state_t;

typedef struct {
  uint32_t handle;
  uint32_t last_seq;
  uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t nonce[REVAC_SYNC_NONCE_BYTES];
} revac_sync_request_t;

typedef struct {
  uint32_t revoked_handle;
  uint8_t pp_hash[REVAC_CONTEXT_HASH_BYTES];
  uint8_t tail_hash[REVAC_CONTEXT_HASH_BYTES];
  /* Keep the public checkpoint in its canonical 12-byte/CRT-coefficient
     representation.  Decoding only on a cache hit halves resident storage
     compared with an acc96_gmat_t (three uint64_t limbs per coefficient). */
  uint8_t *tail4_packed;
  int used;
} revac_q96_tail_cache_entry_t;

typedef struct {
  uint64_t id;
  user_pk_t pk;
  user_sk_t sk;
  uint32_t handle;
  sep_sig_t sig;
  acc_genpp_wit_t witness;
  acc96_relation_witness_t *acc96_witness;
  uint32_t last_seq;
  /* Successful clear MemVer cache.  The key commits to the public
     accumulator and the complete witness, so direct witness mutation cannot
     turn this performance cache into a verification bypass. */
  uint8_t memver_cache_key[REVAC_CONTEXT_HASH_BYTES];
  uint32_t memver_cache_seq;
  int memver_cache_valid;
  /* The accumulator-side exact q-lift depends on exactly the same clear
     MemVer input and is independent of the per-show nonce/link opening.
     Keeping 4*256 int64 coefficients avoids recomputing it for every show. */
  int64_t memver_lift_cache[PARAM_D][PARAM_N];
  uint8_t memver_lift_cache_key[REVAC_CONTEXT_HASH_BYTES];
  int memver_lift_cache_valid;
  /* q88 exact lift is bound to the full validated MemVer digest. */
  acc96_exact_lift_t *memver_q88_lift_cache;
  uint8_t memver_q88_lift_key[REVAC_CONTEXT_HASH_BYTES];
  int memver_q88_lift_valid;
  unsigned memver_q88_lift_bits;
  /* Optional idle-time T4 prefetch cache.  Entries are public update
     acceleration data, not credentials or secret witness material. */
  revac_q96_tail_cache_entry_t q96_tail_cache[ACC96_MAX_MEMBERS];
  int registered;
  int issued;
} revac_uav_t;

typedef struct {
  uint64_t id;
  user_pk_t pk;
  user_sk_t sk;
  uint64_t timestamp;
  int registered;
  int cert_valid;
} revac_gcs_t;

typedef struct {
  show_proof_t sig_show;
  uint8_t nonce[REVAC_SHOW_NONCE_BYTES];
  /* Randomized, hiding commitment to the one handle shared by the signature
     and accumulator proof relations.  Its opening is never serialized. */
  poly_q_vec_d handle_commitment;
  revac_acc_zk_kind_t acc_zk_kind;
  uint8_t *acc_zk_proof;
  size_t acc_zk_proof_len;
} revac_show_proof_t;

typedef struct {
  uint8_t user_seed[SEED_BYTES];
  revac_acc_profile_t acc_profile;
  uint32_t acc_seq;
  uint32_t window_id;
  uint8_t parameter_digest[REVAC_CONTEXT_HASH_BYTES];
  /* Canonical digest of the accumulator value.  This keeps the show context
     constant-size for both q19 and q96 while binding the CRS to the exact
     public accumulator state. */
  uint8_t acc_digest[REVAC_CONTEXT_HASH_BYTES];
  /* Retained for the q19 Labrador bridge.  q96 uses acc_digest and never
     treats an all-zero q19 vector as its accumulator. */
  poly_q_vec_d acc;
} revac_show_context_t;

uint32_t revac_handle_from_identity(const user_pk_t *pk, uint64_t id);
void revac_same_handle_public_matrices(
    poly_q_mat_d_d A_link,
    poly_q_mat_d_d D_link,
    const revac_public_state_t *public_state);
void revac_same_handle_commitment_generate(
    poly_q_vec_d commitment,
    poly_q_vec_d randomness[2],
    const revac_public_state_t *public_state,
    uint32_t handle);
int revac_same_handle_exact_lift(
    int64_t lift[PARAM_D][PARAM_N],
    const revac_public_state_t *public_state,
    uint32_t handle,
    const poly_q_vec_d randomness[2],
    const poly_q_vec_d commitment);
void revac_show_proof_init(revac_show_proof_t *proof);
void revac_show_proof_clear(revac_show_proof_t *proof);
int revac_show_proof_set_acc_zk(revac_show_proof_t *proof,
                                revac_acc_zk_kind_t kind,
                                const uint8_t *bytes,
                                size_t len);
void revac_show_context_init(revac_show_context_t *ctx);
void revac_show_context_clear(revac_show_context_t *ctx);
void revac_show_context_from_ta(revac_show_context_t *ctx,
                                const revac_ta_t *ta,
                                const uint8_t user_seed[SEED_BYTES]);
void revac_show_context_from_public(
    revac_show_context_t *ctx,
    const revac_public_state_t *public_state,
    const uint8_t user_seed[SEED_BYTES]);
int revac_public_accumulator_digest(
    uint8_t digest[REVAC_CONTEXT_HASH_BYTES],
    const revac_public_state_t *public_state);
int revac_show_context_matches_public(
    const revac_show_context_t *ctx,
    const revac_public_state_t *public_state);

void revac_show_lang_init(poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
                          poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K],
                          poly_qshow_mat_k_k A3_embed[PARAM_D],
                          poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D],
                          poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
                          poly_qshow_mat_k_k Dx_embed[PARAM_D],
                          poly_qshow_vec_k u_embed[PARAM_D]);

void revac_show_lang_clear(poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
                           poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K],
                           poly_qshow_mat_k_k A3_embed[PARAM_D],
                           poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D],
                           poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
                           poly_qshow_mat_k_k Dx_embed[PARAM_D],
                           poly_qshow_vec_k u_embed[PARAM_D]);

void revac_show_public_embed(poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
                             poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D * PARAM_K],
                             poly_qshow_mat_k_k A3_embed[PARAM_D],
                             poly_qshow_mat_k_k Ds_embed[PARAM_D][2 * PARAM_D],
                             poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
                             poly_qshow_mat_k_k Dx_embed[PARAM_D],
                             poly_qshow_vec_k u_embed[PARAM_D],
                             const sep_pk_t *ta_pk,
                             const uint8_t user_seed[SEED_BYTES]);

void revac_show_derive_crs_seed(uint8_t crs_seed[CRS_SEED_BYTES],
                                const sep_pk_t *ta_pk,
                                const uint8_t user_seed[SEED_BYTES],
                                const uint8_t nonce[REVAC_SHOW_NONCE_BYTES],
                                revac_acc_profile_t acc_profile,
                                uint32_t acc_seq,
                                uint32_t window_id,
                                const uint8_t parameter_digest[REVAC_CONTEXT_HASH_BYTES],
                                const uint8_t acc_digest[REVAC_CONTEXT_HASH_BYTES]);

#endif
