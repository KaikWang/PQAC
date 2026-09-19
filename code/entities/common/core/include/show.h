#ifndef SHOW_H
#define SHOW_H

#include <stdint.h>
#include "osig.h"
#include "params.h"

// Byte-length of language packing 
#define SHOW_CHALLENGE_BASE_BYTES (1 + CRS_SEED_BYTES + SEED_BYTES) // no need to pack signer public key
// Byte-length of first challenge input
#define CHAL1_SHOW_INPUT_BYTES (SHOW_CHALLENGE_BASE_BYTES + 2*POLYQSHOW_VECD_PACKEDBYTES + POLYQSHOW_VEC256L_PACKEDBYTES)
// Byte-length of second challenge input 
#define CHAL2_SHOW_INPUT_BYTES (CHAL1_SHOW_INPUT_BYTES + PARAM_ARP_SHOW*COEFFQSHOW_PACKEDBYTES)
// Byte-length of third challenge input
#define CHAL3_SHOW_INPUT_BYTES (CHAL2_SHOW_INPUT_BYTES + POLYQSHOW_VECL_PACKEDBYTES)
// Byte-length of fourth challenge input
#define CHAL4_SHOW_INPUT_BYTES (CHAL3_SHOW_INPUT_BYTES + 2*POLYQSHOW_PACKEDBYTES)

// Starting index for v1 in s1
#define IDX_V1_SHOW 0
// Starting index for v12 in s1
#define IDX_V12_SHOW (IDX_V1_SHOW + PARAM_D*PARAM_K_SHOW)
// Starting index for v2 in s1 (TSampler: KH vec_d entries)
#define V2_SHOW_POLYS (PARAM_D*PARAM_KH*PARAM_K_SHOW)
#define IDX_V2_SHOW (IDX_V12_SHOW + PARAM_D*PARAM_K_SHOW + 1) // +1 for four squares polynomial
// Starting index for v3 in s1 (TSampler: scalar poly_q → PARAM_K_SHOW polys)
#define V3_SHOW_POLYS (PARAM_K_SHOW)
#define IDX_V3_SHOW (IDX_V2_SHOW + V2_SHOW_POLYS + 1) // +1 for four squares polynomial
// Starting index for tag in s1
#define IDX_TAG_SHOW (IDX_V3_SHOW + V3_SHOW_POLYS + 1) // +1 for four squares polynomial
// Starting index for usk in s1
#define IDX_USK_SHOW (IDX_TAG_SHOW + PARAM_K_SHOW)
// Starting index for m in s1
#define IDX_M_SHOW (IDX_USK_SHOW + 2*PARAM_D*PARAM_K_SHOW)
// Starting index for hidden revocation handle x in s1
#define IDX_X_SHOW (IDX_M_SHOW + PARAM_M*PARAM_K_SHOW)
#define IDX_LINK_R0_SHOW (IDX_X_SHOW + PARAM_K_SHOW)
#define IDX_LINK_R1_SHOW (IDX_LINK_R0_SHOW + PARAM_D*PARAM_K_SHOW)
/* One past the final witness polynomial.  Keep this equal to PARAM_M1_SHOW. */
#define IDX_SHOW_WITNESS_END (IDX_LINK_R1_SHOW + PARAM_D*PARAM_K_SHOW)

#if PARAM_M1_SHOW != IDX_SHOW_WITNESS_END
#error "PARAM_M1_SHOW does not match the show witness layout"
#endif
#if ((PARAM_Q_SHOW - 1) % PARAM_GAMMA_SHOW) != 0
#error "PARAM_GAMMA_SHOW must divide PARAM_Q_SHOW - 1"
#endif
#if (PARAM_HINT_MOD_SHOW % 2) != 0
#error "Generalized hint modulus must be even"
#endif

typedef struct {
  poly_qshow_vec_d tA;
  poly_qshow_vec_d hints;
  poly_qshow_vec_256_l tB;
  coeff_qshow z3[PARAM_ARP_SHOW];
  poly_qshow_vec_l h;
  poly_qshow t1;
  poly_qshow c;
  uint32_t ctr_c;
  poly_qshow_vec_m1 z1;
  poly_qshow_vec_m2 z2;
} show_proof_t;

void show_power2round_vec_d(poly_qshow_vec_d high, poly_qshow_vec_d low,
                            const poly_qshow_vec_d value);
void show_highbits_vec_d(poly_qshow_vec_d high,
                         const poly_qshow_vec_d value);
void show_make_g_hints_vec_d(poly_qshow_vec_d hints,
                             const poly_qshow_vec_d z,
                             const poly_qshow_vec_d r);
void show_use_g_hints_vec_d(poly_qshow_vec_d high,
                            const poly_qshow_vec_d hints,
                            const poly_qshow_vec_d r);

void show_proof_init(show_proof_t *proof);
void show_proof_clear(show_proof_t *proof);

/* Bound only the current thread's OpenMP teams while the anonymous-signature
   prover overlaps another proof component.  Passing zero restores the normal
   standalone defaults.  This is a scheduling hint and never enters the
   transcript. */
void show_user_prove_set_thread_cap(unsigned int cap);
void show_verify_set_thread_cap(unsigned int cap);

void show_user_embed(
    poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
    poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D*PARAM_K],
    poly_qshow_mat_k_k A3_embed[PARAM_D],
    poly_qshow_mat_k_k Ds_embed[PARAM_D][2*PARAM_D],
    poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
    poly_qshow_vec_k   u_embed[PARAM_D],
    poly_qshow_vec_m1  s1,
    const user_pk_t    *upk,
    const user_sk_t    *usk,
    const sep_pk_t     *pk,
    const sep_sig_t    *sig,
    const uint8_t      msg[PARAM_M*PARAM_N/8]);
void show_user_embed_handle(
    poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D*PARAM_K], 
    poly_qshow_mat_k_k A3_embed[PARAM_D],
    poly_qshow_mat_k_k Ds_embed[PARAM_D][2*PARAM_D], 
    poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], 
    poly_qshow_mat_k_k Dx_embed[PARAM_D],
    poly_qshow_vec_k   u_embed[PARAM_D], 
    poly_qshow_vec_m1  s1, 
    const user_pk_t    *upk, 
    const user_sk_t    *usk, 
    const sep_pk_t     *pk, 
    const sep_sig_t    *sig,
    const uint8_t      msg[PARAM_M*PARAM_N/8],
    uint32_t           x);
void show_user_embed_link_randomness(
    poly_qshow_vec_m1 s1,
    const poly_q_vec_d randomness[2]);
void show_link_public_embed(
    poly_qshow_mat_k_k A_link_embed[PARAM_D][PARAM_D],
    poly_qshow_mat_k_k D_link_embed[PARAM_D],
    poly_qshow_vec_k C_link_embed[PARAM_D],
    const poly_q_mat_d_d A_link,
    const poly_q_mat_d_d D_link,
    const poly_q_vec_d C_link);
void show_user_prove(
    show_proof_t             *proof,
    const poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
    const poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D*PARAM_K],
    const poly_qshow_mat_k_k A3_embed[PARAM_D],
    const poly_qshow_mat_k_k Ds_embed[PARAM_D][2*PARAM_D],
    const poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
    const poly_qshow_vec_m1  s1,
    const uint8_t            crs_seed[CRS_SEED_BYTES],
    const uint8_t            seed[SEED_BYTES]);
void show_user_prove_handle(
    show_proof_t             *proof, 
    const poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D*PARAM_K], 
    const poly_qshow_mat_k_k A3_embed[PARAM_D],
    const poly_qshow_mat_k_k Ds_embed[PARAM_D][2*PARAM_D], 
    const poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], 
    const poly_qshow_mat_k_k Dx_embed[PARAM_D],
    const poly_qshow_mat_k_k A_link_embed[PARAM_D][PARAM_D],
    const poly_qshow_mat_k_k D_link_embed[PARAM_D],
    const poly_qshow_vec_k   C_link_embed[PARAM_D],
    const poly_qshow_vec_m1  s1, 
    const uint8_t            crs_seed[CRS_SEED_BYTES], 
    const uint8_t            seed[SEED_BYTES]);
int show_verify(
    const show_proof_t       *proof,
    const poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D],
    const poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D*PARAM_K],
    const poly_qshow_mat_k_k A3_embed[PARAM_D],
    const poly_qshow_mat_k_k Ds_embed[PARAM_D][2*PARAM_D],
    const poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M],
    const poly_qshow_vec_k   u_embed[PARAM_D],
    const uint8_t            crs_seed[CRS_SEED_BYTES],
    const uint8_t            seed[SEED_BYTES]);
int show_verify_handle(
    const show_proof_t       *proof, 
    const poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D*PARAM_K], 
    const poly_qshow_mat_k_k A3_embed[PARAM_D],
    const poly_qshow_mat_k_k Ds_embed[PARAM_D][2*PARAM_D], 
    const poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], 
    const poly_qshow_mat_k_k Dx_embed[PARAM_D],
    const poly_qshow_mat_k_k A_link_embed[PARAM_D][PARAM_D],
    const poly_qshow_mat_k_k D_link_embed[PARAM_D],
    const poly_qshow_vec_k   C_link_embed[PARAM_D],
    const poly_qshow_vec_k   u_embed[PARAM_D], 
    const uint8_t            crs_seed[CRS_SEED_BYTES], 
    const uint8_t            seed[SEED_BYTES]);

#endif
