#ifndef PARAMS_H
#define PARAMS_H

/*************************************************
* TSampler Signature parameters
* Based on tag-friendly_sampler s4.1 standard model
**************************************************/

/* Domain separators for XOF — TSampler signature */
#define DOMAIN_SEPARATOR_A 0
#define DOMAIN_SEPARATOR_R 1
#define DOMAIN_SEPARATOR_A3 2
#define DOMAIN_SEPARATOR_U 3
#define DOMAIN_SEPARATOR_D 4
/* osig issuance */
#define DOMAIN_SEPARATOR_DS 5
#define DOMAIN_SEPARATOR_S 6
#define DOMAIN_SEPARATOR_RAND 7
#define DOMAIN_SEPARATOR_A1_ISS 8
#define DOMAIN_SEPARATOR_A2_ISS 9
#define DOMAIN_SEPARATOR_BYG_ISS 10
#define DOMAIN_SEPARATOR_B_ISS 11
#define DOMAIN_SEPARATOR_CHAL1_ISS 12
#define DOMAIN_SEPARATOR_CHAL2_ISS 13
#define DOMAIN_SEPARATOR_CHAL3_ISS 14
#define DOMAIN_SEPARATOR_CHAL4_ISS 15
#define DOMAIN_SEPARATOR_RAND_S2_ISS 16
#define DOMAIN_SEPARATOR_RAND_G_ISS 17
/* show */
#define DOMAIN_SEPARATOR_A1_SHOW 18
#define DOMAIN_SEPARATOR_A2_SHOW 19
#define DOMAIN_SEPARATOR_BYG_SHOW 20
#define DOMAIN_SEPARATOR_B_SHOW 21
#define DOMAIN_SEPARATOR_CHAL1_SHOW 22
#define DOMAIN_SEPARATOR_CHAL2_SHOW 23
#define DOMAIN_SEPARATOR_CHAL3_SHOW 24
#define DOMAIN_SEPARATOR_CHAL4_SHOW 25
#define DOMAIN_SEPARATOR_RAND_S2_SHOW 26
#define DOMAIN_SEPARATOR_RAND_G_SHOW 27
#define DOMAIN_SEPARATOR_DX 28
#define DOMAIN_SEPARATOR_REVAC_LINK_A 29
#define DOMAIN_SEPARATOR_REVAC_LINK_D 30
#define DOMAIN_SEPARATOR_REVAC_LINK_R 31

/* Signature parameters */
#define PARAM_N 256
#define PARAM_Q 405769L        /* q_L * q_H = 2401 * 169 */
#define PARAM_QL 2401L         /* b_L^k_L = 49^2 */
#define PARAM_QH 169L          /* b_H^k_H = 13^2 */
#define PARAM_Q_BITLEN 19
#define PARAM_D 4
#define PARAM_KH 2             /* G_H dimension */
#define PARAM_KL 2             /* G_L dimension */
#define PARAM_BH 13            /* G_H base */
#define PARAM_BL 49            /* G_L base */
/* keep PARAM_K for backward compat with proof layer */
#define PARAM_K PARAM_KH
#define PARAM_B PARAM_BL  /* compat: old SEP gadget base */
#define PARAM_M 10
#define PARAM_IT_SPEC_NORM 6
#define PARAM_W 5
#define PARAM_R_MAX_SQ_SPECTRAL_NORM 3898.89415332639373446000

/* Gaussian parameters */
#define PARAM_S1SQ_SLSQ 22816850.56626987829804420471
#define PARAM_S3SQ 22816838.86761695519089698792
#define PARAM_S4 76.49922142570575545051
#define PARAM_SL_DIV_BL 3.38718927843394679300
#define PARAM_SH_DIV_BH 3.38718927843394679300
#define PARAM_SQRT_S2SQ_SLSQ 3.42032935794184789202
#define PARAM_SQRT_S4SQ_SHSQ 62.55545719277940008851

/* Schur complement constants */
#define PARAM_NEG_SHSQS4SQ_DIV_S4SQ_SHSQ -2899.67458311796281122952
#define PARAM_NEG_SHSQ_DIV_S4SQ_SHSQ -0.49549038516081478356

/* Old SEP compat — map to new names where needed */
#define PARAM_S1SQ (PARAM_S1SQ_SLSQ + 2*PARAM_S3SQ) /* approximate compat */
#define PARAM_S2 3.42
#define PARAM_SQRT_S2SQ_SGSQ PARAM_SQRT_S2SQ_SLSQ
#define PARAM_SGINVSQ_S2INVSQ 1.0
#define PARAM_NEGSGSQ_DIV_S2SQ_SGSQ -0.5

/* Verification bounds */
#define PARAM_B11SQ 18149347086L
#define PARAM_B12SQ 4652284329L
#define PARAM_B2SQ 2238819L
#define PARAM_B3SQ 367252L

/* Hiding bounds for show proof */
#define PARAM_B11SQ_HIDING 18157970158L
#define PARAM_B12SQ_HIDING 4656650642L
/* compat: map SEP names to TSampler names */
#define PARAM_B1SQ (PARAM_B11SQ_HIDING + PARAM_B12SQ_HIDING)
#define PARAM_B1SQ_HIDING PARAM_B1SQ

/* Seeds */
#define SEED_BYTES 32
#define CRS_SEED_BYTES 32
#define STATE_BYTES 64

/*************************************************
* Placeholder AC proof params (reused from SEP)
* These need to be recomputed for TSampler.
**************************************************/
#define PARAM_N_ISS 64
#define PARAM_K_ISS 4
#define PARAM_Q_ISS 212704515569L
#define PARAM_Q_ISS_BITLEN 38
#define PARAM_Q1_ISS 524201L
#define PARAM_Q1_ISS_BITLEN 19
#define PARAM_Q2_ISS PARAM_Q
#define PARAM_Q1_INVMOD_Q2_ISS 85651L
#define PARAM_Q2_INVMOD_Q1_ISS 413551L
#define PARAM_D_ISS 20
#define PARAM_M1_ISS 104
#define PARAM_M1_K_ISS 26
#define PARAM_M2_ISS 58
#define PARAM_L_ISS 7
#define PARAM_ARP_ISS 256
#define PARAM_ARP_DIV_N_ISS 4
#define PARAM_ARP_DIV_N_L_ISS 11
#define PARAM_S1_ISS 369050.89730269293067976832
#define PARAM_S1SQ_ISS 136198564799.92280578613281250000
#define PARAM_S2_ISS 275602.77920886297943070531
#define PARAM_S2SQ_ISS 75956891907.64927673339843750000
#define PARAM_S3_ISS 72848.10643310110026504844
#define PARAM_S3SQ_ISS 5306846610.88842582702636718750
#define PARAM_REJ1_ISS 2
#define PARAM_REJ2_ISS 2
#define PARAM_REJ3_ISS 2
#define PARAM_B1SQ_ISS 180657566352359UL
#define PARAM_B2SQ_ISS 60411097594469UL
#define PARAM_B3SQ_ISS 584702794673UL
#define PARAM_RHO_ISS 8
#define PARAM_ETA_ISS 93

#define PARAM_N_SHOW 64
#define PARAM_K_SHOW 4
/*
 * Runtime-compatible 58-bit CRT modulus.  The corrected TSampler parameter
 * model for the pre-link 151-polynomial, ordinary-rejection show relation
 * selects a
 * roughly 75-bit proof modulus.  The current FLINT nmod/uniform-sampling
 * backend is limited to at most 60 bits, so these constants must not be
 * presented as a fresh 128-bit security derivation.  See SHOW_PROOF_SIZE.md.
 */
#define PARAM_Q_SHOW 216093359578546177L
#define PARAM_Q_SHOW_BITLEN 58
#define PARAM_Q1_SHOW 532552658233L
#define PARAM_Q1_SHOW_BITLEN 39
#define PARAM_Q2_SHOW PARAM_Q
#define PARAM_Q1_INVMOD_Q2_SHOW 274618L
#define PARAM_Q2_INVMOD_Q1_SHOW 172129496536L
#define PARAM_D_SHOW 23
#define PARAM_D_ROUND_SHOW 22
#define PARAM_GAMMA_SHOW 1073741824L
#define PARAM_HINT_MOD_SHOW ((PARAM_Q_SHOW - 1) / PARAM_GAMMA_SHOW)
/*
 * Show witness layout (all values are qshow polynomials):
 *   v11 || v12 plus one norm slack, v2 plus one norm slack,
 *   v3 plus one norm slack, tag, user secret key, message, handle x, and
 *   two binary R_q^d opening vectors for the randomized same-handle
 *   commitment shared with the accumulator proof.
 *
 * Keep this derived from the underlying signature parameters.  The previous
 * fixed value (183) left 32 unused witness polynomials after the handle.
 * They now hold the two opening vectors of the shared-handle commitment, so
 * the relation returns to 183 meaningful witness polynomials.
 */
#define PARAM_SHOW_V1_POLYS (2 * PARAM_D * PARAM_K_SHOW + 1)
#define PARAM_SHOW_V2_POLYS (PARAM_D * PARAM_KH * PARAM_K_SHOW + 1)
#define PARAM_SHOW_V3_POLYS (PARAM_K_SHOW + 1)
#define PARAM_SHOW_TAG_POLYS (PARAM_K_SHOW)
#define PARAM_SHOW_USK_POLYS (2 * PARAM_D * PARAM_K_SHOW)
#define PARAM_SHOW_MSG_POLYS (PARAM_M * PARAM_K_SHOW)
#define PARAM_SHOW_HANDLE_POLYS (PARAM_K_SHOW)
#define PARAM_SHOW_LINK_RANDOMNESS_POLYS (2 * PARAM_D * PARAM_K_SHOW)
#define PARAM_M1_SHOW (PARAM_SHOW_V1_POLYS + PARAM_SHOW_V2_POLYS + \
                       PARAM_SHOW_V3_POLYS + PARAM_SHOW_TAG_POLYS + \
                       PARAM_SHOW_USK_POLYS + PARAM_SHOW_MSG_POLYS + \
                       PARAM_SHOW_HANDLE_POLYS + \
                       PARAM_SHOW_LINK_RANDOMNESS_POLYS)
#define PARAM_M2_SHOW 74
#define PARAM_M2_COMMIT_SHOW (PARAM_M2_SHOW - PARAM_D_SHOW)
/* ceil(ceil(128/log2(PARAM_Q2_SHOW))/2) with garbage-term amplification. */
#define PARAM_L_SHOW 4
#define PARAM_ARP_SHOW 256
#define PARAM_ARP_DIV_N_SHOW 4
#define PARAM_ARP_DIV_N_L_SHOW (PARAM_ARP_DIV_N_SHOW + PARAM_L_SHOW)
/* Retuned for m1=183 and for the 2048 additional binary coefficients in the
 * two shared-handle opening vectors.  B1 uses the 2^-(sec+3) Gaussian tail
 * cut over m1*n coefficients; B3 follows the enlarged witness norm. */
#define PARAM_S1_SHOW 582380259.27222955226898193359
#define PARAM_S1SQ_SHOW 339166766389989312.00000000000000000000
#define PARAM_S2_SHOW 311304.54102290823357179761
#define PARAM_S2SQ_SHOW 96910517261.48355102539062500000
#define PARAM_S3_SHOW 114957853.84096613526344299316
#define PARAM_S3SQ_SHOW 13215308159720932.00000000000000000000
#define PARAM_REJ1_SHOW 2
#define PARAM_REJ2_SHOW 2
#define PARAM_REJ3_SHOW 2
#define PARAM_B1SQ_SHOW_LOW64 12307694251824381952UL
#define PARAM_B1SQ_SHOW_HIGH64 40UL
/*
 * Compressed-opening z2 bound:
 * (2*sqrt(B2_uncompressed) +
 *  (2^D*eta + gamma)*sqrt(n*d))^2.
 */
#define PARAM_B2SQ_SHOW_LOW64 1920750681752862720UL
#define PARAM_B2SQ_SHOW_HIGH64 171UL
#define PARAM_B3SQ_SHOW 1456048795079616000UL
#define PARAM_RHO_SHOW 8
#define PARAM_ETA_SHOW 93

/* Testing params */
#define Q_1_MOD 47
#define Q_1_MOD_INV 2
#define Q_MIN_MOD 31
#define Q_MIN_MOD_INV 44
#define Q_HAT_MOD (Q_1_MOD * Q_MIN_MOD)
#define N_HAT_RING 4

#endif /* PARAMS_H */
