#ifndef REVAC_ACC96_EVAL_H
#define REVAC_ACC96_EVAL_H

#include "acc80_eval.h"
#include "acc96_crt.h"

typedef acc80_gmat_t acc96_gmat_t;
typedef acc80_ba_t acc96_ba_t;
typedef acc80_evalfx_delta_t acc96_evalfx_delta_t;

#define ACC96_HANDLE_BITS ACC80_HANDLE_BITS
#define acc96_gadget_identity acc80_gadget_identity
#define acc96_ginv_reconstructs acc80_ginv_reconstructs
#define acc96_evalf acc80_evalf
#define acc96_evalf_with_tail acc80_evalf_with_tail
#define acc96_evalf_apply acc80_evalf_apply
#define acc96_evalfx acc80_evalfx
#define acc96_evalfx_prefix acc80_evalfx_prefix
#define acc96_evalf_evalfx_prefix acc80_evalf_evalfx_prefix
#define acc96_evalf_tail_checkpoint acc80_evalf_tail_checkpoint
#define acc96_evalfx_from_tail acc80_evalfx_from_tail
#define acc96_member_matrix_apply acc80_member_matrix_apply

#endif
