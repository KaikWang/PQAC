#ifndef SAMPLING_H
#define SAMPLING_H

#include <stdint.h>
#include "params.h"
#include "arith.h"

void perturbation_sampler(
    poly_q_vec_d p1_p2[PARAM_KL],
    poly_q_vec_d p3,
    poly_q_vec_d p4[PARAM_KH],
    const poly_q_mat_d_d R[2][PARAM_KH],
    const poly_real_mat_2d_2d S);

void gadget_sampler(
    poly_q_vec_d zL[PARAM_KL],
    poly_q_vec_d zH[PARAM_KH],
    const poly_q_vec_d w,
    const poly_q tag,
    const poly_q taginv);

void sampler(
    poly_q_vec_d v11,
    poly_q_vec_d v12,
    poly_q_vec_d v2[PARAM_KH],
    const poly_q_mat_d_d R[2][PARAM_KH],
    const poly_q_mat_d_d A,
    const poly_q_mat_d_d B[PARAM_KH],
    const poly_q_vec_d u,
    const poly_q tag,
    const poly_q taginv,
    const poly_real_mat_2d_2d S);

#endif /* SAMPLING_H */
