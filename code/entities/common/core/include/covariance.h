#ifndef COVARIANCE_H
#define COVARIANCE_H

#include <stddef.h>
#include <stdint.h>
#include "arith.h"

void fft_precomp_setup(void);
void fft_precomp_teardown(void);

void compute_covariance(poly_real_mat_2d_2d S, const poly_q_mat_d_d RRstar[2][2]);

int sk_sq_spectral_norm(poly_q_mat_d_d RRstar[2][2], const poly_q_mat_d_d R[2][PARAM_KH]);

int sk_sq_spectral_norm_columns(poly_q_mat_d_d RRstar[2][2],
                                size_t columns,
                                const poly_q_mat_d_d R[2][columns],
                                double max_sq_norm);

#endif
