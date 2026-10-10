#ifndef LABRADOS_H
#define LABRADOS_H

#include "proofsystem.h"
#include "dachshund.h"

#define py_init_witness NAMESPACE(py_init_witness)
__attribute__((visibility("default")))
void py_init_witness(witness wt, size_t r, size_t n[]);

#define py_comkey_reset NAMESPACE(py_comkey_reset)
__attribute__((visibility("default")))
void py_comkey_reset(void);

#define py_comkey_nbytes NAMESPACE(py_comkey_nbytes)
__attribute__((visibility("default")))
size_t py_comkey_nbytes(void);

#define py_comkey_serialize NAMESPACE(py_comkey_serialize)
__attribute__((visibility("default")))
int py_comkey_serialize(uint8_t *out, size_t outlen);

#define py_comkey_deserialize NAMESPACE(py_comkey_deserialize)
__attribute__((visibility("default")))
int py_comkey_deserialize(const uint8_t *in, size_t inlen);

#define py_print_witness_vector NAMESPACE(py_print_witness_vector)
__attribute__((visibility("default")))
int py_print_witness_vector(witness wt, size_t idx);

#define py_set_witness_vector NAMESPACE(py_set_witness_vector)
__attribute__((visibility("default")))
int py_set_witness_vector(witness wt, size_t idx, size_t n, size_t deg, 
                          const int64_t s[]);

#define py_init_statement NAMESPACE(py_init_statement)
__attribute__((visibility("default")))
void py_init_statement(statement st, size_t r, size_t n[], uint64_t normsq[], 
                       uint64_t normsq_req[], normtype normty[], 
                       size_t num_rq_cnst, size_t num_zq_cnst, 
                       size_t num_int_cnst);

#define py_statement_hash NAMESPACE(py_statement_hash)
__attribute__((visibility("default")))
void py_statement_hash(uint8_t out[16], const statement st);

#define py_append_constraint NAMESPACE(py_append_constraint)
__attribute__((visibility("default")))
int py_append_constraint(statement st, size_t nvec, const size_t idx[],
                         const size_t n[], size_t deg, int64_t *phi, int64_t *b,
                         int full);

#define py_append_quadratic NAMESPACE(py_append_quadratic)
__attribute__((visibility("default")))
int py_append_quadratic(statement st, size_t nlin, size_t nprod,
                        const size_t idx_lin[], const size_t idx_prod1[],
                        const size_t idx_prod2[], const size_t len_phi[],
                        size_t deg, int64_t *a, int64_t *phi, int64_t *b);

#define py_append_quadratic_monomials NAMESPACE(py_append_quadratic_monomials)
__attribute__((visibility("default")))
int py_append_quadratic_monomials(
    statement st, size_t nlin, size_t nprod, const size_t idx_lin[],
    const size_t idx_prod1[], const size_t idx_prod2[],
    const size_t len_phi[], size_t deg, int64_t *a,
    const uint32_t *monomial_exponents, int64_t *phi, int64_t *b);

#define py_append_deg0_constraint NAMESPACE(py_append_deg0_constraint)
__attribute__((visibility("default")))
int py_append_deg0_constraint(statement st, size_t idx, size_t deg);

#define py_append_deg0_vector_constraints NAMESPACE(py_append_deg0_vector_constraints)
__attribute__((visibility("default")))
int py_append_deg0_vector_constraints(statement st, size_t idx, size_t count);

#define py_append_deg0_logical_vector_constraints NAMESPACE(py_append_deg0_logical_vector_constraints)
__attribute__((visibility("default")))
int py_append_deg0_logical_vector_constraints(statement st, size_t idx,
                                               size_t count, size_t deg);

#define py_seeded_uniform_matrix_bank NAMESPACE(py_seeded_uniform_matrix_bank)
__attribute__((visibility("default")))
int py_seeded_uniform_matrix_bank(int32_t *out, const uint8_t seed[32],
                                  size_t matrices, size_t blocks,
                                  size_t rows, size_t cols, size_t degree,
                                  uint32_t q, uint32_t qbits,
                                  uint32_t domain_separator);

#define py_refresh_set_witness NAMESPACE(py_refresh_set_witness)
__attribute__((visibility("default")))
int py_refresh_set_witness(witness wt, const size_t *chunk_idx,
                           const size_t *lift_idx, const int32_t *bank,
                           const uint8_t *bits, size_t bit_count,
                           size_t blocks, size_t rows, size_t cols,
                           size_t degree, int64_t q, int64_t base,
                           size_t chunks);

#define py_statement_set_norms NAMESPACE(py_statement_set_norms)
__attribute__((visibility("default")))
int py_statement_set_norms(statement st, size_t count,
                           const uint64_t *normsq);

#define py_append_refresh_constraints NAMESPACE(py_append_refresh_constraints)
__attribute__((visibility("default")))
int py_append_refresh_constraints(statement st, size_t r_idx,
                                  const size_t *chunk_idx,
                                  const size_t *lift_idx,
                                  const int32_t *bank, size_t bit_count,
                                  size_t blocks, size_t rows, size_t cols,
                                  size_t degree, int64_t q, int64_t base,
                                  size_t chunks);

#define py_append_refresh_phase_constraints \
  NAMESPACE(py_append_refresh_phase_constraints)
__attribute__((visibility("default")))
int py_append_refresh_phase_constraints(
    statement st, size_t r_idx, const size_t *chunk_idx,
    const size_t *lift_idx, const int32_t *bank, size_t bit_count,
    size_t blocks, size_t lanes, size_t phases, size_t degree,
    const uint32_t *moduli, int64_t base, size_t chunks);

#define py_decompose_set_witness NAMESPACE(py_decompose_set_witness)
__attribute__((visibility("default")))
int py_decompose_set_witness(witness wt, const size_t *indices,
                             const int64_t *source, size_t poly_count,
                             size_t degree, int64_t base, size_t chunks);

#define py_malloc_trim NAMESPACE(py_malloc_trim)
__attribute__((visibility("default")))
int py_malloc_trim(void);

#define py_malloc_configure NAMESPACE(py_malloc_configure)
__attribute__((visibility("default")))
int py_malloc_configure(size_t mmap_threshold, size_t trim_threshold,
                        int arena_max);

#define py_statement_rebind_tail NAMESPACE(py_statement_rebind_tail)
__attribute__((visibility("default")))
int py_statement_rebind_tail(statement st, size_t first, size_t count,
                             const int64_t *rhs,
                             const uint8_t static_hash[16]);

#define py_statement_hash_reset NAMESPACE(py_statement_hash_reset)
__attribute__((visibility("default")))
void py_statement_hash_reset(statement st, const uint8_t hash[16]);

#define py_statement_rebind_constraint_raw NAMESPACE(py_statement_rebind_constraint_raw)
__attribute__((visibility("default")))
int py_statement_rebind_constraint_raw(statement st, size_t index,
                                       size_t nvec, const size_t n[],
                                       size_t deg, const int64_t *phi,
                                       const int64_t *b, size_t nprod,
                                       const int64_t *a);

#define py_gen_params NAMESPACE(py_gen_params)
__attribute__((visibility("default")))
int py_gen_params(dch_pack_params pp, const statement st, int zk, int debug);

#define py_simple_verify NAMESPACE(py_simple_verify)
__attribute__((visibility("default")))
int py_simple_verify(const statement st, const witness wt);

#define py_prove NAMESPACE(py_prove)
__attribute__((visibility("default")))
void py_prove(dch_pack_proof pi, const statement ist, 
              const witness iwt, const dch_pack_params pp);

#define py_last_prove_tail_valid NAMESPACE(py_last_prove_tail_valid)
__attribute__((visibility("default")))
int py_last_prove_tail_valid(void);

#define py_verify NAMESPACE(py_verify)
__attribute__((visibility("default")))
int py_verify(const statement ist, const dch_pack_params pp, 
              const dch_pack_proof pi);

#define py_proof_nbytes NAMESPACE(py_proof_nbytes)
__attribute__((visibility("default")))
size_t py_proof_nbytes(const dch_pack_params pp, const dch_pack_proof pi);

#define py_params_nbytes NAMESPACE(py_params_nbytes)
__attribute__((visibility("default")))
size_t py_params_nbytes(const dch_pack_params pp);

#define py_params_serialize NAMESPACE(py_params_serialize)
__attribute__((visibility("default")))
int py_params_serialize(uint8_t *out, size_t outlen, const dch_pack_params pp);

#define py_params_deserialize NAMESPACE(py_params_deserialize)
__attribute__((visibility("default")))
int py_params_deserialize(dch_pack_params pp, const uint8_t *in, size_t inlen);

#define py_params_is_zk NAMESPACE(py_params_is_zk)
__attribute__((visibility("default")))
int py_params_is_zk(const dch_pack_params pp);

#define py_proof_serialize NAMESPACE(py_proof_serialize)
__attribute__((visibility("default")))
int py_proof_serialize(uint8_t *out, size_t outlen, const dch_pack_params pp,
                       const dch_pack_proof pi);

#define py_proof_deserialize NAMESPACE(py_proof_deserialize)
__attribute__((visibility("default")))
int py_proof_deserialize(dch_pack_proof pi, const dch_pack_params pp,
                         const uint8_t *in, size_t inlen);

#define py_free_witness NAMESPACE(py_free_witness)
__attribute__((visibility("default")))
void py_free_witness(witness wt);

#define py_free_statement NAMESPACE(py_free_statement)
__attribute__((visibility("default")))
void py_free_statement(statement st);

#define py_free_params NAMESPACE(py_free_params)
__attribute__((visibility("default")))
void py_free_params(dch_pack_params pp);

#define py_free_proof NAMESPACE(py_free_proof)
__attribute__((visibility("default")))
void py_free_proof(dch_pack_proof pi);

#endif
