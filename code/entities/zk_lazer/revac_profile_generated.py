"""Generated revocable-credential profile. Do not edit."""

from __future__ import annotations

PROFILE_NAME = 'revac128-candidate-v1'
PROFILE_ID = bytes.fromhex('ef5a80fcf061b3a7c0307d8d046239bac641143e7cef3e6b07c468bf12f2f14a')
PROFILE_ID_HEX = PROFILE_ID.hex()
RELEASE_READY = False
TARGET_COMPUTATIONAL_BITS = 128
TARGET_STATISTICAL_BITS = 128
PARAMETERS = {'accumulator': {'authenticated_random_update_signs': False,
                 'gadget_base': 3,
                 'gadget_dim': 12,
                 'gadget_signed': True,
                 'handle_bits': 32,
                 'max_deletions_storage': 128,
                 'max_dictionary_entries': 128,
                 'randomized_witness_refresh': False,
                 'truncated_sampler_base': 14,
                 'truncated_sampler_k': 5,
                 'truncated_sampler_l': 2,
                 'witness_beta_sq_policy': 9223372036854775807},
 'issue_proof_compiled': {'challenge_eta': 93,
                          'challenge_rho': 8,
                          'module_rank': 20,
                          'n': 64,
                          'q': 212704515569,
                          'q_bitlen': 38,
                          'randomness_polynomials': 58,
                          'soundness_repetitions': 7,
                          'witness_polynomials': 104},
 'labrador_memver': {'linear_decomposition_base': 256,
                     'linear_decomposition_chunks': 4,
                     'proof_modulus': 274877906837,
                     'proof_modulus_bits': 38,
                     'prover_zero_knowledge': True,
                     'right_fresh_group_norm_sq': 8192,
                     'right_group_norm_sq_max': 1073741824,
                     'right_group_norm_sq_per_deletion': 16777216,
                     'right_layout': 'grouped',
                     'small_chunk_norm_sq': 4194304,
                     'tag_norm_sq': 1,
                     'v2_norm_sq': 268435456},
 'ring': {'d': 4, 'n': 256, 'q': 405769, 'q_bitlen': 19},
 'same_handle_commitment': {'matrix_a_domain_separator': 29,
                            'matrix_d_domain_separator': 30,
                            'opening_dimension': 4,
                            'opening_distribution': 'binary',
                            'opening_poly_norm_sq': 256,
                            'opening_vectors': 2,
                            'randomness_domain_separator': 31},
 'show_proof_compiled': {'challenge_eta': 93,
                         'challenge_rho': 8,
                         'compression_gamma': 1073741824,
                         'module_rank': 23,
                         'n': 64,
                         'ordinary_rejection_sampling': True,
                         'q': 216093359578546177,
                         'q_bitlen': 58,
                         'randomness_polynomials': 74,
                         'same_handle_link_randomness_polynomials': 32,
                         'soundness_repetitions': 4,
                         'witness_polynomials': 183},
 'signature': {'base_high': 13,
               'base_low': 49,
               'dim_high': 2,
               'dim_low': 2,
               'message_slots': 10,
               'q_high': 169,
               'q_low': 2401,
               'tag_weight': 5},
 'system': {'bml_capacity': 20,
            'context_wire_version': 3,
            'labrador_compact_encoding_version': 3,
            'labrador_raw_encoding_version': 2,
            'labrador_worker_protocol_version': 3,
            'parameter_hash_bytes': 32,
            'private_delta_layers': 6,
            'public_delta_layers': 4,
            'seed_bytes': 32,
            'show_nonce_bytes': 32,
            'show_packet_encoding_version': 2,
            'sync_nonce_bytes': 32,
            'sync_wire_version': 2}}

RING = PARAMETERS["ring"]
ACCUMULATOR = PARAMETERS["accumulator"]
LABRADOR_MEMVER = PARAMETERS["labrador_memver"]
SYSTEM = PARAMETERS["system"]
SAME_HANDLE_COMMITMENT = PARAMETERS["same_handle_commitment"]
N, D, Q, Q_BITLEN = RING["n"], RING["d"], RING["q"], RING["q_bitlen"]
HANDLE_BITS = ACCUMULATOR["handle_bits"]
GADGET_BASE = ACCUMULATOR["gadget_base"]
GADGET_DIM = ACCUMULATOR["gadget_dim"]
GADGET_SIGNED = ACCUMULATOR["gadget_signed"]
RANDOM_UPDATE_SIGNS = ACCUMULATOR["authenticated_random_update_signs"]
RANDOMIZED_WITNESS_REFRESH = ACCUMULATOR["randomized_witness_refresh"]
TRUNC_K = ACCUMULATOR["truncated_sampler_k"]
TRUNC_L = ACCUMULATOR["truncated_sampler_l"]
TRUNC_K_L = TRUNC_K - TRUNC_L
PUBLIC_DELTA_LAYERS = SYSTEM["public_delta_layers"]
PRIVATE_DELTA_LAYERS = SYSTEM["private_delta_layers"]
BML_CAPACITY = SYSTEM["bml_capacity"]
SEED_BYTES = SYSTEM["seed_bytes"]
HASH_BYTES = SYSTEM["parameter_hash_bytes"]
PROOF_MODULUS = LABRADOR_MEMVER["proof_modulus"]
PROOF_MODULUS_BITS = LABRADOR_MEMVER["proof_modulus_bits"]
LINEAR_DECOMP_BASE = LABRADOR_MEMVER["linear_decomposition_base"]
LINEAR_DECOMP_CHUNKS = LABRADOR_MEMVER["linear_decomposition_chunks"]
SMALL_CHUNK_NORM_SQ = LABRADOR_MEMVER["small_chunk_norm_sq"]
TAG_NORM_SQ = LABRADOR_MEMVER["tag_norm_sq"]
V2_NORM_SQ = LABRADOR_MEMVER["v2_norm_sq"]
RIGHT_FRESH_GROUP_NORM_SQ = LABRADOR_MEMVER["right_fresh_group_norm_sq"]
RIGHT_GROUP_NORM_SQ_PER_DELETION = LABRADOR_MEMVER["right_group_norm_sq_per_deletion"]
RIGHT_GROUP_NORM_SQ_MAX = LABRADOR_MEMVER["right_group_norm_sq_max"]
LINK_OPENING_POLY_NORM_SQ = SAME_HANDLE_COMMITMENT["opening_poly_norm_sq"]
CONTEXT_WIRE_VERSION = SYSTEM["context_wire_version"]
SYNC_WIRE_VERSION = SYSTEM["sync_wire_version"]
SHOW_PACKET_ENCODING_VERSION = SYSTEM["show_packet_encoding_version"]
LABRADOR_RAW_ENCODING_VERSION = SYSTEM["labrador_raw_encoding_version"]
LABRADOR_COMPACT_ENCODING_VERSION = SYSTEM["labrador_compact_encoding_version"]
LABRADOR_WORKER_PROTOCOL_VERSION = SYSTEM["labrador_worker_protocol_version"]
