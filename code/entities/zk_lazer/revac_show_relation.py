#!/usr/bin/env python3
"""Relation and size analysis for a revocable anonymous show proof.

This file intentionally does not import LaZer.  It captures the proof language
we need before generating parameters or writing C bindings.
"""

from dataclasses import dataclass

from revac_profile_generated import (
    D as PARAM_D,
    GADGET_BASE as ACC_GADGET_BASE,
    GADGET_DIM as ACC_GADGET_DIM,
    HANDLE_BITS as ACC_HANDLE_BITS,
    N as PARAM_N,
    PARAMETERS,
    Q_BITLEN as PARAM_Q_BITLEN,
    TRUNC_K_L as ACC_TRUNC_K_L,
)


SHOW_PROFILE = PARAMETERS["show_proof_compiled"]
PARAM_N_SHOW = SHOW_PROFILE["n"]
PARAM_K_SHOW = PARAM_N // PARAM_N_SHOW
POLYQSHOW_PACKED_BYTES = PARAM_N_SHOW * 8

PARAM_M1_SHOW_CURRENT = SHOW_PROFILE["witness_polynomials"]

@dataclass(frozen=True)
class RelationCost:
    name: str
    native_polys: int
    qshow_polys: int

    @property
    def packed_bytes(self) -> int:
        return self.qshow_polys * POLYQSHOW_PACKED_BYTES

    @property
    def packed_kib(self) -> float:
        return self.packed_bytes / 1024.0


def q_polys_for_vec_d(vec_count: int) -> int:
    return vec_count * PARAM_D


def qshow_polys_for_q_polys(poly_count: int) -> int:
    return poly_count * PARAM_K_SHOW


def direct_accumulator_witness_cost() -> RelationCost:
    """Cost of embedding the current full acc_genpp_wit_t directly.

    MemVer currently depends on:
      left, v12                       -> 2 vec_d
      v2[ACC_TRUNC_K_L]               -> 3 vec_d
      tag                             -> 1 scalar poly
      right[32][12]                   -> 384 vec_d

    This is a direct and intentionally pessimistic embedding.  It explains why
    a new proof language should avoid putting the full right side into the show
    proof whenever possible.
    """

    vec_d_count = 2 + ACC_TRUNC_K_L + ACC_HANDLE_BITS * ACC_GADGET_DIM
    native_polys = q_polys_for_vec_d(vec_d_count) + 1
    qshow_polys = qshow_polys_for_q_polys(native_polys)
    return RelationCost("direct full accumulator witness", native_polys, qshow_polys)


def compact_accumulator_witness_cost() -> RelationCost:
    """Cost of only the compact part: left, v12, v2[], tag.

    This is not enough for the current MemVer equation by itself, but it gives a
    target scale for a redesigned accumulator-membership proof.
    """

    vec_d_count = 2 + ACC_TRUNC_K_L
    native_polys = q_polys_for_vec_d(vec_d_count) + 1
    qshow_polys = qshow_polys_for_q_polys(native_polys)
    return RelationCost("compact accumulator witness target", native_polys, qshow_polys)


def handle_cost() -> RelationCost:
    # One degree-64 proof polynomial can hold 32 binary handle bits with room.
    return RelationCost("hidden handle x bits", 1, 1)


def print_cost(cost: RelationCost) -> None:
    print(f"- {cost.name}:")
    print(f"  native q-polys : {cost.native_polys}")
    print(f"  qshow polys    : {cost.qshow_polys}")
    print(f"  packed payload : {cost.packed_kib:.2f} KiB")


def main() -> None:
    direct = direct_accumulator_witness_cost()
    compact = compact_accumulator_witness_cost()
    handle = handle_cost()
    current_show_bytes = PARAM_M1_SHOW_CURRENT * POLYQSHOW_PACKED_BYTES

    print("Revocable anonymous show proof relation")
    print("=======================================")
    print()
    print("Current show witness:")
    print(f"- PARAM_M1_SHOW current : {PARAM_M1_SHOW_CURRENT} qshow polys")
    print(f"- z1 packed payload     : {current_show_bytes / 1024.0:.2f} KiB")
    print()
    print("Additional witness if we extend the show proof:")
    print_cost(handle)
    print_cost(compact)
    print_cost(direct)
    print()
    print("Required proof constraints:")
    print("1. Signature binding:")
    print("   [q1.I | A | -B | A3 | -Ds | -D | -Dx] * witness = u")
    print("   This extends the current show relation by hidden x.")
    print("2. Same-x binding:")
    print("   The x bits used in Dx*x are exactly the x bits used by accumulator MemVer.")
    print("3. Accumulator membership:")
    print("   MemVer(pp_acc, acc, x, s_x) = 1 and ||s_x|| <= beta.")
    print("4. Session binding:")
    print("   GCS nonce/session data must be included in the Fiat-Shamir transcript.")
    print()
    print("Important conclusion:")
    print("A direct embedding of the current full accumulator witness is too large.")
    print("Use LaZer/Labrador to prototype the relation, but redesign the accumulator")
    print("proof language before trying to merge it into the existing C show proof.")


if __name__ == "__main__":
    main()
