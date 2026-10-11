#!/usr/bin/env python3
"""Compile the q88 accumulator MemVer relation into degree-256 LaBRADOR.

It proves all three CRT-lane accumulator relations with one common bounded
witness and one common 24-bit handle.  In production it also proves that this
handle opens the q19 commitment embedded in the anonymous-signature show, so
the two proof components cannot be spliced across different credentials.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import math
import os
import time
from pathlib import Path
from typing import Any

try:
    import numpy as _np
except ImportError:  # Keep the standalone compiler usable without NumPy.
    _np = None

from revac_acc96_wire import (
    GADGET_BASE, GADGET_DIGITS, HANDLE_BITS, LANES, LIMB_RADIX, MODULI,
    InstancePayload, LIFT_LIMBS, LINK_DECOMP_BASE, LINK_DECOMP_LIMBS,
    LinkPayload, PublicLinkPayload, ZeroProverPayload, N,
    NON_TAG_POLYS, PHASE_N, PHASES, ProverPayload, PublicParameters, Q19,
    REFRESH_BITS, REFRESH_BLOCKS, REFRESH_Z_FIRST_SLOT,
    TAG_SLOT, TRUNC_BASE, TRUNC_RIGHT, WITNESS_LIMBS, WITNESS_POLYS,
    _center_lane,
)
from revac_profile_generated import PROOF_MODULUS


VECTOR_COEFF_BOUND = 4095
AGGREGATE_LIMBS = 4
REFRESH_AHAT_LIMBS = 3
try:
    _aggregate_chunk_bits = int(os.environ.get(
        "REVAC_Q96_AGGREGATE_CHUNK_BITS", "0"))
except ValueError as exc:
    raise ValueError(
        "REVAC_Q96_AGGREGATE_CHUNK_BITS must be an integer") from exc
if _aggregate_chunk_bits == 0:
    AGGREGATE_RANKS = (16, 8)
elif (_aggregate_chunk_bits > 0 and
      _aggregate_chunk_bits & (_aggregate_chunk_bits - 1) == 0 and
      HANDLE_BITS % _aggregate_chunk_bits == 0):
    AGGREGATE_RANKS = ((_aggregate_chunk_bits,) *
                       (HANDLE_BITS // _aggregate_chunk_bits))
else:
    raise ValueError(
        "REVAC_Q96_AGGREGATE_CHUNK_BITS must be 0 or a power-of-two "
        "divisor of HANDLE_BITS")
if sum(AGGREGATE_RANKS) != HANDLE_BITS:
    raise ValueError("q88 aggregate extension ranks do not cover the handle")
AGGREGATE_SLICES = tuple(
    (sum(AGGREGATE_RANKS[:index]), sum(AGGREGATE_RANKS[:index + 1]))
    for index in range(len(AGGREGATE_RANKS)))


def _aggregated_quadratics_enabled() -> bool:
    """Select the exact padded extension-ring quadratic aggregation layout.

    The worker profile identifier binds this choice, so proofs from the two
    compilers cannot be mixed.  LaBRADOR extension ranks must be powers of two,
    so q88's 24 handle coordinates are represented by rank-16 and rank-8
    elements.  The direct layout remains a compatibility fallback.
    """
    return os.environ.get("REVAC_Q96_AGGREGATED_QUAD", "1") != "0"


def _full_ring_enabled() -> bool:
    """Use the degree-4096 q88 compiler; zero selects the phased fallback."""
    return os.environ.get("REVAC_Q96_FULL_RING", "1") != "0"


def _aggregation_constraint_count(accumulator_constraints: int) -> int:
    return (accumulator_constraints * len(AGGREGATE_RANKS)
            if _aggregated_quadratics_enabled() else 0)


def _refresh_constraint_count() -> int:
    return LANES * PHASES * REFRESH_BLOCKS

# LaBRADOR stores exact norm bounds in uint64_t, while Dachshund's exact-norm
# path only accepts a much smaller per-vector bound.  The q96 wire digits are
# consequently proved with LaBRADOR's JL approximate-norm path.  This is a
# fixed statement parameter, not a bound selected from the current witness.
#
# Every non-binary backend coefficient is first bounded by 4095.  The public
# JL requirement is derived separately for each logical vector from its exact
# dimension.  We retain an explicit 2^47 squared-norm floor for robust
# parameter generation across the 38-bit deployment backends; applying that
# floor uniformly is still much tighter than the former 2^49 setting.
# Conservatively adding these fixed requirements over the accumulator witness
# slots, and then recomposing the five radix-4095 witness limbs, implies
#
#   ||w||_2 < 2^80,
#
# which is the fail-closed public MemVer gate and lies below the theorem's
# beta_0 (log2(beta_0)=80.097).  beta_1 is the larger extracted M-SIS bound,
# not an admissible MemVer witness bound.  Lift limbs are quotient witnesses
# and are not part of the accumulator short-witness beta claim.
JL_INF_MULT = 9.75
JL_INF_SLACK = JL_INF_MULT / 0.74
ACC96_PUBLIC_WITNESS_NORM_LOG2 = 80.0
# Keeping the floor explicit makes the generated statement deterministic while
# the resulting five-limb recomposition remains strictly below the 2^80
# MemVer gate with the radix-4095 witness limbs.
PACK_APPROX_NORM_SQ_FLOOR = 1 << 47


def _pack_approx_norm_sq_floor() -> int:
    """Return the optional engineering floor above the derived JL bound."""
    value = os.environ.get("REVAC_Q96_PACK_NORM_FLOOR_BITS")
    if value is None or value == "":
        return PACK_APPROX_NORM_SQ_FLOOR
    try:
        bits = int(value, 10)
    except ValueError as exc:
        raise ValueError(
            "REVAC_Q96_PACK_NORM_FLOOR_BITS must be an integer") from exc
    if not 0 <= bits <= 47:
        raise ValueError(
            "REVAC_Q96_PACK_NORM_FLOOR_BITS must be in [0, 47]")
    return 0 if bits == 0 else 1 << bits


class Q96StatementCache:
    """One resident native q96 statement and its generated Pack parameters.

    Only witness slots and the 52 dynamic right-hand-side polynomials vary.
    The 288 A_r selection constraints and 96 rank-16/rank-8 aggregation
    constraints form an immutable prefix.  The 48 accumulator RHS values
    change only with the accumulator and four randomized same-handle
    commitments change for every show.  The native statement owns the
    immutable coefficient graph; keeping it alive avoids rebuilding millions
    of Python/CFFI objects.
    """

    def __init__(self) -> None:
        self.ps: Any | None = None
        self.static_hash: bytes | None = None
        self.constraint_count = 0
        self.dimensions: tuple[int, ...] = ()
        self.params: Any | None = None
        self.params_owned = False
        self._params_bytes: bytes | None = None
        self.profile_id: bytes | None = None
        # Transcript state after the static A_r selection constraints and
        # optional rank-16/rank-8 aggregation constraints.  Dynamic
        # accumulator RHS rebinding starts here.
        self.aggregation_prefix_hash: bytes | None = None
        # Transcript state immediately after the 48 accumulator constraints.
        # It is keyed by the public accumulator value.  The four randomized
        # same-handle commitment constraints can then be rebound without
        # hashing the multi-megabyte accumulator coefficient graph again.
        self.accumulator_prefix_hash: bytes | None = None
        self.accumulator_key: tuple[int, ...] | None = None
        # Logical witness slots for the small per-show q19 same-handle
        # opening.  On an exact private-payload cache hit the 19k-polynomial
        # q96 witness is unchanged; only these slots need refreshing.
        self.link_witness_indices: tuple[int, ...] | None = None
        # Canonical lane/phase expansion of the 129 x 6 seeded Remark-4.3
        # matrix bank.  Shape: matrix, lane, phase, block, coefficient.
        self.refresh_bank: Any | None = None

    def checked_params(self, require_zk: bool = True) -> Any:
        from lazer import lib

        if self.params is None:
            raise ValueError("q96 resident Pack parameters are unavailable")
        if require_zk and lib.labrador38_py_params_is_zk(self.params) != 1:
            raise ValueError("q96 resident Pack parameters are not zk=True")
        return self.params

    def serialized_bytes(self, require_zk: bool = True) -> bytes:
        from revac_acc_memver_import_labrador import (
            serialize_labrador38_params,
        )

        params = self.checked_params(require_zk=require_zk)
        if self._params_bytes is None:
            self._params_bytes = serialize_labrador38_params(
                params, require_zk=require_zk)
        return self._params_bytes

    def bind_profile_id(self) -> bytes:
        if self.profile_id is None:
            # Standalone CLI fallback.  The production worker installs a
            # stronger deployment profile ID derived from the encoded q96 PP
            # and q19 link-matrix seed before constructing this cache.
            self.profile_id = hashlib.shake_256(
                b"REVAC-Q88-LABRADOR38-PACK-PROFILE-V2-REFRESH128\x00" +
                self.serialized_bytes(require_zk=True)
            ).digest(32)
        return self.profile_id

    def __del__(self) -> None:
        # pack_prove() leaves ownership with proof_statement.  The verifier's
        # independently generated parameter tree needs explicit destruction.
        if self.params_owned and self.params is not None:
            try:
                from lazer import lib
                lib.labrador38_py_free_params(self.params)
            except Exception:
                pass
            self.params = None


def _q96_timing(label: str, started: float) -> None:
    value = os.environ.get("REVAC_LABRADOR_Q96_TIMING", "1")
    if value and value[0] != "0":
        print(f"[q96 Labrador timing] {label:<34} "
              f"{time.perf_counter() - started:10.6f} s", flush=True)


def _dynamic_rhs(prover: Any, link: Any | None,
                 constraint_filter: tuple[int, int] | None = None,
                 ) -> list[int]:
    rhs: list[int] = []
    for lane in range(LANES):
        for output_phase in range(PHASES):
            if (constraint_filter is not None and
                    (lane, output_phase) != constraint_filter):
                continue
            rhs.extend(
                _proof_center(_center_lane(
                    prover.accumulator[PHASES * coefficient + output_phase],
                    lane))
                for coefficient in range(PHASE_N)
            )
    if link is not None:
        for row in range(4):
            rhs.extend(
                _proof_center(value - Q19 if value > Q19 // 2 else value)
                for value in link.commitment[row]
            )
    return rhs


def _proof_center(value: int) -> int:
    value %= PROOF_MODULUS
    return value - PROOF_MODULUS if value > PROOF_MODULUS // 2 else value


def _proof_center_scaled_array(values: Any, scale: int = 1) -> Any:
    """Vectorised exact ``scale * values (mod PROOF_MODULUS)``.

    A direct int64 product can be about 70 bits for q96 lane coefficients.
    Split the centered 38-bit scale into two 19-bit limbs, then use
    ``2^38 == 107 (mod PROOF_MODULUS)`` to fold the high product.  The folded
    value stays below 2^63 and needs only one vectorized remainder operation,
    so this is exact integer arithmetic rather than a floating approximation.
    """
    assert _np is not None
    radix_bits = 19
    radix = 1 << radix_bits
    radix_mask = radix - 1
    scale_centered = _proof_center(scale)
    scale_lo = scale_centered & radix_mask
    scale_hi = scale_centered >> radix_bits
    source = _np.asarray(values, dtype=_np.int64)
    low = source * scale_lo
    high = source * scale_hi
    folded = (low + (high & radix_mask) * radix +
              (high >> radix_bits) * ((1 << 38) - PROOF_MODULUS))
    reduced = folded % PROOF_MODULUS
    return _np.where(reduced > PROOF_MODULUS // 2,
                     reduced - PROOF_MODULUS, reduced).astype(
                         _np.int64, copy=False)


def _phase(poly: Any, lane: int,
           output_phase: int, witness_phase: int,
           scale: int = 1) -> Any:
    public_phase = (output_phase - witness_phase) % PHASES
    if _np is not None:
        source = _np.asarray(poly)
        if source.dtype == _np.int64:
            # The q96 compiler pre-centres each public polynomial once per
            # CRT lane.  Reusing that int64 view avoids repeating an 88-bit
            # Python modulo for every output/witness phase and limb.
            values = source[public_phase::PHASES]
        else:
            modulus = MODULI[lane]
            raw = _np.asarray(poly, dtype=object)[public_phase::PHASES]
            canonical = _np.fromiter(
                (int(value) % modulus for value in raw),
                dtype=_np.int64, count=PHASE_N)
            values = _np.where(canonical > modulus // 2,
                               canonical - modulus, canonical)
        if public_phase + witness_phase >= PHASES:
            values = _np.concatenate((-values[-1:], values[:-1]))
        return _proof_center_scaled_array(values, scale)
    values = [_center_lane(poly[PHASES * coefficient + public_phase], lane)
              for coefficient in range(PHASE_N)]
    if public_phase + witness_phase >= PHASES:
        values = [-values[-1]] + values[:-1]
    return [_proof_center(scale * value) for value in values]


def _phase_matrix(centered_polys: Any, output_phase: int,
                  scale: int = 1) -> Any:
    """Vectorize all witness phases for one or more centered polynomials.

    The result is indexed as ``[witness_phase, polynomial, coefficient]``.
    It is byte-for-byte equivalent to calling ``_phase`` separately for each
    polynomial and witness phase, but removes tens of thousands of small
    NumPy allocations from cold q88 statement construction.
    """
    assert _np is not None
    source = _np.asarray(centered_polys, dtype=_np.int64)
    squeeze = source.ndim == 1
    if squeeze:
        source = source.reshape(1, -1)
    if source.ndim != 2 or source.shape[1] != N:
        raise ValueError("q96 public phase matrix shape mismatch")
    phased = source.reshape(source.shape[0], PHASE_N, PHASES)
    witness_phases = _np.arange(PHASES)
    public_phases = (output_phase - witness_phases) % PHASES
    values = _np.transpose(
        phased[:, :, public_phases], (2, 0, 1)).copy()
    wrapped = witness_phases > output_phase
    if _np.any(wrapped):
        wrapped_values = values[wrapped]
        values[wrapped] = _np.concatenate(
            (-wrapped_values[:, :, -1:], wrapped_values[:, :, :-1]),
            axis=2)
    result = _proof_center_scaled_array(values, scale)
    return result[:, 0, :] if squeeze else result


def _phase_matrix_job(args: tuple[Any, int, int]) -> Any:
    return _phase_matrix(*args)


def _selector(position: int, scale: int = 1) -> Any:
    if _np is not None:
        out = _np.zeros(PHASE_N, dtype=_np.int64)
        out[0] = _proof_center(scale)
        return out
    out = [0] * PHASE_N
    out[0] = _proof_center(scale)
    return out


def _vector_bound(dimension: int) -> int:
    return dimension * PHASE_N * VECTOR_COEFF_BOUND * VECTOR_COEFF_BOUND


def _balanced_radix_limbs(values: Any, limbs: int,
                          radix: int = LIMB_RADIX) -> list[Any]:
    """Exact centered odd-radix decomposition over the proof modulus."""
    if radix % 2 != 1:
        raise ValueError("q96 aggregate decomposition needs an odd radix")
    half = radix // 2
    if _np is not None:
        remaining = _np.asarray(values, dtype=_np.int64).copy()
        out = []
        for _ in range(limbs):
            digit = (remaining + half) % radix - half
            out.append(digit)
            remaining = (remaining - digit) // radix
        if _np.any(remaining != 0):
            raise ValueError("q96 aggregate needs more radix limbs")
        return out
    remaining = [[int(value) for value in poly] for poly in values]
    out = []
    for _ in range(limbs):
        digit_polys = []
        next_remaining = []
        for poly in remaining:
            digits = [((value + half) % radix) - half for value in poly]
            digit_polys.append(digits)
            next_remaining.append([
                (value - digit) // radix
                for value, digit in zip(poly, digits)
            ])
        out.append(digit_polys)
        remaining = next_remaining
    if any(value != 0 for poly in remaining for value in poly):
        raise ValueError("q96 aggregate needs more radix limbs")
    return out


def _refresh_lane_phase_bank(public: PublicParameters) -> Any:
    """Expand the seeded q88 A_i bank once into the three CRT phase views."""
    if _np is None:
        raise ValueError("q88 randomized refresh requires NumPy")
    bank = _np.empty(
        (REFRESH_BITS + 1, LANES, PHASES, REFRESH_BLOCKS, PHASE_N),
        dtype=_np.int64)
    for matrix in range(REFRESH_BITS + 1):
        for block in range(REFRESH_BLOCKS):
            canonical = public.refresh(matrix, block)
            for lane, modulus in enumerate(MODULI):
                values = _np.fromiter(
                    (int(value) % modulus for value in canonical),
                    dtype=_np.int64, count=N)
                bank[matrix, lane, :, block, :] = _np.transpose(
                    values.reshape(PHASE_N, PHASES))
    return bank


def _next_power_of_two(value: int) -> int:
    if value <= 0:
        raise ValueError("norm requirement must be positive")
    return 1 << (value - 1).bit_length()


def _approx_norm_requirement(dimension: int) -> int:
    deterministic = _vector_bound(dimension)
    derived = _next_power_of_two(math.ceil(
        deterministic * JL_INF_SLACK * JL_INF_SLACK))
    return max(derived, _pack_approx_norm_sq_floor())


def _audit_norm_profile(
        dimensions: list[int], approximate_bounds: dict[int, int],
        accumulator_witness_indices: set[int],
        *, limb_radix: int = LIMB_RADIX,
        limb_count: int = WITNESS_LIMBS,
        public_norm_sq_bounds: dict[int, int] | None = None) -> None:
    if not accumulator_witness_indices <= approximate_bounds.keys():
        raise ValueError("q96 accumulator witness norm accounting is incomplete")
    for index, requirement in approximate_bounds.items():
        deterministic_bound = (
            public_norm_sq_bounds[index]
            if public_norm_sq_bounds is not None and index in public_norm_sq_bounds
            else _vector_bound(dimensions[index]))
        minimum = math.ceil(
            deterministic_bound * JL_INF_SLACK * JL_INF_SLACK)
        if requirement < minimum or requirement & (requirement - 1):
            raise ValueError("invalid q96 per-vector JL norm requirement")

    # For each accumulator polynomial, Cauchy--Schwarz across the radix
    # limbs gives ||sum_j B^j w_j||^2 <= sum_j B^(2j) sum_j ||w_j||^2.
    # The grouped vectors partition those witness coefficients, so summing
    # only their public JL requirements proves a beta_0 implication.  Exact
    # Lifts and rank-16/rank-8 aggregation auxiliaries are intentionally
    # excluded:
    # they are constrained by the relation but are not part of MemVer's w.
    aggregate_required_norm_sq = sum(
        approximate_bounds[index] for index in accumulator_witness_indices)
    recomposed_log2 = (0.5 * math.log2(aggregate_required_norm_sq) +
                       0.5 * math.log2(sum(
                           limb_radix ** (2 * limb)
                           for limb in range(limb_count))))
    if recomposed_log2 >= ACC96_PUBLIC_WITNESS_NORM_LOG2:
        raise ValueError(
            "q96 Labrador limb norm no longer implies the accumulator "
            f"public witness bound ({recomposed_log2:.3f} >= "
            f"{ACC96_PUBLIC_WITNESS_NORM_LOG2:.3f})")
    print(
        "q96 Labrador 固定范数配置："
        f"per-vector=2^{min(v.bit_length() - 1 for v in approximate_bounds.values())}"
        f"..2^{max(v.bit_length() - 1 for v in approximate_bounds.values())}，"
        f"重组见证 log2(norm)<={recomposed_log2:.3f} (<80.0)",
        flush=True)


def _negacyclic_coefficient(left: list[int] | tuple[int, ...],
                            right: list[int] | tuple[int, ...],
                            output: int) -> int:
    value = 0
    for index, coefficient in enumerate(left):
        right_index = output - index
        if right_index < 0:
            value -= coefficient * int(right[right_index + PHASE_N])
        else:
            value += coefficient * int(right[right_index])
    return value


def _q19_monomial(poly: tuple[int, ...], degree: int) -> list[int]:
    """Return the centered q19 polynomial X^degree*poly in R_q19."""
    out = [0] * PHASE_N
    for index, raw in enumerate(poly):
        destination = index + degree
        value = raw - Q19 if raw > Q19 // 2 else raw
        if destination >= PHASE_N:
            destination -= PHASE_N
            value = -value
        out[destination] = _proof_center(value)
    return out


def audit_compiled_relation(public: PublicParameters, prover: ProverPayload,
                            coefficients: tuple[int, ...] = (0, 1, 2, 127,
                                                                 255),
                            lanes: tuple[int, ...] = (0,),
                            output_phases: tuple[int, ...] = (0,)) -> None:
    """Evaluate the exact phase/limb compiler formula without LaBRADOR.

    This distinguishes an algebra/layout error in this compiler from a native
    LaBRADOR API error.  It intentionally follows ``build_statement`` rather
    than the independent degree-4096 auditor in ``revac_acc96_wire``.
    """
    atrap = public.atrap()
    ba = {(bit, digit): public.ba(bit, digit)
          for bit in range(HANDLE_BITS)
          for digit in range(GADGET_DIGITS)}
    for lane in lanes:
        modulus = MODULI[lane]
        for output_phase in output_phases:
            for coefficient in coefficients:
                expression = 0
                for limb in range(WITNESS_LIMBS):
                    expression += (LIMB_RADIX ** limb) * int(
                        prover.witness_phase(0, limb, output_phase)[coefficient])
                    for witness_phase in range(PHASES):
                        expression += _negacyclic_coefficient(
                            _phase(atrap, lane, output_phase, witness_phase,
                                   LIMB_RADIX ** limb),
                            prover.witness_phase(1, limb, witness_phase),
                            coefficient)

                scale = TRUNC_BASE * TRUNC_BASE
                for source in range(TRUNC_RIGHT):
                    for limb in range(WITNESS_LIMBS):
                        limb_scale = LIMB_RADIX ** limb
                        for witness_phase in range(PHASES):
                            expression -= _negacyclic_coefficient(
                                _phase(public.trunc_b[source], lane,
                                       output_phase, witness_phase,
                                       limb_scale),
                                prover.witness_phase(2 + source, limb,
                                                     witness_phase),
                                coefficient)
                        expression += (_center_lane(scale, lane) * limb_scale *
                                       int(prover.witness_phase(
                                           2 + source, limb,
                                           output_phase)[coefficient]))
                    scale *= TRUNC_BASE

                gadget_scale = 1
                for digit in range(GADGET_DIGITS):
                    centered_gadget = _center_lane(gadget_scale, lane)
                    for bit in range(HANDLE_BITS):
                        slot = TAG_SLOT + 1 + bit * GADGET_DIGITS + digit
                        for limb in range(WITNESS_LIMBS):
                            limb_scale = LIMB_RADIX ** limb
                            for witness_phase in range(PHASES):
                                expression += _negacyclic_coefficient(
                                    _phase(ba[(bit, digit)], lane,
                                           output_phase, witness_phase,
                                           limb_scale),
                                    prover.witness_phase(slot, limb,
                                                         witness_phase),
                                    coefficient)
                            expression -= (centered_gadget * limb_scale *
                                           prover.handle_vector[bit] *
                                           int(prover.witness_phase(
                                               slot, limb,
                                               output_phase)[coefficient]))
                    gadget_scale *= GADGET_BASE

                for limb in range(LIFT_LIMBS):
                    expression += (modulus * LIMB_RADIX ** limb *
                                   int(prover.lift_phase(
                                       lane, limb,
                                       output_phase)[coefficient]))
                rhs = _center_lane(
                    prover.accumulator[PHASES * coefficient + output_phase],
                    lane)
                if (expression - rhs) % PROOF_MODULUS:
                    raise ValueError(
                        "q96 compiled phase relation mismatch at "
                        f"lane={lane}, phase={output_phase}, "
                        f"coefficient={coefficient}: "
                        f"residue={(expression - rhs) % PROOF_MODULUS}")


def debug_native_witness_selector(prover: ProverPayload,
                                  use_polyvec: bool = False) -> None:
    """Check the lowest-level phase-vector loading and selector convention."""
    from lazer import ffi, polyring_t, poly_t, polyvec_t
    from labrados import proof_statement

    ring = polyring_t(PHASE_N, PROOF_MODULUS)
    phase = 0
    polys = [prover.witness_phase(0, 0, item)
             for item in range(PHASES)]
    ps = proof_statement(
        [PHASE_N], [PHASES], [_vector_bound(PHASES)],
        num_constraints=1, primesize="38", zk=True)
    flattened = ffi.new("int64_t[]", [int(value) for poly in polys
                                       for value in poly])
    if use_polyvec:
        index = ps.append_witness(polyvec_t(
            ring, PHASES,
            [poly_t(ring, [int(value) for value in poly])
             for poly in polys]))
    else:
        index = ps.append_witness_pointer(flattened, PHASE_N, PHASES)
    coefficients = [[0] * PHASE_N for _ in range(PHASES)]
    coefficients[phase][0] = 1
    ps.append_statement(
        [polyvec_t(ring, PHASES,
                   [poly_t(ring, item) for item in coefficients])],
        [index], poly_t(ring, [int(value) for value in polys[phase]]))
    if int(ps.simple_verify(ps.statement, ps.witness)) != 1:
        raise ValueError("q96 native phase-vector selector check failed")
    print("q96 native phase-vector selector check passed", flush=True)


def build_statement(public: PublicParameters, prover: Any,
                    link: Any | None,
                    prove: bool,
                    constraint_filter: tuple[int, int] | None = None,
                    fold_handle_for_debug: bool = False,
                    check_witness: bool = True,
                    cache: Q96StatementCache | None = None,
                    reuse_accumulator_witness: bool = False,
                    prepare_prover: bool = False,
                    ) -> tuple[Any, Any, Any]:
    from lazer import ffi, lib, polyring_t, poly_t, polyvec_t
    from labrados import proof_statement

    if _full_ring_enabled():
        from revac_acc96_full_ring import build_statement_full_ring
        return build_statement_full_ring(
            public, prover, link, prove, constraint_filter,
            fold_handle_for_debug, check_witness, cache,
            reuse_accumulator_witness, prepare_prover)

    build_started = time.perf_counter()
    compact_hash = os.environ.get(
        "REVAC_LABRADOR_Q96_COMPACT_HASH", "0") != "0"
    if compact_hash and (cache is None or cache.profile_id is None):
        raise ValueError(
            "q96 compact statement hash requires a bound deployment profile")

    # Exact compact-request reuse is keyed by PP identity, accumulator
    # sequence/value and hidden handle in the C bridge and checked again by
    # the worker.  The accumulator witness already resident in the native
    # statement is consequently identical.  Refresh only the randomized q19
    # link opening and enter Pack Prove without rebuilding aggregate limbs or
    # copying the large witness a second time.
    if (cache is not None and cache.ps is not None and prove and
            reuse_accumulator_witness and link is not None and
            constraint_filter is None):
        ps = cache.ps
        if (cache.link_witness_indices is None or
                cache.accumulator_prefix_hash is None or
                cache.accumulator_key != tuple(
                    int(value) for value in prover.accumulator)):
            raise ValueError("q96 exact witness reuse cache mismatch")
        link_polys = [link.randomness[part][row]
                      for part in range(2) for row in range(4)]
        link_polys.extend(
            link.lift_limbs[row][limb]
            for row in range(4) for limb in range(LINK_DECOMP_LIMBS))
        if len(link_polys) != len(cache.link_witness_indices):
            raise ValueError("q96 link witness reuse layout mismatch")
        stage_started = time.perf_counter()
        for logical_index, poly in zip(
                cache.link_witness_indices, link_polys):
            if _np is not None:
                values = _np.ascontiguousarray(poly, dtype=_np.int64)
                pointer = ffi.from_buffer("int64_t[]", values)
            else:
                pointer = ffi.new("int64_t[]", [int(value) for value in poly])
            error = ps.set_witness_vector(
                ps.witness, ps.map[logical_index], 1, 1, pointer)
            if error != 0:
                raise ValueError(
                    f"q96 link-only witness refresh failed: {error}")
        _q96_timing("refresh link-only witness", stage_started)
        stage_started = time.perf_counter()
        link_rhs = []
        for row in range(4):
            link_rhs.extend(
                _proof_center(value - Q19 if value > Q19 // 2 else value)
                for value in link.commitment[row])
        error = lib.labrador38_py_statement_rebind_tail(
            ps.statement, cache.constraint_count - 4, 4,
            ffi.new("int64_t[]", link_rhs),
            ffi.new("uint8_t[16]", cache.accumulator_prefix_hash))
        if error != 0:
            raise ValueError(f"q96 link-only RHS rebind failed: {error}")
        _q96_timing("rebind link-only RHS", stage_started)
        proof = ffi.new("labrador38_dch_pack_proof")
        stage_started = time.perf_counter()
        ps.prove(proof, ps.statement, ps.witness,
                 cache.checked_params(require_zk=True))
        _q96_timing("native Pack Prove (exact witness reuse)", stage_started)
        return ps, cache.params, proof

    # A resident verifier already owns the complete immutable coefficient
    # graph and never loads a private witness.  Rebuilding tens of thousands
    # of Python zero-polynomial views (including the aggregate digits) used
    # 0.45--0.55 s per show and cannot affect verification.  Rebind only the
    # public accumulator/link RHS transcript in this fast path.
    if (cache is not None and cache.ps is not None and not prove and
            not check_witness):
        ps = cache.ps
        started = time.perf_counter()
        accumulator_constraints = (1 if constraint_filter is not None
                                   else LANES * PHASES)
        aggregation_constraints = _aggregation_constraint_count(
            accumulator_constraints)
        refresh_constraints = _refresh_constraint_count()
        static_constraints = aggregation_constraints + refresh_constraints
        constraint_count = accumulator_constraints + static_constraints
        if link is not None:
            constraint_count += 4
        if cache.constraint_count != constraint_count:
            raise ValueError("q96 resident verifier profile mismatch")
        dynamic_rhs = _dynamic_rhs(prover, link, constraint_filter)
        aggregation_rhs_count = 0
        accumulator_rhs_count = accumulator_constraints * PHASE_N
        accumulator_key = tuple(int(value) for value in prover.accumulator)
        if link is not None and constraint_filter is None:
            if (cache.accumulator_prefix_hash is None or
                    cache.accumulator_key != accumulator_key):
                prefix_hash = cache.aggregation_prefix_hash
                if prefix_hash is None:
                    raise ValueError(
                        "q96 verifier accumulator prefix is unavailable")
                error = lib.labrador38_py_statement_rebind_tail(
                    ps.statement, static_constraints,
                    accumulator_constraints,
                    ffi.new("int64_t[]", dynamic_rhs[
                        aggregation_rhs_count:
                        aggregation_rhs_count + accumulator_rhs_count]),
                    ffi.new("uint8_t[16]", prefix_hash))
                if error != 0:
                    raise ValueError(
                        f"q96 verifier accumulator rebind failed: {error}")
                prefix_out = ffi.new("uint8_t[16]")
                lib.labrador38_py_statement_hash(prefix_out, ps.statement)
                cache.accumulator_prefix_hash = bytes(
                    ffi.buffer(prefix_out, 16))
                cache.accumulator_key = accumulator_key
            error = lib.labrador38_py_statement_rebind_tail(
                ps.statement,
                static_constraints + accumulator_constraints, 4,
                ffi.new("int64_t[]", dynamic_rhs[
                    aggregation_rhs_count + accumulator_rhs_count:]),
                ffi.new("uint8_t[16]", cache.accumulator_prefix_hash))
        else:
            dynamic_hash = cache.aggregation_prefix_hash
            if dynamic_hash is None:
                raise ValueError("q96 verifier transcript prefix is unavailable")
            error = lib.labrador38_py_statement_rebind_tail(
                ps.statement, static_constraints,
                constraint_count - static_constraints,
                ffi.new("int64_t[]", dynamic_rhs[aggregation_rhs_count:]),
                ffi.new("uint8_t[16]", dynamic_hash))
        if error != 0:
            raise ValueError(f"q96 verifier RHS rebind failed: {error}")
        _q96_timing("fast verifier RHS rebind", started)
        return ps, None, None

    ring = polyring_t(PHASE_N, PROOF_MODULUS)
    vectors: list[list[Any]] = []
    dimensions: list[int] = []
    norms: list[int] = []
    approximate_bounds: dict[int, int] = {}
    accumulator_witness_indices: set[int] = set()

    def append_native_witness(polys: list[Any]) -> None:
        """Copy one logical vector without scalar Python/CFFI conversion.

        q96 transports its large private limbs as packed arrays.  The
        former scalar comprehension materialized roughly 3.7 million Python
        integers on every stable proof.  NumPy performs the conversion into
        the established int64 C ABI in native code, preserving the witness
        layout and proof transcript.
        """
        if _np is not None:
            flattened_np = _np.asarray(polys, dtype=_np.int64, order="C")
            if flattened_np.ndim != 2 or flattened_np.shape[1] != PHASE_N:
                raise ValueError("q96 native witness vector shape mismatch")
            ps.append_witness_pointer(
                ffi.from_buffer("int64_t[]", flattened_np),
                PHASE_N, len(polys))
            return
        flattened = ffi.new(
            "int64_t[]", [int(value) for poly in polys for value in poly])
        ps.append_witness_pointer(flattened, PHASE_N, len(polys))

    def append_vector(polys: list[Any], binary: bool = False,
                      exact_norm_sq: int | None = None,
                      accumulator_witness: bool = False,
                      approx_requirement: int | None = None) -> int:
        vectors.append(polys)
        dimensions.append(len(polys))
        index = len(vectors) - 1
        if binary:
            norms.append(0)
        elif exact_norm_sq is not None:
            norms.append(exact_norm_sq)
        else:
            norms.append(_vector_bound(len(polys)))
            approximate_bounds[index] = (
                approx_requirement if approx_requirement is not None else
                _approx_norm_requirement(len(polys)))
            if accumulator_witness:
                accumulator_witness_indices.add(index)
        return index

    native_quadratic_prepare_seconds = 0.0
    native_quadratic_append_seconds = 0.0

    def append_native_quadratic(
            linear_terms: list[Any], linear_witnesses: list[int],
            quadratic_coefficients: list[int],
            quadratic_left: list[int], quadratic_right: list[int],
            rhs: list[int],
            quadratic_exponents: list[int] | None = None) -> None:
        """Append one q96 relation without transient poly_t/polyvec_t trees.

        ``proof_statement.append_quadratic_statement`` first constructs one
        Python/CFFI object for every public polynomial and then flattens the
        complete tree back into an int64 array.  A q96 accumulator relation
        contains thousands of such public polynomials, and the same process
        is repeated for all 48 CRT-lane/phase relations during cold start.
        The native ABI ultimately consumes only the flat coefficient array,
        term lengths and mapped witness indices.  Supplying those values
        directly preserves the exact relation and transcript while avoiding
        the redundant object graph and copy.
        """
        nonlocal native_quadratic_prepare_seconds
        nonlocal native_quadratic_append_seconds
        operation_started = time.perf_counter()
        if len(linear_terms) != len(linear_witnesses):
            raise ValueError("q96 native linear term accounting mismatch")
        if not (len(quadratic_coefficients) == len(quadratic_left) ==
                len(quadratic_right)):
            raise ValueError("q96 native quadratic term accounting mismatch")
        if (quadratic_exponents is not None and
                len(quadratic_exponents) != len(quadratic_coefficients)):
            raise ValueError("q96 quadratic monomial accounting mismatch")

        lengths: list[int] = []
        flat_parts: list[Any] = []
        for term in linear_terms:
            if _np is not None:
                array = _np.asarray(term, dtype=_np.int64, order="C")
                if array.ndim == 1:
                    if array.shape[0] != PHASE_N:
                        raise ValueError("q96 native polynomial shape mismatch")
                    array = array.reshape(1, PHASE_N)
                elif array.ndim != 2 or array.shape[1] != PHASE_N:
                    raise ValueError("q96 native polyvector shape mismatch")
                lengths.append(int(array.shape[0]))
                flat_parts.append(array)
            else:
                if not term:
                    raise ValueError("q96 native linear term is empty")
                polynomials = term if isinstance(term[0], (list, tuple)) else [term]
                if any(len(poly) != PHASE_N for poly in polynomials):
                    raise ValueError("q96 native polynomial shape mismatch")
                lengths.append(len(polynomials))
                flat_parts.extend(polynomials)

        mapped_linear = [ps.map[index] for index in linear_witnesses]
        mapped_left = [ps.map[index] for index in quadratic_left]
        mapped_right = [ps.map[index] for index in quadratic_right]
        if _np is not None:
            flat = _np.ascontiguousarray(
                _np.concatenate(flat_parts, axis=0).reshape(-1),
                dtype=_np.int64)
            rhs_array = _np.ascontiguousarray(rhs, dtype=_np.int64)
            linear_ptr = ffi.from_buffer("int64_t[]", flat)
            rhs_ptr = ffi.from_buffer("int64_t[]", rhs_array)
        else:
            flat_values = [int(value) for poly in flat_parts for value in poly]
            linear_ptr = ffi.new("int64_t[]", flat_values)
            rhs_ptr = ffi.new("int64_t[]", rhs)
        append_started = time.perf_counter()
        native_quadratic_prepare_seconds += append_started - operation_started
        arguments = (
            ps.statement, len(linear_terms), len(quadratic_coefficients),
            ffi.new("size_t[]", mapped_linear),
            ffi.new("size_t[]", mapped_left),
            ffi.new("size_t[]", mapped_right),
            ffi.new("size_t[]", lengths), 1,
            ffi.new("int64_t[]", quadratic_coefficients))
        if quadratic_exponents is None:
            error = ps.append_quadratic(
                *arguments, linear_ptr, rhs_ptr)
        else:
            error = lib.labrador38_py_append_quadratic_monomials(
                *arguments, ffi.new("uint32_t[]", quadratic_exponents),
                linear_ptr, rhs_ptr)
        if error != 0:
            raise ValueError(f"q96 native accumulator constraint failed: {error}")
        native_quadratic_append_seconds += (
            time.perf_counter() - append_started)
        ps.cur_constraint_num += 1

    handle_slices = (AGGREGATE_SLICES if _aggregated_quadratics_enabled()
                     else ((0, HANDLE_BITS),))
    handle_polys = [
        [prover.handle_vector[bit]] + [0] * (PHASE_N - 1)
        for bit in range(HANDLE_BITS)]
    x_idx = [append_vector(handle_polys[start:end], binary=True)
             for start, end in handle_slices]
    refresh_r_idx = append_vector([
        [prover.refresh_r_bits[bit]] + [0] * (PHASE_N - 1)
        for bit in range(REFRESH_BITS)
    ], binary=True)
    refresh_z_idx = [append_vector([
        prover.witness_phase(REFRESH_Z_FIRST_SLOT + block, 0, phase)
        for block in range(REFRESH_BLOCKS)
    ], exact_norm_sq=REFRESH_BLOCKS * PHASE_N)
        for phase in range(PHASES)]

    stage_started = time.perf_counter()
    if cache is not None and cache.refresh_bank is not None:
        refresh_bank = cache.refresh_bank
    else:
        refresh_bank = _refresh_lane_phase_bank(public)
        if cache is not None:
            cache.refresh_bank = refresh_bank
    refresh_raw = refresh_bank[0].copy()
    selected_refresh = _np.flatnonzero(
        _np.asarray(prover.refresh_r_bits, dtype=_np.int64))
    if selected_refresh.size:
        refresh_raw += _np.sum(
            refresh_bank[selected_refresh + 1], axis=0, dtype=_np.int64)
    refresh_ahat_idx: list[list[list[int]]] = []
    refresh_lift_idx: list[list[int]] = []
    refresh_lift_bound = (
        REFRESH_BLOCKS * PHASE_N * (REFRESH_BITS + 1) ** 2)
    for lane, modulus in enumerate(MODULI):
        centered = refresh_raw[lane] % modulus
        centered = _np.where(
            centered > modulus // 2, centered - modulus, centered)
        lane_chunks = _balanced_radix_limbs(centered, REFRESH_AHAT_LIMBS)
        lane_ahat = []
        lane_lifts = []
        lane_lift_values = (refresh_raw[lane] - centered) // modulus
        for phase in range(PHASES):
            lane_ahat.append([
                append_vector(lane_chunks[chunk][phase].tolist())
                for chunk in range(REFRESH_AHAT_LIMBS)
            ])
            lane_lifts.append(append_vector(
                lane_lift_values[phase].tolist(),
                exact_norm_sq=refresh_lift_bound))
        refresh_ahat_idx.append(lane_ahat)
        refresh_lift_idx.append(lane_lifts)
    _q96_timing("prepare Remark-4.3 refresh bank", stage_started)

    left_idx = [append_vector([
        prover.witness_phase(0, limb, phase) for phase in range(PHASES)
    ], accumulator_witness=True) for limb in range(WITNESS_LIMBS)]
    v12_idx = [append_vector([
        prover.witness_phase(1, limb, phase) for phase in range(PHASES)
    ], accumulator_witness=True) for limb in range(WITNESS_LIMBS)]
    v2_idx = [append_vector([
        prover.witness_phase(2 + source, limb, phase)
        for source in range(TRUNC_RIGHT) for phase in range(PHASES)
    ], accumulator_witness=True) for limb in range(WITNESS_LIMBS)]
    right_idx = [[[[append_vector([
        prover.witness_phase(
            TAG_SLOT + 1 + bit * GADGET_DIGITS + digit, limb, phase)
        for bit in range(start, end)
    ], accumulator_witness=True,
        approx_requirement=(_pack_approx_norm_sq_floor() //
                            len(handle_slices)))
        for start, end in handle_slices]
                  for phase in range(PHASES)]
                  for limb in range(WITNESS_LIMBS)]
                 for digit in range(GADGET_DIGITS)]
    aggregate_idx: list[list[list[list[int]]]] | None = None
    if _aggregated_quadratics_enabled():
        if LIMB_RADIX ** AGGREGATE_LIMBS <= PROOF_MODULUS:
            raise ValueError("q96 aggregate radix capacity is too small")
        aggregate_aux_requirement = _next_power_of_two(math.ceil(
            LANES * PHASES * AGGREGATE_LIMBS *
            _vector_bound(HANDLE_BITS) * JL_INF_SLACK * JL_INF_SLACK))
        aggregate_idx = []
        for lane in range(LANES):
            lane_indices = []
            for output_phase in range(PHASES):
                if _np is not None:
                    total = _np.zeros(
                        (HANDLE_BITS, PHASE_N), dtype=_np.int64)
                else:
                    total = [[0] * PHASE_N for _ in range(HANDLE_BITS)]
                gadget_scale = 1
                for digit in range(GADGET_DIGITS):
                    centered_gadget = _center_lane(gadget_scale, lane)
                    for limb in range(WITNESS_LIMBS):
                        coefficient = _proof_center(
                            -centered_gadget * LIMB_RADIX ** limb)
                        for chunk, (start, end) in enumerate(handle_slices):
                            source = vectors[
                                right_idx[digit][limb][output_phase][chunk]]
                            if _np is not None:
                                total[start:end] += coefficient * _np.asarray(
                                    source, dtype=_np.int64)
                            else:
                                for local_bit, bit in enumerate(
                                        range(start, end)):
                                    for coeff in range(PHASE_N):
                                        total[bit][coeff] += (
                                            coefficient *
                                            source[local_bit][coeff])
                    gadget_scale *= GADGET_BASE
                if _np is not None:
                    total %= PROOF_MODULUS
                    total[total > PROOF_MODULUS // 2] -= PROOF_MODULUS
                else:
                    total = [[_proof_center(value) for value in poly]
                             for poly in total]
                digit_indices = []
                for digit_values in _balanced_radix_limbs(
                        total, AGGREGATE_LIMBS):
                    digit_indices.append([
                        append_vector(
                            (digit_values[start:end].tolist()
                             if _np is not None else
                             digit_values[start:end]),
                            approx_requirement=aggregate_aux_requirement)
                        for start, end in handle_slices])
                lane_indices.append(digit_indices)
            aggregate_idx.append(lane_indices)
    lift_idx = [[[append_vector([
        prover.lift_phase(lane, limb, phase) for phase in range(PHASES)
    ]) for limb in range(LIFT_LIMBS)] for lane in range(LANES)]][0]

    link_r_idx: list[list[int]] | None = None
    link_lift_idx: list[list[int]] | None = None
    if link is not None:
        link_r_idx = [[append_vector(
            [link.randomness[part][row]], exact_norm_sq=PHASE_N)
            for row in range(4)] for part in range(2)]
        link_lift_idx = [[append_vector(
            [link.lift_limbs[row][limb]],
            exact_norm_sq=PHASE_N * (LINK_DECOMP_BASE // 2) ** 2)
            for limb in range(LINK_DECOMP_LIMBS)] for row in range(4)]

    base_non_tag_polys = NON_TAG_POLYS - REFRESH_BLOCKS
    expected_backend_polys = (HANDLE_BITS + REFRESH_BITS +
                              base_non_tag_polys * WITNESS_LIMBS * PHASES +
                              REFRESH_BLOCKS * PHASES +
                              LANES * PHASES * REFRESH_BLOCKS *
                              (REFRESH_AHAT_LIMBS + 1) +
                              LANES * LIFT_LIMBS * PHASES)
    if aggregate_idx is not None:
        expected_backend_polys += (
            LANES * PHASES * AGGREGATE_LIMBS * HANDLE_BITS)
    if link is not None:
        expected_backend_polys += 2 * 4 + 4 * LINK_DECOMP_LIMBS
    if sum(dimensions) != expected_backend_polys:
        raise ValueError("q96 Labrador witness layout accounting mismatch")
    _audit_norm_profile(dimensions, approximate_bounds,
                        accumulator_witness_indices)
    approx_norm_list = [
        (index, requirement)
        for index, requirement in approximate_bounds.items()
    ]
    accumulator_constraints = (1 if constraint_filter is not None
                               else LANES * PHASES)
    aggregation_constraints = _aggregation_constraint_count(
        accumulator_constraints)
    refresh_constraints = _refresh_constraint_count()
    static_constraints = aggregation_constraints + refresh_constraints
    constraint_count = accumulator_constraints + static_constraints
    if link is not None:
        constraint_count += 4
    cache_hit = cache is not None and cache.ps is not None
    _q96_timing("prepare witness views/aggregate limbs", build_started)
    if cache_hit:
        assert cache is not None and cache.ps is not None
        if (cache.dimensions != tuple(dimensions) or
                cache.constraint_count != constraint_count):
            raise ValueError("q96 resident statement profile mismatch")
        ps = cache.ps
        stage_started = time.perf_counter()
        if prove or prepare_prover:
            ps.cur_witness_num = 0
            ps.witness_polys = 0
            for polys in vectors:
                append_native_witness(polys)
            _q96_timing("load dynamic witness", stage_started)
        stage_started = time.perf_counter()
        assert cache.static_hash is not None
        dynamic_rhs = _dynamic_rhs(prover, link, constraint_filter)
        accumulator_constraints = (1 if constraint_filter is not None
                                   else LANES * PHASES)
        aggregation_constraints = _aggregation_constraint_count(
            accumulator_constraints)
        refresh_constraints = _refresh_constraint_count()
        static_constraints = aggregation_constraints + refresh_constraints
        aggregation_rhs_count = 0
        accumulator_rhs_count = accumulator_constraints * PHASE_N
        accumulator_key = tuple(int(value) for value in prover.accumulator)
        can_split = link is not None and constraint_filter is None
        if can_split:
            if (cache.accumulator_prefix_hash is None or
                    cache.accumulator_key != accumulator_key):
                prefix_hash = cache.aggregation_prefix_hash
                if prefix_hash is None:
                    raise ValueError(
                        "q96 aggregation transcript prefix is unavailable")
                error = lib.labrador38_py_statement_rebind_tail(
                    ps.statement, static_constraints,
                    accumulator_constraints,
                    ffi.new("int64_t[]",
                            dynamic_rhs[
                                aggregation_rhs_count:
                                aggregation_rhs_count + accumulator_rhs_count]),
                    ffi.new("uint8_t[16]", prefix_hash))
                if error != 0:
                    if static_constraints:
                        raise ValueError(
                            "q96 accumulator-prefix rebind extension is "
                            f"unavailable: {error}")
                    # Compatibility with an unpatched upstream LaZer build:
                    # preserve correctness by rebinding all 52 constraints.
                    # The entities/zk_lazer/lazer_rebind_prefix.patch extension
                    # enables caching the new prefix after an accumulator
                    # transition and is used by the production build.
                    error = lib.labrador38_py_statement_rebind_tail(
                        ps.statement, 0, constraint_count,
                        ffi.new("int64_t[]", dynamic_rhs),
                        ffi.new("uint8_t[16]", cache.static_hash))
                    if error != 0:
                        raise ValueError(
                            f"q96 full RHS compatibility rebind failed: {error}")
                    cache.accumulator_prefix_hash = None
                    cache.accumulator_key = None
                    _q96_timing("rebind full RHS (legacy LaZer)", stage_started)
                    error = None
                else:
                    prefix_hash = ffi.new("uint8_t[16]")
                    lib.labrador38_py_statement_hash(prefix_hash, ps.statement)
                    cache.accumulator_prefix_hash = bytes(
                        ffi.buffer(prefix_hash, 16))
                    cache.accumulator_key = accumulator_key
                    _q96_timing("rebind accumulator RHS prefix", stage_started)
                    stage_started = time.perf_counter()
            else:
                error = 0
            if error is not None:
                error = lib.labrador38_py_statement_rebind_tail(
                    ps.statement,
                    static_constraints + accumulator_constraints, 4,
                    ffi.new("int64_t[]", dynamic_rhs[
                        aggregation_rhs_count + accumulator_rhs_count:]),
                    ffi.new("uint8_t[16]", cache.accumulator_prefix_hash))
                if error != 0:
                    raise ValueError(
                        f"q96 same-handle tail RHS rebind failed: {error}")
                _q96_timing("rebind same-handle RHS tail", stage_started)
        else:
            dynamic_first = static_constraints
            dynamic_count = constraint_count - static_constraints
            dynamic_hash = cache.aggregation_prefix_hash
            if dynamic_hash is None:
                raise ValueError(
                    "q96 dynamic transcript prefix is unavailable")
            error = lib.labrador38_py_statement_rebind_tail(
                ps.statement, dynamic_first, dynamic_count,
                ffi.new("int64_t[]", dynamic_rhs[aggregation_rhs_count:]),
                ffi.new("uint8_t[16]", dynamic_hash))
            if error != 0:
                raise ValueError(
                    f"q96 resident statement RHS rebind failed: {error}")
            _q96_timing("rebind dynamic RHS", stage_started)
        stage_started = time.perf_counter()
        if check_witness and int(ps.simple_verify(
                ps.statement, ps.witness)) != 1:
            raise ValueError("q96 cached Labrador relation check failed")
        if check_witness:
            _q96_timing("clear native relation check", stage_started)
        if not prove:
            return ps, None, None
        if cache.params is None:
            raise ValueError("q96 resident prover has no cached Pack parameters")
        proof = ffi.new("labrador38_dch_pack_proof")
        stage_started = time.perf_counter()
        ps.prove(proof, ps.statement, ps.witness, cache.params)
        _q96_timing("native Pack Prove", stage_started)
        return ps, cache.params, proof

    stage_started = time.perf_counter()
    ps = proof_statement(
        [PHASE_N] * len(vectors), dimensions, norms,
        num_constraints=constraint_count, primesize="38",
        num_deg0_constraints=HANDLE_BITS + REFRESH_BITS, zk=True,
        approx_norm_list=approx_norm_list)
    _q96_timing("allocate cold proof_statement", stage_started)
    stage_started = time.perf_counter()
    for polys in vectors:
        append_native_witness(polys)
    _q96_timing("load cold witness", stage_started)
    stage_started = time.perf_counter()
    for chunk, (start, end) in enumerate(handle_slices):
        error = lib.labrador38_py_append_deg0_vector_constraints(
            ps.statement, ps.map[x_idx[chunk]], end - start)
        if error != 0:
            raise ValueError(
                f"q96 compact handle degree-0 constraints failed: {error}")
        ps.cur_deg0_constraint_num += end - start
    error = lib.labrador38_py_append_deg0_vector_constraints(
        ps.statement, ps.map[refresh_r_idx], REFRESH_BITS)
    if error != 0:
        raise ValueError(f"q88 refresh-r degree-0 constraints failed: {error}")
    ps.cur_deg0_constraint_num += REFRESH_BITS
    _q96_timing("append degree-0 handle/refresh constraints", stage_started)

    stage_started = time.perf_counter()
    refresh_bank_i32 = _np.ascontiguousarray(refresh_bank, dtype=_np.int32)
    refresh_chunk_indices = [
        ps.map[refresh_ahat_idx[lane][phase][chunk]]
        for lane in range(LANES)
        for phase in range(PHASES)
        for chunk in range(REFRESH_AHAT_LIMBS)
    ]
    refresh_lift_indices = [
        ps.map[refresh_lift_idx[lane][phase]]
        for lane in range(LANES)
        for phase in range(PHASES)
    ]
    error = lib.labrador38_py_append_refresh_phase_constraints(
        ps.statement, ps.map[refresh_r_idx],
        ffi.new("size_t[]", refresh_chunk_indices),
        ffi.new("size_t[]", refresh_lift_indices),
        ffi.from_buffer("int32_t[]", refresh_bank_i32), REFRESH_BITS,
        REFRESH_BLOCKS, LANES, PHASES, PHASE_N,
        ffi.new("uint32_t[]", MODULI), LIMB_RADIX, REFRESH_AHAT_LIMBS)
    if error != 0:
        raise ValueError(f"q88 sparse A_r selection constraints failed: {error}")
    ps.cur_constraint_num += refresh_constraints
    if cache is not None:
        static_hash = ffi.new("uint8_t[16]")
        lib.labrador38_py_statement_hash(static_hash, ps.statement)
        cache.static_hash = bytes(ffi.buffer(static_hash, 16))
    _q96_timing("append A_r selection constraints", stage_started)

    if aggregate_idx is not None:
        stage_started = time.perf_counter()
        for lane in range(LANES):
            for output_phase in range(PHASES):
                if (constraint_filter is not None and
                        (lane, output_phase) != constraint_filter):
                    continue
                for chunk, (start, end) in enumerate(handle_slices):
                    rank = end - start
                    aggregate_witnesses: list[int] = []
                    aggregate_coefficients: list[int] = []
                    gadget_scale = 1
                    for digit in range(GADGET_DIGITS):
                        centered_gadget = _center_lane(gadget_scale, lane)
                        for limb in range(WITNESS_LIMBS):
                            aggregate_witnesses.append(
                                right_idx[digit][limb][output_phase][chunk])
                            aggregate_coefficients.append(_proof_center(
                                -centered_gadget * LIMB_RADIX ** limb))
                        gadget_scale *= GADGET_BASE
                    for limb in range(AGGREGATE_LIMBS):
                        aggregate_witnesses.append(
                            aggregate_idx[lane][output_phase][limb][chunk])
                        aggregate_coefficients.append(_proof_center(
                            -(LIMB_RADIX ** limb)))

                    # Each power-of-two chunk is one extension-ring element.
                    # A scalar in coordinate zero multiplies every component,
                    # encoding sum(c_i R_i) - sum(B^j A_j) = 0 exactly.
                    zero_rank_rhs = (
                        _np.zeros((rank, PHASE_N), dtype=_np.int64)
                        if _np is not None else [0] * (rank * PHASE_N))
                    if _np is not None:
                        phi = _np.zeros(
                            (len(aggregate_witnesses), rank, PHASE_N),
                            dtype=_np.int64)
                        phi[:, 0, 0] = aggregate_coefficients
                        phi_ptr = ffi.from_buffer("int64_t[]", phi)
                        rhs_ptr = ffi.from_buffer("int64_t[]", zero_rank_rhs)
                    else:
                        phi_values = []
                        for coefficient in aggregate_coefficients:
                            phi_values.extend(
                                [coefficient] +
                                [0] * (rank * PHASE_N - 1))
                        phi_ptr = ffi.new("int64_t[]", phi_values)
                        rhs_ptr = ffi.new("int64_t[]", zero_rank_rhs)
                    mapped = ffi.new("size_t[]", [
                        ps.map[index] for index in aggregate_witnesses])
                    one_extension_element = ffi.new(
                        "size_t[]", [1] * len(aggregate_witnesses))
                    error = ps.append_constraint(
                        ps.statement, len(aggregate_witnesses), mapped,
                        one_extension_element, rank, phi_ptr, rhs_ptr, 1)
                    if error != 0:
                        raise ValueError(
                            f"q88 rank-{rank} aggregate constraint failed: "
                            f"{error}")
                    ps.cur_constraint_num += 1
        if cache is not None:
            aggregation_hash = ffi.new("uint8_t[16]")
            lib.labrador38_py_statement_hash(
                aggregation_hash, ps.statement)
            cache.aggregation_prefix_hash = bytes(
                ffi.buffer(aggregation_hash, 16))
        _q96_timing(
            "append rank-(16+8) aggregate constraints", stage_started)
    elif cache is not None:
        cache.aggregation_prefix_hash = cache.static_hash

    if compact_hash:
        # The profile binds the canonical q88 PP, q19 matrix context and the
        # versioned compiler layout.  Use it as the static Fiat-Shamir prefix
        # so dynamic RHS rebinding need not re-hash the multi-GiB expanded
        # coefficient graph.  Native append/rebind still absorbs every RHS in
        # constraint order using the established canonical polynomial bytes.
        assert cache is not None and cache.profile_id is not None
        semantic_prefix = hashlib.shake_256(
            b"REVAC-Q88-LABRADOR38-STATEMENT-HASH-V2-REFRESH128\x00" +
            cache.profile_id
        ).digest(16)
        lib.labrador38_py_statement_hash_reset(
            ps.statement, ffi.new("uint8_t[16]", semantic_prefix))
        cache.static_hash = semantic_prefix
        cache.aggregation_prefix_hash = semantic_prefix

    stage_started = time.perf_counter()
    atrap = public.atrap()
    ba = {(bit, digit): public.ba(bit, digit)
          for bit in range(HANDLE_BITS)
          for digit in range(GADGET_DIGITS)}
    if _np is not None:
        def center_public(poly: Any, lane_index: int) -> Any:
            modulus = MODULI[lane_index]
            canonical = _np.fromiter(
                (int(value) % modulus for value in poly),
                dtype=_np.int64, count=N)
            return _np.where(canonical > modulus // 2,
                             canonical - modulus, canonical).astype(
                                 _np.int64, copy=False)

        atrap_lanes = [center_public(atrap, lane)
                       for lane in range(LANES)]
        trunc_lanes = [[center_public(poly, lane)
                        for poly in public.trunc_b]
                       for lane in range(LANES)]
        ba_lanes = [{key: center_public(poly, lane)
                     for key, poly in ba.items()}
                    for lane in range(LANES)]
        ba_lane_batches = [[
            _np.stack([ba_lanes[lane][(bit, digit)]
                       for bit in range(HANDLE_BITS)])
            for digit in range(GADGET_DIGITS)]
            for lane in range(LANES)]
        trunc_lane_batches = [
            _np.stack(trunc_lanes[lane]) for lane in range(LANES)]
    else:
        atrap_lanes = [atrap] * LANES
        trunc_lanes = [public.trunc_b] * LANES
        ba_lanes = [ba] * LANES
    one = [1] + [0] * (PHASE_N - 1)
    phase_executor = None
    if _np is not None:
        try:
            python_statement_threads = int(os.environ.get(
                "REVAC_LABRADOR_PYTHON_STATEMENT_THREADS", "1"))
        except ValueError:
            python_statement_threads = 1
        if python_statement_threads > 1:
            phase_executor = concurrent.futures.ThreadPoolExecutor(
                max_workers=min(python_statement_threads, 32),
                thread_name_prefix="q88-statement")
    for lane, modulus in enumerate(MODULI):
        for output_phase in range(PHASES):
            if (constraint_filter is not None and
                    (lane, output_phase) != constraint_filter):
                continue
            stat_lin: list[Any] = []
            wit_lin: list[int] = []
            if _np is not None:
                limb_scales = [LIMB_RADIX ** limb
                               for limb in range(WITNESS_LIMBS)]
                phase_jobs = [
                    (atrap_lanes[lane], output_phase, limb_scale)
                    for limb_scale in limb_scales]
                phase_jobs.extend(
                    (trunc_lane_batches[lane], output_phase, limb_scale)
                    for limb_scale in limb_scales)
                phase_jobs.extend(
                    (ba_lane_batches[lane][digit], output_phase, limb_scale)
                    for digit in range(GADGET_DIGITS)
                    for limb_scale in limb_scales)
                if phase_executor is not None:
                    phase_results = list(phase_executor.map(
                        _phase_matrix_job, phase_jobs))
                else:
                    phase_results = [
                        _phase_matrix_job(job) for job in phase_jobs]
                phase_offset = 0
                atrap_phase_matrices = phase_results[
                    phase_offset:phase_offset + WITNESS_LIMBS]
                phase_offset += WITNESS_LIMBS
                trunc_phase_matrices = phase_results[
                    phase_offset:phase_offset + WITNESS_LIMBS]
                phase_offset += WITNESS_LIMBS
                ba_phase_matrices = []
                for _digit in range(GADGET_DIGITS):
                    ba_phase_matrices.append(phase_results[
                        phase_offset:phase_offset + WITNESS_LIMBS])
                    phase_offset += WITNESS_LIMBS
                if phase_offset != len(phase_results):
                    raise ValueError("q96 phase job accounting mismatch")
            limb_scale = 1
            for limb in range(WITNESS_LIMBS):
                if _np is not None:
                    coeff = _np.zeros((PHASES, PHASE_N), dtype=_np.int64)
                else:
                    coeff = [[0] * PHASE_N for _ in range(PHASES)]
                coeff[output_phase] = _selector(output_phase, limb_scale)
                stat_lin.append(coeff)
                wit_lin.append(left_idx[limb])
                limb_scale *= LIMB_RADIX

            limb_scale = 1
            for limb in range(WITNESS_LIMBS):
                if _np is not None:
                    stat_lin.append(atrap_phase_matrices[limb])
                else:
                    stat_lin.append([
                        _phase(atrap_lanes[lane], lane, output_phase,
                               witness_phase, limb_scale)
                        for witness_phase in range(PHASES)
                    ])
                wit_lin.append(v12_idx[limb])
                limb_scale *= LIMB_RADIX

            trunc_scale = TRUNC_BASE * TRUNC_BASE
            for limb in range(WITNESS_LIMBS):
                limb_scale = LIMB_RADIX ** limb
                coeffs = []
                scale = trunc_scale
                if _np is not None:
                    trunc_phase_matrix = trunc_phase_matrices[limb]
                for source in range(TRUNC_RIGHT):
                    for witness_phase in range(PHASES):
                        phase_coeff = (
                            trunc_phase_matrix[witness_phase, source]
                            if _np is not None else
                            _phase(trunc_lanes[lane][source], lane,
                                   output_phase, witness_phase, limb_scale))
                        coeff = (-phase_coeff if _np is not None else
                                 [-value for value in phase_coeff])
                        if witness_phase == output_phase:
                            coeff[0] = _proof_center(
                                int(coeff[0]) + limb_scale *
                                _center_lane(scale, lane))
                        coeffs.append(coeff if _np is not None else
                                      [_proof_center(value)
                                       for value in coeff])
                    scale *= TRUNC_BASE
                stat_lin.append(coeffs)
                wit_lin.append(v2_idx[limb])

            for digit in range(GADGET_DIGITS):
                for limb in range(WITNESS_LIMBS):
                    limb_scale = LIMB_RADIX ** limb
                    if _np is not None:
                        ba_phase_matrix = ba_phase_matrices[digit][limb]
                    for witness_phase in range(PHASES):
                        if _np is not None:
                            right_coeffs = ba_phase_matrix[witness_phase]
                            if (fold_handle_for_debug and
                                    witness_phase == output_phase):
                                right_coeffs = right_coeffs.copy()
                                for bit in range(HANDLE_BITS):
                                    if prover.handle_vector[bit]:
                                        right_coeffs[bit, 0] = _proof_center(
                                            int(right_coeffs[bit, 0]) -
                                            _center_lane(
                                                GADGET_BASE ** digit, lane) *
                                            limb_scale)
                        else:
                            right_coeffs = []
                            for bit in range(HANDLE_BITS):
                                coeff = _phase(
                                    ba_lanes[lane][(bit, digit)], lane,
                                    output_phase, witness_phase, limb_scale)
                                if (fold_handle_for_debug and
                                        prover.handle_vector[bit] and
                                        witness_phase == output_phase):
                                    coeff[0] = _proof_center(
                                        coeff[0] - _center_lane(
                                            GADGET_BASE ** digit, lane) *
                                        limb_scale)
                                right_coeffs.append(coeff)
                        for chunk, (start, end) in enumerate(handle_slices):
                            stat_lin.append(right_coeffs[start:end])
                            wit_lin.append(
                                right_idx[digit][limb][witness_phase][chunk])

            for limb in range(LIFT_LIMBS):
                if _np is not None:
                    coeff = _np.zeros((PHASES, PHASE_N), dtype=_np.int64)
                else:
                    coeff = [[0] * PHASE_N for _ in range(PHASES)]
                coeff[output_phase] = _selector(
                    output_phase, modulus * LIMB_RADIX ** limb)
                stat_lin.append(coeff)
                wit_lin.append(lift_idx[lane][limb])

            stat_quad = []
            wit_quad1 = []
            wit_quad2 = []
            quad_exponents = []
            if not fold_handle_for_debug:
                if aggregate_idx is not None:
                    for limb in range(AGGREGATE_LIMBS):
                        for chunk in range(len(handle_slices)):
                            stat_quad.append(_proof_center(
                                LIMB_RADIX ** limb))
                            wit_quad1.append(x_idx[chunk])
                            wit_quad2.append(
                                aggregate_idx[lane][output_phase][limb][chunk])
                            quad_exponents.append(0)
                else:
                    gadget_scale = 1
                    for digit in range(GADGET_DIGITS):
                        centered_gadget = _center_lane(gadget_scale, lane)
                        for limb in range(WITNESS_LIMBS):
                            stat_quad.append(_proof_center(
                                -centered_gadget * LIMB_RADIX ** limb))
                            wit_quad1.append(x_idx[0])
                            wit_quad2.append(
                                right_idx[digit][limb][output_phase][0])
                            quad_exponents.append(0)
                        gadget_scale *= GADGET_BASE
            for chunk in range(REFRESH_AHAT_LIMBS):
                for z_phase in range(PHASES):
                    ahat_phase = (output_phase - z_phase) % PHASES
                    stat_quad.append(_proof_center(LIMB_RADIX ** chunk))
                    wit_quad1.append(
                        refresh_ahat_idx[lane][ahat_phase][chunk])
                    wit_quad2.append(refresh_z_idx[z_phase])
                    quad_exponents.append(1 if z_phase > output_phase else 0)
            rhs = [
                _proof_center(_center_lane(
                    prover.accumulator[PHASES * coefficient + output_phase],
                    lane))
                for coefficient in range(PHASE_N)
            ]
            if cache is not None:
                # A resident worker derives Pack parameters from one
                # canonical zero-RHS statement shared by UAV and GCS.  Build
                # that form directly; constructing the real RHS first forced
                # an otherwise redundant full re-hash of all public matrices.
                rhs = [0] * PHASE_N
            append_native_quadratic(
                stat_lin, wit_lin, stat_quad, wit_quad1, wit_quad2, rhs,
                quad_exponents)
    if phase_executor is not None:
        phase_executor.shutdown(wait=True)
    _q96_timing("append accumulator constraints (native flat)", stage_started)
    if os.environ.get("REVAC_LABRADOR_Q96_PROFILE", "0") != "0":
        print(
            "[q96 Labrador timing]   Python flatten/indices: "
            f"{native_quadratic_prepare_seconds:.6f} s; "
            "native append/hash: "
            f"{native_quadratic_append_seconds:.6f} s",
            flush=True)

    if link is not None:
        stage_started = time.perf_counter()
        assert link_r_idx is not None and link_lift_idx is not None
        link_a, link_d = link.public_matrices()
        one_poly = poly_t(ring, one)
        for row in range(4):
            stat_link: list[Any] = [one_poly]
            wit_link = [link_r_idx[0][row]]
            for col in range(4):
                stat_link.append(poly_t(
                    ring, [_proof_center(value - Q19 if value > Q19 // 2
                                         else value)
                           for value in link_a[row][col]]))
                wit_link.append(link_r_idx[1][col])
            link_handle_coeffs = [
                poly_t(ring, _q19_monomial(link_d[row][0], bit))
                for bit in range(HANDLE_BITS)]
            for chunk, (start, end) in enumerate(handle_slices):
                stat_link.append(polyvec_t(
                    ring, end - start, link_handle_coeffs[start:end]))
                wit_link.append(x_idx[chunk])
            scale = 1
            for limb in range(LINK_DECOMP_LIMBS):
                stat_link.append(poly_t(
                    ring, [_proof_center(Q19 * scale)] +
                    [0] * (PHASE_N - 1)))
                wit_link.append(link_lift_idx[row][limb])
                scale *= LINK_DECOMP_BASE
            rhs = [_proof_center(value - Q19 if value > Q19 // 2 else value)
                   for value in link.commitment[row]]
            if cache is not None:
                rhs = [0] * PHASE_N
            ps.append_statement(stat_link, wit_link, poly_t(ring, rhs))
        _q96_timing("append same-handle link constraints", stage_started)

    if cache is not None and link_r_idx is not None and link_lift_idx is not None:
        cache.link_witness_indices = tuple(
            [link_r_idx[part][row]
             for part in range(2) for row in range(4)] +
            [link_lift_idx[row][limb]
             for row in range(4) for limb in range(LINK_DECOMP_LIMBS)])

    if ps.cur_constraint_num != constraint_count:
        raise ValueError("q96 Labrador constraint accounting mismatch")
    if (check_witness and cache is None and
            int(ps.simple_verify(ps.statement, ps.witness)) != 1):
        raise ValueError("q96 Labrador native simple relation check failed")
    if cache is not None:
        # Parameter generation hashes the statement, including its RHS.  The
        # RHS is dynamic (accumulator and handle commitment), although the
        # parameter dimensions are not.  Canonicalize all RHS polynomials to
        # zero before gen_params so independently provisioned UAV/GCS workers
        # derive one identical public parameter profile.  Restore the real
        # RHS before entering Prove/Verify.
        dynamic_rhs_values = _dynamic_rhs(
            prover, link, constraint_filter)
        dynamic_first = static_constraints
        dynamic_count = constraint_count - static_constraints
        dynamic_prefix_hash = cache.aggregation_prefix_hash
        if dynamic_prefix_hash is None:
            raise ValueError("q96 canonical transcript prefix is unavailable")
        assert cache.static_hash is not None
        params_started = time.perf_counter()
        error = ps.gen_params(ps.params, ps.statement, 1, 0)
        if error != 0:
            raise ValueError(f"q96 Pack parameter generation failed: {error}")
        cache.ps = ps
        cache.constraint_count = constraint_count
        cache.dimensions = tuple(dimensions)
        cache.params = ps.params
        cache.params_owned = True
        cache.bind_profile_id()
        _q96_timing("canonical Pack parameter generation", params_started)
        restore_started = time.perf_counter()
        if link is not None and constraint_filter is None:
            accumulator_rhs_count = accumulator_constraints * PHASE_N
            error = lib.labrador38_py_statement_rebind_tail(
                ps.statement, static_constraints,
                accumulator_constraints,
                ffi.new("int64_t[]",
                        dynamic_rhs_values[:accumulator_rhs_count]),
                ffi.new("uint8_t[16]", dynamic_prefix_hash))
            if error != 0:
                raise ValueError(
                    f"q96 real accumulator RHS restore failed: {error}")
            prefix_hash = ffi.new("uint8_t[16]")
            lib.labrador38_py_statement_hash(prefix_hash, ps.statement)
            cache.accumulator_prefix_hash = bytes(
                ffi.buffer(prefix_hash, 16))
            cache.accumulator_key = tuple(
                int(value) for value in prover.accumulator)
            error = lib.labrador38_py_statement_rebind_tail(
                ps.statement,
                static_constraints + accumulator_constraints, 4,
                ffi.new("int64_t[]",
                        dynamic_rhs_values[accumulator_rhs_count:]),
                ffi.new("uint8_t[16]", cache.accumulator_prefix_hash))
        else:
            error = lib.labrador38_py_statement_rebind_tail(
                ps.statement, dynamic_first, dynamic_count,
                ffi.new("int64_t[]", dynamic_rhs_values),
                ffi.new("uint8_t[16]", dynamic_prefix_hash))
        if error != 0:
            raise ValueError(f"q96 real Pack RHS restore failed: {error}")
        _q96_timing("restore real dynamic RHS once", restore_started)
        if check_witness and int(ps.simple_verify(
                ps.statement, ps.witness)) != 1:
            raise ValueError("q96 restored Labrador relation check failed")
        if not prove:
            return ps, None, None
        proof = ffi.new("labrador38_dch_pack_proof")
        prove_started = time.perf_counter()
        ps.prove(proof, ps.statement, ps.witness, cache.params)
        _q96_timing("native Pack Prove (cold statement)", prove_started)
        return ps, cache.params, proof
    if not prove:
        return ps, None, None
    error, params, proof = ps.pack_prove()
    if error != 0:
        raise ValueError(f"q96 Labrador parameter generation failed: {error}")
    return ps, params, proof


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--public", required=True, type=Path)
    parser.add_argument("--prover", type=Path)
    parser.add_argument(
        "--instance", type=Path,
        help="public GCS verifier instance (never contains a witness)")
    parser.add_argument(
        "--link", type=Path,
        help="q19 same-handle commitment/opening bridge payload")
    parser.add_argument("--prove", action="store_true")
    parser.add_argument("--proof-out", type=Path)
    parser.add_argument("--proof-in", type=Path)
    parser.add_argument(
        "--skip-verify", action="store_true",
        help="bring-up only: do not self-verify a freshly generated proof")
    parser.add_argument("--audit-compiler", action="store_true")
    parser.add_argument(
        "--debug-constraint", metavar="LANE:PHASE",
        help="construct and check only one CRT lane/phase constraint")
    parser.add_argument(
        "--debug-fold-handle", action="store_true",
        help="fold the private handle into linear coefficients for diagnosis")
    parser.add_argument("--debug-witness-selector", action="store_true")
    parser.add_argument("--debug-witness-polyvec", action="store_true")
    args = parser.parse_args()
    public = PublicParameters.decode(args.public.read_bytes())
    if args.prove:
        if args.prover is None or args.link is None or args.instance is not None:
            parser.error("--prove requires --prover and --link, not --instance")
        prover: Any = ProverPayload(args.prover.read_bytes())
        link: Any = LinkPayload.decode(args.link.read_bytes())
        link.audit_relation(prover.handle_vector)
        print("q88/q19 同句柄 exact-link 审计通过", flush=True)
        check_witness = True
    elif args.proof_in is not None:
        if args.instance is None or args.prover is not None or args.link is not None:
            parser.error("--proof-in requires only --instance as relation input")
        instance = InstancePayload.decode(args.instance.read_bytes())
        prover = ZeroProverPayload(instance)
        link = PublicLinkPayload(instance)
        check_witness = False
    else:
        if args.prover is None:
            parser.error("relation audit requires --prover")
        prover = ProverPayload(args.prover.read_bytes())
        link = LinkPayload.decode(args.link.read_bytes()) if args.link else None
        if link is not None:
            link.audit_relation(prover.handle_vector)
            print("q88/q19 同句柄 exact-link 审计通过", flush=True)
        check_witness = True
    if args.audit_compiler:
        audit_compiled_relation(public, prover)
        print("q88 phase/limb compiler relation audit passed", flush=True)
        if not args.debug_constraint and not args.prove:
            return 0
    if args.debug_witness_selector:
        debug_native_witness_selector(prover, args.debug_witness_polyvec)
        return 0
    constraint_filter = None
    if args.debug_constraint:
        fields = args.debug_constraint.split(":")
        if len(fields) != 2:
            parser.error("--debug-constraint must be LANE:PHASE")
        constraint_filter = (int(fields[0]), int(fields[1]))
        if not (0 <= constraint_filter[0] < LANES and
                0 <= constraint_filter[1] < PHASES):
            parser.error("q96 lane/phase is out of range")
    started = time.perf_counter()
    ps, params, proof = build_statement(
        public, prover, link, args.prove, constraint_filter,
        args.debug_fold_handle, check_witness=check_witness)
    elapsed = time.perf_counter() - started
    if args.prove:
        from lazer import lib
        from labrados import pack_verify
        from revac_acc_memver_import_labrador import serialize_labrador38_proof
        proof_bytes = int(lib.labrador38_py_proof_nbytes(params, proof))
        params_bytes = int(lib.labrador38_py_params_nbytes(params))
        proof_scope = "q88 MemVer + q19 same-handle" if link else "q88 accumulator-only"
        print(f"{proof_scope} zk=True proof generated: {elapsed:.3f} s, "
              f"proof={proof_bytes / 1024:.2f} KiB, "
              f"params={params_bytes / 1024:.2f} KiB")
        if not args.skip_verify:
            verify_started = time.perf_counter()
            if not pack_verify((params, proof), ps.output_statement(), "38"):
                raise ValueError(f"{proof_scope} Pack Verify failed")
            print(f"{proof_scope} Pack Verify passed: "
                  f"{time.perf_counter() - verify_started:.3f} s")
        if args.proof_out is not None:
            args.proof_out.write_bytes(serialize_labrador38_proof(
                params, proof, require_zk=True))
            print(f"q88 proof envelope written: {args.proof_out} "
                  f"({args.proof_out.stat().st_size / 1024:.2f} KiB)")
    elif args.proof_in is not None:
        from labrados import pack_verify
        from revac_acc_memver_import_labrador import (
            deserialize_labrador38_proof,
        )
        params_in, proof_in, refs = deserialize_labrador38_proof(
            args.proof_in.read_bytes(), require_zk=True)
        verify_started = time.perf_counter()
        try:
            if not pack_verify((params_in, proof_in),
                               ps.output_statement(), "38"):
                raise ValueError("q96 public-instance Pack Verify failed")
        finally:
            from lazer import lib
            lib.labrador38_py_free_proof(proof_in)
            lib.labrador38_py_free_params(params_in)
            refs.clear()
        print("q88 public-instance Pack Verify passed: "
              f"{time.perf_counter() - verify_started:.3f} s")
    else:
        print(f"q88 accumulator-only native relation passed: {elapsed:.3f} s; "
              f"witness backend polys={ps.witness_polys}, "
              f"constraints={ps.num_constraints}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
