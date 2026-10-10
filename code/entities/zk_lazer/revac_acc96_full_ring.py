"""Experimental degree-4096 compiler for the q88 MemVer relation.

The stable compiler splits every q88 polynomial into sixteen degree-256 phase
polynomials.  This compiler presents the same coefficients to LaZer as native
degree-4096 extension-ring elements, reducing the number of quadratic
constraints and terms without changing the mathematical relation or wire ABI.
"""

from __future__ import annotations

import hashlib
import math
import os
import time
from typing import Any

import numpy as np

from revac_acc96_wire import (
    GADGET_BASE, GADGET_DIGITS, HANDLE_BITS, LANES, LIMB_RADIX, MODULI,
    LIFT_LIMBS, LINK_DECOMP_BASE, LINK_DECOMP_LIMBS, N, NON_TAG_POLYS,
    PHASE_N, PHASES, Q19, REFRESH_BITS, REFRESH_BLOCKS,
    REFRESH_Z_FIRST_SLOT, TAG_SLOT, TRUNC_BASE, TRUNC_RIGHT, WITNESS_LIMBS,
    _center_lane,
)
from revac_profile_generated import PROOF_MODULUS

from revac_acc96_memver_labrador import (
    AGGREGATE_SLICES, JL_INF_SLACK,
    REFRESH_AHAT_LIMBS, VECTOR_COEFF_BOUND,
    Q96StatementCache, _aggregated_quadratics_enabled, _audit_norm_profile,
    _balanced_radix_limbs, _next_power_of_two, _proof_center,
    _proof_center_scaled_array, _q19_monomial, _q96_timing,
    _refresh_lane_phase_bank, _pack_approx_norm_sq_floor,
)


DEGREE_MULTIPLE = N // PHASE_N


def _direct_right_quadratics_enabled() -> bool:
    """Use the original bounded right limbs directly in the q88 products.

    The former layout first copied each lane's public linear combination of
    the right limbs into four auxiliary radix limbs and then multiplied those
    limbs by the handle bits.  Substituting that defining equality into the
    quadratic constraint gives the same relation while removing 4,608 backend
    witness polynomials and six large linear constraints.
    """
    return os.environ.get("REVAC_Q96_DIRECT_RIGHT_QUAD", "0") != "0"


def _coefficient_bound(backend_polynomials: int) -> int:
    return (backend_polynomials * PHASE_N * VECTOR_COEFF_BOUND *
            VECTOR_COEFF_BOUND)


def _approx_requirement(backend_polynomials: int) -> int:
    deterministic = _coefficient_bound(backend_polynomials)
    derived = _next_power_of_two(math.ceil(
        deterministic * JL_INF_SLACK * JL_INF_SLACK))
    return max(derived, _pack_approx_norm_sq_floor())


def _from_phases(phases: Any) -> np.ndarray:
    array = np.asarray(phases, dtype=np.int64)
    if array.shape[-2:] != (PHASES, PHASE_N):
        raise ValueError("q88 full-ring phase shape mismatch")
    return np.ascontiguousarray(
        np.swapaxes(array, -2, -1).reshape(array.shape[:-2] + (N,)),
        dtype=np.int64)


def _witness_poly(prover: Any, slot: int, limb: int) -> np.ndarray:
    return _from_phases([
        prover.witness_phase(slot, limb, phase) for phase in range(PHASES)
    ])


def _compact_witness_enabled() -> bool:
    """Use a compact internal encoding while retaining the q88 wire ABI."""
    return os.environ.get("REVAC_Q96_COMPACT_WITNESS_LIMBS", "0") != "0"


def _compact_witness_limb_count() -> int:
    count = int(os.environ.get("REVAC_Q96_COMPACT_WITNESS_COUNT", "4"))
    if count not in (3, 4):
        raise ValueError("q88 compact witness count must be 3 or 4")
    return count


def _compact_witness_radix(limb_count: int) -> int:
    """Return a safe compact proof radix for the 50000-delete window.

    Four radix-8191 limbs cover only about 52 signed bits and therefore do
    not cover the deterministic q88 window coefficient gate.  Radix 32767
    keeps every digit below 2^14 while its four-limb signed range covers the
    complete gate.  The three-limb experiment uses radix 2^20-1.  Reject
    smaller experimental radices rather than silently proving only the
    particular witness sampled by a benchmark run.
    """
    default_radix = 1048575 if limb_count == 3 else 32767
    radix = int(os.environ.get(
        "REVAC_Q96_COMPACT_WITNESS_RADIX", str(default_radix)))
    if radix < 3 or radix % 2 != 1 or radix > 4194303:
        raise ValueError(
            "q88 compact witness radix must be odd in [3, 4194303]")
    window_coefficient_bound = (
        50000 * 6 * N * 13003 * 13003 + ((1 << 30) - 1))
    signed_capacity = (radix ** limb_count - 1) // 2
    if signed_capacity < window_coefficient_bound:
        raise ValueError(
            "q88 compact witness radix does not cover the 50000-delete "
            f"coefficient gate ({signed_capacity} < "
            f"{window_coefficient_bound})")
    return radix


def _aggregate_limb_count(radix: int) -> int:
    """Smallest balanced-radix width covering every proof-field residue."""
    limbs = 1
    capacity = radix
    while capacity < PROOF_MODULUS:
        limbs += 1
        capacity *= radix
    return limbs


def _scalar_refresh_bits_enabled() -> bool:
    """Keep the hidden Remark-4.3 selector in the degree-256 base ring."""
    return os.environ.get("REVAC_Q96_SCALAR_REFRESH_BITS", "0") != "0"


def _proof_witness_limbs(prover: Any, slot: int) -> list[np.ndarray]:
    """Return the accumulator witness in the selected proof-only radix.

    The C wire format always carries five centered radix-4095 limbs.  The
    compact path reconstructs that exact integer coefficient-wise and
    decomposes it into four centered radix-8191 limbs before appending the
    native statement, so no protocol payload or accumulator arithmetic is
    changed.
    """
    source = [_witness_poly(prover, slot, limb)
              for limb in range(WITNESS_LIMBS)]
    if not _compact_witness_enabled():
        return source
    limb_count = _compact_witness_limb_count()
    radix = _compact_witness_radix(limb_count)
    value = np.zeros_like(source[0], dtype=np.int64)
    scale = 1
    for limb in source:
        value += scale * limb
        scale *= LIMB_RADIX
    return [np.asarray(limb, dtype=np.int64) for limb in
            _balanced_radix_limbs(value, limb_count, radix=radix)]


def _lift_poly(prover: Any, lane: int, limb: int) -> np.ndarray:
    return _from_phases([
        prover.lift_phase(lane, limb, phase) for phase in range(PHASES)
    ])


def _scalar_polynomials(values: Any, degree: int) -> np.ndarray:
    source = np.asarray(values, dtype=np.int64)
    result = np.zeros((source.size, degree), dtype=np.int64)
    result[:, 0] = source.reshape(-1)
    return result


def _dynamic_rhs(prover: Any, link: Any | None) -> np.ndarray:
    values = []
    accumulator = np.asarray(prover.accumulator, dtype=object)
    for lane in range(LANES):
        modulus = MODULI[lane]
        canonical = np.fromiter(
            (int(value) % modulus for value in accumulator),
            dtype=np.int64, count=N)
        centered = np.where(canonical > modulus // 2,
                            canonical - modulus, canonical)
        values.append(np.asarray(
            [_proof_center(int(value)) for value in centered],
            dtype=np.int64))
    if link is not None:
        for row in range(4):
            values.append(np.asarray([
                _proof_center(value - Q19 if value > Q19 // 2 else value)
                for value in link.commitment[row]
            ], dtype=np.int64))
    return np.ascontiguousarray(np.concatenate(values), dtype=np.int64)


def build_statement_full_ring(
        public: Any, prover: Any, link: Any | None, prove: bool,
        constraint_filter: tuple[int, int] | None = None,
        fold_handle_for_debug: bool = False,
        check_witness: bool = True,
        cache: Q96StatementCache | None = None,
        reuse_accumulator_witness: bool = False,
        prepare_prover: bool = False) -> tuple[Any, Any, Any]:
    from lazer import ffi, lib, polyring_t, poly_t, polyvec_t
    from labrados import proof_statement

    if constraint_filter is not None:
        raise ValueError("degree-4096 compiler does not support phase filtering")
    if fold_handle_for_debug:
        raise ValueError("degree-4096 compiler does not support folded handles")
    if not _aggregated_quadratics_enabled():
        raise ValueError("degree-4096 compiler requires aggregated quadratics")

    build_started = time.perf_counter()
    compact_hash = os.environ.get(
        "REVAC_LABRADOR_Q96_COMPACT_HASH", "0") != "0"
    if compact_hash and (cache is None or cache.profile_id is None):
        raise ValueError(
            "q88 full-ring compact hash requires a bound deployment profile")

    scalar_refresh_bits = _scalar_refresh_bits_enabled()
    refresh_constraints = (LANES * PHASES * REFRESH_BLOCKS
                           if scalar_refresh_bits else
                           LANES * REFRESH_BLOCKS)
    direct_right_quadratics = _direct_right_quadratics_enabled()
    compact_witness = _compact_witness_enabled()
    proof_limb_count = (_compact_witness_limb_count()
                        if compact_witness else WITNESS_LIMBS)
    proof_limb_radix = (_compact_witness_radix(proof_limb_count)
                        if compact_witness else LIMB_RADIX)
    proof_limb_bound = proof_limb_radix // 2
    aggregate_limb_count = _aggregate_limb_count(proof_limb_radix)
    aggregation_constraints = (0 if direct_right_quadratics else
                               LANES * len(AGGREGATE_SLICES))
    accumulator_constraints = LANES
    static_constraints = refresh_constraints + aggregation_constraints
    constraint_count = static_constraints + accumulator_constraints
    if link is not None:
        constraint_count += 4

    def rebind_dynamic(ps: Any) -> None:
        assert cache is not None
        dynamic = _dynamic_rhs(prover, link)
        accumulator_coefficients = LANES * N
        accumulator_key = tuple(int(value) for value in prover.accumulator)
        if link is not None:
            if (cache.accumulator_prefix_hash is None or
                    cache.accumulator_key != accumulator_key):
                prefix = cache.aggregation_prefix_hash
                if prefix is None:
                    raise ValueError("q88 full-ring static prefix is unavailable")
                error = lib.labrador38_py_statement_rebind_tail(
                    ps.statement, static_constraints, accumulator_constraints,
                    ffi.from_buffer("int64_t[]", dynamic[:accumulator_coefficients]),
                    ffi.new("uint8_t[16]", prefix))
                if error != 0:
                    raise ValueError(
                        f"q88 full-ring accumulator RHS rebind failed: {error}")
                prefix_out = ffi.new("uint8_t[16]")
                lib.labrador38_py_statement_hash(prefix_out, ps.statement)
                cache.accumulator_prefix_hash = bytes(
                    ffi.buffer(prefix_out, 16))
                cache.accumulator_key = accumulator_key
            error = lib.labrador38_py_statement_rebind_tail(
                ps.statement, static_constraints + accumulator_constraints, 4,
                ffi.from_buffer("int64_t[]", dynamic[accumulator_coefficients:]),
                ffi.new("uint8_t[16]", cache.accumulator_prefix_hash))
        else:
            prefix = cache.aggregation_prefix_hash
            if prefix is None:
                raise ValueError("q88 full-ring static prefix is unavailable")
            error = lib.labrador38_py_statement_rebind_tail(
                ps.statement, static_constraints, accumulator_constraints,
                ffi.from_buffer("int64_t[]", dynamic),
                ffi.new("uint8_t[16]", prefix))
        if error != 0:
            raise ValueError(f"q88 full-ring dynamic RHS rebind failed: {error}")

    if (cache is not None and cache.ps is not None and prove and
            reuse_accumulator_witness and link is not None):
        ps = cache.ps
        if (cache.link_witness_indices is None or
                cache.accumulator_prefix_hash is None or
                cache.accumulator_key != tuple(
                    int(value) for value in prover.accumulator)):
            raise ValueError("q88 full-ring exact witness reuse cache mismatch")
        link_polys = [link.randomness[part][row]
                      for part in range(2) for row in range(4)]
        link_polys.extend(
            link.lift_limbs[row][limb]
            for row in range(4) for limb in range(LINK_DECOMP_LIMBS))
        stage_started = time.perf_counter()
        for logical_index, poly in zip(cache.link_witness_indices, link_polys):
            values = np.ascontiguousarray(poly, dtype=np.int64)
            error = ps.set_witness_vector(
                ps.witness, ps.map[logical_index], 1, 1,
                ffi.from_buffer("int64_t[]", values))
            if error != 0:
                raise ValueError(
                    f"q88 full-ring link witness refresh failed: {error}")
        _q96_timing("full-ring refresh link witness", stage_started)
        stage_started = time.perf_counter()
        link_rhs = _dynamic_rhs(prover, link)[LANES * N:]
        error = lib.labrador38_py_statement_rebind_tail(
            ps.statement, constraint_count - 4, 4,
            ffi.from_buffer("int64_t[]", link_rhs),
            ffi.new("uint8_t[16]", cache.accumulator_prefix_hash))
        if error != 0:
            raise ValueError(
                f"q88 full-ring link-only RHS rebind failed: {error}")
        _q96_timing("full-ring rebind link RHS", stage_started)
        proof = ffi.new("labrador38_dch_pack_proof")
        stage_started = time.perf_counter()
        ps.prove(proof, ps.statement, ps.witness,
                 cache.checked_params(require_zk=True))
        _q96_timing("full-ring native Pack Prove (reuse)", stage_started)
        return ps, cache.params, proof

    if (cache is not None and cache.ps is not None and not prove and
            not check_witness):
        if cache.constraint_count != constraint_count:
            raise ValueError("q88 full-ring verifier profile mismatch")
        stage_started = time.perf_counter()
        rebind_dynamic(cache.ps)
        _q96_timing("full-ring verifier RHS rebind", stage_started)
        return cache.ps, None, None

    vectors: list[tuple[int, np.ndarray]] = []
    degrees: list[int] = []
    dimensions: list[int] = []
    backend_dimensions: list[int] = []
    norms: list[int] = []
    approximate_bounds: dict[int, int] = {}
    accumulator_witness_indices: set[int] = set()

    def append_vector(polys: Any, degree: int = N, binary: bool = False,
                      exact_norm_sq: int | None = None,
                      accumulator_witness: bool = False,
                      approx_requirement: int | None = None) -> int:
        array = np.ascontiguousarray(polys, dtype=np.int64)
        if array.ndim == 1:
            array = array.reshape(1, -1)
        if array.ndim != 2 or array.shape[1] != degree:
            raise ValueError("q88 full-ring witness shape mismatch")
        index = len(vectors)
        vectors.append((degree, array))
        degrees.append(degree)
        dimensions.append(int(array.shape[0]))
        backend_dimension = int(array.shape[0]) * degree // PHASE_N
        backend_dimensions.append(backend_dimension)
        if binary:
            norms.append(0)
        elif exact_norm_sq is not None:
            norms.append(exact_norm_sq)
        else:
            coeff_bound = (proof_limb_bound if compact_witness and
                           accumulator_witness else VECTOR_COEFF_BOUND)
            coefficient_bound = (backend_dimension * PHASE_N *
                                 coeff_bound * coeff_bound)
            if coefficient_bound >= 1 << 64:
                raise ValueError(
                    "q88 compact witness squared-norm exceeds uint64_t")
            norms.append(coefficient_bound)
            derived_requirement = _next_power_of_two(math.ceil(
                coefficient_bound * JL_INF_SLACK * JL_INF_SLACK))
            approximate_bounds[index] = (
                max(approx_requirement, derived_requirement)
                if approx_requirement is not None else derived_requirement)
            if accumulator_witness:
                accumulator_witness_indices.add(index)
        return index

    handle_slices = AGGREGATE_SLICES
    x_idx = [append_vector(
        _scalar_polynomials(prover.handle_vector[start:end], N), binary=True)
        for start, end in handle_slices]
    x_link_idx = [append_vector(
        _scalar_polynomials(prover.handle_vector[start:end], PHASE_N),
        degree=PHASE_N, binary=True)
        for start, end in handle_slices]
    refresh_r_degree = PHASE_N if scalar_refresh_bits else N
    refresh_r_idx = append_vector(
        _scalar_polynomials(prover.refresh_r_bits, refresh_r_degree),
        degree=refresh_r_degree, binary=True)
    refresh_z_idx = append_vector([
        _witness_poly(prover, REFRESH_Z_FIRST_SLOT + block, 0)
        for block in range(REFRESH_BLOCKS)
    ], exact_norm_sq=REFRESH_BLOCKS * N)

    stage_started = time.perf_counter()
    if cache is not None and cache.refresh_bank is not None:
        refresh_bank = cache.refresh_bank
    else:
        refresh_bank = _refresh_lane_phase_bank(public)
        if cache is not None:
            cache.refresh_bank = refresh_bank
    refresh_raw = refresh_bank[0].copy()
    selected = np.flatnonzero(
        np.asarray(prover.refresh_r_bits, dtype=np.int64))
    if selected.size:
        refresh_raw += np.sum(
            refresh_bank[selected + 1], axis=0, dtype=np.int64)
    refresh_ahat_idx: list[list[int]] = []
    refresh_lift_idx: list[int] = []
    refresh_lift_bound = REFRESH_BLOCKS * N * (REFRESH_BITS + 1) ** 2
    for lane, modulus in enumerate(MODULI):
        centered = refresh_raw[lane] % modulus
        centered = np.where(centered > modulus // 2,
                            centered - modulus, centered)
        chunks = _balanced_radix_limbs(centered, REFRESH_AHAT_LIMBS)
        refresh_ahat_idx.append([
            append_vector(_from_phases(
                np.swapaxes(np.asarray(chunk), 0, 1)))
            for chunk in chunks
        ])
        lift_values = (refresh_raw[lane] - centered) // modulus
        refresh_lift_idx.append(append_vector(
            _from_phases(np.swapaxes(lift_values, 0, 1)),
            exact_norm_sq=refresh_lift_bound))
    _q96_timing("full-ring prepare Remark-4.3 bank", stage_started)

    left_idx = [append_vector(poly, accumulator_witness=True)
                for poly in _proof_witness_limbs(prover, 0)]
    v12_idx = [append_vector(poly, accumulator_witness=True)
               for poly in _proof_witness_limbs(prover, 1)]
    v2_limbs = [_proof_witness_limbs(prover, 2 + source)
                for source in range(TRUNC_RIGHT)]
    v2_idx = [append_vector(
        [v2_limbs[source][limb] for source in range(TRUNC_RIGHT)],
        accumulator_witness=True)
        for limb in range(proof_limb_count)]
    right_limbs = [
        [_proof_witness_limbs(
            prover, TAG_SLOT + 1 + bit * GADGET_DIGITS + digit)
         for bit in range(HANDLE_BITS)]
        for digit in range(GADGET_DIGITS)]
    right_idx = [[[ 
        append_vector([
            right_limbs[digit][bit][limb]
            for bit in range(start, end)
        ], accumulator_witness=True,
            approx_requirement=(
                _pack_approx_norm_sq_floor() // len(handle_slices)))
        for start, end in handle_slices]
        for limb in range(proof_limb_count)]
        for digit in range(GADGET_DIGITS)]

    if os.environ.get("REVAC_Q96_WITNESS_LIMB_AUDIT", "0") != "0":
        limb_max = []
        audit_slots = (
            [0, 1] + list(range(2, 2 + TRUNC_RIGHT)) +
            [TAG_SLOT + 1 + bit * GADGET_DIGITS + digit
             for bit in range(HANDLE_BITS)
             for digit in range(GADGET_DIGITS)] +
            [REFRESH_Z_FIRST_SLOT + block
             for block in range(REFRESH_BLOCKS)])
        for limb in range(WITNESS_LIMBS):
            maximum = 0
            for slot in audit_slots:
                for phase in range(PHASES):
                    values = prover.witness_phase(slot, limb, phase)
                    maximum = max(maximum, max(abs(int(value))
                                               for value in values))
            limb_max.append(maximum)
        print(f"[q88 limb audit] witness limb maxima={limb_max}", flush=True)

    aggregate_idx: list[list[list[int]]] = []
    if not direct_right_quadratics:
        aggregate_aux_requirement = _next_power_of_two(math.ceil(
            LANES * aggregate_limb_count *
                    _coefficient_bound(HANDLE_BITS * DEGREE_MULTIPLE) *
            JL_INF_SLACK * JL_INF_SLACK))
        for lane in range(LANES):
            total = np.zeros((HANDLE_BITS, N), dtype=np.int64)
            gadget_scale = 1
            for digit in range(GADGET_DIGITS):
                centered_gadget = _center_lane(gadget_scale, lane)
                for limb in range(proof_limb_count):
                    coefficient = _proof_center(
                        -centered_gadget * proof_limb_radix ** limb)
                    for chunk, (start, end) in enumerate(handle_slices):
                        total[start:end] += coefficient * vectors[
                            right_idx[digit][limb][chunk]][1]
                gadget_scale *= GADGET_BASE
            total %= PROOF_MODULUS
            total[total > PROOF_MODULUS // 2] -= PROOF_MODULUS
            lane_indices = []
            for digit_values in _balanced_radix_limbs(
                    total, aggregate_limb_count, radix=proof_limb_radix):
                lane_indices.append([
                    append_vector(
                        digit_values[start:end],
                        approx_requirement=aggregate_aux_requirement)
                    for start, end in handle_slices])
            aggregate_idx.append(lane_indices)

    lift_idx = [[append_vector(_lift_poly(prover, lane, limb))
                 for limb in range(LIFT_LIMBS)]
                for lane in range(LANES)]

    link_r_idx: list[list[int]] | None = None
    link_lift_idx: list[list[int]] | None = None
    if link is not None:
        link_r_idx = [[append_vector(
            link.randomness[part][row], degree=PHASE_N,
            exact_norm_sq=PHASE_N)
            for row in range(4)] for part in range(2)]
        link_lift_idx = [[append_vector(
            link.lift_limbs[row][limb], degree=PHASE_N,
            exact_norm_sq=PHASE_N * (LINK_DECOMP_BASE // 2) ** 2)
            for limb in range(LINK_DECOMP_LIMBS)] for row in range(4)]

    base_non_tag_polys = NON_TAG_POLYS - REFRESH_BLOCKS
    expected_backend = (
        HANDLE_BITS * DEGREE_MULTIPLE + HANDLE_BITS +
        REFRESH_BITS * (1 if scalar_refresh_bits else DEGREE_MULTIPLE) +
        base_non_tag_polys * proof_limb_count * PHASES +
        REFRESH_BLOCKS * PHASES +
        LANES * REFRESH_BLOCKS * PHASES * (REFRESH_AHAT_LIMBS + 1) +
        LANES * LIFT_LIMBS * PHASES)
    if not direct_right_quadratics:
        expected_backend += (
            LANES * aggregate_limb_count * HANDLE_BITS * PHASES)
    if link is not None:
        expected_backend += 2 * 4 + 4 * LINK_DECOMP_LIMBS
    if sum(backend_dimensions) != expected_backend:
        raise ValueError("q88 full-ring witness layout accounting mismatch")
    _audit_norm_profile(
        backend_dimensions, approximate_bounds, accumulator_witness_indices,
        limb_radix=proof_limb_radix, limb_count=proof_limb_count)
    approx_norm_list = list(approximate_bounds.items())
    _q96_timing("full-ring prepare witness views", build_started)

    cache_hit = cache is not None and cache.ps is not None
    if cache_hit:
        assert cache is not None and cache.ps is not None
        if (cache.dimensions != tuple(backend_dimensions) or
                cache.constraint_count != constraint_count):
            raise ValueError("q88 full-ring resident profile mismatch")
        ps = cache.ps
        if prove or prepare_prover:
            stage_started = time.perf_counter()
            ps.cur_witness_num = 0
            ps.witness_polys = 0
            for degree, polys in vectors:
                ps.append_witness_pointer(
                    ffi.from_buffer("int64_t[]", polys), degree, len(polys))
            _q96_timing("full-ring load dynamic witness", stage_started)
        stage_started = time.perf_counter()
        rebind_dynamic(ps)
        _q96_timing("full-ring rebind dynamic RHS", stage_started)
        if check_witness and int(ps.simple_verify(
                ps.statement, ps.witness)) != 1:
            raise ValueError("q88 full-ring cached relation check failed")
        if not prove:
            return ps, None, None
        if cache.params is None:
            raise ValueError("q88 full-ring resident Pack parameters missing")
        proof = ffi.new("labrador38_dch_pack_proof")
        stage_started = time.perf_counter()
        ps.prove(proof, ps.statement, ps.witness, cache.params)
        _q96_timing("full-ring native Pack Prove", stage_started)
        return ps, cache.params, proof

    ps = proof_statement(
        degrees, dimensions, norms, num_constraints=constraint_count,
        primesize="38", num_int_constraints=HANDLE_BITS,
        num_deg0_constraints=2 * HANDLE_BITS + REFRESH_BITS, zk=True,
        approx_norm_list=approx_norm_list)
    stage_started = time.perf_counter()
    for degree, polys in vectors:
        ps.append_witness_pointer(
            ffi.from_buffer("int64_t[]", polys), degree, len(polys))
    _q96_timing("full-ring load cold witness", stage_started)

    stage_started = time.perf_counter()
    for chunk, (start, end) in enumerate(handle_slices):
        count = end - start
        error = lib.labrador38_py_append_deg0_logical_vector_constraints(
            ps.statement, ps.map[x_idx[chunk]], count, DEGREE_MULTIPLE)
        if error != 0:
            raise ValueError(f"q88 full-ring handle degree-zero failed: {error}")
        ps.cur_deg0_constraint_num += count
    if scalar_refresh_bits:
        error = lib.labrador38_py_append_deg0_vector_constraints(
            ps.statement, ps.map[refresh_r_idx], REFRESH_BITS)
    else:
        error = lib.labrador38_py_append_deg0_logical_vector_constraints(
            ps.statement, ps.map[refresh_r_idx], REFRESH_BITS,
            DEGREE_MULTIPLE)
    if error != 0:
        raise ValueError(f"q88 full-ring refresh-r degree-zero failed: {error}")
    ps.cur_deg0_constraint_num += REFRESH_BITS
    for chunk, (start, end) in enumerate(handle_slices):
        count = end - start
        error = lib.labrador38_py_append_deg0_vector_constraints(
            ps.statement, ps.map[x_link_idx[chunk]], count)
        if error != 0:
            raise ValueError(f"q88 q19 handle degree-zero failed: {error}")
        ps.cur_deg0_constraint_num += count
        for local_bit in range(count):
            full_coefficients = [0] * (count * N)
            link_coefficients = [0] * (count * PHASE_N)
            full_coefficients[local_bit * N] = 1
            link_coefficients[local_bit * PHASE_N] = -1
            ps.append_integer_statement(
                full_coefficients + link_coefficients,
                [x_idx[chunk], x_link_idx[chunk]], 0)
    _q96_timing("full-ring append scalar bindings", stage_started)

    def append_native_constraint(
            linear_terms: list[np.ndarray], witnesses: list[int],
            rhs: np.ndarray, degree_multiple: int) -> None:
        lengths = [int(term.shape[0]) for term in linear_terms]
        flat = np.ascontiguousarray(
            np.concatenate(linear_terms, axis=0).reshape(-1), dtype=np.int64)
        rhs_values = np.ascontiguousarray(rhs, dtype=np.int64).reshape(-1)
        error = ps.append_constraint(
            ps.statement, len(witnesses),
            ffi.new("size_t[]", [ps.map[index] for index in witnesses]),
            ffi.new("size_t[]", lengths), degree_multiple,
            ffi.from_buffer("int64_t[]", flat),
            ffi.from_buffer("int64_t[]", rhs_values), 1)
        if error != 0:
            raise ValueError(f"q88 full-ring linear constraint failed: {error}")
        ps.cur_constraint_num += 1

    stage_started = time.perf_counter()
    if scalar_refresh_bits:
        refresh_bank_i32 = np.ascontiguousarray(refresh_bank, dtype=np.int32)
        refresh_chunk_indices = [
            ps.map[refresh_ahat_idx[lane][chunk]]
            for lane in range(LANES)
            for _phase in range(PHASES)
            for chunk in range(REFRESH_AHAT_LIMBS)
        ]
        refresh_lift_indices = [
            ps.map[refresh_lift_idx[lane]]
            for lane in range(LANES)
            for _phase in range(PHASES)
        ]
        error = lib.labrador38_py_append_refresh_phase_constraints(
            ps.statement, ps.map[refresh_r_idx],
            ffi.new("size_t[]", refresh_chunk_indices),
            ffi.new("size_t[]", refresh_lift_indices),
            ffi.from_buffer("int32_t[]", refresh_bank_i32), REFRESH_BITS,
            REFRESH_BLOCKS, LANES, PHASES, PHASE_N,
            ffi.new("uint32_t[]", MODULI), LIMB_RADIX,
            REFRESH_AHAT_LIMBS)
        if error != 0:
            raise ValueError(
                f"q88 scalar refresh constraints failed: {error}")
        ps.cur_constraint_num += refresh_constraints
    else:
        for lane, modulus in enumerate(MODULI):
            for block in range(REFRESH_BLOCKS):
                r_coefficients = _from_phases(
                    refresh_bank[1:, lane, :, block, :])
                terms = [r_coefficients]
                witnesses = [refresh_r_idx]
                scale = 1
                for chunk in range(REFRESH_AHAT_LIMBS):
                    selector = np.zeros((REFRESH_BLOCKS, N), dtype=np.int64)
                    selector[block, 0] = _proof_center(-scale)
                    terms.append(selector)
                    witnesses.append(refresh_ahat_idx[lane][chunk])
                    scale *= LIMB_RADIX
                selector = np.zeros((REFRESH_BLOCKS, N), dtype=np.int64)
                selector[block, 0] = _proof_center(-modulus)
                terms.append(selector)
                witnesses.append(refresh_lift_idx[lane])
                rhs = -_from_phases(refresh_bank[0, lane, :, block, :])
                rhs = np.asarray([_proof_center(int(value)) for value in rhs],
                                 dtype=np.int64)
                append_native_constraint(
                    terms, witnesses, rhs, DEGREE_MULTIPLE)
    if cache is not None:
        hash_out = ffi.new("uint8_t[16]")
        lib.labrador38_py_statement_hash(hash_out, ps.statement)
        cache.static_hash = bytes(ffi.buffer(hash_out, 16))
    _q96_timing("full-ring append refresh constraints", stage_started)

    stage_started = time.perf_counter()
    if not direct_right_quadratics:
        for lane in range(LANES):
            for chunk, (start, end) in enumerate(handle_slices):
                rank = end - start
                witnesses = []
                coefficients = []
                gadget_scale = 1
                for digit in range(GADGET_DIGITS):
                    centered_gadget = _center_lane(gadget_scale, lane)
                    for limb in range(proof_limb_count):
                        witnesses.append(right_idx[digit][limb][chunk])
                        coefficients.append(_proof_center(
                            -centered_gadget * proof_limb_radix ** limb))
                    gadget_scale *= GADGET_BASE
                for limb in range(aggregate_limb_count):
                    witnesses.append(aggregate_idx[lane][limb][chunk])
                    coefficients.append(
                        _proof_center(-(proof_limb_radix ** limb)))
                terms = []
                extension_degree = rank * DEGREE_MULTIPLE
                for coefficient in coefficients:
                    term = np.zeros((1, rank * N), dtype=np.int64)
                    term[0, 0] = coefficient
                    terms.append(term)
                append_native_constraint(
                    terms, witnesses, np.zeros(rank * N, dtype=np.int64),
                    extension_degree)
    if cache is not None:
        hash_out = ffi.new("uint8_t[16]")
        lib.labrador38_py_statement_hash(hash_out, ps.statement)
        cache.aggregation_prefix_hash = bytes(ffi.buffer(hash_out, 16))
    _q96_timing("full-ring append aggregate constraints", stage_started)

    if compact_hash:
        assert cache is not None and cache.profile_id is not None
        semantic_prefix = hashlib.shake_256(
            b"REVAC-Q88-LABRADOR38-STATEMENT-HASH-V3-FULL-RING\x00" +
            cache.profile_id).digest(16)
        lib.labrador38_py_statement_hash_reset(
            ps.statement, ffi.new("uint8_t[16]", semantic_prefix))
        cache.static_hash = semantic_prefix
        cache.aggregation_prefix_hash = semantic_prefix

    def center_public(poly: Any, lane: int) -> np.ndarray:
        modulus = MODULI[lane]
        canonical = np.fromiter(
            (int(value) % modulus for value in poly),
            dtype=np.int64, count=N)
        return np.where(canonical > modulus // 2,
                        canonical - modulus, canonical).astype(
                            np.int64, copy=False)

    atrap = public.atrap()
    ba = {(bit, digit): public.ba(bit, digit)
          for bit in range(HANDLE_BITS)
          for digit in range(GADGET_DIGITS)}
    one = np.zeros(N, dtype=np.int64)
    one[0] = 1

    def append_native_quadratic(
            linear_terms: list[np.ndarray], linear_witnesses: list[int],
            quadratic_coefficients: list[int], quadratic_left: list[int],
            quadratic_right: list[int], rhs: np.ndarray) -> None:
        lengths = [int(term.shape[0]) for term in linear_terms]
        flat = np.ascontiguousarray(
            np.concatenate(linear_terms, axis=0).reshape(-1), dtype=np.int64)
        rhs_values = np.ascontiguousarray(rhs, dtype=np.int64)
        error = ps.append_quadratic(
            ps.statement, len(linear_terms), len(quadratic_coefficients),
            ffi.new("size_t[]", [ps.map[i] for i in linear_witnesses]),
            ffi.new("size_t[]", [ps.map[i] for i in quadratic_left]),
            ffi.new("size_t[]", [ps.map[i] for i in quadratic_right]),
            ffi.new("size_t[]", lengths), DEGREE_MULTIPLE,
            ffi.new("int64_t[]", quadratic_coefficients),
            ffi.from_buffer("int64_t[]", flat),
            ffi.from_buffer("int64_t[]", rhs_values))
        if error != 0:
            raise ValueError(f"q88 full-ring quadratic failed: {error}")
        ps.cur_constraint_num += 1

    stage_started = time.perf_counter()
    accumulator = np.asarray(prover.accumulator, dtype=object)
    for lane, modulus in enumerate(MODULI):
        terms: list[np.ndarray] = []
        witnesses: list[int] = []
        limb_scale = 1
        for limb in range(proof_limb_count):
            terms.append((_proof_center_scaled_array(one, limb_scale)).reshape(1, N))
            witnesses.append(left_idx[limb])
            limb_scale *= proof_limb_radix
        atrap_lane = center_public(atrap, lane)
        limb_scale = 1
        for limb in range(proof_limb_count):
            terms.append(_proof_center_scaled_array(
                atrap_lane, limb_scale).reshape(1, N))
            witnesses.append(v12_idx[limb])
            limb_scale *= proof_limb_radix
        trunc_lanes = [center_public(poly, lane) for poly in public.trunc_b]
        for limb in range(proof_limb_count):
            limb_scale = proof_limb_radix ** limb
            scale = TRUNC_BASE * TRUNC_BASE
            coefficients = []
            for source in range(TRUNC_RIGHT):
                coefficient = -_proof_center_scaled_array(
                    trunc_lanes[source], limb_scale)
                coefficient = coefficient.copy()
                coefficient[0] = _proof_center(
                    int(coefficient[0]) + limb_scale *
                    _center_lane(scale, lane))
                coefficients.append(coefficient)
                scale *= TRUNC_BASE
            terms.append(np.stack(coefficients))
            witnesses.append(v2_idx[limb])
        for digit in range(GADGET_DIGITS):
            ba_lane = [center_public(ba[(bit, digit)], lane)
                       for bit in range(HANDLE_BITS)]
            for limb in range(proof_limb_count):
                scaled = np.stack([
                    _proof_center_scaled_array(poly, proof_limb_radix ** limb)
                    for poly in ba_lane
                ])
                for chunk, (start, end) in enumerate(handle_slices):
                    terms.append(scaled[start:end])
                    witnesses.append(right_idx[digit][limb][chunk])
        for limb in range(LIFT_LIMBS):
            coefficient = np.zeros((1, N), dtype=np.int64)
            coefficient[0, 0] = _proof_center(
                modulus * LIMB_RADIX ** limb)
            terms.append(coefficient)
            witnesses.append(lift_idx[lane][limb])

        quad_coefficients = []
        quad_left = []
        quad_right = []
        if direct_right_quadratics:
            gadget_scale = 1
            for digit in range(GADGET_DIGITS):
                centered_gadget = _center_lane(gadget_scale, lane)
                for limb in range(proof_limb_count):
                    coefficient = _proof_center(
                        -centered_gadget * proof_limb_radix ** limb)
                    for chunk in range(len(handle_slices)):
                        quad_coefficients.append(coefficient)
                        quad_left.append(x_idx[chunk])
                        quad_right.append(right_idx[digit][limb][chunk])
                gadget_scale *= GADGET_BASE
        else:
            for limb in range(aggregate_limb_count):
                for chunk in range(len(handle_slices)):
                    quad_coefficients.append(
                        _proof_center(proof_limb_radix ** limb))
                    quad_left.append(x_idx[chunk])
                    quad_right.append(aggregate_idx[lane][limb][chunk])
        for chunk in range(REFRESH_AHAT_LIMBS):
            quad_coefficients.append(_proof_center(LIMB_RADIX ** chunk))
            quad_left.append(refresh_ahat_idx[lane][chunk])
            quad_right.append(refresh_z_idx)

        canonical = np.fromiter(
            (int(value) % modulus for value in accumulator),
            dtype=np.int64, count=N)
        centered = np.where(canonical > modulus // 2,
                            canonical - modulus, canonical)
        rhs = np.asarray([_proof_center(int(value)) for value in centered],
                         dtype=np.int64)
        if cache is not None:
            rhs = np.zeros(N, dtype=np.int64)
        append_native_quadratic(
            terms, witnesses, quad_coefficients, quad_left, quad_right, rhs)
    _q96_timing("full-ring append accumulator constraints", stage_started)

    if link is not None:
        assert link_r_idx is not None and link_lift_idx is not None
        stage_started = time.perf_counter()
        ring = polyring_t(PHASE_N, PROOF_MODULUS)
        one_poly = poly_t(ring, [1] + [0] * (PHASE_N - 1))
        link_a, link_d = link.public_matrices()
        for row in range(4):
            stat_link: list[Any] = [one_poly]
            wit_link = [link_r_idx[0][row]]
            for col in range(4):
                stat_link.append(poly_t(ring, [
                    _proof_center(value - Q19 if value > Q19 // 2 else value)
                    for value in link_a[row][col]
                ]))
                wit_link.append(link_r_idx[1][col])
            handle_coefficients = [
                poly_t(ring, _q19_monomial(link_d[row][0], bit))
                for bit in range(HANDLE_BITS)]
            for chunk, (start, end) in enumerate(handle_slices):
                stat_link.append(polyvec_t(
                    ring, end - start, handle_coefficients[start:end]))
                wit_link.append(x_link_idx[chunk])
            scale = 1
            for limb in range(LINK_DECOMP_LIMBS):
                stat_link.append(poly_t(
                    ring, [_proof_center(Q19 * scale)] +
                    [0] * (PHASE_N - 1)))
                wit_link.append(link_lift_idx[row][limb])
                scale *= LINK_DECOMP_BASE
            rhs = [_proof_center(
                value - Q19 if value > Q19 // 2 else value)
                for value in link.commitment[row]]
            if cache is not None:
                rhs = [0] * PHASE_N
            ps.append_statement(stat_link, wit_link, poly_t(ring, rhs))
        cache_link_indices = tuple(
            [link_r_idx[part][row] for part in range(2) for row in range(4)] +
            [link_lift_idx[row][limb]
             for row in range(4) for limb in range(LINK_DECOMP_LIMBS)])
        if cache is not None:
            cache.link_witness_indices = cache_link_indices
        _q96_timing("full-ring append q19 link", stage_started)

    if ps.cur_constraint_num != constraint_count:
        raise ValueError("q88 full-ring constraint accounting mismatch")
    if ps.cur_int_constraint_num != HANDLE_BITS:
        raise ValueError("q88 full-ring handle-binding accounting mismatch")
    if check_witness and cache is None and int(ps.simple_verify(
            ps.statement, ps.witness)) != 1:
        raise ValueError("q88 full-ring native relation check failed")

    if cache is not None:
        if cache.aggregation_prefix_hash is None:
            raise ValueError("q88 full-ring transcript prefix missing")
        dynamic = _dynamic_rhs(prover, link)
        params_started = time.perf_counter()
        error = ps.gen_params(ps.params, ps.statement, 1, 0)
        if error != 0:
            raise ValueError(f"q88 full-ring parameter generation failed: {error}")
        cache.ps = ps
        cache.constraint_count = constraint_count
        cache.dimensions = tuple(backend_dimensions)
        cache.params = ps.params
        cache.params_owned = True
        cache.bind_profile_id()
        _q96_timing("full-ring Pack parameter generation", params_started)
        stage_started = time.perf_counter()
        accumulator_coefficients = LANES * N
        error = lib.labrador38_py_statement_rebind_tail(
            ps.statement, static_constraints, accumulator_constraints,
            ffi.from_buffer("int64_t[]", dynamic[:accumulator_coefficients]),
            ffi.new("uint8_t[16]", cache.aggregation_prefix_hash))
        if error != 0:
            raise ValueError(
                f"q88 full-ring accumulator restore failed: {error}")
        prefix_out = ffi.new("uint8_t[16]")
        lib.labrador38_py_statement_hash(prefix_out, ps.statement)
        cache.accumulator_prefix_hash = bytes(ffi.buffer(prefix_out, 16))
        cache.accumulator_key = tuple(int(value) for value in prover.accumulator)
        if link is not None:
            error = lib.labrador38_py_statement_rebind_tail(
                ps.statement, static_constraints + accumulator_constraints, 4,
                ffi.from_buffer("int64_t[]", dynamic[accumulator_coefficients:]),
                ffi.new("uint8_t[16]", cache.accumulator_prefix_hash))
            if error != 0:
                raise ValueError(f"q88 full-ring link restore failed: {error}")
        _q96_timing("full-ring restore dynamic RHS", stage_started)
        if check_witness and int(ps.simple_verify(
                ps.statement, ps.witness)) != 1:
            raise ValueError("q88 full-ring restored relation check failed")
        if not prove:
            return ps, None, None
        proof = ffi.new("labrador38_dch_pack_proof")
        stage_started = time.perf_counter()
        ps.prove(proof, ps.statement, ps.witness, cache.params)
        _q96_timing("full-ring native Pack Prove (cold)", stage_started)
        return ps, cache.params, proof

    if not prove:
        return ps, None, None
    error, params, proof = ps.pack_prove()
    if error != 0:
        raise ValueError(f"q88 full-ring Pack parameter generation failed: {error}")
    return ps, params, proof
