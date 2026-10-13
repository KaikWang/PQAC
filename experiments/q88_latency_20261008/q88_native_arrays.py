"""Experiment-only array assembler using the production native Pack API.

Only host-side representation changes. Vector boundaries, norms, permutation,
extension degrees and native constraint calls stay identical to labrados.
"""
import numpy as np
import hashlib
import os
import struct
from labrados import ffi, lib


class ArrayPolynomial:
    def __init__(self, ring, values):
        self.ring = ring
        self.dim = 1
        self.values = np.ascontiguousarray(values, dtype=np.int64).reshape(1, ring.deg)


class ArrayVector:
    def __init__(self, ring, dim):
        self.ring, self.dim = ring, dim
        self.values = np.empty((dim, ring.deg), dtype=np.int64)

    def __setitem__(self, index, poly):
        if poly.ring.deg != self.ring.deg:
            raise ValueError("native array polynomial degree mismatch")
        self.values[index] = poly.values[0]


def install_array_assembler(ps, *, resident=False):
    """Install the baseline's buffer path on this experimental instance only."""
    if os.environ.get("REVAC_LABRADOR_Q96_COMPACT_HASH", "0")[:1] not in ("", "0"):
        raise ValueError("isolated profile requires full coefficient hashing; compact hash is not configured")
    position = 0
    if not resident:
        ps.array_constraint_keys = []
    ps.array_resident = resident
    def append_witness(vector):
        return ps.append_witness_pointer(
            ffi.from_buffer("int64_t[]", vector.values), vector.ring.deg, vector.dim)

    def append(linear, witnesses, rhs, quadratic=None):
        nonlocal position
        if len(linear) != len(witnesses) or ps.cur_constraint_num >= ps.num_constraints:
            raise ValueError("native array constraint accounting mismatch")
        if any(term.ring.deg != rhs.ring.deg for term in linear):
            raise ValueError("native array constraint degree mismatch")
        if len(set(witnesses)) != len(witnesses):
            merged, order = {}, []
            for term, index in zip(linear, witnesses):
                if index not in merged:
                    vector = ArrayVector(term.ring, term.dim)
                    vector.values[:] = term.values
                    merged[index] = vector
                    order.append(index)
                else:
                    merged[index].values += term.values
            linear, witnesses = [merged[i] for i in order], order
        coefficients = (np.ascontiguousarray(np.concatenate([term.values for term in linear]))
                        if linear else np.empty((0, rhs.ring.deg), dtype=np.int64))
        indices = ffi.new("size_t[]", [ps.map[i] for i in witnesses])
        lengths = ffi.new("size_t[]", [term.dim for term in linear])
        coeff_ptr = ffi.from_buffer("int64_t[]", coefficients)
        rhs_ptr = ffi.from_buffer("int64_t[]", rhs.values)
        degree_multiple = rhs.ring.deg // 256
        identity = hashlib.sha256()
        identity.update(struct.pack("<QQ", rhs.ring.deg, len(linear)))
        identity.update(np.asarray(witnesses, dtype="<i8").tobytes())
        identity.update(np.asarray([term.dim for term in linear], dtype="<i8").tobytes())
        identity.update(coefficients.astype("<i8", copy=False).tobytes())
        if quadratic is not None:
            for sequence in quadratic:
                identity.update(struct.pack("<Q", len(sequence)))
                identity.update(np.asarray(sequence, dtype="<i8").tobytes())
        key = identity.digest()
        if resident:
            if position >= len(ps.array_constraint_keys) or key != ps.array_constraint_keys[position]:
                raise ValueError("resident statement coefficients/layout changed; rebuild required")
            if position == 0:
                lib.labrador38_py_statement_hash_reset(
                    ps.statement, ffi.new("uint8_t[16]", ps.array_static_prefix))
            scales = quadratic[0] if quadratic is not None else []
            # Hash the original signed RHS buffer, as cold assembly does.
            # RHS-only rebind canonicalizes through polx and can change the
            # hash for negative/noncanonical inputs despite equal equations.
            error = lib.labrador38_py_statement_rebind_constraint_raw(
                ps.statement, position, len(linear), lengths, degree_multiple,
                coeff_ptr, rhs_ptr, len(scales), ffi.new("int64_t[]", scales))
            if error:
                raise ValueError(f"resident raw RHS rebind failed: {error}")
            position += 1
            ps.cur_constraint_num += 1
            return
        if position == 0:
            prefix = ffi.new("uint8_t[16]")
            lib.labrador38_py_statement_hash(prefix, ps.statement)
            ps.array_static_prefix = bytes(ffi.buffer(prefix, 16))
        if quadratic is None:
            error = ps.append_constraint(ps.statement, len(linear), indices, lengths,
                                         degree_multiple, coeff_ptr, rhs_ptr, 1)
        else:
            scales, left, right = quadratic
            if not len(scales) == len(left) == len(right):
                raise ValueError("native array quadratic shape mismatch")
            error = ps.append_quadratic(
                ps.statement, len(linear), len(scales), indices,
                ffi.new("size_t[]", [ps.map[i] for i in left]),
                ffi.new("size_t[]", [ps.map[i] for i in right]), lengths,
                degree_multiple, ffi.new("int64_t[]", scales), coeff_ptr, rhs_ptr)
        if error:
            raise ValueError(f"native array constraint failed: {error}")
        ps.array_constraint_keys.append(key)
        position += 1
        ps.cur_constraint_num += 1

    def append_quadratic(linear, witnesses, scales, left, right, rhs):
        append(linear, witnesses, rhs, (scales, left, right))

    ps.append_witness = append_witness
    ps.append_statement = append
    ps.append_quadratic_statement = append_quadratic
    def finish():
        if resident:
            if position != len(ps.array_constraint_keys):
                raise ValueError("resident constraint inventory changed")
    ps.finish_array_assembly = finish


def append_scalar_refresh(ps, terms, witnesses, rhs, scalar_ring):
    """Exact phase expansion for scalar-selector refresh bridge equations.

    Full-ring operands here have scalar coefficients. Consequently each
    extension phase is independent; selectors multiply the corresponding
    public phase polynomial. No dense full-ring product is approximated.
    """
    phases = rhs.ring.deg // scalar_ring.deg
    expanded = []
    for term, index in zip(terms, witnesses):
        witness_degree = ps.deg_list[ps.map[index]]
        if witness_degree == scalar_ring.deg:
            if term.dim != 1:
                raise ValueError("scalar refresh selector must be singleton")
            expanded.append((True, term.values.reshape(scalar_ring.deg, phases).T,
                             getattr(term, "selector_offset", 0)))
        else:
            if witness_degree != rhs.ring.deg or np.any(term.values[:, 1:]):
                raise ValueError("phase expansion requires scalar full-ring coefficients")
            expanded.append((False, term.values[:, 0], 0))
    right_sides = rhs.values.reshape(scalar_ring.deg, phases).T
    for phase in range(phases):
        phase_terms, phase_witnesses = {}, []
        for index, (selector, values, offset) in zip(witnesses, expanded):
            if index not in phase_terms:
                dim = ps.num_pols_list[ps.map[index]] * ps.deg_list[ps.map[index]] // scalar_ring.deg
                vector = ArrayVector(scalar_ring, dim)
                vector.values.fill(0)
                phase_terms[index] = vector
                phase_witnesses.append(index)
            vector = phase_terms[index]
            if selector:
                vector.values[offset] += values[phase]
            else:
                vector.values[phase::phases, 0] += values
        ps.append_statement([phase_terms[i] for i in phase_witnesses], phase_witnesses,
                            ArrayPolynomial(scalar_ring, right_sides[phase]))
