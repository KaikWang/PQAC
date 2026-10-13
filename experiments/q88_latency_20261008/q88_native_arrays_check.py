"""Check native/object transcript equality and resident RHS rebinding."""
import random
from lazer import polyring_t, poly_t, polyvec_t
from labrados import proof_statement, ffi, lib
from q88_native_arrays import ArrayPolynomial, ArrayVector, install_array_assembler


def transcript(ps):
    out = ffi.new("uint8_t[16]")
    lib.labrador38_py_statement_hash(out, ps.statement)
    return bytes(ffi.buffer(out, 16))


def assemble(values, *, arrays=False, resident=None, change_coefficient=False):
    n = len(values)
    ring = polyring_t(n, (1 << 38)-107)
    scalar_ring = polyring_t(256, (1 << 38)-107)
    ps = resident or proof_statement([n, n, 256], [1, 1, 1],
                                    [n*64**2, n*64**2, 0], 3, primesize="38",
                                    num_int_constraints=1, num_deg0_constraints=1, zk=True)
    if resident:
        ps.cur_witness_num = ps.witness_polys = 0
        ps.cur_constraint_num = ps.cur_deg0_constraint_num = ps.cur_int_constraint_num = 0
    if arrays:
        install_array_assembler(ps, resident=resident is not None)
    polynomial = ArrayPolynomial if arrays else poly_t
    vector = ArrayVector if arrays else polyvec_t
    for r, coefficients in ((ring, values), (ring, values), (scalar_ring, [1]+[0]*255)):
        v = vector(r, 1)
        v[0] = polynomial(r, coefficients)
        ps.append_witness(v)
    if resident:
        ps.cur_int_constraint_num += 1
        ps.cur_deg0_constraint_num += 1
    else:
        ps.append_integer_statement([1]+[0]*(n-1), [0], 1)
        ps.append_deg0_statement(2)
    one = [1+int(change_coefficient)]+[0]*(n-1)
    zero = polynomial(ring, [0]*n)
    ps.append_statement([polynomial(ring, one)], [0], polynomial(ring, values))
    # Large public scalars exercise native modular conversion as well.
    huge = [405769*256**3]+[0]*(n-1)
    ps.append_quadratic_statement(
        [polynomial(ring, huge), polynomial(ring, [-v for v in huge])],
        [0, 1], [1, -1], [0, 1], [0, 1], zero)
    ps.append_quadratic_statement([polynomial(scalar_ring, [-1]+[0]*255)],
                                  [2], [1], [2], [2], polynomial(scalar_ring, [0]*256))
    if arrays:
        ps.finish_array_assembly()
    if not change_coefficient:
        assert ps.simple_verify(ps.statement, ps.witness) == 1
    return ps


def main():
    rng = random.Random(20261013)
    for n in (256, 4096):
        first = [1]+[rng.randrange(-63, 64) for _ in range(n-1)]
        second = [1]+[rng.randrange(-63, 64) for _ in range(n-1)]
        old = assemble(first)
        resident = assemble(first, arrays=True)
        assert transcript(old) == transcript(resident), "object/native transcript mismatch"
        expected = assemble(second)
        assemble(second, arrays=True, resident=resident)
        assert transcript(expected) == transcript(resident), "resident RHS transcript mismatch"
        try:
            assemble(second, arrays=True, resident=resident, change_coefficient=True)
        except ValueError as error:
            assert "coefficients/layout changed" in str(error)
        else:
            raise AssertionError("resident cache accepted changed coefficients")
    print("object/native transcript, changed RHS and cache invalidation: passed")


if __name__ == "__main__":
    main()
