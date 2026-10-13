#!/usr/bin/env python3
"""One isolated single-prime Pack containing main, handle, q19 and refresh.

The captured CRT witness provides only a bounded sample. A new main target is
derived in the independent prime field. This is an integration/performance
experiment; it does not certify the new modulus's scheme security or preserve
the production beta bound automatically.
"""
import argparse
from pathlib import Path
import time

from q88_bridge_core import (SINGLE88_PRIME, compile_single88_bridge_fast,
                             build_multi_bridge_pack_statement)
from q88_bridge_relation import build_single88_handle_bound_relation
from q88_refresh_six_block import make_bridge
from q88_refresh_products import refresh_product_constraints
from q88_single88_profile import single88_binding_constraints
from revac_acc96_wire import (PublicParameters, ProverPayload, LinkPayload,
                              PROVER_BYTES, REFRESH_BLOCKS)


def build_joint_profile(public, prover, link):
    started = time.perf_counter()
    refresh = [make_bridge(public, prover, b) for b in range(REFRESH_BLOCKS)]
    ew, eq = refresh_product_constraints(refresh, prover)
    product_values = {item["ref"]: item["values"][0] for item in ew}
    matrix, witness, meta = build_single88_handle_bound_relation(
        public, prover, modulus=SINGLE88_PRIME)
    # Drop the six effective refresh coefficients. Reconstruct the refresh
    # products using the same native vectors already constrained above.
    matrix, witness = list(matrix[:-6]), list(witness[:-6])
    n = len(witness[0])
    aliases = {}
    main_block = len(refresh)
    # The prefix needs three digits; truncation needs two; ternary handle
    # right witnesses and bounded six handle sums need one. Capacity checks
    # in the native assembler reject samples that do not fit this profile.
    counts = [3, 3] + [2]*22 + [1]*144 + [1]*6
    for block in range(REFRESH_BLOCKS):
        for digit in range(13):
            for layer in range(3):
                ref = ("refresh-product", block, digit, layer)
                col = len(witness)
                witness.append(product_values[ref])
                scalar = pow(128, digit+layer, SINGLE88_PRIME)
                matrix.append(tuple([scalar] + [0]*(n-1)))
                counts.append(1)
                aliases[(main_block, "w", 0, col)] = ref
    main = compile_single88_bridge_fast(matrix, witness, witness_digits=3,
                                        quotient_digits=8)
    bridge_seconds = time.perf_counter()-started
    ewh, eqh, deg0, linear, integers = single88_binding_constraints(witness, meta, link)
    def remap(ref):
        return (main_block, *ref) if isinstance(ref, tuple) and ref[0] in ("w", "q", "c") else ref
    for item in eqh:
        item["linear"] = [(c, remap(r)) for c, r in item.get("linear", ())]
        item["quadratic"] = [(c, remap(l), remap(r))
                             for c, l, r in item.get("quadratic", ())]
    statement_started = time.perf_counter()
    statement, refs = build_multi_bridge_pack_statement(
        refresh, extra_bridges=[main], public_radix_norms=True,
        witness_layer_counts=[[1]*len(b.w_digits) for b in refresh] + [counts],
        witness_aliases=aliases, extra_witness=ew+ewh,
        extra_quadratics=eq+eqh, extra_deg0=deg0,
        extra_linear=linear, extra_integers=integers)
    for alias, target in aliases.items():
        if refs[alias] != refs[target]:
            raise RuntimeError("main relation failed to share refresh product vector")
    audit_values = {item["ref"]: item["values"][0] for item in ew+ewh}
    audit_values[(0, "w", 0, 0)] = refresh[0].w_digits[0][0]
    return statement, refs, audit_values, bridge_seconds, time.perf_counter()-statement_started


def check_tampering(statement, refs, values):
    from labrados import ffi
    # Change one native coefficient while preserving the original norm
    # budget. These checks exercise native relations, not Python prechecks.
    cases = [("selector", (0, "w", 0, 0), True),
             ("refresh-product", ("refresh-product", 0, 0, 0), False),
             ("shared-handle", ("q19bit", 0), True)]
    for label, ref, binary in cases:
        original = list(values[ref])
        modified = original[:]
        modified[0] = 1-original[0] if binary else original[0]+(1 if original[0] < 63 else -1)
        def set_values(poly):
            error = statement.set_witness_vector(
                statement.witness, statement.map[refs[ref]], 1, len(poly)//256,
                ffi.new("int64_t[]", poly))
            if error:
                raise RuntimeError(f"tampering witness write failed: {error}")
        set_values(modified)
        try:
            if int(statement.simple_verify(statement.statement, statement.witness)) == 1:
                raise RuntimeError(f"native Clear Verify accepted tampered {label}")
        finally:
            set_values(original)
        if int(statement.simple_verify(statement.statement, statement.witness)) != 1:
            raise RuntimeError(f"Clear Verify failed after restoring {label}")
        print(f"negative_check={label} rejected restored=passed")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--negative-check", action="store_true")
    args = parser.parse_args()
    public = PublicParameters.decode(Path("/tmp/q88bridge_real/public.bin").read_bytes())
    raw = Path("/tmp/q88bridge_real/prover_link.bin").read_bytes()
    prover = ProverPayload(raw[:PROVER_BYTES])
    link = LinkPayload.decode(raw[PROVER_BYTES:])
    link.audit_relation(prover.handle_vector)
    statement, refs, values, bridge_seconds, statement_seconds = build_joint_profile(public, prover, link)
    if args.negative_check:
        check_tampering(statement, refs, values)
    from labrados import ffi, lib, pack_verify
    started = time.perf_counter()
    error, params, proof = statement.pack_prove()
    prove_seconds = time.perf_counter()-started
    if error:
        raise RuntimeError(f"joint single88 Pack Prove failed: {error}")
    started = time.perf_counter()
    if not pack_verify((params, proof), statement.output_statement(), "38"):
        raise RuntimeError("joint single88 Pack Verify failed")
    verify_seconds = time.perf_counter()-started
    size = int(lib.labrador38_py_proof_nbytes(params, proof))
    serialized = ffi.new("uint8_t[]", size)
    if lib.labrador38_py_proof_serialize(serialized, size, params, proof):
        raise RuntimeError("joint single88 serialization failed")
    decoded = ffi.new("labrador38_dch_pack_proof")
    if lib.labrador38_py_proof_deserialize(decoded, params, serialized, size):
        raise RuntimeError("joint single88 deserialization failed")
    if not pack_verify((params, decoded), statement.output_statement(), "38"):
        raise RuntimeError("joint single88 reload Verify failed")
    print(f"single88_joint_profile modulus={SINGLE88_PRIME}")
    print("main=bound handle=bound q19=bound refresh_selectors=bound refresh_products=bound")
    print(f"native_witness_vectors={statement.num_witness_vectors} "
          f"constraints={statement.num_constraints} "
          f"deg0_constraints={statement.num_deg0_constraints} "
          f"integer_constraints={statement.num_int_constraints}")
    print(f"bridge_coefficient_bound={statement.bridge_coefficient_bound} "
          "terminal_carry=constrained-zero")
    print(f"bridge={bridge_seconds:.3f}s statement={statement_seconds:.3f}s")
    print(f"clear=passed prove={prove_seconds:.3f}s verify={verify_seconds:.3f}s")
    print(f"serialized={size} bytes ({size/1024:.3f} KiB) reload=passed")


if __name__ == "__main__":
    main()
