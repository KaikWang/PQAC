#!/usr/bin/env python3
"""One isolated single-prime Pack containing main, handle, q19 and refresh.

The captured CRT witness provides only a bounded sample. A new main target is
derived in the independent prime field. This is an integration/performance
experiment; it does not certify the new modulus's scheme security or preserve
the production beta bound automatically.
"""
import argparse
import hashlib
import json
import os
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


def build_joint_profile(public, prover, link, *, reference_preparation=False, progress=None):
    convolution = None
    if not reference_preparation:
        from q88_profile_fft import digit_sum_negacyclic
        convolution = digit_sum_negacyclic
    started = time.perf_counter()
    stages, previous = {}, started
    def checkpoint(label):
        nonlocal previous
        now = time.perf_counter()
        stages[label] = now-previous
        previous = now
        if progress:
            progress(f"preparation_stage={label} seconds={stages[label]:.3f}")
    refresh = [make_bridge(public, prover, b, convolution=convolution)
               for b in range(REFRESH_BLOCKS)]
    checkpoint("refresh_bridges")
    ew, eq = refresh_product_constraints(refresh, prover)
    checkpoint("refresh_products")
    product_values = {item["ref"]: item["values"][0] for item in ew}
    matrix, witness, meta = build_single88_handle_bound_relation(
        public, prover, modulus=SINGLE88_PRIME, include_refresh=False)
    checkpoint("handle_relation")
    # Request the handle prefix directly, without constructing effective
    # refresh columns that would be discarded. Register the bound refresh
    # product vectors below as aliases, rather than duplicate witnesses.
    matrix, witness = list(matrix), list(witness)
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
                                        quotient_digits=8, convolution=convolution)
    checkpoint("main_bridge")
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
    statement.preparation_stages = stages
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


def measure_round(statement, cached_params=None):
    from labrados import ffi, lib, pack_verify
    from revac_acc_memver_worker import (native_thread_scope, polxvec_parallel_scope,
                                        verify_aggregation_parallel_scope)
    started = time.perf_counter()
    with native_thread_scope("prove"), polxvec_parallel_scope("prove"), verify_aggregation_parallel_scope("prove"):
        if cached_params is None:
            error, params, proof = statement.pack_prove()
        else:
            params = cached_params
            proof = ffi.new("labrador38_dch_pack_proof")
            statement.prove(proof, statement.statement, statement.witness, params)
            error = 0
    prove_seconds = time.perf_counter()-started
    if error:
        raise RuntimeError(f"joint single88 Pack Prove failed: {error}")
    started = time.perf_counter()
    with native_thread_scope("verify"), polxvec_parallel_scope("verify"), verify_aggregation_parallel_scope("verify"):
        verified = pack_verify((params, proof), statement.output_statement(), "38")
    if not verified:
        raise RuntimeError("joint single88 Pack Verify failed")
    verify_seconds = time.perf_counter()-started
    size = int(lib.labrador38_py_proof_nbytes(params, proof))
    serialized = ffi.new("uint8_t[]", size)
    if lib.labrador38_py_proof_serialize(serialized, size, params, proof):
        raise RuntimeError("joint single88 serialization failed")
    decoded = ffi.new("labrador38_dch_pack_proof")
    if lib.labrador38_py_proof_deserialize(decoded, params, serialized, size):
        raise RuntimeError("joint single88 deserialization failed")
    with native_thread_scope("verify"), polxvec_parallel_scope("verify"), verify_aggregation_parallel_scope("verify"):
        verified = pack_verify((params, decoded), statement.output_statement(), "38")
    if not verified:
        raise RuntimeError("joint single88 reload Verify failed")
    digest = hashlib.sha256(bytes(ffi.buffer(serialized, size))).hexdigest()
    statement.free_proof(decoded)
    if cached_params is not None:
        statement.free_proof(proof)
    return params, {"prove": prove_seconds, "verify": verify_seconds,
                    "bytes": size, "proof_sha256": digest}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--negative-check", action="store_true")
    parser.add_argument("--reference-preparation", action="store_true")
    parser.add_argument("--runs", type=int, default=1)
    parser.add_argument("--result-json")
    parser.add_argument("--parallel-kernels", action="store_true",
                        help="opt in to baseline-derived Gaussian/projection kernel budgets")
    args = parser.parse_args()
    if args.runs < 1:
        parser.error("--runs must be positive")
    threads = str(min(8, max(1, os.cpu_count() or 1)))
    if args.parallel_kernels:
        for variable in ("REVAC_LABRADOR_GAUSSIAN_THREADS",
                         "REVAC_LABRADOR_FROMPOLY_THREADS",
                         "REVAC_LABRADOR_LDD_PROJECTION_THREADS"):
            os.environ.setdefault(variable, threads)
        os.environ.setdefault("REVAC_LABRADOR_DCH_ADDCHECK_PARALLEL", "1")
        os.environ.setdefault("REVAC_LABRADOR_LDD_PROJECTION_PARALLEL", "1")
    public = PublicParameters.decode(Path("/tmp/q88bridge_real/public.bin").read_bytes())
    raw = Path("/tmp/q88bridge_real/prover_link.bin").read_bytes()
    prover = ProverPayload(raw[:PROVER_BYTES])
    link = LinkPayload.decode(raw[PROVER_BYTES:])
    link.audit_relation(prover.handle_vector)
    statement, refs, values, bridge_seconds, statement_seconds = build_joint_profile(
        public, prover, link, reference_preparation=args.reference_preparation, progress=print)
    if args.negative_check:
        check_tampering(statement, refs, values)
    params = None
    rounds = []
    for iteration in range(args.runs):
        params, result = measure_round(statement, params)
        if any(old["proof_sha256"] == result["proof_sha256"] for old in rounds):
            raise RuntimeError("resident run repeated a previous proof byte-for-byte")
        rounds.append(result)
        print(f"round={iteration+1} mode={'cold' if iteration == 0 else 'resident-fixed-instance'} "
              f"prove={result['prove']:.3f}s verify={result['verify']:.3f}s "
              f"serialized={result['bytes']} bytes reload=passed")
    prove_seconds, verify_seconds, size = rounds[-1]["prove"], rounds[-1]["verify"], rounds[-1]["bytes"]
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
    if args.result_json:
        Path(args.result_json).write_text(json.dumps({
            "bridge_seconds": bridge_seconds, "statement_seconds": statement_seconds,
            "preparation_seconds": bridge_seconds+statement_seconds,
            "preparation_stages": statement.preparation_stages,
            "witness_vectors": statement.num_witness_vectors,
            "constraints": statement.num_constraints,
            "bridge_coefficient_bound": statement.bridge_coefficient_bound,
            "fixed_instance_repeats": True, "rounds": rounds,
            "runtime_settings": {key: os.environ.get(key) for key in (
                "REVAC_LABRADOR_THREADS", "REVAC_LABRADOR_THREADS_PROVE",
                "REVAC_LABRADOR_THREADS_VERIFY", "REVAC_LABRADOR_GAUSSIAN_THREADS",
                "REVAC_LABRADOR_FROMPOLY_THREADS", "REVAC_LABRADOR_LDD_PROJECTION_THREADS")},
        }, indent=2)+"\n")


if __name__ == "__main__":
    main()
