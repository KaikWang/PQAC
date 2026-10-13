#!/usr/bin/env python3
"""One isolated single-prime Pack containing main, handle, q19 and refresh.

The captured CRT witness provides only a bounded sample. A new main target is
derived in the independent prime field. This is an integration/performance
experiment; it does not certify the new modulus's scheme security or preserve
the production beta bound automatically.
"""
import argparse
import copy
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


class PreparedPublic:
    """Keep the baseline's immutable public bank reuse in the calling layer."""
    def __init__(self, public):
        self.seed, self.trunc_b = public.seed, public.trunc_b
        self.public, self.bank = public, {}

    def get(self, method, *indices):
        key = (method, *indices)
        if key not in self.bank:
            self.bank[key] = tuple(getattr(self.public, method)(*indices))
        return self.bank[key]

    def atrap(self):
        return self.get("atrap")

    def ba(self, bit, digit):
        return self.get("ba", bit, digit)

    def refresh(self, matrix, block):
        return self.get("refresh", matrix, block)


def build_joint_profile(public, prover, link, *, reference_preparation=False,
                        native_arrays=True, linear_group_size=1, progress=None,
                        convolution_cache=None, resident_statement=None,
                        scalar_refresh=False, group_refresh=False, group_handles=False):
    convolution = None
    if not reference_preparation:
        from q88_profile_fft import digit_sum_negacyclic
        convolution = convolution_cache or digit_sum_negacyclic
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
    vector_groups = []
    if group_handles:
        if not scalar_refresh or not group_refresh:
            raise ValueError("handle grouping requires grouped scalar refresh")
        # Native BIN compilation already proves every coefficient binary.
        # Keep all degree-zero and shared-handle predicates explicitly.
        for bit, item in enumerate(eqh[:24]):
            if item["quadratic"] != [(1, ("bit", bit), ("bit", bit))]:
                raise ValueError("unexpected handle binary predicate layout")
        eqh = eqh[24:]
        vector_groups = [[(0, "w", 0, bit) for bit in range(128)],
                         [("bit", bit) for bit in range(24)],
                         [("q19bit", bit) for bit in range(24)]]
        for digit, item in enumerate(eqh):
            right = [(main_block, "w", 0, col) for col in meta["right_refs"][digit]]
            vector_groups.append(right)
            # A single vector dot product is exactly the old sum of 24
            # bit/right products, rather than 24 copies of that dot product.
            item["quadratic"] = [(1, ("bit", 0), right[0])]
    statement_started = time.perf_counter()
    statement, refs = build_multi_bridge_pack_statement(
        refresh, extra_bridges=[main], public_radix_norms=True,
        witness_layer_counts=[[1]*len(b.w_digits) for b in refresh] + [counts],
        witness_aliases=aliases, extra_witness=ew+ewh,
        extra_quadratics=eq+eqh, extra_deg0=deg0,
        extra_linear=linear, extra_integers=integers, native_arrays=native_arrays,
        linear_group_size=linear_group_size, resident_statement=resident_statement,
        scalar_refresh=scalar_refresh, group_refresh=group_refresh and not group_handles,
        vector_groups=vector_groups)
    for alias, target in aliases.items():
        if refs[alias] != refs[target]:
            raise RuntimeError("main relation failed to share refresh product vector")
    audit_values = {item["ref"]: item["values"][0] for item in ew+ewh}
    audit_values[(0, "w", 0, 0)] = refresh[0].w_digits[0][0][:256 if scalar_refresh else n]
    if group_refresh:
        audit_values[(0, "w", 0, 0)] = [v for bit in range(128)
                                       for v in refresh[0].w_digits[bit][0][:256]]
    statement.audit_shapes = {(0, "w", 0, 0): (128, 256)} if group_refresh else {}
    if group_handles:
        audit_values[("q19bit", 0)] = [v for bit in range(24)
                                      for v in audit_values[("q19bit", bit)]]
        statement.audit_shapes[("q19bit", 0)] = (24, 256)
    statement.preparation_stages = stages
    return statement, refs, audit_values, bridge_seconds, time.perf_counter()-statement_started


def check_tampering(statement, refs, values):
    from labrados import ffi
    # Change one native coefficient while preserving the original norm
    # budget. These checks exercise native relations, not Python prechecks.
    cases = [("selector", (0, "w", 0, 0), True, 0),
             ("refresh-product", ("refresh-product", 0, 0, 0), False, 0),
             ("shared-handle", ("q19bit", 0), True, 0),
             ("binary-selector", (0, "w", 0, 0), False, 0)]
    if ("q19bit", 0) in statement.audit_shapes:
        cases.append(("last-shared-handle", ("q19bit", 0), True, 23*256))
    for label, ref, binary, position in cases:
        original = list(values[ref])
        modified = original[:]
        modified[position] = (2 if label == "binary-selector" else
                              1-original[position] if binary else
                              original[position]+(1 if original[position] < 63 else -1))
        def set_values(poly):
            dim, degree = statement.audit_shapes.get(ref, (1, len(poly)))
            error = statement.set_witness_vector(
                statement.witness, statement.map[refs[ref]], dim, degree//256,
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
    parser.add_argument("--object-assembler", action="store_true",
                        help="retain the old poly_t conversion path for comparison")
    parser.add_argument("--linear-group-size", type=int, default=1,
                        help="fixed quotient/carry groups with the original single-vector norm budget")
    parser.add_argument("--runs", type=int, default=1)
    parser.add_argument("--prepare-runs", type=int, default=1,
                        help="rebuild bridges for changed selectors through a resident statement")
    parser.add_argument("--public-cache-mib", type=int, default=768)
    parser.add_argument("--check-resident-transcript", action="store_true")
    parser.add_argument("--scalar-refresh", action="store_true",
                        help="baseline-style degree-256 scalar selectors and exact phase constraints")
    parser.add_argument("--group-refresh", action="store_true",
                        help="one native BIN selector vector; retain all 128 degree-zero predicates")
    parser.add_argument("--group-handles", action="store_true",
                        help="baseline-style handle dot products with the original single-vector norm budget")
    parser.add_argument("--baseline-layout", action="store_true",
                        help="enable scalar refresh, grouped BIN selectors and grouped handle dot products")
    parser.add_argument("--result-json")
    parser.add_argument("--parallel-kernels", action="store_true",
                        help="opt in to baseline-derived Gaussian/projection kernel budgets")
    args = parser.parse_args()
    if args.baseline_layout:
        args.scalar_refresh = args.group_refresh = args.group_handles = True
    if args.object_assembler and args.scalar_refresh:
        parser.error("scalar refresh requires native arrays")
    if args.group_refresh and not args.scalar_refresh:
        parser.error("group refresh requires scalar refresh")
    if args.group_handles and (not args.group_refresh or args.linear_group_size != 1):
        parser.error("handle grouping requires group refresh and excludes carry grouping")
    if args.runs < 1 or args.prepare_runs < 1 or args.public_cache_mib < 0:
        parser.error("run counts must be positive; cache size must be nonnegative")
    if args.prepare_runs > 1 and args.object_assembler:
        parser.error("resident preparation requires the native array assembler")
    threads = str(min(8, max(1, os.cpu_count() or 1)))
    if args.parallel_kernels:
        for variable in ("REVAC_LABRADOR_GAUSSIAN_THREADS",
                         "REVAC_LABRADOR_FROMPOLY_THREADS",
                         "REVAC_LABRADOR_LDD_PROJECTION_THREADS"):
            os.environ.setdefault(variable, threads)
        os.environ.setdefault("REVAC_LABRADOR_DCH_ADDCHECK_PARALLEL", "1")
        os.environ.setdefault("REVAC_LABRADOR_LDD_PROJECTION_PARALLEL", "1")
    public = PreparedPublic(PublicParameters.decode(Path("/tmp/q88bridge_real/public.bin").read_bytes()))
    from q88_profile_fft import PublicDigitCache
    convolution_cache = PublicDigitCache(args.public_cache_mib << 20) if args.prepare_runs > 1 else None
    raw = Path("/tmp/q88bridge_real/prover_link.bin").read_bytes()
    prover = ProverPayload(raw[:PROVER_BYTES])
    link = LinkPayload.decode(raw[PROVER_BYTES:])
    link.audit_relation(prover.handle_vector)
    statement, refs, values, bridge_seconds, statement_seconds = build_joint_profile(
        public, prover, link, reference_preparation=args.reference_preparation,
        native_arrays=not args.object_assembler, linear_group_size=args.linear_group_size,
        convolution_cache=convolution_cache, scalar_refresh=args.scalar_refresh,
        group_refresh=args.group_refresh, group_handles=args.group_handles, progress=print)
    cold_stages = dict(statement.preparation_stages)
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
    preparation_rounds = []
    for iteration in range(1, args.prepare_runs):
        changed = copy.copy(prover)
        bits = list(prover.refresh_r_bits)
        bits[(iteration-1) % len(bits)] ^= 1
        changed.refresh_r_bits = tuple(bits)
        statement, refs, values, bridge_time, statement_time = build_joint_profile(
            public, changed, link, native_arrays=True, linear_group_size=args.linear_group_size,
            convolution_cache=convolution_cache, resident_statement=statement,
            scalar_refresh=args.scalar_refresh, group_refresh=args.group_refresh,
            group_handles=args.group_handles, progress=print)
        if args.negative_check:
            check_tampering(statement, refs, values)
        if args.check_resident_transcript:
            from q88_native_arrays_check import transcript
            expected, _, _, _, _ = build_joint_profile(
                public, changed, link, native_arrays=True, linear_group_size=args.linear_group_size,
                convolution_cache=convolution_cache, scalar_refresh=args.scalar_refresh,
                group_refresh=args.group_refresh, group_handles=args.group_handles)
            if transcript(expected) != transcript(statement):
                raise RuntimeError("resident/fresh full statement transcript mismatch")
            del expected
            print("resident_transcript=fresh-rebuild-identical")
        params, result = measure_round(statement, params)
        result.update(bridge_seconds=bridge_time, statement_seconds=statement_time,
                      preparation_seconds=bridge_time+statement_time,
                      changed_selector=(iteration-1) % len(bits))
        preparation_rounds.append(result)
        print(f"changed_instance={iteration} preparation={bridge_time+statement_time:.3f}s "
              f"prove={result['prove']:.3f}s verify={result['verify']:.3f}s "
              f"serialized={result['bytes']} bytes reload=passed")
    prove_seconds, verify_seconds, size = rounds[-1]["prove"], rounds[-1]["verify"], rounds[-1]["bytes"]
    print(f"single88_joint_profile modulus={SINGLE88_PRIME}")
    print("main=bound handle=bound q19=bound refresh_selectors=bound refresh_products=bound")
    print(f"native_witness_vectors={statement.num_witness_vectors} "
          f"native_witness_polynomials={statement.witness_polys} "
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
            "preparation_stages": cold_stages,
            "witness_vectors": statement.num_witness_vectors,
            "witness_polynomials": statement.witness_polys,
            "constraints": statement.num_constraints,
            "bridge_coefficient_bound": statement.bridge_coefficient_bound,
            "fixed_instance_repeats": True, "rounds": rounds,
            "native_array_assembler": not args.object_assembler,
            "linear_group_size": args.linear_group_size,
            "scalar_refresh": args.scalar_refresh,
            "group_refresh": args.group_refresh,
            "group_handles": args.group_handles,
            "changed_instance_rounds": preparation_rounds,
            "public_cache_bytes": convolution_cache.bytes if convolution_cache else 0,
            "public_cache_hits": convolution_cache.hits if convolution_cache else 0,
            "resident_transcript_checked": args.check_resident_transcript,
            "runtime_settings": {key: os.environ.get(key) for key in (
                "REVAC_LABRADOR_THREADS", "REVAC_LABRADOR_THREADS_PROVE",
                "REVAC_LABRADOR_THREADS_VERIFY", "REVAC_LABRADOR_GAUSSIAN_THREADS",
                "REVAC_LABRADOR_FROMPOLY_THREADS", "REVAC_LABRADOR_LDD_PROJECTION_THREADS")},
        }, indent=2)+"\n")


if __name__ == "__main__":
    main()
