#!/usr/bin/env python3
"""Fail closed when the implemented q88 accumulator drifts from candidate D."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
import sys
from pathlib import Path

from audit_accumulator_security import audit


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = ROOT / "params/revac-q88-window10000-candidate-v1.json"
CRT_HEADER = ROOT / "code/entities/common/core/include/acc80_crt.h"
EVAL_HEADER = ROOT / "code/entities/common/core/include/acc80_eval.h"
RELATION_HEADER = ROOT / "code/entities/common/core/include/acc96_relation.h"
PROOF_HEADER = ROOT / "code/entities/common/core/include/acc96_proof_layout.h"
TRUNC_HEADER = (
    ROOT / "code/entities/common/core/include/acc96_truncated_params_generated.h"
)
WIRE_MODULE = ROOT / "code/entities/zk_lazer/revac_acc96_wire.py"


class CheckError(RuntimeError):
    pass


def literal_define(path: Path, name: str) -> int:
    text = path.read_text(encoding="utf-8")
    match = re.search(
        rf"^\s*#define\s+{re.escape(name)}\s+"
        rf"(?:(?:UINT64_C|INT32_C)\()?([0-9]+)\)?U?\s*$",
        text,
        re.MULTILINE,
    )
    if match is None:
        raise CheckError(f"cannot parse {name} from {path}")
    return int(match.group(1))


def close(actual: float, expected: float, tolerance: float = 5e-3) -> bool:
    return abs(actual - expected) <= tolerance


def load_wire():
    spec = importlib.util.spec_from_file_location("revac_q88_wire_check", WIRE_MODULE)
    if spec is None or spec.loader is None:
        raise CheckError("cannot load q88 wire module")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args()
    manifest_path = args.manifest
    if not manifest_path.is_absolute():
        manifest_path = ROOT / manifest_path
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    acc = manifest["accumulator"]
    gadget = acc["gadget"]
    sampler = acc["truncated_sampler"]
    window = manifest["window"]
    layout = manifest["proof_layout"]
    sizes = manifest["size_projection"]

    factors = tuple(int(value) for value in acc["crt_factors"])
    q = math.prod(factors)
    if q != int(acc["modulus"]):
        raise CheckError("CRT factors do not reconstruct the manifest modulus")
    if q.bit_length() != acc["modulus_bits"]:
        raise CheckError("manifest modulus bit length is inconsistent")
    if any((factor - 1) % (2 * acc["ring_degree"]) for factor in factors):
        raise CheckError("a CRT factor is not NTT-compatible")

    c_factors = tuple(literal_define(CRT_HEADER, f"ACC80_P{i}") for i in range(3))
    if c_factors != factors:
        raise CheckError("C CRT factors drifted from candidate D")
    expected_literals = {
        (CRT_HEADER, "ACC80_MODULUS_BITS"): acc["modulus_bits"],
        (CRT_HEADER, "ACC80_GADGET_BASE"): gadget["base"],
        (CRT_HEADER, "ACC80_GADGET_DIGITS"): gadget["digits"],
        (CRT_HEADER, "ACC80_GADGET_BOUND"): gadget["digit_bound"],
        (CRT_HEADER, "ACC80_COEFF_PACKED_BYTES"): acc["canonical_coefficient_bytes"],
        (EVAL_HEADER, "ACC80_HANDLE_BITS"): acc["handle_bits"],
        (RELATION_HEADER, "ACC96_TRUNC_K"): sampler["k"],
        (RELATION_HEADER, "ACC96_TRUNC_L"): sampler["l"],
        (RELATION_HEADER, "ACC96_MAX_DELETIONS_PER_WINDOW"): window["runtime_deletion_limit"],
        (PROOF_HEADER, "ACC96_WINDOW_PROOF_LIMBS"): layout["witness_limbs"],
        (PROOF_HEADER, "ACC96_LIFT_PROOF_LIMBS"): layout["exact_lift_limbs"],
    }
    for (path, name), expected in expected_literals.items():
        if literal_define(path, name) != expected:
            raise CheckError(f"{name} drifted from candidate D")

    trunc_text = TRUNC_HEADER.read_text(encoding="utf-8")
    norm_match = re.search(
        r"^#define ACC96_PUBLIC_WITNESS_NORM_SQ_BOUND 0x1p([0-9]+)L$",
        trunc_text,
        re.MULTILINE,
    )
    if norm_match is None or int(norm_match.group(1)) != \
            window["public_witness_norm_sq_bound_log2"]:
        raise CheckError("public witness norm gate drifted from candidate D")

    if gadget["base"] ** gadget["digits"] < q:
        raise CheckError("balanced gadget does not cover the modulus")
    if gadget["digit_bound"] != (gadget["base"] - 1) // 2:
        raise CheckError("balanced gadget digit bound is inconsistent")

    item = audit(
        q,
        ring_degree=acc["ring_degree"],
        module_rows=1,
        deletion_bound=window["runtime_deletion_limit"],
        handle_bits=acc["handle_bits"],
        eval_base=gadget["base"],
        eval_digit_bound=gadget["digit_bound"],
        trunc_base=sampler["base"],
        trunc_l=sampler["l"],
        statistical_bits=sampler["statistical_bits"],
        spectral_slack=sampler["spectral_slack"],
        target_quantum_bits=0,
        q_factors=factors,
    )
    for actual, key in (
        (item.beta0_log2, "beta0_log2_at_10000"),
        (item.beta1_log2, "beta1_log2_at_10000"),
        (item.beta1_margin_log2, "log2_q_over_beta1_at_10000"),
    ):
        if not close(actual, float(window[key])):
            raise CheckError(f"{key} drifted from the arithmetic audit")
    conditional_limit = math.floor(
        window["runtime_deletion_limit"] * 2 ** item.beta1_margin_log2
    )
    if conditional_limit != window["conditional_arithmetic_limit_beta1_lt_q"]:
        raise CheckError("conditional beta1<q deletion limit drifted")

    coefficient_bytes = acc["ring_degree"] * acc["canonical_coefficient_bytes"]
    non_tag_polys = (
        2 + sampler["right_columns"] +
        acc["handle_bits"] * gadget["digits"]
    )
    expected_sizes = {
        "accumulator_kib": coefficient_bytes / 1024,
        "transmitted_full_witness_mib": non_tag_polys * coefficient_bytes / 2**20,
        "ba_mib": acc["handle_bits"] * gadget["digits"] * coefficient_bytes / 2**20,
    }
    for key, actual in expected_sizes.items():
        if not close(actual, float(sizes[key]), 1e-9):
            raise CheckError(f"{key} does not follow from the active layout")
    if non_tag_polys != layout["transmitted_non_tag_polynomials"]:
        raise CheckError("transmitted witness polynomial count drifted")
    if non_tag_polys + 1 != layout["witness_polynomials_including_fixed_tag"]:
        raise CheckError("fixed-tag witness polynomial count drifted")

    def delta4_kib(accumulated_updates: int) -> float:
        coefficient_bound = (
            accumulated_updates * gadget["digits"] * acc["ring_degree"] *
            gadget["digit_bound"] ** 2
        )
        zigzag_bits = (2 * coefficient_bound).bit_length()
        return (4 * gadget["digits"] * acc["ring_degree"] *
                zigzag_bits / 8 / 1024)

    expected_delta_sizes = {
        "delta4_single_worst_case_kib": delta4_kib(1),
        "delta4_accumulated_10000_worst_case_kib": delta4_kib(
            window["runtime_deletion_limit"]),
    }
    for key, actual in expected_delta_sizes.items():
        if not close(actual, float(sizes[key]), 1e-9):
            raise CheckError(f"{key} does not follow from the active layout")

    wire = load_wire()
    expected_wire = (
        wire.Q96 == q and tuple(wire.MODULI) == factors and
        wire.HANDLE_BITS == acc["handle_bits"] and
        wire.GADGET_BASE == gadget["base"] and
        wire.GADGET_DIGITS == gadget["digits"] and
        wire.TRUNC_RIGHT == sampler["right_columns"] and
        wire.NON_TAG_POLYS == non_tag_polys and
        wire.COEFF_BYTES == acc["canonical_coefficient_bytes"] and
        wire.PUBLIC_MAGIC == b"RV88PP02" and wire.PROVER_MAGIC == b"RV88PV02"
    )
    if not expected_wire:
        raise CheckError("Python LaBRADOR wire ABI drifted from candidate D")

    print(f"q88 candidate consistency: passed ({manifest['profile']})")
    print(
        f"window={window['runtime_deletion_limit']}, "
        f"conditional beta1<q limit={conditional_limit}, "
        f"log2(q/beta1)={item.beta1_margin_log2:.3f}"
    )
    print(
        f"accumulator={sizes['accumulator_kib']:.0f} KiB, "
        f"witness={sizes['transmitted_full_witness_mib']:.3f} MiB, "
        f"Ba={sizes['ba_mib']:.3f} MiB"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CheckError as exc:
        print(f"q88 candidate check failed: {exc}")
        raise SystemExit(1)
