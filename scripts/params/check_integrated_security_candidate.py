#!/usr/bin/env python3
"""Fail-closed consistency checker for the integrated paper candidate."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "params/revac112pq-q96-candidate-v1.json"
ACC96_RELATION_HEADER = (
    ROOT / "code/entities/common/core/include/acc96_relation.h"
)
ACC96_PROOF_HEADER = (
    ROOT / "code/entities/common/core/include/acc96_proof_layout.h"
)
ACC96_TRUNCATED_HEADER = (
    ROOT / "code/entities/common/core/include/acc96_truncated_params_generated.h"
)


class CheckError(RuntimeError):
    pass


def read_json(relative: str | Path) -> dict:
    path = Path(relative)
    if not path.is_absolute():
        path = ROOT / path
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise CheckError(f"missing evidence file: {path}") from exc


def close(actual: float, expected: float, tolerance: float = 5e-3) -> bool:
    return abs(actual - expected) <= tolerance


def read_c_uint_define(path: Path, name: str) -> int:
    text = path.read_text(encoding="utf-8")
    match = re.search(
        rf"^\s*#define\s+{re.escape(name)}\s+([0-9]+)U?\s*$",
        text,
        re.MULTILINE,
    )
    if match is None:
        raise CheckError(f"missing literal C parameter {name} in {path}")
    return int(match.group(1))


def selected_accumulator(screen: dict, modulus: int) -> dict:
    rows = [row for row in screen["candidates"] if row["q"] == modulus]
    if len(rows) != 1:
        raise CheckError("selected accumulator row is missing or ambiguous")
    return rows[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT)
    parser.add_argument("--require-release", action="store_true")
    args = parser.parse_args()

    manifest = read_json(args.manifest)
    target = manifest["security_claim_target"]
    sig = manifest["signature_and_show"]
    acc = manifest["accumulator"]
    evidence = manifest["evidence_status"]

    target_window = target["maximum_deletions_per_window"]
    code_window = read_c_uint_define(
        ACC96_RELATION_HEADER, "ACC96_MAX_DELETIONS_PER_WINDOW"
    )
    proof_header = ACC96_PROOF_HEADER.read_text(encoding="utf-8")
    runtime_norm_sq_log2 = int(acc["runtime_witness_norm_sq_log2"])
    truncated_header = ACC96_TRUNCATED_HEADER.read_text(encoding="utf-8")
    norm_match = re.search(
        r"^\s*#define\s+ACC96_PUBLIC_WITNESS_NORM_SQ_BOUND\s+0x1p(\d+)L\s*$",
        truncated_header,
        re.MULTILINE,
    )
    if norm_match is None:
        raise CheckError("cannot parse q96 runtime witness norm gate")
    if int(norm_match.group(1)) != runtime_norm_sq_log2:
        raise CheckError("runtime witness norm gate drifted from candidate manifest")
    if not close(float(acc["runtime_witness_norm_log2"]),
                 runtime_norm_sq_log2 / 2.0):
        raise CheckError("runtime witness norm and squared-norm gates disagree")
    if target_window != 1000:
        raise CheckError("the selected candidate is not the required 1000-deletion window")
    if code_window != target_window:
        raise CheckError("runtime deletion limit drifted from the candidate manifest")
    if acc["window_proof_encoding_max_deletions"] != target_window:
        raise CheckError("proof encoding window drifted from the candidate manifest")
    if not re.search(
        r"^\s*#define\s+ACC96_WINDOW_MAX_DELETIONS\s+"
        r"ACC96_MAX_DELETIONS_PER_WINDOW\s*$",
        proof_header,
        re.MULTILINE,
    ):
        raise CheckError("proof layout does not reuse the runtime deletion limit")
    expected_gate = (
        target_window * acc["window_per_update_coefficient_bound"]
        + acc["window_fresh_witness_coefficient_gate"]
    )
    if acc["window_1000_update_coefficient_gate"] != expected_gate:
        raise CheckError("stored 1000-update coefficient gate drifted")
    if evidence.get("q96_window_1000_runtime_gate") is not True:
        raise CheckError("1000-deletion runtime gate lacks recorded regression evidence")

    tag = read_json(sig["screen_file"])
    selected = [
        row for row in tag["composition_candidates"]
        if row["issue_m2"] == sig["issue_m2"]
        and row["show_m2"] == sig["show_m2"]
    ]
    if len(selected) != 1:
        raise CheckError("selected Tag-Friendly composition row is missing")
    composition = selected[0]
    classical_min = min(composition["classical"].values())
    quantum_min = min(composition["quantum"].values())
    if classical_min < target["classical_bits"]:
        raise CheckError("Tag-Friendly classical target is not met")
    if quantum_min < target["quantum_bits"]:
        raise CheckError("Tag-Friendly quantum target is not met")
    stored_composition = (
        sig["screened_classical_unforgeability_bits"],
        sig["screened_classical_anonymity_bits"],
        sig["screened_quantum_unforgeability_bits"],
        sig["screened_quantum_anonymity_bits"],
    )
    actual_composition = (
        composition["classical"]["unforgeability_bits"],
        composition["classical"]["anonymity_bits"],
        composition["quantum"]["unforgeability_bits"],
        composition["quantum"]["anonymity_bits"],
    )
    if any(not close(actual, stored) for actual, stored in zip(
        actual_composition, stored_composition
    )):
        raise CheckError("stored Tag-Friendly composition values drifted")

    factors = acc["crt_factors"]
    if math.prod(factors) != acc["modulus"]:
        raise CheckError("accumulator CRT factors do not reconstruct q_acc")
    for gamma, file_key, beta_key, qsec_key in (
        (500, "screen_500_file", "beta1_log2_at_500", "quantum_core_svp_bits_at_500"),
        (1000, "screen_1000_file", "beta1_log2_at_1000", "quantum_core_svp_bits_at_1000"),
    ):
        row = selected_accumulator(read_json(acc[file_key]), acc["modulus"])
        if row["deletion_bound"] != gamma or not row["theorem_screen_pass"]:
            raise CheckError(f"accumulator gamma={gamma} theorem screen failed")
        if not close(row["beta1_log2"], acc[beta_key]):
            raise CheckError(f"accumulator gamma={gamma} beta1 value drifted")
        if gamma == target_window and not (
            float(acc["runtime_witness_norm_log2"]) < row["beta0_log2"]
        ):
            raise CheckError("runtime MemVer gate is not strictly below beta0")
        if not close(row["quantum_core_svp_bits"], acc[qsec_key], 0.11):
            raise CheckError(f"accumulator gamma={gamma} M-SIS value drifted")
        if row["quantum_core_svp_bits"] < target["quantum_bits"]:
            raise CheckError(f"accumulator gamma={gamma} quantum target failed")

    mlwe = read_json(acc["mlwe_screen_file"])
    if mlwe.get("q") != acc["modulus"]:
        raise CheckError("accumulator M-LWE evidence uses a different modulus")
    if min(
        item["best_rop_log2"] for item in mlwe["instances"].values()
    ) < target["quantum_bits"]:
        raise CheckError("accumulator M-LWE rough screen misses target")

    blockers = manifest.get("release_blockers", [])
    unfinished = sorted(key for key, value in evidence.items() if value is not True)
    if args.require_release and (
        manifest.get("release_status") != "frozen" or blockers or unfinished
    ):
        raise CheckError(
            "release check failed: candidate status, blockers, or unfinished evidence remain"
        )

    print(f"candidate consistency: passed ({manifest['profile_name']})")
    print(f"enforced deletion window: {target_window}")
    print(
        f"screened lower bound: {classical_min:.2f} classical bits, "
        f"{quantum_min:.2f} quantum bits"
    )
    print(f"release ready: {'yes' if not blockers and not unfinished else 'no'}")
    if unfinished:
        print("unfinished evidence: " + ", ".join(unfinished))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CheckError as exc:
        print(f"candidate check failed: {exc}")
        raise SystemExit(1)
