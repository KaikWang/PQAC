#!/usr/bin/env python3
"""Search auditable Tag-Friendly Issue/Show parameter candidates.

This wrapper deliberately reuses the upstream parameter classes instead of
copying their formulas.  The upstream file currently executes one default
profile at import time, so we load only the definitions above that profile.
Run it with Sage/Python from an environment containing the upstream estimator.

The result is a *component-security screen*, not an end-to-end theorem claim:
the anonymous-credential reduction and the accumulator composition still have
to account for their query and hybrid losses explicitly.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import math
import sys
import time
from pathlib import Path


def load_upstream(path: Path, quantum: bool) -> dict[str, object]:
    source = path.read_text(encoding="utf-8")
    marker = "\nZK_OPTIMIZATIONS = True"
    if marker not in source:
        raise RuntimeError("upstream profile marker was not found")
    sys.path.insert(0, str(path.parent))
    namespace: dict[str, object] = {
        "__file__": str(path),
        "__name__": "tag_friendly_parameter_definitions",
    }
    exec(compile(source.split(marker, 1)[0], str(path), "exec"), namespace)
    namespace["QUANTUM"] = quantum
    return namespace


def construct(cls, **kwargs):
    # Estimator helpers can be verbose.  Keep machine-readable output stable.
    with contextlib.redirect_stdout(io.StringIO()):
        return cls(**kwargs)


def zk_row(kind: str, obj) -> dict[str, object]:
    return {
        "kind": kind,
        "n": obj.n,
        "d": obj.d,
        "m1": obj.m_1,
        "m2": obj.m_2,
        "q": int(obj.q),
        "q_bits": int(obj.q).bit_length(),
        "q_min": int(obj.q_min),
        "l": obj.l,
        "mlwe_classical_bits": float(obj.mlwe_coresvp_c),
        "mlwe_quantum_bits": float(obj.mlwe_coresvp_q),
        "msis_classical_bits": float(obj.msis_coresvp_c),
        "msis_quantum_bits": float(obj.msis_coresvp_q),
        "soundness_bits": -math.log2(float(obj.soundness_error)),
        "proof_kib": float(obj.proof_bitsize) / 8192.0,
        "crs_kib": float(obj.crs_size) / 8192.0,
    }


def zk_error(obj, quantum: bool) -> float:
    hardness = obj.mlwe_coresvp_q if quantum else obj.mlwe_coresvp_c
    return (
        obj.eps_1 / obj.M_1
        + obj.eps_2 / obj.M_2
        + obj.eps_3 / obj.M_3
        + 2.0 ** (-hardness)
    )


def soundness_error(kind: str, obj, quantum: bool) -> float:
    hardness = obj.msis_coresvp_q if quantum else obj.msis_coresvp_c
    first = obj.q_min ** (-obj.n / (2 if kind == "show" else obj.kappa))
    result = first + 2 / obj.challenge_space_size + 2.0 ** (-hardness)
    result += obj.q_min ** (-2 * obj.l)
    return float(result)


def composition(sig, issue, show, quantum: bool) -> dict[str, float]:
    e_zk_issue = zk_error(issue, quantum)
    e_zk_show = zk_error(show, quantum)
    e_anonymity = 2 * e_zk_show
    sig_user = sig.mlwe_U_coresvp_q if quantum else sig.mlwe_U_coresvp_c
    sig_misis = sig.misis_coresvp_q if quantum else sig.misis_coresvp_c
    e_unforgeability = soundness_error("show", show, quantum) + 2 * (
        e_zk_issue
        + e_zk_show
        + 2 * soundness_error("issue", issue, quantum)
        + 2.0 ** (-sig_user)
        + sig.tag_space_size * 2.0 ** (-sig_misis)
        + sig.forgery_advantage
    )
    return {
        "anonymity_bits": -math.log2(e_anonymity),
        "unforgeability_bits": -math.log2(e_unforgeability),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--upstream",
        type=Path,
        default=Path("/home/wkk/tag-friendly_sampler/scripts/parameters_anonymous_credentials.py"),
    )
    parser.add_argument("--target", type=float, default=128.0)
    parser.add_argument(
        "--protocol-sec",
        type=int,
        default=128,
        help="statistical/challenge security parameter passed to the upstream classes",
    )
    parser.add_argument("--issue-m2", default="56,60,64,68,72,76,80")
    parser.add_argument("--show-m2", default="74,78,82,86,90,94,98,102,106")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    ns = load_upstream(args.upstream.resolve(), quantum=True)
    SIG = ns["SIG_Parameters"]
    Issue = ns["Issue_ZKP_Parameters"]
    Show = ns["Show_ZKP_Parameters"]

    started = time.monotonic()
    ns["QUANTUM"] = False
    sig_classical = construct(
        SIG,
        target_bitsec=args.protocol_sec,
        n=256,
        d=4,
        b_H=13,
        k_H=2,
        b_L=49,
        k_L=2,
        m=10,
        no_guessing=False,
    )
    ns["QUANTUM"] = True
    sig = construct(
        SIG,
        target_bitsec=args.protocol_sec,
        n=256,
        d=4,
        b_H=13,
        k_H=2,
        b_L=49,
        k_L=2,
        m=10,
        no_guessing=False,
    )
    signature = {
        "n": sig.n,
        "d": sig.d,
        "q": int(sig.q),
        "signature_bytes": float(sig.sig_bitsize) / 8.0,
        "signer_mlwe_classical_bits": float(sig.mlwe_S_coresvp_c),
        "signer_mlwe_quantum_bits": float(sig.mlwe_S_coresvp_q),
        "user_mlwe_classical_bits": float(sig.mlwe_U_coresvp_c),
        "user_mlwe_quantum_bits": float(sig.mlwe_U_coresvp_q),
        "forgery_bits_with_classical_costs": float(sig_classical.security),
        "forgery_bits_with_quantum_costs": float(sig.security),
    }
    print(
        f"[search] signature finished in {time.monotonic() - started:.1f}s",
        file=sys.stderr,
        flush=True,
    )

    rows: list[dict[str, object]] = []
    objects: dict[str, list[tuple[object, dict[str, object]]]] = {
        "issue": [],
        "show": [],
    }
    for m2 in (int(x) for x in args.issue_m2.split(",") if x):
        point_started = time.monotonic()
        try:
            obj = construct(
                Issue,
                target_bitsec=args.protocol_sec,
                n=64,
                d=16,
                m_2=m2,
                q_start=math.ceil(2**47.7),
                n_attr=0,
                gamma=2**19,
                D=11,
                sig=sig,
                bimodal=True,
                compression=True,
                garbage=True,
            )
            row = zk_row("issue", obj)
            row["component_screen_pass"] = min(
                row["mlwe_quantum_bits"],
                row["msis_quantum_bits"],
                row["soundness_bits"],
            ) >= args.target
            rows.append(row)
            objects["issue"].append((obj, row))
        except Exception as exc:  # Preserve rejected points in the audit trail.
            rows.append({"kind": "issue", "m2": m2, "error": str(exc)})
        print(
            f"[search] issue m2={m2} finished in "
            f"{time.monotonic() - point_started:.1f}s",
            file=sys.stderr,
            flush=True,
        )

    for m2 in (int(x) for x in args.show_m2.split(",") if x):
        point_started = time.monotonic()
        try:
            obj = construct(
                Show,
                target_bitsec=args.protocol_sec,
                n=64,
                d=23,
                m_2=m2,
                q_start=2**75,
                n_attr=0,
                gamma=2**30,
                D=22,
                sig=sig,
                handle_bits=32,
                bimodal=False,
                compression=True,
                garbage=True,
                crt_lift=True,
            )
            row = zk_row("show", obj)
            row["component_screen_pass"] = min(
                row["mlwe_quantum_bits"],
                row["msis_quantum_bits"],
                row["soundness_bits"],
            ) >= args.target
            rows.append(row)
            objects["show"].append((obj, row))
        except Exception as exc:
            rows.append({"kind": "show", "m2": m2, "error": str(exc)})
        print(
            f"[search] show m2={m2} finished in "
            f"{time.monotonic() - point_started:.1f}s",
            file=sys.stderr,
            flush=True,
        )

    compositions: list[dict[str, object]] = []
    for issue_obj, issue_row in objects["issue"]:
        for show_obj, show_row in objects["show"]:
            classical = composition(
                sig_classical, issue_obj, show_obj, quantum=False
            )
            quantum = composition(sig, issue_obj, show_obj, quantum=True)
            compositions.append(
                {
                    "issue_m2": issue_row["m2"],
                    "show_m2": show_row["m2"],
                    "classical": classical,
                    "quantum": quantum,
                    "classical_screen_pass": min(classical.values()) >= args.target,
                    "quantum_screen_pass": min(quantum.values()) >= args.target,
                }
            )

    report = {
        "schema": "pqac-tag-friendly-component-screen-v1",
        "target_quantum_bits": args.target,
        "protocol_security_parameter": args.protocol_sec,
        "cost_model": ns["COST_MODEL"],
        "rough": bool(ns["ROUGH"]),
        "quantum_losses_enabled": True,
        "upstream": str(args.upstream.resolve()),
        "signature": signature,
        "zk_candidates": rows,
        "composition_candidates": compositions,
        "release_ready": False,
        "blockers": [
            "component estimates are not an end-to-end composition bound",
            "the signature profile must meet the selected end-to-end target",
            "the accumulator and same-handle reductions must be composed",
            "candidate constants must be generated, implemented, and regression-tested",
        ],
    }
    rendered = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
