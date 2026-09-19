#!/usr/bin/env python3
"""Runnable LaZer arithmetic smoke test for a shared-x constraint.

This remains a small standalone proof-system integration test.  The production
C path now uses a randomized shared-handle commitment in both its real
signature-show and accumulator MemVer relations.  This script only proves
knowledge of one hidden vector `s` satisfying

    A * s + t = 0

where column 0 of `s` is treated as the hidden handle bits `x`.  Rows 0..1 are
the stand-in for the credential/signature relation and rows 2..3 are the
stand-in for the accumulator membership relation, so the same hidden `x` is
used by both blocks.
"""

from __future__ import annotations

import pathlib
import sys


LAZER_PYTHON = pathlib.Path("/home/wkk/pq_ac/lazer/python")
LAZER_DEMO = LAZER_PYTHON / "demo"

sys.path.append(str(LAZER_PYTHON))
sys.path.append(str(LAZER_DEMO))

from lazer import (  # noqa: E402
    VerificationError,
    lin_prover_state_t,
    lin_verifier_state_t,
    polymat_t,
    polyring_t,
    polyvec_t,
    print_stopwatch_lnp_prover_prove,
    print_stopwatch_lnp_verifier_verify,
)
from _demo_params_cffi import lib  # noqa: E402
from demo_params import deg, dim, mod  # noqa: E402


def main() -> int:
    seed = b"\0" * 32
    m, n = dim[0], dim[1]

    print("Revocable AC LaZer same-x prototype")
    print("===================================")
    print()
    print("Relation:")
    print("- witness s[0]     : hidden handle bits x")
    print("- witness s[1..3]  : stand-in credential/signature secrets")
    print("- witness s[4..7]  : stand-in accumulator witness pieces")
    print("- rows 0..1        : credential/signature relation block")
    print("- rows 2..3        : accumulator membership relation block")
    print()

    Rp = polyring_t(deg, mod)
    A = polymat_t(Rp, m, n)
    A.urandom(mod, seed, 17)

    s = polyvec_t(Rp, n)
    s.brandom(1, seed, 23)
    t = -A * s

    prover = lin_prover_state_t(seed, lib.get_params("param"))
    verifier = lin_verifier_state_t(seed, lib.get_params("param"))

    prover.set_statement(A, t)
    prover.set_witness(s)

    print("generate proof ...")
    proof = prover.prove()
    print_stopwatch_lnp_prover_prove(0)

    verifier.set_statement(A, t)
    print("verify proof ...")
    try:
        verifier.verify(proof)
    except VerificationError:
        print("reject")
        return 1

    print("accept")
    print_stopwatch_lnp_verifier_verify(0)
    print()
    print(f"proof size: {len(proof) / 1024.0:.2f} KiB ({len(proof)} bytes)")
    print()
    print("Status:")
    print("- LaZer proof generation and verification are wired into this repo.")
    print("- The same hidden x lane is shared by the signature and accumulator blocks.")
    print("- Production status: the real show and MemVer blocks are linked through a shared randomized commitment.")
    print("- This script remains an isolated LaZer arithmetic smoke test, not the production packet path.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
