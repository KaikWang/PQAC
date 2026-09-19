#!/usr/bin/env python3
"""Labrador prototype for the core revocable anonymous show relation.

This demo models the important algebraic shape of the final proof:

  1. the credential/signature side contains the hidden handle bits x;
  2. the accumulator membership side uses the same x;
  3. membership includes the quadratic selector terms x_i * right_i;
  4. x_i is constrained to be a binary degree-0 witness;
  5. the remaining witness pieces have norm bounds.

It is still a small proof-language prototype.  The public coefficients are toy
constants, not the full C matrices from AccGenpp.  The shape of the second
constraint is the one we need for MemVer:

    left + B0*r0 + B1*r1 - G0*(x0*r0) - G1*(x1*r1) = acc.

The next step is mechanical but large: emit these constraints for all
PARAM_D * PARAM_N coefficients and all ACC_GENPP_HANDLE_BITS *
ACC_GENPP_GADGET_DIM right vectors, with lifted modulus handling.
"""

from __future__ import annotations

import pathlib
import sys
import argparse


LAZER_PYTHON = pathlib.Path("/home/wkk/pq_ac/lazer/python")
sys.path.append(str(LAZER_PYTHON))

from lazer import polyring_t, poly_t  # noqa: E402
from labrados import pack_verify, proof_statement  # noqa: E402
from revac_profile_generated import (  # noqa: E402
    GADGET_BASE,
    N,
    PROOF_MODULUS,
)


RINGS = {
    "32": polyring_t(N, 2**32 - 99),
    "36": polyring_t(N, 2**36 - 243),
    "38": polyring_t(N, PROOF_MODULUS),
}


def cpoly(ring: polyring_t, value: int) -> poly_t:
    return poly_t(ring, [value] + [0] * (ring.deg - 1))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primesize", choices=sorted(RINGS), default="38")
    parser.add_argument("--no-zk", action="store_true", help="debug relation without zero-knowledge packing")
    args = parser.parse_args()

    ring = RINGS[args.primesize]
    one = cpoly(ring, 1)

    # Hidden handle bits.  In the real proof this is 32 degree-0 binary polys.
    x0 = cpoly(ring, 1)
    x1 = cpoly(ring, 0)

    # Toy signature and accumulator witnesses.
    sig_secret = cpoly(ring, 7)
    left = cpoly(ring, 5)
    right0 = cpoly(ring, 3)
    right1 = cpoly(ring, -2)

    # Public toy matrices/gadget constants.
    dx0 = cpoly(ring, 11)
    dx1 = cpoly(ring, 17)
    b0 = cpoly(ring, 19)
    b1 = cpoly(ring, -23)
    # First two powers of the production balanced base-3 gadget.
    g0 = 1
    g1 = GADGET_BASE

    # Public statement values.
    sig_public = sig_secret + dx0 * x0 + dx1 * x1
    acc_public = left + b0 * right0 + b1 * right1 - cpoly(ring, g0) * x0 * right0 - cpoly(ring, g1) * x1 * right1

    print("Revocable AC Labrador quadratic prototype")
    print("=========================================")
    print()
    print("Constraints:")
    print("- signature block: sig_secret + Dx0*x0 + Dx1*x1 = sig_public")
    print("- accumulator block: left + B0*r0 + B1*r1 - G0*(x0*r0) - G1*(x1*r1) = acc")
    print("- x0,x1 are binary and degree-0")
    print(f"- Labrador prime size: {args.primesize}")
    print(f"- zero knowledge packing: {'off' if args.no_zk else 'on'}")
    print()

    deg_list = [ring.deg] * 6
    num_pols = [1] * 6
    # 0 means Labrador's binary norm type.  The other bounds are deliberately
    # loose for this proof-language integration test.
    norm_bounds = [0, 0, 2**20, 2**20, 2**20, 2**20]

    ps = proof_statement(
        deg_list,
        num_pols,
        norm_bounds,
        num_constraints=2,
        primesize=args.primesize,
        num_deg0_constraints=2,
        zk=not args.no_zk,
    )

    idx_x0 = ps.append_witness(x0)
    idx_x1 = ps.append_witness(x1)
    idx_sig_secret = ps.append_witness(sig_secret)
    idx_left = ps.append_witness(left)
    idx_right0 = ps.append_witness(right0)
    idx_right1 = ps.append_witness(right1)

    ps.append_deg0_statement(idx_x0)
    ps.append_deg0_statement(idx_x1)

    ps.append_statement(
        [dx0, dx1, one],
        [idx_x0, idx_x1, idx_sig_secret],
        sig_public,
    )
    ps.append_quadratic_statement(
        [one, b0, b1],
        [idx_left, idx_right0, idx_right1],
        [-g0, -g1],
        [idx_x0, idx_x1],
        [idx_right0, idx_right1],
        acc_public,
    )

    statement = ps.output_statement()
    ps.smpl_verify()
    proof = ps.pack_prove()
    if proof[0] != 0:
        print(f"pack prove failed with error {proof[0]}")
        return 1
    if not pack_verify(proof[1:3], statement, args.primesize):
        return 1

    print()
    print("Status:")
    print("- Labrador quadratic proof generation and verification succeeded.")
    print("- The same hidden binary handle bits are used by both blocks.")
    print("- This is the exact proof shape needed for the x_i * right_i terms in MemVer.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
