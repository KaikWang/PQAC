# Complete CRT single-Pack timing research

## Findings

The LaBRADOR/LaZer implementation is already using the main practical
optimizations from the baseline: resident parameters, parallel witness and LDD
preparation, aggregated quadratic constraints, and a serial transcript-sensitive
Pack reducer. The resident-cache implementation in `q88_bridge_core.py` removes
statement construction and parameter generation, but the second run still
measures about 6.14 s prove and 3.32 s verify. The remaining cost is therefore
inside Pack reduction and verification.

The LaZer paper presents composable lattice proofs and the public artifact
documents that its benchmark timings are for a single-core implementation;
this supports treating thread tuning as an engineering factor rather than a new
proof-size or asymptotic improvement:

- https://doi.org/10.1145/3658644.3690330
- https://github.com/lazer-crypto/lazer

The accumulator paper uses binary gadget decompositions and powers-of-two
packing. Its useful transferable idea is to make one bounded gadget relation
carry many coefficients, but replacing the current q88 quadratic relation with
that construction would require a new security reduction:

- https://link.springer.com/article/10.1007/s00145-023-09470-6

Recent lattice folding work reports that packing many field constraints into a
single ring constraint is the route to lower verification cost, but it changes
the proof system and parameter selection. It is not a drop-in optimization for
the current LaBRADOR statement:

- https://eprint.iacr.org/2025/247.pdf

## Implementation decision

The safe immediate optimization is resident statement/parameter reuse. The
next viable change is a quotient/carry bridge that shares the existing handle
gadget witness across the three CRT lanes. The current full-ring compiler
duplicates lane-specific aggregate products; direct substitution removes those
auxiliary vectors but was measured slower and has a larger quadratic term.
Combining the lanes requires new exact carry witnesses and a norm proof, so it
must be implemented as a separate statement profile and accepted only after
clear verification, Pack verification, serialization, deserialization, and
parameter-bound checks all pass.

Changing only thread counts, Pack thresholds, or disabling splitting has not
reached 3 s / 1.5 s and will not be enabled as a claimed optimization.
