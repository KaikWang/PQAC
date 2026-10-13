# Isolated single-88-bit profile

This profile is independent of the production three-prime CRT parameters. It
uses the certified 88-bit prime

`309333912899981927159480321`, with `p = 1 mod 8192`, and keeps the q88
degree-4096 witness inventory: 174 public/witness columns. The carry bridge is
constructed over the integers with radix 128; native Pack uses the existing
38-bit proof field.

## Measured on `autodl-pqa11`

The source witness is a captured CRT show only for the bounded witness shape.
The single-prime target is derived afresh in the new field, so this is not a
conversion proof for the existing CRT accumulator.

| Pack grouping | bridge | statement | Prove | Verify | serialized proof |
|---:|---:|---:|---:|---:|---:|
| 12 | 13.115 s | 31.611 s | 7.551 s | 5.152 s | 163459 B (159.628 KiB) |
| 16 | 13.307 s | 31.803 s | 6.656 s | 4.999 s | 163235 B (159.409 KiB) |
| 24 | 14.138 s | 32.065 s | 6.749 s | 4.802 s | 163651 B (159.815 KiB) |

All rows passed native Clear Verify, Pack Verify, proof serialization,
deserialization, and Pack Verify after reload. Grouping 16 is the default in
`q88_single88_profile.py` because it was the fastest proving row while still
passing the exact norm gate.

## Security boundary

The current adapter still uses the effective hidden-handle and refresh
coefficients to define the independent profile target. It is therefore a
full-shape performance baseline, not the requested security-complete single
88-bit q88 proof. The next implementation gate is to replace those effective
coefficients with auxiliary product witnesses and append their binary,
quotient/carry, CRT-lift, and q19 same-handle bindings to the same native
statement.

## Handle-bound milestone

The `--handle-bound` mode now removes the effective hidden-handle shortcut.
It registers six unscaled handle product witnesses, adds 24 binary
constraints, and adds six quadratic product bindings in the same statement.
The refresh selector is still effective in this milestone.

Measured on `autodl-pqa11` with the exact serialized reload gate:

`bridge=13.656 s`, `statement=34.291 s`, `Prove=13.576 s`,
`Verify=9.651 s`, `serialized=169228 B = 165.262 KiB`.

Clear Verify, Pack Verify, serialization, deserialization, and Pack Verify
after reload all passed. The q19 part is now also available in `--handle-bound`:
it adds degree-256 randomness and lift witnesses, four exact q19 linear rows,
24 degree-4096 to degree-256 handle integer links, and binary/degree-zero
checks for the small handle representation. A fresh run passed all gates with:

`bridge=14.011 s`, `statement=39.761 s`, `Prove=12.246 s`,
`Verify=9.427 s`, `serialized=169196 B = 165.230 KiB`.

Refresh selector binding remains the next missing security component, so this
is still an intermediate rather than the final secure single-prime profile.

## Refresh quotient/carry KAT

`q88_refresh_binding_kat.py` now proves one complete refresh block using the
production structure: 128 selector witnesses, 13 radix-128 centered refresh
digits, one 88-bit modulus lift, and 128 binary/degree-zero constraints. On
`autodl-pqa11` it passed all native and serialization gates with:

`bridge=8.824 s`, `statement=26.293 s`, `Prove=9.462 s`,
`Verify=8.153 s`, `serialized=169772 B = 165.793 KiB`.

This is a standalone one-block KAT. The six blocks still need to share their
witness allocation and be appended to the handle plus q19 statement before
the profile can be called complete.

## Shared-selector six-block Pack milestone

The six refresh quotient/carry bridges now use one shared native selector
allocation. The prior implementation duplicated 128 selector vectors per
block and tied copies with integer equations. The new implementation aliases
all six blocks to the same 128 degree-4096 binary vectors, so the selector
binary/degree-zero checks are emitted once while each block retains its exact
quotient/carry equations.

A fresh run on `autodl-pqa11` passed Clear Verify, Pack Verify, serialization,
deserialization, and Pack Verify after reload:

`bridge=54.605 s`, `statement=158.216 s`, `Prove=41.790 s`,
`Verify=33.440 s`, `serialized=169548 B = 165.574 KiB`.

Compared with the duplicated-selector run (`174741 B`, `46.801 s`,
`35.791 s`), this saves 5193 bytes and about 5.0 s proving / 2.4 s verifying.
This is still a refresh binding milestone: the main accumulator relation's
refresh product terms have not yet been replaced by these selector bridges.

## Compact layout and joint statement (2026-10-13)

The previous refresh adapter decomposed binary selectors and already
radix-128-decomposed refresh digits into nine further layers. The compact
adapter uses one layer for these values. It removes the zero-modulo-p lift
column; the compiler's exact p quotient and carry chain remain. Three signed
quotient digits suffice for this selector relation: with scalar binary
selectors, public coefficients < p and 13 bounded refresh digits, each
coefficient of the quotient is bounded by 128 + 13*64 + 1 = 961.
The legacy layout remains available with `--legacy-layout`.

The following intermediate measurements (before explicit terminal-carry
hardening) passed native Clear Verify, Pack Verify,
serialization, deserialization and Pack Verify after reload on
`autodl-pqa11`:

| Profile | witness vectors | polynomial constraints | Prove | Verify | serialized proof |
|---|---:|---:|---:|---:|---:|
| Shared selectors, legacy repeat | 7472 | 302 | 42.229 s | 36.050 s | 169196 B (165.230 KiB) |
| Compact six refresh bridges | 416 | 224 | 6.086 s | 4.978 s | 163907 B (160.065 KiB) |
| Compact refresh plus 78 product bindings | 656 | 302 | 7.237 s | 5.589 s | 169804 B (165.824 KiB) |
| Joint main + handle + q19 + refresh | 1020 | 357 | 9.007 s | 7.282 s | 175701 B (171.583 KiB) |

The product and joint profiles use fixed public radix norm budgets
`n*64^2`, a refresh-z norm budget `n`, and fixed carry capacities (two
radix layers for refresh bridges; four in the initial main bridge,
three in the hardened version below). The main raw
witness uses three layers for its first two columns, two for truncation,
and one for handle-right and handle-product columns. Capacity checks reject
samples that do not fit this independent profile.

`q88_refresh_products.py` binds each of the 13 refresh digits times z using
a degree-4096 native quadratic equation and a three-layer product
reconstruction. The native equality cannot wrap in the Pack field under
these norm budgets: the product coefficient is bounded by `4096*64`, and
the reconstructed coefficient by `64*64*(1+128+128^2)`; their sum is
`67899392 < 2^38-107`.

`q88_single88_joint_profile.py` removes the six effective refresh columns
from the main relation and replaces them with public radix powers times
these bound product layers. Those main columns alias the already registered
product vectors directly. There are 176 scalar/degree-zero constraints,
24 q19 shared-handle integer links, and one Pack proof, rather than separate
main and refresh proofs. Bridge construction took 59.178 s and statement
construction (including Clear Verify) took 28.268 s; Prove/Verify timings
above exclude those preparation costs.

### Hardened joint gate

The final run reduced the main carry capacity to three radix layers, added
one explicit terminal-carry-zero equation per bridge, and checked a public
Cauchy-Schwarz bound for each integer bridge equation before native
construction. Any bridge whose equation coefficients could wrap in the
38-bit field is rejected. These checks cover the bridge equations; they do
not replace a full audit of q19 parameters or the production beta bound.

`--negative-check` alters one native coefficient of the selector, refresh
product and q19 handle witness in turn. Native Clear Verify rejected all
three cases and passed after each restoration. The restored witness then
passed Pack Prove, Pack Verify, serialization, deserialization and reload
Pack Verify. This was measured on `autodl-pqa11` as:

| witness vectors | polynomial constraints | degree-zero | integer links | Prove | Verify | serialized proof |
|---:|---:|---:|---:|---:|---:|---:|
| 999 | 364 | 176 | 24 | 9.863 s | 7.477 s | 175701 B (171.583 KiB) |

Preparation remains separate: bridge construction took 59.507 s, and
statement construction including Clear Verify took 32.316 s. Negative checks,
serialization and reload verification are excluded from Prove/Verify timing.
The 3 s / 1.5 s latency target has not been reached.

### Remaining security and comparison boundary

This is a joined independent single-prime profile with verified bindings,
not a certified replacement for the CRT production scheme. Its target is
derived afresh from a captured CRT witness in the new field; the original
production beta bound is not automatically implied by per-digit budgets.
The new field's accumulator/trapdoor parameters, aggregate norm bound,
and public-input/target distribution still require a separate audit before
a production claim or a like-for-like security comparison. The measured
171.583 KiB body is smaller than the historical 193.75-194.22 KiB CRT body,
but this does not establish equivalent security or improved latency.

Reproduce the joined profile with:

```sh
cd /root/autodl-tmp/pq_ac
PYTHONPATH=OURS/experiments/q88_latency_20261008:OURS/code/entities/zk_lazer:lazer/python:Lazer-Acc/python \
LD_LIBRARY_PATH=lazer /root/miniconda3/envs/sage/bin/python -u \
OURS/experiments/q88_latency_20261008/q88_single88_joint_profile.py --negative-check
```
