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
