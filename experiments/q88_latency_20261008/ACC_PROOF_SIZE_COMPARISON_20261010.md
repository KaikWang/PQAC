# q88 accumulator proof size versus the Lazer-Acc grouped bridge

## Measured sizes

The current q88 production envelope is **193.88--194.53 KiB** in the
five-show/default measurements.  The complete compressed show packet is about
213--219 KiB.  The envelope includes the compact profile framing and the
serialized native LaBRADOR Pack proof; it is the object sent as the accumulator
proof by the current worker.

The `Latice-Accumulator` grouped bridge reports **126.44 KB**, but that number
is an in-process Pack-size estimate.  Its Python binding does not serialize
that proof.  It is also a different 100-bit accumulator relation, with radix
128 quotient/carry digits and 340 witness polynomials grouped into 12-polynomial
vectors.  The 34-block LNP proof in the same repository is about 27.42 MiB and
must not be compared with the 126 KB single-Pack estimate.

## Why q88 is larger

The q88 degree-4096 compiler still accounts for the Pack witness in degree-256
backend units.  With the current four-limb proof-only encoding, the default
profile has approximately:

| source | degree-256-equivalent vectors |
| --- | ---: |
| 24-bit handle binding | 384 + 24 scalar vectors |
| 128-bit refresh selector | 2,048 |
| 168 non-refresh witness polynomials × 4 limbs × 16 phases | 10,752 |
| refresh z and A/carry vectors | 1,248 |
| quotient/lift vectors | 336 |
| right-side quadratic aggregation auxiliaries | 4,608 |
| **total** | **19,400** |

The exact Pack proof size is not a linear function of this count, because
Dachshund parameter selection, norm bounds, challenge fields, and proof
encoding also contribute.  Therefore removing vectors alone does not imply a
126 KiB proof.

## Experiments

`REVAC_Q96_SCALAR_REFRESH_BITS=1` is relation-preserving: it keeps the hidden
refresh selector as degree-zero data.  One complete remote show passed with
4.716 s prove, 2.645 s verify, and **193.75 KiB** envelope.  This is a small
size change, so it is retained only as a valid optional optimization.

`REVAC_Q96_DIRECT_RIGHT_QUAD=1` substitutes the right-side aggregation equality
directly into the quadratic relation and removes 4,608 auxiliary vectors in
the compiler model.  The current q88 security/norm profile cannot generate a
valid Dachshund parameter set for this layout (`no valid split/merge of the
witness was found`), so it is **not** enabled.  Enabling it without deriving a
new valid norm/split bound would be unsafe.

Enabling both switches fails at parameter generation for the same reason.

The existing three-limb compact experiment was also tested with its exact
capacity setting (`radix=1048575`, which covers the complete window gate).  It
reached statement construction but Dachshund rejected the resulting norm
profile with `no valid split/merge of the witness was found`.  The four-limb
profile therefore remains the smallest currently valid q88 witness encoding.

The native proof profile shows that the largest serialized component is a
Labrador subproof entry (`lab10`, about 65.7 KiB in the current profile),
followed by the other Lab entries and Pack challenge vectors. Changing only
the outer envelope or removing metadata therefore cannot provide a meaningful
reduction; the next reduction must lower valid Pack parameter dimensions or
replace the current CRT/lift layout with a compact quotient/carry statement.

## Conclusion and next implementation boundary

The Lazer-Acc ideas are applicable, but the safe transplant is not a simple
radix change.  The next useful prototype is a q88-specific full-ring bridge
that keeps the native q88 88-bit CRT relation, introduces exact quotient/carry
digits only where the lift crosses the Pack proof prime, and then measures a
real serialized proof.  It must pass: clear relation verification, Pack
verification, exact signed-lift reconstruction, the existing norm gate, and a
wire-size regression check.  Until that adapter exists, **194 KiB is the
verified q88 accumulator-proof size**; 126 KB is a promising target estimate,
not a current q88 result.
