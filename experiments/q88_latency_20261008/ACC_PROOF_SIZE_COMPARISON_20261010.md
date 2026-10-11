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
size change, so it is retained as the default optional optimization.  A second
two-show run produced 193.78 KiB for both proofs, with prove times 5.469/4.653 s
and verify times 2.562/2.555 s; both Pack verifications passed.

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

Changing the handle aggregation chunks from the default 16+8 split to three
rank-8 chunks did pass clear relation checking and Pack Verify, but produced the
same **194.03 KiB** envelope while increasing the one-show measurement to
6.236 s prove and 3.128 s verify. It is therefore rejected as a size
optimization.

A grouped-refresh prototype was then tested with six block equations packed
into one padded rank-8 extension relation per CRT lane. The existing Python
constraint API rejects this mixed degree-256 selector/degree-4096 witness shape
before relation checking. The prototype was removed from the active source;
the next implementation must add a native mixed-degree sparse helper in the
LaBRADOR bridge before this optimization can be evaluated.

## q88 bridge shape audit

The new `q88_bridge_shape.py` audit applies the Lazer-Acc layout to the actual
q88 profile without changing production code. It combines the three CRT
moduli into an 88-bit integer modulus, keeps the existing 50,000-deletion
witness gate, and uses radix 128 with groups of 12 polynomials. The initial
shape estimate is:

```text
public digits       13
witness digits       9  (58-bit conservative gate)
quotient digits     12  (79-bit preliminary bound)
carry layers        25
grouped vectors    168  (before refresh, handle, and link constraints)
```

This is the first q88-specific estimate in the same direction as the 126 KiB
Lazer-Acc construction. The quotient and carry bounds are deliberately
conservative estimates until the CRT-combined public polynomial compiler is
implemented; this is not yet a proof-size claim.

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

## 2026-10-10: exact CRT bridge KAT

`q88_bridge_core.py` now contains an engine-independent compiler for the first
real Lazer-Acc step: it CRT-lifts the three lane public polynomials, checks the
single 88-bit congruence, expands it in radix 128, and verifies every integer
quotient/carry layer before a Pack adapter can consume it.  The degree-8 KAT is
runnable with:

```bash
PYTHONPATH=experiments/q88_latency_20261008 \
  python3 experiments/q88_latency_20261008/q88_bridge_kat.py
```

This KAT passes and is not a size claim.  It does not yet include q88's hidden
handle products, refresh product, or q19 link constraints.  Those products
must be represented by bound auxiliary witnesses and their gadget relations;
otherwise a bridge would be a fixed-handle shortcut and would not preserve the
security statement.  The next integration gate is therefore a full q88
relation dump into this compiler, followed by a native single-statement Pack
proof and serialized-byte measurement.

### q88 relation decomposition used by the bridge

The wire relation has 22 truncation-right polynomials (`TRUNC_RIGHT=22`),
168 non-tag witness slots, and six refresh `z` polynomials.  A safe single
relation must preserve the following terms:

| term | bridge treatment |
| --- | --- |
| `w0 + A_trap*w1` | direct common witness columns after CRT lift |
| 22 truncation terms | direct columns plus their public radix scales |
| `B_{bit,digit}*right_{bit,digit}` | auxiliary aggregate witness plus a bound gadget equation |
| refresh `A_hat(r)*z_block` | auxiliary composed public polynomial plus selector binding |
| accumulator RHS and lane lifts | one CRT target and quotient/carry digits |

The last two rows cannot be replaced by prover-selected effective public
coefficients in a production proof.  That would reveal or fix the hidden
handle/selector.  The current KAT deliberately stops before those rows; its
purpose is to validate the exact integer bridge and quotient/carry arithmetic,
not to claim a secure q88 proof yet.

## 2026-10-10: serialized single-Pack measurement

The degree-4096 single-statement Pack KAT was run on `autodl-pqac8` with the
native serializer enabled.  It passed clear verification and Pack Verify:

```text
backend witness polynomials: 48
Pack Prove: 2.865 s
Pack Verify: 1.477 s
serialized single Pack proof: 157370 bytes = 153.68 KiB
```

This is a real serialized single Pack proof, but it is the existing arithmetic
KAT relation, not the complete q88 accumulator relation.  It is therefore a
lower-bound reference for the adapter and must not be reported as the q88
result.  The remaining integration work is to feed the 174 q88 witness columns
plus the bound handle/refresh auxiliary equations into the same native
statement, then rerun this exact serialization measurement.

## 2026-10-10: refreshed complete q88 baseline

The complete production q88 path was rerun on `autodl-pqac8` with one show:

```text
accumulator proof envelope: 194.53 KiB
complete show packet: 213.56 KiB
prove: 4.960 s
verify: 2.607 s
```

This is the comparison baseline for the full bridge.  The earlier 193.78 KiB
value is a lower sample from the same profile; proof framing is challenge and
parameter dependent, so the comparison must use the same serialized envelope
format and repeated-show protocol.

A two-show scalar-refresh rerun produced envelopes of **193.75 KiB** and
**193.78 KiB** (complete show packets 218.97 and 219.00 KiB), with prove times
4.820/4.626 s and verify times 2.574/2.429 s.  These are the current like-for-
like production comparison values.

## Full-security acceptance gate

A reduced q88 proof is accepted only if one serialized Pack proof contains all
of the following and passes both clear and Pack verification:

1. the three-lane accumulator relation after exact 88-bit CRT lifting;
2. radix-128 quotient and carry equations with pinned norm bounds;
3. handle gadget auxiliary witnesses and their hidden-bit binding equations;
4. refresh composed-polynomial auxiliary witnesses and selector binding;
5. accumulator lift reconstruction;
6. the q19 link relation and shared-handle binding;
7. native proof serialization and deserialization with the same statement hash.

Any proof that substitutes effective public handle/refresh coefficients, drops
one lane, or omits the link relation is an audit prototype only and must not be
compared with the 193.75/193.78 KiB production envelopes.

The isolated worker dump also measured the complete serialized production
q88 envelope directly: **198882 bytes = 194.220703 KiB**.  Pack Verify passed;
prove was 4.673 s and verify was 2.466 s.  This byte count is the exact
serialized baseline for the next full bridge proof.

## 2026-10-10: grouped native Pack bridge KAT

The ported grouped statement builder was executed on `autodl-pqac8` with a
degree-256 exact bridge KAT.  It produced one native Pack proof, passed Pack
Verify, and serialized successfully:

```text
serialized proof: 151066 bytes = 147.525 KiB
```

This validates the Lazer-Acc grouped statement construction and serializer
path.  It is not the complete q88 accumulator relation; the real q88 bridge
still needs the 4096-degree witness input and all auxiliary binding equations.

## 2026-10-11: real q88 main relation on autodl-pqa11

The new host alias `autodl-pqa11` was configured for the user's port 10657
endpoint and the authoritative local `OURS` plus modified lowercase `lazer`
trees were synchronized one way from WSL.  The real degree-4096 CRT bridge
initially failed native clear verification because raw carry coefficients can
reach about 41047, beyond the signed 16-bit witness limb accepted by
LaBRADOR.  Splitting each carry into three signed radix-128 witness layers
fixed that failure.

The resulting single Pack **main accumulator relation** passed on the new
host:

```text
statement vectors: 63648
constraints:       36
Pack Verify:       passed
serialized proof:   174869 bytes = 170.7705078125 KiB
```

A timed rerun on the same host measured Pack Prove **12.393 s**, Pack Verify
**8.736 s**, and serialization **174805 bytes = 170.7080078125 KiB**.  The
small byte variation is from randomized proof generation.  These timings are
for the main relation only and are not comparable as a complete q88 proof
until all binding constraints are included.

This is not yet the full q88 proof.  Handle, refresh, lift-reconstruction and
q19-link binding constraints still need to be appended to the same statement
before comparing it with the 193.75--194.22 KiB production envelope.

## 2026-10-11: complete CRT single-Pack acceptance run on autodl-pqa11

The complete dynamic worker payload (`ProverPayload || LinkPayload`) was decoded
and cross-checked before statement construction.  The production full-ring
quadratic builder then included the hidden-handle binding, refresh selector and
product binding, all three CRT lift reconstruction constraints, and the four
q19 link rows with shared-handle binding in one native statement.  Clear
witness checks passed.  Pack Prove, Pack Verify, serialization, deserialization,
and Pack Verify of the deserialized proof all passed:

```text
Pack Prove:             12.2195369601 s
Pack Verify:             7.6178674698 s
serialized proof:       180542 bytes = 176.310546875 KiB
deserialized Pack Verify: passed
```

The complete proof is therefore smaller than the 193.75--194.22 KiB production
envelope by 17.44--17.91 KiB (about 9.0--9.2%).  It does not yet meet the
4 s / 2 s latency target.  The 170.7 KiB result above remains a main-relation
bridge measurement only; the 176.31 KiB result is the first complete single-Pack
measurement with the full security constraint set.

## 2026-10-11: baseline time optimizations applied to complete Pack

The complete acceptance function now inherits the production worker's
relation-preserving tuning when called standalone: twelve OpenMP workers,
parallel witness/statement/LDD preparation, Gaussian and refresh worker pools,
witness/dynamic binding caches, lazy rotation scratch, and the production
`PACK_IMPROVE_THRESHOLD=0.50`. The transcript-sensitive Pack reducer remains
serial. On `autodl-pqa11`, the optimized run passed clear verification, both
Pack verifications, serialization, and deserialization:

```text
Pack Prove:             5.533 s
Pack Verify:             2.947 s
serialized proof:       198338 bytes = 193.69 KiB
deserialized Pack Verify: passed
```

This is a substantial improvement over the uncached direct invocation
(12.220/7.618 s), but it does not reach 3/1.5 s. Setting reducer threads to 2
was slower (5.805/3.106 s); disabling Pack splitting failed with “No LNP params
found, ZK not possible”. Direct right-side quadratic substitution remained
valid but regressed to 14.880/7.702 s. Four-limb compact witness mode failed
the native norm check, so none of these unsafe or slower variants is enabled.

Against the production baseline, this should not be reported as a performance
win: the 193.69 KiB sample is effectively the same size as 193.75--194.22 KiB,
while 5.533/2.947 s is slower than the baseline 4.626--4.960/2.429--2.782 s.
The current result is an integration acceptance point only. A real speed or
size claim now requires a lower-dimensional valid witness/Pack layout, or the
resident online worker path; changing thread counts and fold thresholds alone
is insufficient.
