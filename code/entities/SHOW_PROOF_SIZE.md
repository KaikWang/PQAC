# Show proof size and parameter audit

Before the same-handle binding was added, this tree produced a compact
anonymous-signature show proof of about 58,820 bytes (57.44 KiB).  The current
linked relation uses the previously removed 32 response polynomials to open a
randomized handle commitment.  One local run measured 66,624 bytes
(65.06 KiB); `test/revac_proof_codec_demo.c` enforces a separate 70 KiB
post-link regression budget and prints a per-field breakdown.

## Fixed causes

- The pre-link show witness dimension was hard-coded to 183 although that
  relation ended after 151 polynomials.  The 32 unused Gaussian responses were
  removed at that point.  They are now deliberately restored as two binary
  `R_q^d` openings for the shared-handle commitment.
- The verifier reconstructed the Fiat--Shamir challenge but did not compare it
  with the proof challenge.
- The verifier hashed a prover-supplied `t0_stored` instead of the `t0`
  reconstructed from the statement and responses.
- The redundant `t1_gadget` compensation did not belong in the quadratic
  relation.
- `z2` norms were truncated to 64 bits, and unsigned norm subtraction in the
  rejection probability could underflow.
- Uniform and Huffman decoders accepted non-canonical or malformed encodings.
- The proof serialized a full ABDLOP commitment/opening.  It now uses
  `Power2Round(tA)`, generalized high-bit hints, and transmits only the
  `m2-d` opening head.  `A2` is sampled as `[A2' | I_d]`; the unused tail
  columns of `Byg` and `b` are zero.

Typical current signature field sizes are:

| Field | Bytes |
|---|---:|
| `tA` high part | 6,624 |
| generalized hints | 368 |
| `tB` | 3,712 |
| `z3` | about 905 |
| `h` | 1,856 |
| `t1` | 464 |
| challenge and counter | 44 |
| `z1` | about 44,620 |
| transmitted `z2` head | about 8,010 |
| total | about 66,624 |

Gaussian fields vary slightly between runs.

## Full revocable show packet

The outer `RVACP2` packet also transmits the canonical 2,432-byte randomized
handle commitment.  Therefore the signature-only packet is currently about
69,110 bytes (67.49 KiB), including nonce and length/type fields.

Before the same-handle relation was added, one online Labrador run produced:

| Component | Bytes |
|---|---:|
| compact anonymous-signature show proof | 58,826 |
| raw serialized Labrador proof | 342,774 |
| deflated Labrador payload in the wire packet | 181,933 |
| outer nonce, lengths, kinds, and magic | 54 |
| complete revocable show wire packet | 240,813 |

That 235.17 KiB value is a historical pre-link baseline, not a current result.
The current Labrador witness adds 24 fixed witness vectors/polynomials and four
linear constraints, so the complete packet must be remeasured on the native
server.  Both Gaussian Huffman output and generated Labrador proof parameters
may vary between runs.

## Why this is not the paper's approximately 50 KiB profile

Running
`scripts/parameters_anonymous_credentials.py` in the tag-friendly-sampler
repository with the pre-link dimensions (`n=64`, `d=23`, `m1=151`, `m2=74`)
and the implementation's CRT lift gives about 61.37 KiB for ordinary
rejection sampling.  The smaller paper profile uses bimodal rejection sampling
and a proved hidden sign bit.  With the same implementation dimensions and CRT
lift, the corrected script reports 51.49 KiB for that profile (`m1=152`, about
a 65-bit modulus).  Both profiles use `l=4`, derived from the fixed
`q_signature=405769` CRT factor and the implemented garbage-term optimization;
the previous `l=7` was another unrecomputed placeholder.

Those estimates do not yet cover the current `m1=183` linked relation.  Its
Gaussian widths, response bounds, modulus, and concrete security must be
rerun together before the new relation is used for a paper security claim.
The corrected pre-link ordinary-rejection model selects a roughly 75-bit proof modulus.
The present arithmetic uses FLINT `nmod_poly`, signed 64-bit coefficients, and
uniform samplers that explicitly reject moduli above 60 bits.  Consequently,
both the approximately 65-bit bimodal profile and the larger
ordinary-rejection profile require a different backend; changing only codec
constants or Gaussian widths cannot safely produce either profile.

The estimator currently reports about 100.2 bits for the ordinary CRT profile.
For the 51.49 KiB bimodal CRT profile it reports about 118.6 bits.  These
figures are useful for comparing profiles but are not a 128-bit security claim;
the final parameter search must also retune the dimensions until every required
concrete-security bound is met.

Reaching that profile requires a protocol/backend migration:

1. implement bimodal rejection sampling and add its hidden sign witness and
   quadratic constraint;
2. replace the `nmod_poly` show arithmetic with a greater-than-64-bit or RNS
   backend, including canonical serialization and sampling;
3. regenerate `q`, `gamma`, `D`, `d`, `m2`, `l`, Gaussian widths, norm bounds,
   and Huffman tables together;
4. rerun concrete M-SIS/M-LWE estimates and all transcript, tamper, codec, and
   end-to-end tests before claiming the target security level or size.

For comparison, set `TSAMPLER_SHOW_BIMODAL=1` when running the parameter script
to report the protocol variant rather than the current C sampler.
