# Revocable Anonymous Credential Entity Layout

This directory separates the revocable anonymous credential workflow by entity.

- `common/`: shared public types and helpers, including `handle = H(PK, ID)`.
- `ta/`: trusted authority logic: setup, UAV/GCS registration, credential issuance, revocation, and BML maintenance.
- `uav/`: UAV-side logic: key generation, update application, and missed-broadcast recovery.
- `gcs/`: ground control station logic: key generation and credential/witness verification.
- `zk_lazer/`: LaZer/Labrador prototype for the future revocable anonymous show proof.

For deployment experiments, copy the entity directory you need plus `common/`.
For example, the UAV-side portable package is:

- `uav/`: UAV-side wrapper and demo build file.
- `common/`: shared revocable-credential helpers.
- `common/core/include/`: copied cryptographic headers.
- `common/core/src/`: copied cryptographic core sources.
- `common/core/flint_headers/`: copied FLINT headers.
- `common/core/lib/libflint.so`: bundled FLINT dynamic library.

Build the UAV package after copying `uav/` and `common/` together:

```sh
cmake -S uav -B /tmp/revac_uav_build
cmake --build /tmp/revac_uav_build -j2
/tmp/revac_uav_build/uav_demo
```

Build the full standalone entity package after copying `entities/`:

```sh
cmake -S entities -B /tmp/revac_entities_build
cmake --build /tmp/revac_entities_build -j2
/tmp/revac_entities_build/revocable_ac_demo
```

The production regression in `test/revocable_ac_demo.c` uses the q88 candidate
D accumulator and covers
Setup, Issue, full anonymous show, GCS verification, signed BML deletion
updates, and authenticated TA-side MemWitSync recovery.  Production
applications call `revac_ta_setup()`, which selects q88 candidate D.  Historical
`acc96_*`, `REVAC_ACC_PROFILE_Q96`, and q96 environment-variable names remain
as compatibility identifiers; the active modulus, dimensions, encodings, and
`RV88*` wire magics are q88.  The older `Plan_process.pdf` q19 regression has moved to
`test/revocable_ac_q19_legacy_demo.c`; q19 experiments must request
`revac_ta_setup_q19()` explicitly.  `revac_q96_entities_demo` remains as a
temporary executable alias for scripts written during the earlier q96
integration.
TA setup, UAV/GCS registration, handle derivation, UAV credential issuance,
GCS certificate gating, GCS verification, revocation by `Delete`, direct
broadcast of an authenticated tail-context profile, context-assisted
`MemWitUp`, offline context recovery through BML, fallback resynchronization
by `MemWitSync`, and final verification.

## q88 candidate D accumulator

The accumulator uses the three NTT-prime CRT factors
`1073692673`, `536813569`, and `536690689`, whose product is the 88-bit modulus
`309333912899981927159406593`.  Its inverse maps coefficients to six balanced
base-26007 digits in `[-13003, 13003]`:

```text
v = sum(j=0..5) digit[j] * 26007^j  (mod q)
```

The revocation handle is 24 bits, the truncated SampleLeft gadget uses
`base=13, k=24, l=2`, and a window accepts at most 10000 deletions before a
fresh setup/refresh.  The parameter audit gives the conditional arithmetic
screen `beta_1 < q` up to 56245 deletions and margin
`log2(q/beta_1)=2.492` at 10000 deletions.  This is an arithmetic condition,
not by itself a complete security proof.

Canonical compact size projections are 44 KiB for the accumulator,
7.219 MiB for the transmitted full witness, and 6.188 MiB for `Ba`.  Under the
conservative zigzag bit-width projection, one four-layer DELTA4 is 516 KiB and
an accumulated 10000-update DELTA4 is 684 KiB.  These DELTA4 projections are
not the same object as the signed T4-prefetch packet reported by the current
demo.

Existing APIs remain unchanged, but persisted q19/q96 public parameters,
witnesses, contexts, and proofs must be regenerated.  The compatibility names
do not imply wire compatibility.

### Lossless DELTA4 block codec

The q88 implementation also contains an optional lossless DELTA4 codec in
`common/core/include/acc96_delta4_codec.h`.  It does not change the algebraic
update: each centered coefficient is ZigZag encoded and packed with the
minimum width of its independent 256-coefficient block.  The decoder checks
the modulus range, padding bits, descriptor count and total length, and only
installs the result after the complete payload has decoded successfully.

For the measured q88 candidate-D trajectory, the four-layer DELTA4 changes
from 444.00 KiB with a single fixed 37-bit width to 432.66 KiB with block
widths, a 2.56% reduction.  Prefix packets are 108.16, 216.34, 324.53 and
432.66 KiB for one through four layers.  A receiver whose first differing
handle bit is `j` needs only the first `j+1` layers; under uniform 24-bit
handles the expected prefix payload is 175.77 KiB over all survivors (187.49
KiB conditioned on the DELTA4-covered 15/16 of survivors).  This is a
receiver-specific/unicast or layered-multicast optimization, not a claim that
one single full broadcast is 175.77 KiB.

The codec was tested against the exact relation update, including all prefix
lengths and rejection of a forged width descriptor without mutating the output.
The DELTA4 sample had no zero coefficients, and DEFLATE reduced the packed
payload by effectively 0%; generic compression is therefore not enabled.

### Authenticated selective DELTA6

The production q88 entity path extends the same exact prefix relation to six
layers.  TA first sends a signed 8.62 KiB manifest containing the deletion
sequence, revoked handle, post-delete accumulator digest, parameter digest,
and the hash and length of each independently encoded layer.  Only after
verifying this manifest does a UAV compute
`j = ctz(member_handle XOR revoked_handle)`:

- `j < 6`: fetch packets `0..j`, verify their sequence, order, hash and
  canonical codec, then atomically add the decoded prefix;
- `j >= 6`: fetch no layer packet and run the exact LOCAL `EvalFx` recurrence.

The six measured layer packets are 108.19--108.25 KiB each.  For independent
uniform 24-bit handles, DELTA6 therefore covers `63/64 = 98.4375%` of
survivors, versus `15/16 = 93.75%` for DELTA4.  The expected number of layer
packets per survivor is

```text
sum(k=1..6) k / 2^k = 1.875,
```

or about 202.9 KiB of layer data plus the 8.62 KiB manifest.  Conditioned on a
DELTA6 hit, the corresponding layer expectation is about 206.1 KiB.  A UAV on
the 1/64 LOCAL path receives only the manifest.

These are per-receiver values for on-demand transfer, six multicast groups,
or another layer-aware transport.  On a single radio channel where TA always
broadcasts all six packets, the airtime remains about 658 KiB per update;
discarding unused packets at the UAV saves decoding but not channel capacity.
The implementation retains the authenticated T4 prefetch API as a backwards-
compatible alternative.  Tests cover every first-difference bit 0--23 and
reject a missing, reordered, cross-sequence or modified layer without changing
the witness or its sequence number.

## Authenticated tail-context deployment

The original local `MemWitUp` mode is still the default and remains available
through `revac_uav_apply_update()` or `revac_uav_apply_update_local()`.  A TA can
now call `revac_ta_revoke_uav_context()` instead.  That call captures `T1` to
`T12` while the normal `Delete/EvalF` chain is already running; it does not run
a second EvalF.  The TA retains a 12-tail catalog for the same 20-update horizon
as the BML and exports a selected network packet with
`revac_ta_export_context()`.

Each packet uses canonical 19-bit coefficient packing and carries:

- the deletion sequence and revoked handle;
- the accumulator-parameter hash and hashes of all 12 catalog tails;
- a TA lattice signature over that manifest;
- only the tail chunks selected by the network profile.

The UAV verifies the version/base/dimension, sequence, handle, public-parameter
hash, TA signature, and every received chunk hash before changing its witness.
`revac_uav_apply_update_context()` can reject a bad/missing packet or fall back
to local computation.  `revac_uav_apply_update_auto()` minimizes estimated
download plus local recurrence time and obtains the chosen packet through a
transport callback.  A real deployment can therefore plug in HTTP, QUIC,
broadcast storage, or another receiver without linking TA code into the UAV.
Likewise, `revac_uav_recover_updates_auto()` applies the same policy to each
missing BML entry. If the BML horizon is exhausted, the UAV no longer calls
`MemWitSync` locally: it sends a nonce-bearing request to a TA endpoint, the TA
uses its private trapdoor to compute the current witness, and the UAV verifies
and installs the signed response.

## Cross-device role boundary

The entity APIs now enforce the following ownership boundary:

| Role | Locally owns/executes |
|---|---|
| TA | signing secret key, accumulator trapdoor/state, member database, `Setup`, `Issue`, `Delete`, `MemWitSync`, context cache |
| UAV | its credential/signing secrets and membership witness; `MemWitUp`, show-proof generation, response verification |
| GCS | its own keys/certificate and public state; show/credential verification |

UAV and GCS functions accept `revac_public_state_t`, not `revac_ta_t`. The
former contains only the TA public key, accumulator public parameters/current
value/sequence, the public-parameter hash, and the public BML window. In the
single-process tests it is a non-owning view created by
`revac_ta_public_state()`; across devices those same public fields must be
provisioned or decoded into device-local storage and refreshed atomically.

The remote recovery protocol is split at the transport callback:

1. UAV creates `(handle, last_seq, pp_hash, random_nonce)`.
2. The authenticated TA endpoint binds the session identity to `handle`,
   rejects unissued/revoked/mismatched members, and alone runs `MemWitSync`.
3. TA returns a canonical full witness plus current sequence, accumulator hash,
   parameter hash, nonce and a TA lattice signature over all of them.
4. UAV checks version/parameters, nonce, sequence, hashes, TA signature and
   `MemVer` using only public state. It writes the witness only after every
   check succeeds.

The TA response signature provides application-layer authenticity and replay
binding, but the request identity and witness confidentiality remain transport
requirements. A deployment must use a mutually authenticated, encrypted
channel (for example mTLS/QUIC) and pass the authenticated device handle to
`revac_ta_export_memwitsync()`; it must not expose that function as an
unauthenticated handle-only endpoint. Revocation broadcasts/public-state
snapshots likewise need an authenticated TA channel.

With the current base-3 parameters, one tail is exactly 114 KiB.  The signed
manifest/signature overhead is about 8.76 KiB.  The available public profiles
are:

| Profile | Packet | Expected local composites |
|---|---:|---:|
| `T4` | 122.76 KiB | 4.625 |
| `T2,T6` | 236.76 KiB | 2.125 |
| `T1,T3,T7` | 350.76 KiB | 1.0547 |
| `T1,T2,T4,T8` | 464.76 KiB | 0.5234 |
| `T1,T2,T3,T5,T9` | 578.76 KiB | 0.2598 |
| `T1,T2,T3,T4,T6,T10` | 692.76 KiB | 0.1289 |
| `T1,T2,T3,T4,T5,T7,T11` | 806.76 KiB | 0.0640 |
| `T1..T6,T8,T12` | 920.76 KiB | 0.0317 |

Run the context integration test to verify all profiles against the original
local witness update, packet tamper rejection, automatic selection, and BML
recovery:

```sh
cd entities/build
./revac_context_demo
```

By default, the full demo keeps the expensive Labrador AccGenpp MemVer proof
disabled so the normal entity test remains quick.  In that mode the GCS verifies
the anonymous signature-show proof, then correctly rejects the full revocable
show because no accumulator-membership ZK proof bytes are attached.

To run the online Labrador-backed show path, enable:

```sh
cd entities/build
REVAC_LABRADOR_ONLINE=1 ./revocable_ac_demo
```

With `REVAC_LABRADOR_ONLINE=1`, the UAV serializes a Labrador proof into
`revac_show_proof_t.acc_zk_proof`, and `revac_gcs_verify_show()` verifies those
proof bytes through `zk_lazer/revac_acc_memver_import_labrador.py`.

### q96 worker prewarm and reproducible LaZer build

The q96 compiler has a large immutable statement.  Production code may move
that construction out of the first latency-sensitive request by calling
`revac_uav_memver_cache_prewarm()` and
`revac_labrador_q96_prover_prewarm()` on the UAV after Issue and
`revac_labrador_q96_verifier_prewarm()` on the GCS after provisioning the
public accumulator parameters.  The UAV cache prewarm performs the same clear
MemVer check required by a show and binds its digest to the complete witness.
The prover prewarm then installs the current private witness only in the UAV's
mode-0600 local worker; neither operation is a network message.  The first and
every subsequent show still use fresh link randomness and generate a fresh
zero-knowledge proof.

`revocable_ac_demo` performs both prewarms by default and reports their time
separately from online show latency.  Set
`REVAC_LABRADOR_Q96_PREWARM=0` to measure a true cold first show.  In a real
cross-device deployment, UAV and GCS prewarm independently and never share a
worker, statement object, or private witness.

When the guarded GCS path verifies the anonymous-signature and accumulator
components concurrently, it caps the signature verifier at four OpenMP
threads by default.  The q96 worker currently serializes transcript-sensitive
native Pack reducers because their shared scratch state is not reentrant;
`OMP_NUM_THREADS=12` still accelerates independent statement-import work.
`REVAC_SHOW_PARALLEL_VERIFY_THREADS` can override the signature cap for a
measured deployment.  This scheduling policy does not change the statement,
transcript, proof encoding, or security parameters.

All LaZer changes required by q96 are stored under `entities/zk_lazer`, so a
clean LaZer checkout can be prepared without relying on an old server image:

```sh
export LAZER_ROOT=/path/to/lazer
code/entities/tools/apply_q96_lazer_patches.sh "$LAZER_ROOT"
make -C "$LAZER_ROOT" -B liblabrador38.so
make -C "$LAZER_ROOT/python" -B
```

Multi-threaded Pack proving is diagnostic-only: repeated pqac8 tests with 12
threads exposed proofs that failed a later tail-reduction self-check, including
when individual optional kernels were disabled.  Production therefore leaves
`REVAC_LABRADOR_THREADS=1` until all reducer scratch becomes request-local.
Host-specific PGO changes code layout only; it does not alter proof parameters
or bytes.  Build it from a representative complete-show workload:

```sh
export LAZER_ROOT=/path/to/lazer
export REVAC_Q96_PGO_DEMO=/path/to/revocable_ac_demo
code/entities/tools/build_q96_labrador_pgo.sh
```

Always run a separate multi-show prove/verify stress test after PGO.  The PGO
script intentionally leaves the optimized library at
`$LAZER_ROOT/liblabrador38.so` and prints its SHA-256 digest, so experiment
records can identify the exact native binary.

The compact transport codec is in `common/revac_proof_codec.c`.  It does two
things:

- packs the signature show proof with field-aware bit packing plus the existing
  Huffman tables for Gaussian response vectors;
- optionally deflate-compresses the Labrador accumulator proof bytes and lets
  the GCS verifier inflate them before calling the normal Labrador verifier.

Run the codec demo with:

```sh
cd entities/build
./revac_proof_codec_demo
REVAC_LABRADOR_ONLINE=1 ./revac_proof_codec_demo
```

The online run verifies the compressed accumulator proof bytes end to end.  The
codec demo first revokes one UAV, applies `MemWitUp` to another UAV, and then
proves and verifies the surviving UAV's updated witness.  The full
`test/revocable_ac_demo.c` flow also uses this compact UAV-to-GCS wire format
through `revac_uav_show_prove_wire()` and `revac_gcs_verify_show_wire()`, and
prints the compressed show packet size in KB during phase 3b.  The older
struct-level APIs remain available for local tests, but the entity demo models
the compressed transport path.

The bridge exports `gadget_base`, `gadget_dim`, `gadget_signed`, and the
public accumulator deletion counter.  The verifier uses the same values when
rebuilding the MemVer statement.  For compatibility, the Python importer also
accepts old JSON without the new gadget fields and interprets it as the former
base-2 gadget.

The Labrador relation groups the 32 handle-indexed `right` polynomials for
each `(gadget digit, output row)` pair.  This reduces 1536 separate right-side
norm/witness parts to 48 vector parts.  Its current post-update norm schedule is
a prototype engineering policy selected from the public deletion counter; it
still needs to be replaced by a bound from the final accumulator parameter and
security analysis before deployment.

The production show path now binds the handle in the signature relation to the
handle in the accumulator MemVer relation through one randomized lattice
commitment `C_x = r0 + A_link*r1 + D_link*x`.  The signature proof proves a
short binary opening using its hidden handle, while the Labrador proof proves a
short opening of the same `C_x` using its own hidden handle.  The GCS verifies
both openings, so mixing independently generated signature and accumulator
proof components is rejected.  These remain two proof systems connected by a
binding commitment, rather than one monolithic proof; the final paper still
needs to state the corresponding Module-SIS/Module-LWE binding and hiding
assumptions and retune the enlarged signature relation parameters.

The UAV show packet intentionally does not carry `user_seed`, `acc_seq`,
`acc`, or `crs_seed`.  The GCS-side verifier receives those from a
`revac_show_context_t`, modeling data synchronized from TA before verification:
`user_seed` for rebuilding the public `Dx` matrix, and `acc_seq/acc` for the
current accumulator state.  The UAV packet contains only the nonce, randomized
same-handle commitment, signature show proof, accumulator proof kind, and
accumulator proof bytes.

### q96 authenticated T4 witness updates

The q96 production revocation wire uses an exact intermediate EvalF tail rather
than transmitting four already-applied DELTA vectors.  During Delete, TA's
single high-to-low traversal both produces the final indicator and captures
`T4`.  The BML signature binds the checkpoint hash, sequence number, revoked
handle, parameter digest and post-deletion accumulator digest.  A surviving
UAV resumes the recurrence from `T4`, computes only `T3,T2,T1`, and applies the
required prefix through its first differing low-order handle bit.  Handles
whose low four bits equal the revoked handle use the original LOCAL path.
Consequently the exact T4 path covers 15/16 = 93.75% of uniformly distributed
survivors without changing MemVer, witness values, or the later ZK proof.

Wire version 5 contains five full q96 polynomials: `5 * 4096 * 96` bits, or
240 KiB, plus the existing TA signature and manifest.  The measured complete
packet is 248.44 KiB.  On the server used for the q96 experiments, T4 resume
took 0.222 s versus 1.105 s for LOCAL, while the earlier applied-DELTA4 packet
was 468.44 KiB.  The core regression checks all first-difference positions
0--3 against LOCAL byte for byte and rejects position 4 without mutation.

T2 and T8 have the same payload size because each is one matrix checkpoint.
Measured worst-covered-handle resume times were 0.101 s and 0.465 s, with
coverage 75% and 99.609%, respectively.  T4 is the default compromise: compared
with T2 it substantially reduces LOCAL fallback, while T8 spends over twice as
much UAV time for only 5.859 percentage points of additional coverage.

T4 also supports an idle-time prefetch profile.  The checkpoint is a
deterministic function of the public parameters and the candidate revoked
handle, so it may be authenticated and cached before the revocation enters the
latency-critical path.  At revocation time the TA sends only the signed
manifest, which binds the sequence, revoked handle, public-parameter digest,
post-delete accumulator digest, T4 digest, codec and payload length.  A cache
hit executes the exact T4 recurrence; a cache miss uses the same authenticated
handle to run the exact LOCAL recurrence, without another download.  Thus the
current real-time broadcast is always the 8.44 KiB manifest.  Each prefetched
packet is 248.41 KiB and its canonical cache payload is 240 KiB; caching 128
candidates costs 30 MiB per UAV (or 240 MiB for the implementation maximum of
1024 candidates).  Deployments may cache only high-risk candidates and accept
LOCAL latency on a miss.  This profile moves optional acceleration data from
the real-time path to registration or idle periods; it does not claim to
compress away the checkpoint entropy.

The q19 DELTA4 result is much smaller for two independent reasons.  Its ring
degree is 256 rather than 4096, and its applied increments are usually small
enough for the ZigZag-9 codec.  Four q19 layers contain
`4 * 12 * 4 * 256 = 49,152` coefficients, whose 9-bit fast representation is
54 KiB before escapes and authentication.  The measured signed packet is
about 62.47 KiB.  In q96, a DELTA4 contains
`4 * 5 * 4096 = 81,920` coefficients and the observed exact representation
needs roughly 46 bits per coefficient, producing the earlier 468.44 KiB
packet.  The q96 T4 checkpoint removes the four-layer factor but must retain
`5 * 4096` full 96-bit CRT coefficients, hence its irreducible 240 KiB payload.
These figures are not directly attributable to base 3 alone: q96 deliberately
uses a much larger ring and modulus to obtain the target security and deletion
capacity.

The checkpoint is pseudorandom-looking under the public matrix chain, so a
generic compressor cannot encode it in tens of KiB.  The truly small endpoint
is the authenticated handle-only update, but decoding that token is exactly
the full local EvalFx computation.  A smaller broadcast with simultaneously
negligible UAV work therefore requires a protocol-level resource such as
pre-distributed checkpoint state, a trusted edge updater, or a newly structured
seed-expandable accumulator matrix, together with the corresponding storage,
trust, or security-proof cost.

The default paths are:

- `REVAC_LABRADOR_PYTHON=/home/wkk/miniconda3/envs/sage/bin/python`
- `REVAC_LABRADOR_SDE=/home/wkk/tools/sde-external-10.8.0-2026-03-15-lin/sde64`
- `REVAC_LABRADOR_SCRIPT=/home/wkk/pq_ac/lattice-anonymous-credentials/code/entities/zk_lazer/revac_acc_memver_import_labrador.py`

Set `REVAC_LABRADOR_USE_SDE=0` only on hardware that can run the LaZer AVX-512
code directly.  The `zk_lazer/` directory also keeps standalone proof-layer
demos and import checks:

```sh
cd entities/zk_lazer
make describe
make same-x-demo
```

The bundled FLINT library may still need compatible system GMP/MPFR runtime
libraries on the target. The current random backend uses AES-NI intrinsics, so
non-x86 UAV targets will need a portable random backend before deployment.
