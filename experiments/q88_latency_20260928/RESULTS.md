# q88 latency optimization log

All measurements below are complete online shows on `autodl-pqac8`, with the
resident LaBRADOR worker and four or three repetitions.  The credential wire
format and q88 parameters were unchanged.

| case | generation (s) | verification (s) | envelope | CUDA audit |
|---|---:|---:|---:|---|
| final CPU profile, 4 shows | 5.331, 5.276, 6.326, 5.907 | 3.010, 2.871, 2.950, 2.977 | 193.97--194.53 KiB | 54,844 calls; max length 3033; 0 GPU dispatches |
| `OMP_WAIT_POLICY=PASSIVE`, 3 shows | 5.669, 6.062, 6.393 | 3.062, 3.146, 3.294 | 194.00--194.53 KiB | 0 GPU dispatches |
| first-round `nmax=2000`, 3 shows | 5.607, 5.530, 5.969 | 2.922, 2.900, 2.973 | 194.41--194.53 KiB | 0 GPU dispatches |
| exact CRT disabled, 3 shows | 6.041, 5.667, 5.324 | 3.308, 3.174, 2.973 | 193.97--194.10 KiB | 0 GPU dispatches |
| prove rotation aggregation enabled, 5 shows | 5.461, 4.965, 6.223, 5.556, 5.448 | 2.923, 3.076, 2.917, 3.082, 3.079 | 194.00--194.53 KiB | 0 GPU dispatches |
| default comparison, 5 shows | 5.358, 5.404, 5.978, 7.071, 5.277 | 3.007, 2.930, 2.932, 2.834, 2.884 | 194.00--194.53 KiB | 0 GPU dispatches |
| fused CRT scale factors, 3 shows | 5.571, 5.832, 5.096 | 2.898, 2.765, 2.827 | 194.03--194.53 KiB | CRT diagnostic passed; 0 GPU dispatches |
| actual q88 verifier phase budget=16, 3 shows | 5.339, 5.067, 5.735 | 2.909, 2.961, 2.969 | 194.00--194.28 KiB | all shows passed |
| batched rotation experiment, 3 shows | 5.664, 5.048, 7.270 | 2.648, 2.780, 3.055 | 194.00--194.06 KiB | coefficient audit passed |
| LNP norm threads=8, 3 shows | 5.641, 4.806, 5.112 | 2.812, 2.867, 2.881 | 194.03--194.53 KiB | all shows passed |
| LNP sections parallel, 3 shows | 5.621, 4.849, 5.465 | 2.981, 2.913, 2.750 | 194.00--194.53 KiB | all shows passed |
| LNP z1 sigma ceil, 3 shows | 5.754, 5.261, 5.315 | 3.483, 2.974, 2.895 | 205.43 KiB | rejected: size +5.6% |
| direct LNP response addition, 3 shows | 4.928, 5.475, 5.234 | 2.802, 2.810, 2.744 | 194.03--194.53 KiB | coefficient audit passed |
| bounded CRT reduction, 3 shows | 5.685, 5.088, 5.552 | 2.904, 2.771, 2.824 | 194.35--194.41 KiB | 100026 arithmetic tests; full audit passed |
| fused LNP rejection products, 3 shows | 5.303, 5.612, 5.309 | 2.840, 2.772, 2.897 | 194.06--194.53 KiB | old/new products audited |

The conservative CUDA threshold is now opt-in.  The complete q88 call graph
only creates sprod calls up to about 3033 coefficients, while the previous
threshold was 32768; therefore the earlier apparent GPU gain was measurement
noise, not GPU work.  The audited default is CPU and the CUDA code remains
available for a future batched implementation.

The exact CRT conversion now precomputes `inv_scale[k] * crt_inv[k]` once and
uses one modular reduction per coefficient instead of two.  The optional
`REVAC_LABRADOR_CRT_SCALE_CHECK=1` path compared the fused and unfused formulas
inside a complete show and passed.  The fused path is retained because it
changes only arithmetic scheduling, not the reconstructed residue.

Stage timing on one profiled show identifies the main cost as the first two
Labradoodle folds (about 1.84 s and 1.20 s in prove; 0.93 s and 0.69 s in
verify).  LNP is only about 0.29 s in prove and 0.09 s in verify, so reducing
its repetitions alone cannot reach the target while preserving security.

Correction (2026-10-01): the phase-specific scope in
`revac_acc_memver_worker.py` is not called by `revac_acc96_memver_worker.py`.
The rows labelled "prove rotation aggregation enabled" therefore do not
measure an activated q88 optimization.  The unsupported default change was
reverted.  These measurements are retained as raw observations only.

The actual q88 worker now accepts phase-specific thread variables such as
`REVAC_LABRADOR_LDD_PROJECTION_THREADS_VERIFY=16`.  This is an opt-in tuning
interface; stable defaults remain unchanged because this sample did not
establish a consistent sub-2-second verifier result.

The first phase-budget sample was discarded because the initial implementation
allowed verify settings to leak into the next prove request.  The worker now
restores the startup budget at every request boundary; the table contains only
the corrected run.

Research note (2026-10-01): the current LaZer/Labrador source already uses
AVX-512 pointwise NTT-domain kernels.  Orthus is a separate
batch-verification construction in the upstream Labrador repository, not a
drop-in replacement for this q88 credential transcript.  CUDA batching only
becomes useful after host/device transfers and allocations are amortized; the
current q88 calls are many short sprod operations, so the existing per-call
CUDA experiment remains disabled.  A GPU path must batch calls and keep device
buffers resident before it can be judged.

The rotation batching path is controlled by
`REVAC_LABRADOR_ROTATION_BATCH=1`; its optional
`REVAC_LABRADOR_ROTATION_CHECK=1` mode evaluates both implementations and
aborts on any coefficient or width mismatch.  It remains opt-in because the
current sample has not shown a stable latency improvement.

The q88 worker now defaults to `REVAC_LABRADOR_CRT_BOUNDED_REDUCE=1`.
Exact CRT accumulation has six nonnegative residues, so the centered sum is
bounded below `6P`; a fixed sequence of conditional subtractions replaces the
128-bit division.  Setting the variable to `2` enables a reference comparison
on every coefficient for regression testing.

The LNP response path now has an opt-in
`REVAC_LABRADOR_LNP_RESPONSE_ADD=1` mode that adds the coefficient-domain
Gaussian response directly.  Mode `2` compares it against the original NTT
round-trip path coefficient by coefficient and passed complete shows.  It is
not yet the default because the short benchmark did not establish a stable
latency gain.

`REVAC_LABRADOR_LNP_REJECTION_PRODUCTS=1` fuses the two response inner
products into one coefficient pass; value `2` compares both products with the
old two-pass implementation. Complete-show audits passed, but the short
sample is not yet enough to make it the default.

The `REVAC_LNP_Z1_SIGMA=ceil` sampler-grid experiment was rejected: it did not
reduce latency and increased the proof envelope from about 194 KiB to 205 KiB.
