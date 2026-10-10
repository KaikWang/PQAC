# LNP outer retry sampling (2026-10-08)

The outer LNP retry loop already keeps the fixed folded relation, commitments,
and g0/g1 inputs outside the inner retry.  A fresh outer mask and each inner
challenge change the responses, so those responses cannot be cached across
their respective retries.  Profiling found discrete-Gaussian sampling and
coefficient materialization to be the dominant outer-loop operation.

`polzvec_gaussian_scratch()` already had a fixed-chunk, domain-separated
sampler, but its team size was tied to `REVAC_LABRADOR_THREADS=1` in the q88
worker.  `REVAC_LABRADOR_GAUSSIAN_THREADS` now sets only that sampler and its
coefficient conversion team.  The q88 worker defaults to 8, capped by the
available OpenMP threads and work.  Setting it to 1 restores serial sampling.
Chunk boundaries (65,536 coefficients), the discrete-Gaussian sampler,
standard deviation, rejection checks, proof format, and verifier are unchanged.

| Mode | Outer attempts | Gaussian and conversion total | Per attempt |
| --- | --- | --- | --- |
| 8 threads | 8, 1, 9 | 0.0609, 0.0070, 0.0716 s | 7.0-8.0 ms |
| 1 thread | 1, 7, 5 | 0.0317, 0.2161, 0.1505 s | 30.1-31.7 ms |

All six profiled full q88 shows verified.  Whole-show prove/verify times were
5.091/2.764, 4.628/2.645, and 4.924/2.849 s with 8 threads; with one
thread they were 4.871/2.795, 5.139/2.687, and 4.829/2.760 s.  Different
rejection counts make these small samples unsuitable for estimating a stable
whole-show speedup.  The kernel improvement is much clearer and reduces
long-tail cost roughly in proportion to the number of outer attempts.

Security rationale: each fixed chunk uses a key derived from seed, nonce,
chunk index, and a domain tag via SHAKE128, then invokes the same AES-based
Gaussian sampler.  Under the existing SHAKE128/AES pseudorandomness
assumptions, the chunk streams are computationally indistinguishable from
independent random streams, so the Gaussian mask distribution used by LNP is
preserved up to that assumption.  This is an implementation argument, not a
new formal proof or a change to LNP security parameters.

The 4 s / 2 s accumulator target remains unmet.

Five additional unprofiled shows with the new q88 default all verified.
Full-show prove times were 5.339, 5.128, 4.778, 5.085, and 4.724 s;
verification times were 2.734, 2.689, 2.690, 2.737, and 2.772 s.
Accumulator proof envelopes ranged from 193.88 to 194.44 KiB.  The native
timing log is `lnp_gauss_default.worker.log`.
