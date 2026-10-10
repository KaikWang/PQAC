# Rotation scratch allocation (2026-10-08)

`rotation_aggregate_internal()` used to allocate two temporary polynomial
vectors on every call, even when every block took the scalar fast path and
never used them.  The current implementation allocates those vectors only
upon the first nonscalar block.  The arithmetic and proof format are
unchanged.  `REVAC_LABRADOR_ROTATION_LAZY_SCRATCH=0` restores eager allocation
for comparisons.

In two profiled q88 full shows per mode, the first 8,670-constraint sparse phi
aggregation took 0.235 and 0.238 s with lazy allocation versus 0.280 and
0.278 s with eager allocation during proving.  The corresponding verify
measurements were 0.231, 0.259 s versus 0.240, 0.248 s; that difference is
within run variation.  All four Pack proofs verified.  Whole-show times were
4.788/2.732 and 4.718/2.686 s for lazy, versus 5.056/2.981 and
5.351/2.851 s for eager (prove/verify); the sample is too small to assign
all of the whole-show difference to this change.

Raw timing logs: `rotation_lazy.worker.log`, `rotation_eager.worker.log`.
The 4 s / 2 s target remains unmet.

After keeping both this change and the Dachshund scheduling change enabled by
default, two unprofiled complete shows passed.  Full show prove/verify was
4.762/2.782 and 5.489/2.799 s; Pack prove request was 4.693 and 5.419 s,
and Pack Verify call was 2.755 and 2.773 s.  Proof envelopes were 194.53
and 194.03 KiB.  Logs: `q88_current.show.log`, `q88_current.worker.log`.
