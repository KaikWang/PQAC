# Lifecycle measurement status (2026-10-13)

Production q88 optimized compiler/native retained. Added phase timing records
and a post-revocation complete Show generation/verification to
`code/entities/test/revocable_ac_demo.c`. Local and remote C builds passed.

Remote command executed three independent runs, PQAC_REPETITIONS=1, explicit
PQAC_DEMO=/root/autodl-tmp/pq_ac/OURS/code/entities/build-baseline-20261013/revocable_ac_demo.
The SSH execution session completed with exit status 0. Result files have not
been downloaded: subsequent SSH calls to autodl-pqa11 (port10657) repeatedly
returned Connection refused. No new measured phase values are available locally.

Remote results:
/root/autodl-tmp/pq_ac/OURS/experiments/q88_latency_20261008/results/lifecycle_v3_20261013/run1
/root/autodl-tmp/pq_ac/OURS/experiments/q88_latency_20261008/results/lifecycle_v3_20261013/run2
/root/autodl-tmp/pq_ac/OURS/experiments/q88_latency_20261008/results/lifecycle_v3_20261013/run3

When server access resumes, download only logs/JSON/public artifacts from
these directories, then run tools/summarize_q88_lifecycle.py locally against
this results directory and write REPORT.md with actual stage means/ranges.
Do not rerun the completed three samples unless their downloaded evidence
shows an incomplete run. Private prover fixtures remain on server.

An initial attempt selected an obsolete remotely copied build-ours binary.
Its separate old_binary_diagnostic directory is excluded from lifecycle data.
The runner now picks the newest build if PQAC_DEMO is not explicitly supplied.

Report full test wall including checks and cleanup separately from protocol
operation times; prewarm separately; registration is 3 UAVs +1 GCS; Issue
is measured individually for each of the three UAVs. Sync recovery uses a
local authenticated TA callback and includes no real network latency.

## Superseded by fresh successful rerun (2026-10-14)

At the user explicit request after restoring the server, three new complete
lifecycle runs finished successfully and their logs/JSON/public artifacts were
downloaded into ../lifecycle_v3_20261014. See that directory REPORT.md and
lifecycle_measurement.json for verified current timings. The connection
block above is resolved. No old uncollected result was used in those timings.
