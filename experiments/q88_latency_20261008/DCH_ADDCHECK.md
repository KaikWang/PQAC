# q88 Dachshund quadratic addchecks (2026-10-08)

This experiment changes the scheduling of three independent Dachshund
quadratic constraint builders.  The quadratic-left builder writes a separate
output range and now runs on a pthread while the main thread constructs the
quadratic-right and aggregate relations.  The main thread retains the
existing 8-thread OpenMP exact CRT refresh.  Set
`REVAC_LABRADOR_DCH_ADDCHECK_PARALLEL=0` to disable the new scheduling.
The proof format, relation, challenges and parameters are unchanged.

An OpenMP `parallel sections` trial was discarded: its nested exact CRT
refresh ran serially and increased the quadratic aggregate from about 0.17 s
to about 1.0 s.  The pthread version passed two complete q88 online shows;
the serial control used the same compiled native library.

| Mode | Dachshund prove (s) | Dachshund verify (s) | Pack Verify call (s) | Full show prove/verify (s) |
| --- | --- | --- | --- | --- |
| pthread enabled | 0.342, 0.361 | 0.234, 0.259 | 2.699, 2.699 | 4.743/2.723, 5.342/2.721 |
| pthread disabled | 0.514, 0.466 | 0.327, 0.375 | 2.733, 2.928 | 6.562/2.758, 5.773/2.951 |

Both modes produced valid Pack proofs.  Accumulator proof envelopes across
these four randomized runs were 194.00-194.44 KiB; the scheduling change
cannot change the proof format or length for a fixed proof transcript.
Two runs per mode establish the stage improvement, but are too few to infer a
stable whole-show speedup.  The 4 s / 2 s goal remains unmet.

Raw native timing logs: `dch_pthread.worker.log`, `dch_serial.worker.log`.
