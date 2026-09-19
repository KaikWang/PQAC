# PQAC q88 Source Bundle

This directory is the local, curated source bundle for the current revocable
post-quantum anonymous credential implementation.

## Contents

- `code/entities/`: protocol, q88 accumulator, anonymous signature, revocation,
  DELTA6 synchronization, proof bridge, tests, and the Python q96/q88 worker.
- `params/`: parameter manifests, generated profiles, and audit results.
- `scripts/params/`: parameter generation and consistency checks used by the
  CMake build and q88 regression tests.
- `code/entities/common/core/flint_headers/`: bundled FLINT headers.
- `code/entities/common/core/lib/libflint.so`: bundled FLINT runtime required by
  the entity binaries.
- `code/entities/common/core/lib/libflint.so.18`: SONAME compatibility symlink
  required by the bundled FLINT runtime loader.

## External dependency

The online Labrador accumulator proof uses the modified LaZer source tree:

`/home/wkk/pq_ac/lazer`

`OURS` contains the bridge, worker scripts, wire format, and patch files, but
does not duplicate LaZer. The clean upstream reference is
`/home/wkk/pq_ac/LAZER` and must remain read-only.

## Build from this bundle

```sh
cd /home/wkk/pq_ac/OURS
cmake -S code/entities -B code/entities/build-ours \
  -DCMAKE_BUILD_TYPE=Release \
  -DREVAC_LAZER_PYTHON=/home/wkk/miniconda3/envs/sage/bin/python \
  -DCMAKE_C_COMPILER=/usr/bin/cc
cmake --build code/entities/build-ours -j12
ctest --test-dir code/entities/build-ours --output-on-failure \
  -R '^(revac_params_check|revac_q88_params_check|revocable_ac_demo|revocable_ac_q19_legacy_demo|revac_acc96_relation_demo)$'
```

For the online proof path, build the modified LaZer tree first and set
`LD_LIBRARY_PATH`, `PYTHONPATH`, `REVAC_LABRADOR_PYTHON`, and
`REVAC_LABRADOR_ONLINE=1` as documented in `code/entities/README.md`.

## Source-of-truth and experiment workflow

`OURS` is the authoritative editable source tree for the anonymous credential
scheme. The complete system consists of this repository plus the modified
LaZer dependency at `/home/wkk/pq_ac/lazer`.

Source changes are made only in the local WSL repositories. For an experiment,
the local sources are synchronized in one direction to the server selected for
that experiment, then built and tested on that server. Server copies must not
overwrite the local source, and Intel SDE is not used.

Build directories, Python caches, Git metadata, and server artifacts are
intentionally excluded from version control.
