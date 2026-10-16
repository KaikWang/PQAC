#!/usr/bin/env bash
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root=$(cd -- "$here/../.." && pwd)
entities="$root/code/entities"
build=${PQAC_ENTITY_BUILD:-$entities/build-ours}
core="$entities/common/core"
output=${PQAC_REFRESH_AUDIT_OUTPUT:-$here/clear_audit}
cc -std=gnu11 -O2 -fopenmp \
  -I"$core/include" -I"$core/flint_headers" \
  "$here/fresh_witness_audit.c" "$build/librevac_entities_core.a" \
  -L"$core/lib" -Wl,-rpath,"$core/lib" -lflint -lz -lmpfr -lgmp -lm \
  -o "$output"
"$output"
