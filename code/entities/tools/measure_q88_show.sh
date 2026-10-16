#!/usr/bin/env bash
# Measure the production complete q88 Show packet, including serialization.
set -euo pipefail
entities=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
workspace=$(cd -- "$entities/../../.." && pwd)
python_bin=${PQAC_PYTHON:-/root/miniconda3/envs/sage/bin/python}
# Pick the newest available build unless the caller specifies a binary.
default_demo="$entities/build-ours/revocable_ac_demo"
for candidate in "$entities/build-ours/revocable_ac_demo" \
                 "$entities/build-baseline-20261013/revocable_ac_demo"; do
    if [[ -x "$candidate" && ( ! -x "$default_demo" || "$candidate" -nt "$default_demo" ) ]]; then
        default_demo="$candidate"
    fi
done
demo=${PQAC_DEMO:-$default_demo}
results=${PQAC_RESULTS_DIR:-$workspace/OURS/experiments/q88_latency_20261008/results/production_v3_$$}
socket=/tmp/pqac_production_v3_$$.sock
mkdir -p -- "$results"
export PYTHONPATH="$entities/zk_lazer:$workspace/lazer/python"
export LD_LIBRARY_PATH="$entities/common/core/lib:$workspace/lazer"
unset LD_PRELOAD
export REVAC_LABRADOR_ONLINE=1 REVAC_LABRADOR_USE_SDE=0
export REVAC_LABRADOR_Q96_SOCKET="$socket"
export REVAC_Q96_SHOW_REPETITIONS=${PQAC_REPETITIONS:-5}
export REVAC_Q88_SHOW_PACKET_DUMP="$results/last_show_packet.bin"
export REVAC_Q88_ACCUMULATOR_DUMP="$results/last_accumulator_envelope.bin"
"$python_bin" -u "$entities/zk_lazer/revac_acc96_memver_worker.py" --socket "$socket" > "$results/worker.log" 2>&1 &
worker_pid=$!
cleanup() {
    kill "$worker_pid" 2>/dev/null || true
    wait "$worker_pid" 2>/dev/null || true
    rm -f -- "$socket"
}
trap cleanup EXIT
for ((attempt=0; attempt<100; attempt++)); do
    if [[ -S "$socket" ]]; then break; fi
    if ! kill -0 "$worker_pid" 2>/dev/null; then
        cat -- "$results/worker.log" >&2
        exit 1
    fi
    sleep 0.1
done
[[ -S "$socket" ]] || { printf 'Worker startup timed out\n' >&2; exit 1; }
cd -- "$(dirname -- "$demo")"
"$demo" > "$results/show.log" 2>&1
"$python_bin" "$entities/tools/summarize_q88_show.py" "$results"
printf 'Results: %s\n' "$results"
