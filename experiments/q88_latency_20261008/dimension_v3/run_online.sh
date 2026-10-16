#!/usr/bin/env bash
# Run the independent compiler with the existing production binary and engine.
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
workspace=$(cd -- "$script_dir/../../../.." && pwd)
python_bin=${PQAC_PYTHON:-/root/miniconda3/envs/sage/bin/python}
demo=${PQAC_DEMO:-$workspace/OURS/code/entities/build-baseline-20261013/revocable_ac_demo}
results=${PQAC_RESULTS_DIR:-$workspace/OURS/experiments/q88_latency_20261008/results/dimension_v3_run_$$}
socket=/tmp/pqac_dimension_v3_$$.sock
mkdir -p -- "$results"
export PYTHONPATH="$script_dir:$workspace/OURS/code/entities/zk_lazer:$workspace/lazer/python"
export LD_LIBRARY_PATH="$workspace/OURS/code/entities/common/core/lib:$workspace/lazer"
source "$script_dir/profile.sh"
export REVAC_LABRADOR_ONLINE=1 REVAC_LABRADOR_USE_SDE=0
export REVAC_LABRADOR_Q96_SOCKET="$socket"
export REVAC_Q96_SHOW_REPETITIONS=${PQAC_REPETITIONS:-5}
export REVAC_Q96_EXPERIMENT_PROOF_DUMP="$results/last_accumulator_envelope.bin"
"$python_bin" -u "$script_dir/revac_acc96_memver_worker.py" --socket "$socket" > "$results/worker.log" 2>&1 &
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
printf 'Results: %s\n' "$results"
