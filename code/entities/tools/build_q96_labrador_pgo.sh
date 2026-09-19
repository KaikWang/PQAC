#!/usr/bin/env bash
set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
entities_dir=$(CDPATH= cd -- "$script_dir/.." && pwd)
lazer_root=${LAZER_ROOT:-${1:-}}
demo=${REVAC_Q96_PGO_DEMO:-$entities_dir/build-autodl/revocable_ac_demo}
train_repetitions=${REVAC_Q96_PGO_REPETITIONS:-6}

if [[ -z "$lazer_root" ]]; then
  echo "set LAZER_ROOT or pass the LaZer checkout as argument 1" >&2
  exit 2
fi
if [[ ! -x "$demo" ]]; then
  echo "q96 demo is not executable: $demo" >&2
  exit 2
fi
if ! [[ "$train_repetitions" =~ ^[1-9][0-9]*$ ]] ||
   (( train_repetitions > 20 )); then
  echo "REVAC_Q96_PGO_REPETITIONS must be in [1, 20]" >&2
  exit 2
fi

if [[ ${REVAC_Q96_PGO_APPLY_PATCHES:-1} != 0 ]]; then
  "$script_dir/apply_q96_lazer_patches.sh" "$lazer_root"
else
  echo "q96 LaZer patch application skipped; using the audited source snapshot"
fi

work_dir=$(mktemp -d "${TMPDIR:-/tmp}/revac-q96-pgo.XXXXXX")
profile_dir="$work_dir/profile"
backup="$work_dir/liblabrador38.before.so"
socket_path="$work_dir/worker.sock"
worker_log="$work_dir/worker.log"
training_log="$work_dir/training.log"
mkdir -p "$profile_dir"

library="$lazer_root/liblabrador38.so"
had_library=0
if [[ -f "$library" ]]; then
  cp -p "$library" "$backup"
  had_library=1
fi
completed=0
cleanup() {
  if (( completed == 0 && had_library == 1 )); then
    cp -p "$backup" "$library"
  fi
  rm -rf "$work_dir"
}
trap cleanup EXIT INT TERM

echo "[1/3] building instrumented liblabrador38.so"
make -C "$lazer_root" -B liblabrador38.so \
  LABRADOR_EXTRA_CFLAGS="-fprofile-generate=$profile_dir"

# U/G prewarm plus one prove and one verify request per repetition.  A finite
# worker lifetime is required because GCC writes .gcda counters on clean exit.
max_requests=$((2 + 2 * train_repetitions))
echo "[2/3] training with $train_repetitions complete q96 shows"
env \
  LD_LIBRARY_PATH="$lazer_root${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
  PYTHONPATH="$lazer_root/python${PYTHONPATH:+:$PYTHONPATH}" \
  REVAC_LABRADOR_ONLINE=1 \
  REVAC_LABRADOR_USE_SDE=0 \
  REVAC_LABRADOR_Q96_SOCKET="$socket_path" \
  REVAC_LABRADOR_Q96_WORKER_LOG="$worker_log" \
  REVAC_LABRADOR_Q96_MAX_REQUESTS="$max_requests" \
  REVAC_LABRADOR_Q96_PREWARM=1 \
  REVAC_Q96_SHOW_REPETITIONS="$train_repetitions" \
  "$demo" >"$training_log" 2>&1

if ! find "$profile_dir" -type f -name '*.gcda' -print -quit | grep -q .; then
  echo "training completed but no GCC profile counters were written" >&2
  exit 1
fi

echo "[3/3] building profile-guided liblabrador38.so"
make -C "$lazer_root" -B liblabrador38.so \
  LABRADOR_EXTRA_CFLAGS="-fprofile-use=$profile_dir -fprofile-correction -Wno-missing-profile"

completed=1
echo "PGO library: $library"
echo "SHA256: $(sha256sum "$library" | awk '{print $1}')"
echo "Run a separate correctness/stress pass before deploying this host-specific library."
