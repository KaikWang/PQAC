#!/usr/bin/env bash
set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
entities_dir=$(CDPATH= cd -- "$script_dir/.." && pwd)
lazer_root=${1:-${LAZER_ROOT:-}}

if [[ -z "$lazer_root" ]]; then
  echo "usage: $0 /path/to/lazer" >&2
  echo "or set LAZER_ROOT" >&2
  exit 2
fi
if [[ ! -e "$lazer_root/.git" || ! -e "$lazer_root/src/labrados/.git" ]]; then
  echo "not a LaZer checkout with the labrados submodule: $lazer_root" >&2
  exit 2
fi

apply_once() {
  local repository=$1
  local patch_file=$2

  if git -C "$repository" apply --check "$patch_file" 2>/dev/null; then
    git -C "$repository" apply "$patch_file"
    echo "applied: $patch_file"
  elif git -C "$repository" apply --check --reverse "$patch_file" 2>/dev/null; then
    echo "already applied: $patch_file"
  else
    echo "cannot apply cleanly and is not already applied: $patch_file" >&2
    return 1
  fi
}

apply_once "$lazer_root" \
  "$entities_dir/zk_lazer/lazer_q96_outer_full.patch"
apply_once "$lazer_root/src/labrados" \
  "$entities_dir/zk_lazer/lazer_q96_labrados_full.patch"

profile_source="$entities_dir/zk_lazer/labrados_overlay/revac_profile.h"
profile_target="$lazer_root/src/labrados/revac_profile.h"
if [[ ! -e "$profile_target" ]]; then
  install -m 0644 "$profile_source" "$profile_target"
  echo "installed: $profile_target"
elif cmp -s "$profile_source" "$profile_target"; then
  echo "already installed: $profile_target"
else
  echo "refusing to overwrite a different $profile_target" >&2
  exit 1
fi

echo "q96 LaZer source patches are ready."
