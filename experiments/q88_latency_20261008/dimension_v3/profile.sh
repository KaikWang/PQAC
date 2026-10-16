#!/usr/bin/env bash
# Independent complete production CRT experiment. Source before benchmarking.
profile_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
workspace=$(cd -- "$profile_dir/../../../.." && pwd)
export PYTHONPATH="$profile_dir:$workspace/OURS/code/entities/zk_lazer:$workspace/lazer/python"
export LD_LIBRARY_PATH="$workspace/OURS/code/entities/common/core/lib:$workspace/lazer"
export LD_PRELOAD="${PQAC_NATIVE_LIBRARY:-$profile_dir/liblabrador38_sprod.so}"
export REVAC_Q88_TIGHT_WITNESS_TOP=1 REVAC_Q96_AGGREGATE_TIGHT_NORM=1
export REVAC_Q88_NORM_GROUP_SIZE=16 REVAC_Q88_APPROX_ORDER=limb
export REVAC_Q88_SCALAR_HANDLE_NORM=1 REVAC_Q88_TIGHT_REFRESH_NORM=1
export REVAC_Q88_AGGREGATE_RADIX=3779 REVAC_Q88_PRESERVE_AUX_REQUIREMENT=1
export REVAC_Q88_AUX_TOP_UNBALANCED=1 REVAC_Q88_SCALED_REFRESH_Z=1
export REVAC_Q88_BATCH_SPLIT=1 REVAC_Q88_FUSED_LING=1 REVAC_Q88_FUSED_LING_THREADS=4
export REVAC_LABRADOR_LINEAR_THREADS=4 REVAC_LABRADOR_SCRATCH_POOL_MIB=1024
export REVAC_LABRADOR_SCRATCH_POOL_AUDIT=1

# Selected kernel; default to independently scheduled six-prime projections.
export REVAC_Q88_SPROD_THREADS=${REVAC_Q88_SPROD_THREADS:-6}
export REVAC_LABRADOR_LDD_ADDCHECK_THREADS=1 REVAC_Q88_DCH_COST=0
