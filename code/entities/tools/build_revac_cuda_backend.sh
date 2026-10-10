#!/usr/bin/env bash
set -euo pipefail

# Build the optional runtime CUDA backend.  The default PTX target keeps this
# usable on later NVIDIA architectures; set REVAC_CUDA_ARCH to a local
# compute capability (for example, 89) when a deployment wants native SASS.
root_dir=$(cd "$(dirname "$0")/../../../.." && pwd)
cuda_src="$root_dir/lazer/cuda/revac_labrador_cuda_sprod.cu"
output=${REVAC_CUDA_OUTPUT:-$root_dir/lazer/librevac_labrador_cuda.so}
nvcc_bin=${NVCC:-nvcc}

if ! command -v "$nvcc_bin" >/dev/null 2>&1; then
  echo "nvcc not found; leave the CPU backend enabled" >&2
  exit 2
fi

if [[ -n "${REVAC_CUDA_ARCH:-}" ]]; then
  arch_flags=("-arch=sm_${REVAC_CUDA_ARCH}")
else
  arch_flags=("-gencode" "arch=compute_75,code=compute_75")
fi

"$nvcc_bin" -O3 -std=c++14 -Xcompiler -fPIC "${arch_flags[@]}" \
  -shared "$cuda_src" -o "$output"
echo "built $output"
