#include <cuda_runtime.h>
#include <stdint.h>
#include <stddef.h>

/* Optional, generic CUDA backend for one large CRT vector inner product.
 * The host library loads this symbol dynamically, so this file is never a
 * build requirement for CPU-only deployments. */
__global__ static void revac_sprod_kernel(const int16_t *const *a,
                                          const int16_t *const *b,
                                          int16_t *const *r,
                                          size_t len, size_t n, size_t k,
                                          const int16_t *p,
                                          const int16_t *inv_r) {
  const size_t index = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
  const size_t total = k * n;
  if(index >= total) return;
  const size_t prime = index / n;
  const size_t coefficient = index % n;
  const int64_t modulus = p[prime];
  int64_t sum = 0;
  for(size_t i = 0; i < len; i++) {
    sum += (int64_t)a[prime][i * n + coefficient] *
           (int64_t)b[prime][i * n + coefficient];
  }
  sum %= modulus;
  sum = (sum * inv_r[prime]) % modulus;
  if(sum > modulus / 2) sum -= modulus;
  if(sum < -modulus / 2) sum += modulus;
  r[prime][coefficient] = (int16_t)sum;
}

extern "C" int revac_cuda_sprod(const int16_t *const *a,
                                 const int16_t *const *b,
                                 int16_t *const *r, size_t len, size_t n,
                                 size_t k, const int16_t *p,
                                 const int16_t *inv_r) {
  if(a == nullptr || b == nullptr || r == nullptr || p == nullptr ||
     inv_r == nullptr || len == 0 || n == 0 || k == 0) return 0;
  int device_count = 0;
  if(cudaGetDeviceCount(&device_count) != cudaSuccess || device_count < 1)
    return 0;

  const size_t input_bytes = len * n * sizeof(int16_t);
  const size_t output_bytes = n * sizeof(int16_t);
  const size_t pointer_bytes = k * sizeof(void *);
  int16_t **da = nullptr, **db = nullptr, **dr = nullptr;
  int16_t *dpa = nullptr, *dpb = nullptr, *dpr = nullptr;
  int16_t *dp = nullptr, *dinv = nullptr;
  int16_t **ha = new int16_t *[k];
  int16_t **hb = new int16_t *[k];
  int16_t **hr = new int16_t *[k];
  int ok = 0;

  if(cudaMalloc(&dpa, input_bytes * k) != cudaSuccess ||
     cudaMalloc(&dpb, input_bytes * k) != cudaSuccess ||
     cudaMalloc(&dpr, output_bytes * k) != cudaSuccess ||
     cudaMalloc(&dp, k * sizeof(int16_t)) != cudaSuccess ||
     cudaMalloc(&dinv, k * sizeof(int16_t)) != cudaSuccess ||
     cudaMalloc(&da, pointer_bytes) != cudaSuccess ||
     cudaMalloc(&db, pointer_bytes) != cudaSuccess ||
     cudaMalloc(&dr, pointer_bytes) != cudaSuccess) goto cleanup;

  for(size_t i = 0; i < k; i++) {
    ha[i] = dpa + i * len * n;
    hb[i] = dpb + i * len * n;
    hr[i] = dpr + i * n;
    if(cudaMemcpy(ha[i], a[i], input_bytes, cudaMemcpyHostToDevice) != cudaSuccess ||
       cudaMemcpy(hb[i], b[i], input_bytes, cudaMemcpyHostToDevice) != cudaSuccess)
      goto cleanup;
  }
  if(cudaMemcpy(dp, p, k * sizeof(int16_t), cudaMemcpyHostToDevice) != cudaSuccess ||
     cudaMemcpy(dinv, inv_r, k * sizeof(int16_t), cudaMemcpyHostToDevice) != cudaSuccess ||
     cudaMemcpy(da, ha, pointer_bytes, cudaMemcpyHostToDevice) != cudaSuccess ||
     cudaMemcpy(db, hb, pointer_bytes, cudaMemcpyHostToDevice) != cudaSuccess ||
     cudaMemcpy(dr, hr, pointer_bytes, cudaMemcpyHostToDevice) != cudaSuccess)
    goto cleanup;

  revac_sprod_kernel<<<(unsigned)((k * n + 255) / 256), 256>>>(
      (const int16_t *const *)da, (const int16_t *const *)db, dr,
      len, n, k, dp, dinv);
  if(cudaGetLastError() != cudaSuccess || cudaDeviceSynchronize() != cudaSuccess)
    goto cleanup;
  for(size_t i = 0; i < k; i++)
    if(cudaMemcpy(r[i], hr[i], output_bytes, cudaMemcpyDeviceToHost) != cudaSuccess)
      goto cleanup;
  ok = 1;

cleanup:
  cudaFree(da); cudaFree(db); cudaFree(dr); cudaFree(dpa); cudaFree(dpb);
  cudaFree(dpr); cudaFree(dp); cudaFree(dinv);
  delete[] ha; delete[] hb; delete[] hr;
  return ok;
}
