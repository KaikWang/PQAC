#ifndef Q88_FAST_MOD_H
#define Q88_FAST_MOD_H
#include <stdint.h>
#include <limits.h>
/* floor(2^32 / p) gives a quotient underestimated by at most one for
 * 0 <= x < 2^32. One correction yields exactly x mod p. */
static inline int64_t q88_fast_small_mod(int64_t x, uint32_t p, uint64_t reciprocal) {
  uint64_t magnitude = x < 0 ? (uint64_t)(-x) : (uint64_t)x;
  uint64_t quotient = (magnitude * reciprocal) >> 32;
  uint64_t residue = magnitude - quotient * p;
  if (residue >= p) residue -= p;
  return x < 0 && residue ? (int64_t)(p-residue) : (int64_t)residue;
}
/* Exact reduction for q=2^bits-offset. Replacing 2^bits by offset
 * preserves the residue. The guard prevents truncation/overflow even for
 * signed 128-bit extremes; those inputs use the reference division. */
static inline int64_t q88_fast_center(__int128 value, unsigned bits, uint64_t offset) {
  const uint64_t radix = (uint64_t)1 << bits;
  const uint64_t mask = radix-1, q = radix-offset;
  const unsigned __int128 magnitude = value < 0
      ? -(unsigned __int128)value : (unsigned __int128)value;
  const unsigned __int128 high = magnitude >> bits;
  if (high > (UINT64_MAX-mask)/offset) {
    __int128 r = value % q;
    if (r > (__int128)(q/2)) r -= q;
    else if (r < -(__int128)(q/2)) r += q;
    return (int64_t)r;
  }
  uint64_t r = (uint64_t)(magnitude & mask) + (uint64_t)high*offset;
  /* For the supported (32,99), (36,243), (38,107) moduli,
   * two further folds of a uint64 leave r < 2q. */
  r = (r & mask) + (r >> bits)*offset;
  r = (r & mask) + (r >> bits)*offset;
  if (r >= q) r -= q;
  int64_t centered = r > q/2 ? (int64_t)r-(int64_t)q : (int64_t)r;
  return value < 0 ? -centered : centered;
}
#endif
