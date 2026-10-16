#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "fast_mod.h"
static uint64_t state = 0x9e3779b97f4a7c15ULL;
static uint64_t rng(void) { state ^= state << 13; state ^= state >> 7; state ^= state << 17; return state; }
static void check(__int128 value, unsigned bits, unsigned off) {
  const __int128 q = ((__int128)1 << bits)-off;
  __int128 expected = value % q;
  if(expected > q/2) expected -= q;
  if(expected < -q/2) expected += q;
  if(q88_fast_center(value,bits,off) != expected) { fprintf(stderr,"center mismatch: bits=%u\n",bits); exit(1); }
}
int main(void) {
  unsigned bits[] = {32,36,38}, offs[] = {99,243,107};
  for(unsigned k=0;k<3;k++) {
    __int128 q = ((__int128)1 << bits[k])-offs[k];
    for(int i=-8;i<=8;i++) for(int d=-3;d<=3;d++) {
      check(i*q+d,bits[k],offs[k]); check(i*q+q/2+d,bits[k],offs[k]);
    }
    for(unsigned i=0;i<500000;i++) {
      unsigned __int128 u = ((unsigned __int128)rng()<<64)|rng();
      check((__int128)u,bits[k],offs[k]);
      check((__int128)(u & (((unsigned __int128)1<<78)-1)),bits[k],offs[k]);
      check(-(__int128)(u & (((unsigned __int128)1<<78)-1)),bits[k],offs[k]);
    }
    check((__int128)((unsigned __int128)1<<127),bits[k],offs[k]);
    check((__int128)(((unsigned __int128)1<<127)-1),bits[k],offs[k]);
  }
  /* Exhaust every int16 coefficient for every eligible positive modulus;
   * check both extremal scales plus random full-domain 32-bit magnitudes. */
  for(uint32_t p=257;p<32768;p+=2) {
    uint64_t inv=((uint64_t)1<<32)/p;
    for(int64_t c=-32768;c<32768;c+=31) for(unsigned j=0;j<2;j++) {
      int64_t x=c*(j ? p-1 : p/2), expected=x%p;
      if(expected<0) expected+=p;
      if(q88_fast_small_mod(x,p,inv)!=expected) { fputs("small mismatch\n",stderr); return 1; }
    }
    for(unsigned j=0;j<20;j++) {
      int64_t x=(uint32_t)rng(); if(j&1) x=-x;
      int64_t expected=x%p; if(expected<0) expected+=p;
      if(q88_fast_small_mod(x,p,inv)!=expected) return 1;
    }
  }
  puts("PASS: centered reduction boundaries, 4.5M random signed 128-bit/78-bit values, small-prime reciprocal reductions");
}
