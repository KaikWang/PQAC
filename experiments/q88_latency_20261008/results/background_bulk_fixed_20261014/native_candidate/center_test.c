#include <stdint.h>
#include <stdio.h>
static int64_t q88_center_fold(__int128 value) {
  const unsigned __int128 mask=((unsigned __int128)1 << 38)-1;
  const uint64_t q=((uint64_t)1 << 38)-107;
  const int negative=value<0;
  unsigned __int128 u=(unsigned __int128)value;
  if(negative) u=0-u;
  while(u>mask) u=(u & mask)+(u >> 38)*107;
  uint64_t residue=(uint64_t)u;
  if(residue>=q) residue-=q;
  int64_t centered=residue>q/2 ? (int64_t)residue-(int64_t)q : (int64_t)residue;
  return negative ? -centered : centered;
}


static uint64_t step(uint64_t *s) {*s ^= *s<<13;*s ^= *s>>7;*s ^= *s<<17;return *s;}
static int check(__int128 v) {
 const int64_t q=((int64_t)1<<38)-107;
 __int128 r=v%q;if(r>q/2)r-=q;else if(r< -q/2)r+=q;
 return q88_center_fold(v)==r;
}
int main(void){
 uint64_t seed=0x9173468021ULL;
 for(unsigned i=0;i<1000000;i++){
  uint64_t hi=step(&seed),lo=step(&seed);
  __int128 v=(__int128)(((unsigned __int128)hi<<64)|lo);
  if(!check(v))return 1;
 }
 const unsigned __int128 top=(unsigned __int128)1<<127;
 __int128 edges[]={0,1,-1,(__int128)top,(__int128)(top-1),((__int128)1<<90),-((__int128)1<<90)};
 for(unsigned i=0;i<sizeof(edges)/sizeof(*edges);i++)if(!check(edges[i]))return 2;
 puts("Q88_FOLDED_CENTER_1000007_EXACT_MATCH_PASS");return 0;
}
