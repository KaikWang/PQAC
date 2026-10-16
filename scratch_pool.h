#ifndef Q88_SCRATCH_POOL_H
#define Q88_SCRATCH_POOL_H
#include <pthread.h>
#include <stdint.h>
/* Reuse allocations only, never witness values, masks, or challenges.
 * Exact-size matching preserves the existing polxvec allocation contract.
 * The cache has a fixed global byte budget and a hard entry-count limit. */
#define Q88_POOL_SLOTS 128
struct q88_pool_slot { void *data; size_t bytes; };
static struct q88_pool_slot q88_pool[Q88_POOL_SLOTS];
static size_t q88_pool_bytes, q88_pool_limit;
static unsigned long q88_pool_hits, q88_pool_misses;
static pthread_mutex_t q88_pool_lock = PTHREAD_MUTEX_INITIALIZER;
static pthread_once_t q88_pool_once = PTHREAD_ONCE_INIT;
static void q88_pool_init(void) {
  const char *env = getenv("REVAC_LABRADOR_SCRATCH_POOL_MIB");
  unsigned long mib = env ? strtoul(env, NULL, 10) : 0;
  if (mib > 4096) mib = 4096;
  q88_pool_limit = (size_t)mib * 1024 * 1024;
}
static void *q88_pool_take(size_t bytes) {
  pthread_once(&q88_pool_once, q88_pool_init);
  if (bytes < 65536 || bytes > q88_pool_limit) return _aligned_alloc(64,bytes);
  pthread_mutex_lock(&q88_pool_lock);
  for (size_t i=0;i<Q88_POOL_SLOTS;i++) if (q88_pool[i].data && q88_pool[i].bytes==bytes) {
    void *data=q88_pool[i].data; q88_pool[i].data=NULL;
    q88_pool_bytes-=bytes; q88_pool_hits++;
    pthread_mutex_unlock(&q88_pool_lock); return data;
  }
  q88_pool_misses++; pthread_mutex_unlock(&q88_pool_lock);
  return _aligned_alloc(64,bytes);
}
static void q88_pool_give(void *data,size_t bytes) {
  if (!data) return;
  pthread_once(&q88_pool_once,q88_pool_init);
  if(bytes<65536 || bytes>q88_pool_limit) { free(data);return; }
  pthread_mutex_lock(&q88_pool_lock);
  if (bytes<=q88_pool_limit-q88_pool_bytes) for(size_t i=0;i<Q88_POOL_SLOTS;i++) if(!q88_pool[i].data) {
    q88_pool[i].data=data; q88_pool[i].bytes=bytes; q88_pool_bytes+=bytes;
    pthread_mutex_unlock(&q88_pool_lock);return;
  }
  pthread_mutex_unlock(&q88_pool_lock);free(data);
}
__attribute__((destructor)) static void q88_pool_cleanup(void) {
  for(size_t i=0;i<Q88_POOL_SLOTS;i++) free(q88_pool[i].data);
  if(getenv("REVAC_LABRADOR_SCRATCH_POOL_AUDIT")) fprintf(stderr,
      "[q88 scratch pool] retained_bytes=%zu hits=%lu misses=%lu limit=%zu\n",
      q88_pool_bytes,q88_pool_hits,q88_pool_misses,q88_pool_limit);
}
#endif
