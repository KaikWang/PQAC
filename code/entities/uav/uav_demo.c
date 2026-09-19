#include <inttypes.h>
#include <stdio.h>

#include "arith.h"
#include "random.h"
#include "randombytes.h"
#include "revac_uav.h"

int main(void) {
  revac_uav_t uav;
  uint8_t seed[SEED_BYTES];
  uint64_t id = UINT64_C(0x1000000000000001);

  arith_setup();
  random_init();
  revac_uav_init(&uav);

  randombytes(seed, sizeof(seed));
  revac_uav_keygen(&uav, id, seed);
  uav.handle = revac_handle_from_identity(&uav.pk, uav.id);

  printf("UAV initialized\n");
  printf("id: 0x%016" PRIx64 "\n", uav.id);
  printf("handle: 0x%08" PRIx32 "\n", uav.handle);

  revac_uav_clear(&uav);
  arith_teardown();
  return 0;
}
