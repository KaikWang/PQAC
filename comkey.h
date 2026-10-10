#ifndef COMKEY_H
#define COMKEY_H

#include "polx.h"

extern polxvec comkey;
extern size_t comkey_nonce;

void comkey_init(size_t len);
void comkey_free();
void comkey_reset();

#endif
