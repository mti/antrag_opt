#include "randombytes.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

prng p;

#ifdef CST_TIME
int poisoned_rng = 0;
#endif

void randombytes(uint8_t *buf, size_t n){
  Zf(prng_get_bytes)(&p, buf, n);
  #ifdef CST_TIME
  if(poisoned_rng) {
    poison(buf, n);
  }
  #endif
}


uint64_t get64(){
  uint64_t v = prng_get_u64(&p);
  #ifdef CST_TIME
  if(poisoned_rng) {
    poison(&v, sizeof(uint64_t));
  }
  #endif
  return v;
} 

uint8_t get8(){
  uint8_t v = prng_get_u8(&p);
  #ifdef CST_TIME
  if(poisoned_rng) {
    poison(&v, sizeof(uint8_t));
  }
  #endif
  return v;
}

void seed_rng(void){
  inner_shake256_context sc;
  struct timespec t;

  clock_gettime(CLOCK_MONOTONIC, &t);

  inner_shake256_init(&sc);
  inner_shake256_inject(&sc, (uint8_t*)&t, sizeof(struct timespec));
  inner_shake256_flip(&sc);

  Zf(prng_init)(&p, &sc);
}

