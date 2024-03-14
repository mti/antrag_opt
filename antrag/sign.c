#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <math.h> 
#include "api.h"
#include "poly.h"
#include "param.h"
#include "fips202.h"
#include "randombytes.h"
#include "normaldist.h"
#include "samplerZ.h"

#include "../gen/const.h"

void H(size_t m_len, const uint8_t m[m_len], const uint8_t salt[SALT_BYTES], u16_poly* c1){
  uint64_t state[25];
  shake128_absorb(state, salt, SALT_BYTES);
  shake128_absorb(state, m, m_len);

  unsigned int out_idx = 0;
  while(out_idx < ANTRAG_D) 
  {
    uint8_t buf[SHAKE128_RATE];
    shake128_squeezeblocks(buf, 1, state);
    for(int i = 0; i < SHAKE128_RATE && out_idx < ANTRAG_D; i += 2)
    {
      uint16_t val = (buf[i] | ((uint16_t) buf[i+1] << 8));
      if(val < Q_MULT16)
      {
        c1->coeffs[out_idx] = val % ANTRAG_Q;
        out_idx++;
      }
    }
  }
}

static inline uint16_t decenter(int16_t x) {
  uint16_t xu = x;
  return xu + (ANTRAG_Q & -(xu >> 15));
}
static inline int16_t recenter(uint16_t x) {
  uint16_t t = x - ANTRAG_Q/2;
  return x - (ANTRAG_Q & ((t >> 15) - 1));
}

int check_norm(const u16_poly* p1, const u16_poly* p2){
  int64_t s = 0;
  for(int i = 0; i < ANTRAG_D; ++i) {
    int16_t x1 = recenter(p1->coeffs[i]);
    int16_t x2 = recenter(p2->coeffs[i]);
    s += x1*x1 + x2*x2;
  }
  #if ANTRAG_LOG3_D
  for(int i = 0; i < ANTRAG_D/2; ++i) {
    s += recenter(p1->coeffs[i]) * recenter(p1->coeffs[i + ANTRAG_D/2]);
    s += recenter(p2->coeffs[i]) * recenter(p2->coeffs[i + ANTRAG_D/2]);
  }
  #endif
  //printf("s=%ld / %ld\n", s, (int64_t) GAMMA_SQUARE);
  return s <= (int64_t) GAMMA_SQUARE;
}

void sampler(const secret_key* sk, const u16_poly* c2i, u16_poly* v1i, u16_poly* v2i) {
  poly c2, v1, v2;
  for(int i = 0; i < ANTRAG_D; i++) {
    c2.coeffs[i].v = recenter(c2i->coeffs[i]);
  }
  
  // offline
  poly y1, y2;
  normaldist(&y1); normaldist(&y2);
  //FFT(&y1); FFT(&y2); /* FFT unnecessary: just a scaling in normaldist */
  pointwise_mul(&y1, &(sk->sigma1)); pointwise_mul(&y2, &(sk->sigma2));

  poly d, temp_c1, temp_c2;
  poly* x;

  // online
  // first nearest plane
  for(int i=0; i < ANTRAG_D; ++i) temp_c1.coeffs[i].v = 0;
  set_poly(&temp_c2, &c2);
  FFT(&temp_c2);
  set_poly(&d, &temp_c2);
  pointwise_mul(&d, &(sk->beta21)); //d2 fft form

  x = &d; //x2
  poly_sub(x, &y2); // x2 = d2 - y2

  invFFT(x);
  sample_discrete_gauss(x);
  
  //second nearest plane

  FFT(x);
  set_poly(&v1, x); set_poly(&v2, x);
  pointwise_mul(&v1, &(sk->b20)); pointwise_mul(&v2, &(sk->b21));
  poly_sub(&temp_c1, &v1); poly_sub(&temp_c2, &v2); // Alg 13 line 14 c1 in the paper is [temp_c1,temp_c2] in the code

  poly temp;
  set_poly(&d, &temp_c1); set_poly(&temp, &temp_c2);
  pointwise_mul(&d, &(sk->beta10));
  pointwise_mul(&temp, &(sk->beta11));
  poly_add(&d, &temp); //d1 fft form

  x = &d; //x1
  poly_sub(x, &y1); // x1 = d1 - y1
  invFFT(x);
  sample_discrete_gauss(x);
  FFT(x);

  set_poly(&temp, x);
  pointwise_mul(&temp, &(sk->b10)); pointwise_mul(x, &(sk->b11));
  poly_add(&v1, &temp); poly_add(&v2, x); // v1 = v2 + x1*b1

  invFFT(&v1); invFFT(&v2);
  
  for(int i = 0; i < ANTRAG_D; i++) {
    v1i->coeffs[i] = decenter(round(v1.coeffs[i].v));
    v2i->coeffs[i] = decenter(round(v2.coeffs[i].v));
  }
}

void sign(size_t m_len, const uint8_t m[m_len], const secret_key* sk, signature* s){
  uint8_t salt[SALT_BYTES];
  u16_poly c2, v1, v2;
  
  do{
    randombytes(salt, SALT_BYTES);
    H(m_len, m, salt, &c2);
    
    sampler(sk, &c2, &v1, &v2);
    
    u16_poly_sub(&c2, &v2);
  } while (declassify(!check_norm(&v1, &c2)));
  
  //printf("# Signature generated after %d trials\n\n", trials);
  u16_poly_neg(&s->s1, &v1);
  for(int i=0; i < SALT_BYTES; ++i) (s->r)[i] = salt[i];
}

int verify(size_t m_len, const uint8_t m[m_len], const public_key* pk, const signature* s){
  u16_poly c1;
  H(m_len, m, s->r, &c1);

  u16_poly s2;
  memcpy(&s2, &s->s1, sizeof(u16_poly));
  u16_poly_mul(&s2, &pk->h_ntt);
  u16_poly_add(&s2, &c1);
  
  return check_norm(&s->s1, &s2);
}
