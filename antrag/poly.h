#ifndef POLY_H
#define POLY_H

#include "param.h"
#include "../falcon/fpr.h"

typedef struct{
  fpr coeffs[ANTRAG_D];
} poly;

typedef struct {
  uint16_t coeffs[ANTRAG_D];
} u16_poly;

typedef struct{
  int8_t f[ANTRAG_D];
  int8_t g[ANTRAG_D];
  int8_t F[ANTRAG_D];
  int8_t G[ANTRAG_D];
  poly b10;
  poly b11;
  poly b20;
  poly b21;
  poly GSO_b10; //~b1[0]/<~b1, ~b1>
  poly GSO_b11; //~b1[1]/<~b1, ~b1>
  poly GSO_b20; //~b2[0]/<~b2, ~b2>
  poly GSO_b21; //~b2[1]/<~b2, ~b2>
  poly beta10;
  poly beta11;
  poly beta20;
  poly beta21;
  poly sigma1;
  poly sigma2;
} secret_key;

typedef struct{
  u16_poly h_ntt;
} public_key;

typedef struct{
  u16_poly s1;
  uint8_t r[SALT_BYTES];
} signature;


void print_poly(const poly* p);
void FFT(poly* p);
void invFFT(poly* p);
void pointwise_mul(poly* p1, const poly* p2);
void FFT_adj(poly* p);
void FFT_mul_adj(poly* p1, const poly* p2);
void FFT_mul_selfadj(poly *p1);
void poly_div_FFT(poly* p1, const poly* p2);
 

void set_poly(poly* p1, const poly* p2);
void poly_add(poly* p1, const poly* p2);
void poly_sub(poly* p1, const poly* p2);

void poly_recenter(poly* p);

void u16_poly_add(u16_poly* p1, const u16_poly* p2);
void u16_poly_sub(u16_poly* p1, const u16_poly* p2);
void u16_poly_neg(u16_poly* p1, const u16_poly* p2);
void u16_poly_mul(u16_poly* p1, const u16_poly* p2_ntt);

#endif
