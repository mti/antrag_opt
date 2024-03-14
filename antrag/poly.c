#include "poly.h"
#include "../falcon/inner.h"
#include "samplerZ.h"
#include <stdio.h>

void print_poly(const poly* p){
  int PRINT_SIZE = ANTRAG_D-1;
  for(int i=0; i < PRINT_SIZE; ++i) printf("%.1f, ", p->coeffs[i].v);
  printf("%.1f \n", p->coeffs[PRINT_SIZE].v);
}


void FFT(poly* p){
  Zf(FFT)(p->coeffs);

}

void invFFT(poly* p){
  Zf(iFFT)(p->coeffs);
}


void pointwise_mul(poly* p1, const poly* p2){
  /* Multiplication in FFT form, result in p1*/
  Zf(poly_mul_fft)(p1->coeffs, p2->coeffs);

}

void poly_add(poly* p1, const poly* p2){
  /* Addition result in p1*/
  Zf(poly_add)(p1->coeffs, p2->coeffs);
}

void poly_sub(poly* p1, const poly* p2){
  /* Subtraction p1-p2 result in p1*/
  Zf(poly_sub)(p1->coeffs, p2->coeffs);
}

void poly_div_FFT(poly* p1, const poly* p2){
  Zf(poly_div_fft)(p1->coeffs, p2->coeffs);
}

void FFT_adj(poly* p){
  Zf(poly_adj_fft)(p->coeffs);
}

void FFT_mul_adj(poly* p1, const poly* p2){
  Zf(poly_muladj_fft)(p1->coeffs, p2->coeffs);
}

void FFT_mul_selfadj(poly *p1){
  Zf(poly_mulselfadj_fft)(p1->coeffs);
}



void set_poly(poly* p1, const poly* p2){
  for(int i=0; i < ANTRAG_D; ++i)
    p1->coeffs[i].v = p2->coeffs[i].v;
}

void scalar_mul_FFT_form(double sigma2, poly* temp){
for(int i=0; i < ANTRAG_D; ++i)
    temp->coeffs[i].v = sigma2*(temp->coeffs[i].v);
}

void poly_recenter(poly* p) {
  for(int i=0; i < ANTRAG_D; ++i) {
    p->coeffs[i].v = fmod(p->coeffs[i].v, (double) ANTRAG_Q);
    if(p->coeffs[i].v > ANTRAG_Q/2)
      p->coeffs[i].v -= (double) ANTRAG_Q;
    else if(p->coeffs[i].v < -ANTRAG_Q/2)
	    p->coeffs[i].v += (double) ANTRAG_Q;
  }
}

void u16_poly_add(u16_poly* p1, const u16_poly* p2) {
  Zf(u16_poly_add)(p1->coeffs, p2->coeffs);
}
void u16_poly_sub(u16_poly* p1, const u16_poly* p2) {
  Zf(u16_poly_sub)(p1->coeffs, p2->coeffs);
}
void u16_poly_neg(u16_poly* p1, const u16_poly* p2) {
  Zf(u16_poly_neg)(p1->coeffs, p2->coeffs);
}
void u16_poly_mul(u16_poly* p1, const u16_poly* p2_ntt) {
  Zf(u16_poly_mul)(p1->coeffs, p2_ntt->coeffs);
}
