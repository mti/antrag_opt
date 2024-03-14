#include "../gen/const.h"
#include "normaldist.h"
#include "api.h"
#include "poly.h"
#include "randombytes.h"
#include <math.h>
#include <stdint.h>
#include <x86intrin.h>
#include <string.h>
#include <stdio.h>

#include "../ntru/ntru.h"

static void compute_GSO(secret_key* sk)
{
  poly temp1, temp2, temp3;

  for(int i=0; i < ANTRAG_D; ++i){
    (sk->GSO_b10).coeffs[i].v = (sk->b10).coeffs[i].v;
    (sk->GSO_b11).coeffs[i].v = (sk->b11).coeffs[i].v;
    (sk->GSO_b20).coeffs[i].v = (sk->b20).coeffs[i].v;
    (sk->GSO_b21).coeffs[i].v = (sk->b21).coeffs[i].v;
  }

  
  set_poly(&temp1, &(sk->b20)); set_poly(&temp2, &(sk->b21));

  FFT_mul_adj(&temp1, &(sk->b10)); FFT_mul_adj(&temp2, &(sk->b11));
  poly_add(&temp1, &temp2); // temp1 = <b1, b2> 
  
  set_poly(&temp2, &(sk->b10)); set_poly(&temp3, &(sk->b11));

  FFT_mul_selfadj(&temp2); FFT_mul_selfadj(&temp3);
  poly_add(&temp2, &temp3); // temp2 = <b1, b1>
  
  poly_div_FFT(&temp1, &temp2);
  set_poly(&temp2, &temp1); // temp2 = temp1 = <b1, b2>/<b1, b1>

  pointwise_mul(&temp1, &(sk->GSO_b10));
  pointwise_mul(&temp2, &(sk->GSO_b11)); // [temp1,temp2] = (<b1, b2>/<b1, b1>) * ~b1

  poly_sub(&(sk->GSO_b20), &temp1); poly_sub(&(sk->GSO_b21), &temp2);
}

static void compute_sigma(secret_key* sk)
{
  poly temp1, temp2, poly_r_square;

  for(int i=0; i < ANTRAG_D/2; ++i) {
    (sk->sigma1).coeffs[i].v = SIGMA_SQUARE;
    (sk->sigma2).coeffs[i].v = SIGMA_SQUARE;
    poly_r_square.coeffs[i].v = R_SQUARE;

  }
  for(int i=ANTRAG_D/2; i < ANTRAG_D; ++i) {
    poly_r_square.coeffs[i].v = 0;
    (sk->sigma1).coeffs[i].v = 0;
    (sk->sigma2).coeffs[i].v = 0;
  }


  set_poly(&temp1, &(sk->GSO_b10)); set_poly(&temp2, &(sk->GSO_b11));

  FFT_mul_selfadj(&temp1); FFT_mul_selfadj(&temp2);
  poly_add(&temp1, &temp2);

  poly_div_FFT(&(sk->sigma1), &temp1);
  poly_sub(&(sk->sigma1), &poly_r_square);

  for(int i=0; i < ANTRAG_D; ++i)
      (sk->sigma1).coeffs[i].v = sqrt((sk->sigma1).coeffs[i].v);

  set_poly(&temp1, &(sk->GSO_b20)); set_poly(&temp2, &(sk->GSO_b21));

  FFT_mul_selfadj(&temp1); FFT_mul_selfadj(&temp2);
  poly_add(&temp1, &temp2);
  poly_div_FFT(&(sk->sigma2), &temp1);
  poly_sub(&(sk->sigma2), &poly_r_square);

  for(int i=0; i < ANTRAG_D; ++i)
      (sk->sigma2).coeffs[i].v = sqrt((sk->sigma2).coeffs[i].v);
}

static void compute_beta_hat(secret_key* sk)
{
  poly temp1, temp2;
  // beta1
  set_poly(&temp1, &(sk->GSO_b10)); set_poly(&temp2, &(sk->GSO_b11));

  FFT_mul_selfadj(&temp1); FFT_mul_selfadj(&temp2);
  poly_add(&temp1, &temp2);

  set_poly(&(sk->beta10), &(sk->GSO_b10));
  set_poly(&(sk->beta11), &(sk->GSO_b11));
  FFT_adj(&(sk->beta10)); FFT_adj(&(sk->beta11));
  poly_div_FFT(&(sk->beta10), &temp1); poly_div_FFT(&(sk->beta11), &temp1);

  // beta2
  set_poly(&temp1, &(sk->GSO_b20)); set_poly(&temp2, &(sk->GSO_b21));

  FFT_mul_selfadj(&temp1); FFT_mul_selfadj(&temp2);
  poly_add(&temp1, &temp2);

  set_poly(&(sk->beta20), &(sk->GSO_b20));
  set_poly(&(sk->beta21), &(sk->GSO_b21));
  FFT_adj(&(sk->beta20)); FFT_adj(&(sk->beta21));
  poly_div_FFT(&(sk->beta20), &temp1); poly_div_FFT(&(sk->beta21), &temp1);
}

static void simple_frand(double *r, uint64_t *buf, size_t n) {
  static const double pow2m64 = 5.421010862427522e-20; // pow(2,-64);
  randombytes((uint8_t*)buf, n*sizeof(uint64_t));
  for(size_t i=0; i<n; i++) {
    r[i] = ((double)buf[i]) * pow2m64;
  }
}

/* 
 * Decode a real vector as an integer vector with odd sum.
 *
 * In other words: solve CVP in the non trivial coset of the D_n lattice.
 *
 * This is done using the decoder described in Conway & Sloane 20.2:
 * round all coefficients to the nearest integer, and if the sum is even,
 * round the worst coefficient (the one farthest from Z) in the other
 * direction.
 *
 * This implementation below is deliberately not constant time (since we
 * are not aiming for a constant-time keygen), but this is
 * straightforward to fix if deemed necessary.
 */
static void decode_odd(int8_t u[ANTRAG_D], const poly *utilde)
{
  uint8_t umod2 = 0;
  int8_t ui, wi = 0;
  size_t worst_coeff = 0;
  double maxdiff = -1, uitilde, diff;

  for(size_t i=0; i<ANTRAG_D; i++) {
    uitilde = utilde->coeffs[i].v;
    ui      = lrint(uitilde);
    umod2  ^= ui;
    diff    = fabs(uitilde - (double)ui);
    if(diff > maxdiff) {
      worst_coeff = i;
      maxdiff = diff;
      if(uitilde > (double)ui)
        wi = ui + 1;
      else
        wi = ui - 1;
    }
    u[i] = ui;
  }
  if((umod2 & 1) == 0)
    u[worst_coeff] = wi;
}


int keygen_fg(secret_key *sk)
{
  double af[ANTRAG_D/2], ag[ANTRAG_D/2], f[ANTRAG_D], g[ANTRAG_D];
  const double rad = sqrt((double)ANTRAG_Q) * ANTRAG_RADIUS,
    qlow  = ((double)ANTRAG_Q)/(ANTRAG_ALPHA*ANTRAG_ALPHA),
    qhigh = ((double)ANTRAG_Q)*ANTRAG_ALPHA*ANTRAG_ALPHA;

  double r[2*ANTRAG_D];
  uint64_t rint[2*ANTRAG_D];

  int trials = 0;
  int check = 1;

  do {
    trials++;
    simple_frand(r, rint, 3*ANTRAG_D/2);

    for(size_t i=0; i<ANTRAG_D/2; i++) {
      af[i]            = rad   * cos(M_PI/2*r[i]);
      ag[i]            = rad   * sin(M_PI/2*r[i]);
      f [i]            = af[i] * cos(2*M_PI*r[i+  ANTRAG_D/2]);
      f [i+ANTRAG_D/2] = af[i] * sin(2*M_PI*r[i+  ANTRAG_D/2]);
      g [i]            = ag[i] * cos(2*M_PI*r[i+2*ANTRAG_D/2]);
      g [i+ANTRAG_D/2] = ag[i] * sin(2*M_PI*r[i+2*ANTRAG_D/2]);
    }
    for(size_t i=0; i<ANTRAG_D; i++) {
      sk->b10.coeffs[i].v = f[i];
      sk->b11.coeffs[i].v = g[i];
    }
    invFFT(&sk->b10);
    invFFT(&sk->b11);

    decode_odd(sk->f, &(sk->b10));
    decode_odd(sk->g, &(sk->b11));

    for(size_t i=0; i<ANTRAG_D; i++) {
      sk->b10.coeffs[i].v = (double)sk->f[i];
      sk->b11.coeffs[i].v = (double)sk->g[i];
    }

    FFT(&sk->b10);
    FFT(&sk->b11);
    check = 0;
    for(size_t i=0; i<ANTRAG_D/2; i++) {
      double zi =
        fpr_sqr(sk->b10.coeffs[i]).v +
        fpr_sqr(sk->b10.coeffs[i+ANTRAG_D/2]).v +
        fpr_sqr(sk->b11.coeffs[i]).v +
        fpr_sqr(sk->b11.coeffs[i+ANTRAG_D/2]).v;
      if(zi < qlow || zi > qhigh) {
        check = 1;
        break;
      }
    }
  } while(check);
  return trials;
}

int keygen_full(secret_key *sk, public_key *pk)
{
  int trials = 0;

  uint8_t tmp[2*ANTRAG_D];
  
  while(1) {
    trials += keygen_fg(sk);

    if (!Zf(compute_public)(pk->h_ntt.coeffs, sk->f, sk->g, tmp))
      continue;
    
    if(solve_ntru(sk->f, sk->g, sk->F, sk->G))
      break;
  }

  for(int i=0; i<ANTRAG_D; i++) {
    sk->b20.coeffs[i].v = (double)sk->F[i];
    sk->b21.coeffs[i].v = (double)sk->G[i];
  }
  FFT(&sk->b20);
  FFT(&sk->b21);

  Zf(NTT)(pk->h_ntt.coeffs);

  compute_GSO(sk);
  compute_sigma(sk);
  compute_beta_hat(sk);

  return trials;
}
