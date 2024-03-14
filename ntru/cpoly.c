#include "cpoly.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "ipoly.h"
#include "../gen/const.h"
#include "../gen/ntru.h"

void cpoly_new(CPoly* p, int degree) {
	p->size  = degree/2;
	p->coefs = malloc(degree * sizeof(double));
}

void cpoly_log(CPoly* p) {
	for(int i = 0; i < p->size; i++) {
		printf(i == 0 ? "[" : ", ");
		printf("%.17g + %.17gi", p->coefs[i], p->coefs[i + p->size]);
	}
	printf("]\n");
}

// computes FFT of p, dividing coefficients by 2^exp
void fft(CPoly* cp, const IPoly* p, size_t exp) {
	assert(cp->size == p->degree/2);
	
	if(exp == 0) {
		for(int i = 0; i < p->degree; i++) {
			cp->coefs[i] = mpz_get_d(p->coefs[i]);
		}
	} else {
		mpz_t top;
		mpz_init(top);
		for(int i = 0; i < p->degree; i++) {
			mpz_tdiv_q_2exp(top, p->coefs[i], exp);
			cp->coefs[i] = mpz_get_d(top);
		}
		mpz_clear(top);
	}
	
	int max_d = cp->size*2;
	int d = 1;
	int s = max_d;
	int l = ANTRAG_D / max_d;
	while(d < max_d) {
		if(s % 2 == 0) {
			d *= 2;
			s /= 2;
		} else {
			d *= 3;
			s /= 3;
		}
		for(int o = 0; o < s; o++) {
			if(d == 2) {
				if(ANTRAG_LOG3_D) { // m = 6
					double a = cp->coefs[o], b = cp->coefs[o+s];
					double complex w = MTH_ROOTS[1*s*l];
					cpoly_set(cp, o, a + w*b);
				} else { // m = 4
					// evaluate in z=i -> no-op in this representation
				}

			} else if(d % 3 == 0) {
				double complex j = MTH_ROOTS[ANTRAG_M/3];
				for (int i = 0; i*6 < d; i++) {
					double complex w = MTH_ROOTS[ROOT_IDX[i*s*l*3]*s*l];
					double complex w2 = w * w;
					int k = o + 3*s*i;
					double complex a = cpoly_get(cp, k);
					double complex b = cpoly_get(cp, k+s) * w;
					double complex c = cpoly_get(cp, k+2*s) * w2;
					double complex bj = b*j, cj = c*j;
					double complex bj2 = bj*j, cj2 = cj*j;
					cpoly_set(cp, k,     a + b   + c);
					cpoly_set(cp, k+s,   a + bj  + cj2);
					cpoly_set(cp, k+2*s, a + bj2 + cj);
				}

			} else if(d % 2 == 0) {
				for (int i = 0; i*4 < d; i++) {
					double complex w = MTH_ROOTS[ROOT_IDX[i*s*l*2]*s*l];
					int k = o + 2*s*i;
					double complex a = cpoly_get(cp, k);
					double complex b = cpoly_get(cp, k+s) * w;
					cpoly_set(cp, k,   a + b);
					cpoly_set(cp, k+s, a - b);
				}
			}
		}
	}
}

void cpoly_free(CPoly *p) {
	free(p->coefs);
}

// computes inverse FFT of cp, multiplying coefficients by 2^exp and rounding the result
void ifft_rnd(IPoly* p, CPoly* cp, int64_t exp) {
	int max_d = cp->size*2;
	int d = max_d;
	int s = 1;
	int l = ANTRAG_D / max_d;
	while(d > 1) {
		for(int o = 0; o < s; o++) {
			if(d == 2) {
				if(ANTRAG_LOG3_D) { // m = 6
					double complex z = cpoly_get(cp, o);
					double complex w = MTH_ROOTS[1*s*l];
					cp->coefs[o+s] = cimag(z) / cimag(w);
					cp->coefs[o] = creal(z) - cp->coefs[o+s] * creal(w);
				} else { // m = 4
					// no-op
				}

			} else if(d % 3 == 0) {
				double complex j = MTH_ROOTS[ANTRAG_M/3];
				for (int i = 0; i*6 < d; i++) {
					double complex w = MTH_ROOTS[ROOT_IDX[i*s*l*3]*s*l];
					double complex w2 = w * w;
					int k = o + 3*s*i;
					double complex a = cpoly_get(cp, k);
					double complex b = cpoly_get(cp, k+s);
					double complex c = cpoly_get(cp, k+2*s);
					double complex bj = b*j, cj = c*j;
					double complex bj2 = bj*j, cj2 = cj*j;
					cpoly_set(cp, k,     (a + b   + c));
					cpoly_set(cp, k+s,   (a + bj2 + cj)  / w);
					cpoly_set(cp, k+2*s, (a + bj  + cj2) / w2);
				}

			} else if(d % 2 == 0) {
				for (int i = 0; i*4 < d; i++) {
					double complex w = MTH_ROOTS[ROOT_IDX[i*s*l*2]*s*l];
					int k = o + 2*s*i;
					double complex a = cpoly_get(cp, k);
					double complex b = cpoly_get(cp, k+s);
					cpoly_set(cp, k,   (a + b));
					cpoly_set(cp, k+s, (a - b) / w);
				}
			}
		}
		
		if(d % 3 == 0) {
			d /= 3;
			s *= 3;
		} else {
			d /= 2;
			s *= 2;
		}
	}

	ipoly_set_degree(p, cp->size * 2);
	for(int i = 0; i < p->degree; i++) {
		double res = cp->coefs[i] / cp->size;
		int exp2 = exp;
		if(exp > 0) {
			// Make sure we don't exceed DBL_MAX when left-shifting
			int my_exp;
			frexp(res, &my_exp);
			if(my_exp + exp2 > 1024) { // 1024 = Maximum possible exponent for double
				exp2 = 1024 - my_exp;
			}
		}
		if(exp != 0) {
			res = ldexp(res, exp2);
		}
		mpz_set_d(p->coefs[i], round(res));
		if(exp != 0 && exp > exp2) {
			mpz_mul_2exp(p->coefs[i], p->coefs[i], exp - exp2);
		}
	}
}

void adjoint(CPoly* p2, const CPoly* p) {
	assert(p2->size == p->size);
	for(int i = 0; i < p2->size; i++) {
		cpoly_set(p2, i, conj(cpoly_get(p,i)));
	}
}

// computes a*b + c*d
// the result is stored in a and returned for convenience
void cpoly_sum_pairs(CPoly* a, const CPoly* b, const CPoly* c, const CPoly* d) {
	assert(a->size == b->size && a->size == c->size && a->size == d->size);
	for(int i = 0; i < a->size; i++) {
		cpoly_set(a, i, cpoly_get(a,i)*cpoly_get(b,i) + cpoly_get(c,i)*cpoly_get(d,i));
	}
}

void cpoly_div(CPoly* f, const CPoly* g) {
	for(int i = 0; i < f->size; i++) {
		cpoly_set(f, i, cpoly_get(f,i) / cpoly_get(g,i));
	}
}
