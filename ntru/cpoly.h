#ifndef NTRU_CPOLY_H
#define NTRU_CPOLY_H

#include "stdint.h"
#include "complex.h"

#include "ipoly.h"

// real polynomial in FFT representation stored as doubles
typedef struct {
	int size; // number of complex evaluations / dimension of original ring
	double* coefs; // first half is real parts, second half is complex parts
} CPoly;

void cpoly_new(CPoly* p, int degree);
void cpoly_log(CPoly* p);

static inline double complex cpoly_get(const CPoly* p, int i) {
	return p->coefs[i] + I*p->coefs[i + p->size];
}
static inline void cpoly_set(CPoly* p, int i, double complex z) {
	p->coefs[i] = creal(z);
	p->coefs[i + p->size] = cimag(z);
}

// computes FFT of p, dividing coefficients by 2^exp
void fft(CPoly* cp, const IPoly* p, size_t exp);

void cpoly_free(CPoly *p);

// computes inverse FFT of cp, multiplying coefficients by 2^exp and rounding the result
void ifft_rnd(IPoly* p, CPoly* cp, int64_t exp);

void adjoint(CPoly* p2, const CPoly* p);

// computes a*b + c*d
// the result is stored in a
void cpoly_sum_pairs(CPoly* a, const CPoly* b, const CPoly* c, const CPoly* d);

void cpoly_div(CPoly* f, const CPoly* g);


#endif
