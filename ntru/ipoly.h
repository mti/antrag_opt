#ifndef NTRU_IPOLY_H
#define NTRU_IPOLY_H

#include <stdbool.h>

#include <stdio.h>
#include <gmp.h>

// polynomial with integer coefficients
typedef struct {
	int degree; // number of coefficients / degree of ring modulus
	int capacity; // size of allocated array
	mpz_t* coefs;
} IPoly;


void ipoly_new(IPoly* p);
void ipoly_zeros(IPoly* p, int degree);
void ipoly_set_degree(IPoly* p, int degree);
void ipoly_shrink(IPoly* p);
void ipoly_free(IPoly* p);
void ipoly_copy(IPoly* p2, const IPoly* p);
bool ipoly_is_zero(const IPoly* p);
void ipoly_log(const IPoly* p);
void ipoly_add(IPoly* f, const IPoly* g);
void ipoly_sub(IPoly* f, const IPoly* g);
void ipoly_mul(IPoly* h, const IPoly* f, const IPoly* g);

// computes product of conjugates of p (dimension d)
// relative to subfield of dimension d/2
// (there is only one conjugate in this case)
void ipoly_conj2(IPoly* pc, const IPoly* p);

// same thing, for subfield of dimension d/3
// (two conjugates in this case)
void ipoly_conj3(IPoly* pc, const IPoly* p);

// reduce p (dimension d) to subfield of dimension d/ext
void ipoly_reduce(IPoly* p, int ext);

// expand p (dimension d) to extension field of dimension d*ext
void ipoly_expand(IPoly* p, int ext);

size_t ipoly_size(const IPoly* f);
size_t extra_bits(const IPoly* f, const IPoly* g, size_t max_bits);

bool ipoly_eq(const IPoly* f, const IPoly* g);

#endif
