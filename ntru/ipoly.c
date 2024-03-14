#include "ipoly.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <complex.h>

#include <gmp.h>

#include "../gen/const.h"
#include "rns.h"

void ipoly_new(IPoly* p) {
	p->degree = 0;
	p->capacity = 0;
	p->coefs = NULL;
}

void ipoly_zeros(IPoly* p, int degree) {
	p->degree = degree;
	p->capacity = degree;
	p->coefs = malloc(sizeof(mpz_t) * degree);
	for(int i = 0; i < degree; i++) {
		mpz_init(p->coefs[i]);
	}
}

void ipoly_set_degree(IPoly* p, int degree) {
	if(degree > p->degree) {
		if(degree > p->capacity) {
			p->coefs = realloc(p->coefs, degree * sizeof(mpz_t));
			assert(p->coefs);
			p->capacity = degree;
		}
		for(int i = p->degree; i < degree; i++) {
			mpz_init(p->coefs[i]);
		}
	} else if(degree < p->degree) {
		for(int i = degree; i < p->degree; i++) {
			mpz_clear(p->coefs[i]);
		}
	}
	p->degree = degree;
}

void ipoly_shrink(IPoly* p) {
	if(p->degree < p->capacity) {
		p->capacity = p->degree;
		p->coefs = realloc(p->coefs, p->capacity * sizeof(mpz_t));
		assert(p->coefs);
	}
}

void ipoly_free(IPoly* p) {
	for(int i = 0; i < p->degree; i++) {
		mpz_clear(p->coefs[i]);
	}
	free(p->coefs);
}

void ipoly_copy(IPoly* p2, const IPoly* p) {
	ipoly_set_degree(p2, p->degree);
	for(int i = 0; i < p->degree; i++) {
		mpz_set(p2->coefs[i], p->coefs[i]);
	}
}

bool ipoly_is_zero(const IPoly* p) {
	for(int i = 0; i < p->degree; i++) {
		if(mpz_cmp_si(p->coefs[i], 0) != 0) {
			return false;
		}
	}
	return true;
}

void ipoly_log(const IPoly* p) {
	for(int i = 0; i < p->degree; i++) {
		printf(i == 0 ? "[" : ", ");
		mpz_out_str(stdout, 10, p->coefs[i]);
	}
	printf("]\n");
}

void ipoly_add(IPoly* f, const IPoly* g) {
	assert(f->degree == g->degree);
	for(int i = 0; i < f->degree; i++) {
		mpz_add(f->coefs[i], f->coefs[i], g->coefs[i]);
	}
}

void ipoly_sub(IPoly* f, const IPoly* g) {
	assert(f->degree == g->degree);
	for(int i = 0; i < f->degree; i++) {
		mpz_sub(f->coefs[i], f->coefs[i], g->coefs[i]);
	}
}

void ipoly_mul_naive(IPoly* h, const IPoly* f, const IPoly* g) {
	int n = f->degree;
	if(n == 1) {
		mpz_mul(h->coefs[0], f->coefs[0], g->coefs[0]);
		return;
	}
	
	for(int i = 0; i < h->degree; i++) {
		mpz_set_si(h->coefs[i], 0);
	}
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			mpz_addmul(h->coefs[i+j], f->coefs[i], g->coefs[j]);
		}
	}

	for(int i = h->degree-1; i >= n; i--) {
		mpz_sub(h->coefs[i-n], h->coefs[i-n], h->coefs[i]);
#if ANTRAG_LOG3_D
		mpz_add(h->coefs[i-n/2], h->coefs[i-n/2], h->coefs[i]);
#endif
	}
}

void ipoly_mul_rns(IPoly* h, const IPoly* f, const IPoly* g) {
	RnsPoly f_rns, g_rns;
	int mod_cnt = (ipoly_size(f) + ipoly_size(g)) / 30 + 2;
	rnspoly_from_ipoly(&f_rns, f, mod_cnt);
	rnspoly_from_ipoly(&g_rns, g, mod_cnt);
	uint32_t tmp[RNS_BUF_SIZE];
	rns_ntt(&f_rns, tmp);
	rns_ntt(&g_rns, tmp);
	rns_elem_mul(&f_rns, &g_rns);
	rns_intt(&f_rns, tmp);
	rnspoly_to_ipoly(h, &f_rns);
	rnspoly_free(&g_rns);
}

void ipoly_mul(IPoly* h, const IPoly* f, const IPoly* g) {
	assert(f->degree == g->degree);
	int n = f->degree;
	if(n <= 64) {
		ipoly_set_degree(h, 2*n);
		ipoly_mul_naive(h, f, g);
		ipoly_set_degree(h, n);
	} else {
		ipoly_set_degree(h, n);
		ipoly_mul_rns(h, f, g);
	}
}

// computes product of conjugates of p (degree d)
// relative to subfield of degree d/2
// (there is only one conjugate in this case)
void ipoly_conj2(IPoly* pc, const IPoly* p) {
	ipoly_copy(pc, p);
	for(int i = 1; i < p->degree; i += 2) {
		mpz_neg(pc->coefs[i], pc->coefs[i]);
	}
	if(p->degree == 2 && ANTRAG_LOG3_D) { // m = 6
		for(int i = 0; i < p->degree; i += 2) {
			mpz_add(pc->coefs[i], pc->coefs[i], p->coefs[i+1]);
		}
	}
}

// same thing, for subfield of degree d/3
// (two conjugates in this case)
void ipoly_conj3(IPoly* pc, const IPoly* p) {
	// j = x^(d/2) - 1 → j^2 = -x^(d/2)
	// pc1 = p0 - p1*x  + (p1*x - p2*x2)*x^(d/2)
	// pc2 = p0 - p2*x2 + (p2*x2 - p1*x)*x^(d/2)
	// (f + g*x^(d/2)) x^(d/2) = -g + (f+g)*x^(d/2)
	int n = p->degree;

	IPoly pc1, pc2;
	ipoly_zeros(&pc1, n);
	ipoly_zeros(&pc2, n);

	for(int i = 0; i < n; i += 3) {
		mpz_set(pc1.coefs[i], p->coefs[i]);
	}
	for(int i = 1; i < n/2; i += 3) {
		mpz_add(pc1.coefs[i], p->coefs[i], p->coefs[i + n/2]);
		mpz_neg(pc1.coefs[i], pc1.coefs[i]);
	}
	for(int i = n/2 + 1; i < n; i += 3) {
		mpz_set(pc1.coefs[i], p->coefs[i - n/2]);
	}
	for(int i = 2; i < n/2; i += 3) {
		mpz_set(pc1.coefs[i], p->coefs[i + n/2]);
	}
	for(int i = n/2 + 2; i < n; i += 3) {
		mpz_add(pc1.coefs[i], p->coefs[i - n/2], p->coefs[i]);
		mpz_neg(pc1.coefs[i], pc1.coefs[i]);
	}
	
	for(int i = 0; i < n; i += 3) {
		mpz_set(pc2.coefs[i], p->coefs[i]);
	}
	for(int i = 2; i < n/2; i += 3) {
		mpz_add(pc2.coefs[i], p->coefs[i], p->coefs[i + n/2]);
		mpz_neg(pc2.coefs[i], pc2.coefs[i]);
	}
	for(int i = n/2 + 2; i < n; i += 3) {
		mpz_set(pc2.coefs[i], p->coefs[i - n/2]);
	}
	for(int i = 1; i < n/2; i += 3) {
		mpz_set(pc2.coefs[i], p->coefs[i + n/2]);
	}
	for(int i = n/2 + 1; i < n; i += 3) {
		mpz_add(pc2.coefs[i], p->coefs[i - n/2], p->coefs[i]);
		mpz_neg(pc2.coefs[i], pc2.coefs[i]);
	}
	
	/*
	printf("calling ipoly_mul from ipoly_conj3 (degs: %d, %d, %d)\n",
			pc->degree, pc1.degree, pc2.degree);
	*/
	ipoly_mul(pc, &pc1, &pc2);
	ipoly_free(&pc1);
	ipoly_free(&pc2);
}

// reduce p (degree d) to subfield of degree d/ext
void ipoly_reduce(IPoly* p, int ext) {
	for(int i = 1; i < p->degree; i++) {
		if(i % ext == 0) {
			*p->coefs[i / ext] = *p->coefs[i];
		} else {
			assert(mpz_cmp_si(p->coefs[i], 0) == 0);
			mpz_clear(p->coefs[i]);
		}
	}
	p->degree /= ext;
	// not necessary, but helps debugging with address sanitizer
	ipoly_shrink(p);
}

// expand p (degree d) to extension field of degree d*ext
void ipoly_expand(IPoly* p, int ext) {
	p->coefs = realloc(p->coefs, p->degree * ext * sizeof(mpz_t));
	for(int i = p->degree; i < p->degree*ext; i++) {
		mpz_init(p->coefs[i]);
	}
	for(int i = p->degree-1; i >= 1; i--) {
		mpz_set(p->coefs[i*ext], p->coefs[i]);
		mpz_set_si(p->coefs[i], 0);
	}
	p->degree *= ext;
}

size_t ipoly_size(const IPoly* f) {
	size_t max = 0;
	for(int i = 0; i < f->degree; i++) {
		size_t bits = mpz_sizeinbase(f->coefs[i], 2);
		if(bits > max) max = bits;
	}
	return max;
}

size_t extra_bits(const IPoly* f, const IPoly* g, size_t max_bits) {
	size_t max = max_bits;
	size_t bits = ipoly_size(f);
	if(bits > max) max = bits;
	bits = ipoly_size(g);
	if(bits > max) max = bits;
	return max - max_bits;
}

bool ipoly_eq(const IPoly* f, const IPoly* g) {
	assert(f->degree == g->degree);
	for(int i = 0; i < f->degree; i++) {
		if(mpz_cmp(f->coefs[i], g->coefs[i]) != 0)
			return false;
	}
	return true;
}
