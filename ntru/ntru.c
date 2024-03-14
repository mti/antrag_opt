#include "ntru.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <math.h>

#include <gmp.h>

#include "../gen/ntru.h"

#include "ipoly.h"
#include "cpoly.h"

bool reduce_babai(const IPoly* f, const IPoly* g, IPoly* F, IPoly* G) {
	/*
	printf("reduce_babai start. sizes: |f|=%zu, |g|=%zu, |F|=%zu, |G|=%zu.\n",
			ipoly_size(f), ipoly_size(g), ipoly_size(F), ipoly_size(G));
	*/
	// compute constants
	size_t exp1 = extra_bits(f, g, 500);
	CPoly den, cg, cfa, cga;

	cpoly_new(&den, f->degree);
	cpoly_new(&cfa, f->degree);
	cpoly_new(&cg,  g->degree);
	cpoly_new(&cga, g->degree);

	fft(&den, f, exp1);
	fft(&cg,  g, exp1);
	adjoint(&cfa, &den);
	adjoint(&cga, &cg);
	cpoly_sum_pairs(&den, &cfa, &cg, &cga);
	cpoly_free(&cg);

	CPoly num, cG;
	cpoly_new(&num, F->degree);
	cpoly_new(&cG,  G->degree);

	IPoly k, kf, kg;
	ipoly_new(&k);
	ipoly_new(&kf);
	ipoly_new(&kg);

	int initial_F_size = ipoly_size(F);
	
	bool success = true;
	for(int iter = 0;; iter++) {
		size_t exp2 = extra_bits(F, G, 500);
		fft(&num, F, exp2);
		fft(&cG,  G, exp2);

		cpoly_sum_pairs(&num, &cfa, &cG, &cga);
		cpoly_div(&num, &den);

		ifft_rnd(&k, &num, (ssize_t) exp2 - (ssize_t) exp1);
		
		if(ipoly_is_zero(&k)) {
			break;
		}

		ipoly_mul(&kf, &k, f);
		ipoly_mul(&kg, &k, g);

		ipoly_sub(F, &kf);
		ipoly_sub(G, &kg);

		if(iter > initial_F_size*2) {
			#ifdef DEBUG_NTRU
				printf("babai reduction does not converge\n");
			#endif
			success = false;
			break;
		}
	}
	
	cpoly_free(&cfa); cpoly_free(&cga);
	cpoly_free(&den);
	
	cpoly_free(&num); cpoly_free(&cG);
	ipoly_free(&kf); ipoly_free(&kg); ipoly_free(&k);

	return success;
}


bool tower_solver(const IPoly* f, const IPoly* g, int d, int m, IPoly* F, IPoly* G) {
	if(d == 1) {
		mpz_ptr a = f->coefs[0];
		mpz_t b;
		mpz_init(b);
		mpz_neg(b, g->coefs[0]);
		
		mpz_t gcd, u, v;
		mpz_inits(gcd, u, v, NULL);
		mpz_gcdext(gcd, u, v, a, b);
		long gcd2 = mpz_get_ui(gcd);
		if(!mpz_fits_ulong_p(gcd) || gcd2 != 1) {
			mpz_clears(gcd, u, v, b, NULL);
			#ifdef DEBUG_NTRU
				printf("gcd(f,g) != 1\n");
			#endif
			return false;
		}
		mpz_clears(gcd, b, NULL);
		
		long k = ANTRAG_Q / gcd2;

		ipoly_set_degree(F, 1);
		ipoly_set_degree(G, 1);
		mpz_mul_si(F->coefs[0], v, k);
		mpz_mul_si(G->coefs[0], u, k);
		
		mpz_clears(u, v, NULL);

		return true;
		
	} else {
		int ext, d2, m2;
		if(m == 6) {
			ext = 2;
			d2 = d / 2;
			m2 = m / 3;
		} else {
			ext = d % 3 == 0 ? 3 : 2;
			d2 = d / ext;
			m2 = m / ext;
		}
		
		// compute product of conjugates
		IPoly fc, gc;
		ipoly_new(&fc);
		ipoly_new(&gc);
		if(ext == 2) {
			ipoly_conj2(&fc, f);
			ipoly_conj2(&gc, g);
		} else if(ext == 3) {
			ipoly_conj3(&fc, f);
			ipoly_conj3(&gc, g);
		}
		
		// compute algebraic norm of f and g
		IPoly fn, gn;
		ipoly_new(&fn);
		ipoly_new(&gn);
		/*
		printf("calling ipoly_mul from tower_solver (degs: %d, %d, %d)\n",
				fn.degree, f->degree, fc.degree);
		*/
		ipoly_mul(&fn, f, &fc);
		ipoly_mul(&gn, g, &gc);
		ipoly_reduce(&fn, ext);
		ipoly_reduce(&gn, ext);
		
		// recurse into subfield
		bool res = tower_solver(&fn, &gn, d2, m2, F, G);
		if(!res) {
			ipoly_free(&fn);
			ipoly_free(&gn);
			ipoly_free(&fc);
			ipoly_free(&gc);
			return false;
		}
		
		// lift solutions
		ipoly_expand(F, ext);
		ipoly_expand(G, ext);
		
		/*
		printf("calling ipoly_mul from tower_solver (degs: %d, %d, %d)\n",
				fn.degree, F->degree, gc.degree);
		*/
		ipoly_mul(&fn, F, &gc);
		ipoly_mul(&gn, G, &fc);

		ipoly_copy(F, &fn);
		ipoly_copy(G, &gn);

		ipoly_free(&fn);
		ipoly_free(&gn);
		ipoly_free(&fc);
		ipoly_free(&gc);
		
		if(!reduce_babai(f, g, F, G)) {
			return false;
		}
		
		return true;
	}
}

bool solve_ntru(const int8_t fa[ANTRAG_D], const int8_t ga[ANTRAG_D], int8_t Fa[ANTRAG_D], int8_t Ga[ANTRAG_D]) {
	IPoly f, g, F, G;
	ipoly_zeros(&f, ANTRAG_D);
	ipoly_zeros(&g, ANTRAG_D);
	ipoly_new(&F);
	ipoly_new(&G);

	for(int i = 0; i < ANTRAG_D; i++) {
		mpz_set_si(f.coefs[i], fa[i]);
		mpz_set_si(g.coefs[i], ga[i]);
	}
	
	if(!tower_solver(&f, &g, ANTRAG_D, ANTRAG_M, &F, &G)) {
		ipoly_free(&f);
		ipoly_free(&g);
		ipoly_free(&F);
		ipoly_free(&G);
		return false;
	}
	
	#ifndef NDEBUG
		IPoly res, gF;
		ipoly_new(&res);
		ipoly_new(&gF);
		ipoly_mul(&res, &f, &G);
		ipoly_mul(&gF,  &g, &F);
		ipoly_sub(&res, &gF);
		ipoly_free(&gF);

		for(int i = 0; i < ANTRAG_D; i++) {
			assert(mpz_cmp_si(res.coefs[i], i == 0 ? ANTRAG_Q : 0) == 0);
		}
		ipoly_free(&res);
	#endif

	ipoly_free(&f);
	ipoly_free(&g);
	
	for(int i = 0; i < ANTRAG_D; i++) {
		assert(mpz_fits_slong_p(F.coefs[i]) && mpz_fits_slong_p(G.coefs[i]));
		long Fi = mpz_get_si(F.coefs[i]);
		long Gi = mpz_get_si(G.coefs[i]);
		if(Fi < INT8_MIN || Fi >= INT8_MAX || Gi < INT8_MIN || Gi >= INT8_MAX) {
			// F/G coefficients are too large
			// Happens rarely for standard Q values
			return false;
		}
		Fa[i] = Fi;
		Ga[i] = Gi;
	}
	ipoly_free(&F);
	ipoly_free(&G);
	
	return true;
}

/*
vim: ts=4
*/
