#include "../ntru/rns.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gmp.h>

#include "../gen/const.h"

void ipoly_mul_naive(IPoly* h, const IPoly* f, const IPoly* g);
void ipoly_mul_rns(IPoly* h, const IPoly* f, const IPoly* g);

int main() {
	struct timespec ts;
	timespec_get(&ts, TIME_UTC);
	srand(ts.tv_sec * 1000 + ts.tv_nsec / 1000000);
	
	int mod_cnt = MAX_MODULI;
	
	RNS(n_rns);
	for(int qi = 0; qi < mod_cnt; qi++) {
		n_rns[qi] = rand() % MODULI[qi].q;
	}
	
	mpz_t n;
	mpz_init(n);
	rns_to_mpz(1, &n, mod_cnt, n_rns);
	
	RNS(n_rns2);
	rns_from_mpz(mod_cnt, n_rns2, n);
	
	if(rns_eq(mod_cnt, n_rns, n_rns2)) {
		printf("RNS implementation roundtrips.\n");
	} else {
		printf("RNS implementation does not roundtrip!\n");
		return 1;
	}
	
	IPoly p, p2;
	int deg = ANTRAG_D;
	ipoly_zeros(&p, deg);
	ipoly_zeros(&p2, deg);
	for(int i = 0; i < deg; i++) {
		mpz_set_ui(p.coefs[i], rand() % 1000000000);
		mpz_set_ui(p2.coefs[i], rand() % 1000000000);
	}
	
	RnsPoly rp;
	rnspoly_from_ipoly(&rp, &p, mod_cnt);
	
	uint32_t tmp[RNS_BUF_SIZE];
	rns_ntt(&rp, tmp);
	rns_intt(&rp, tmp);
	
	IPoly p_rt;
	ipoly_zeros(&p_rt, deg);
	rnspoly_to_ipoly(&p_rt, &rp);
	
	if(ipoly_eq(&p, &p_rt)) {
		printf("NTT implementation roundtrips.\n");
	} else {
		printf("NTT implementation does not roundtrip!\n");
		return 1;
	}
	ipoly_free(&p_rt);
	
	IPoly p_prod_rns;
	ipoly_zeros(&p_prod_rns, deg);
	ipoly_mul_rns(&p_prod_rns, &p, &p2);
	
	IPoly p_prod_naive;
	ipoly_zeros(&p_prod_naive, 2*deg);
	ipoly_mul_naive(&p_prod_naive, &p, &p2);
	ipoly_set_degree(&p_prod_naive, deg);
	
	if(ipoly_eq(&p_prod_rns, &p_prod_naive)) {
		printf("NTT multiplication works.\n");
	} else {
		printf("NTT multiplication does not work!\n");
		return 1;
	}
	
	ipoly_free(&p);
	ipoly_free(&p2);
	ipoly_free(&p_prod_rns);
	ipoly_free(&p_prod_naive);
}