#include "rns.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include <gmp.h>

#include "../gen/const.h"

static inline uint32_t mq_montymul(uint32_t x, uint32_t y, int qi) {
	uint64_t q = MODULI[qi].q;
	uint64_t z = (uint64_t) x * (uint64_t) y;
	z += ((z * MODULI[qi].m_inv_q) & 0xFFFFFFFF) * q; // z += z * (-q^-1) % 2^32 * q
	z = z >> 32; // z is now 0 mod 2^32
	z -= q;
	z += q & -(z >> 63);
	return z;
}

static inline uint32_t to_monty(uint32_t x, int qi) {
	return mq_montymul(x, MODULI[qi].r2, qi);
}

static inline uint32_t from_monty(uint32_t x, int qi) {
	return mq_montymul(x, 1, qi);
}

static inline uint32_t mq_add(uint32_t x, uint32_t y, int qi) {
	uint64_t q = MODULI[qi].q;
	uint64_t d = (uint64_t) x + (uint64_t) y - q;
	d += MODULI[qi].q & -(d >> 63);
	return d;
}
static inline uint32_t mq_sub(uint32_t x, uint32_t y, int qi) {
	uint64_t d = (uint64_t) x - (uint64_t) y;
	d += MODULI[qi].q & -(d >> 63);
	return d;
}

void rns_from_mpz(int mod_cnt, RNS(ri), const mpz_t mi) {
	for(int qi = 0; qi < mod_cnt; qi++) {
		uint32_t r = mpz_fdiv_ui(mi, MODULI[qi].q);
		ri[qi] = to_monty(r, qi);
	}
}

void rns_to_mpz(int cnt, mpz_t mi[cnt], int mod_cnt, const uint32_t ri[mod_cnt*cnt]) {
	mpz_t n1;
	mpz_init_set_ui(n1, MODULI[0].q);
	for(int i = 0; i < cnt; i++) {
		mpz_set_ui(mi[i], from_monty(ri[i*mod_cnt], 0));
	}
	
	for(int qi = 1; qi < mod_cnt; qi++) {
		uint32_t n2 = MODULI[qi].q;
		uint32_t m1 = MODULI[qi].m1; // n1^-1 mod n2
		
		for(int i = 0; i < cnt; i++) {
			uint32_t a2 = ri[i*mod_cnt + qi]; // (montgomery)
			uint32_t a1mn2 = to_monty(mpz_fdiv_ui(mi[i], n2), qi); // a1 mod n2 (montgomery)
			uint32_t x = mq_montymul(m1, mq_sub(a2, a1mn2, qi), qi); // x = m1 (a2 - a1) mod n2
			mpz_addmul_ui(mi[i], n1, x); // a1 <- a1 + n1 x
			// = a1 (1 - n1 m1) + a2 n1 m1 = a1 n2 m2 + a2 n1 m1
		}
		
		mpz_mul_ui(n1, n1, n2);
	}
	
	mpz_t hn;
	mpz_init(hn);
	for(int i = 0; i < cnt; i++) {
		mpz_fdiv_r(mi[i], mi[i], n1); // mi = mi % n1
		
		mpz_fdiv_q_2exp(hn, n1, 1); // hn = n1 >> 1
		if(mpz_cmp(mi[i], hn) >= 0) {
			mpz_sub(mi[i], mi[i], n1);
		}
	}
	mpz_clears(hn, n1, NULL);
}

void rns_cpy(int mod_cnt, RNS(dst), const RNS(src)) {
	for(int qi = 0; qi < mod_cnt; qi++) {
		dst[qi] = src[qi];
	}
}
void rns_add(int mod_cnt, RNS(dst), const RNS(a), const RNS(b)) {
	for(int qi = 0; qi < mod_cnt; qi++) {
		dst[qi] = mq_add(a[qi], b[qi], qi);
	}
}
void rns_sub(int mod_cnt, RNS(dst), const RNS(a), const RNS(b)) {
	for(int qi = 0; qi < mod_cnt; qi++) {
		dst[qi] = mq_sub(a[qi], b[qi], qi);
	}
}
void rns_mul(int mod_cnt, RNS(dst), const RNS(a), const RNS(b)) {
	for(int qi = 0; qi < mod_cnt; qi++) {
		dst[qi] = mq_montymul(a[qi], b[qi], qi);
	}
}

void rns_root_pow(int mod_cnt, RNS(out), int pow, RNS(tmp)) {
	assert(pow >= 0);
	if(pow == 0) {
		for(int qi = 0; qi < mod_cnt; qi++) {
			out[qi] = 1;
		}
	} else {
		for(int qi = 0; qi < mod_cnt; qi++) {
			tmp[qi] = to_monty(1, qi);
			out[qi] = RNS_ROOT[qi];
		}
		while(pow > 1) {
			if(pow % 2 == 1) {
				rns_mul(mod_cnt, tmp, tmp, out);
			}
			rns_mul(mod_cnt, out, out, out);
			pow = pow / 2;
		}
		rns_mul(mod_cnt, out, out, tmp);
	}
}

void rns_log(int mod_cnt, const RNS(a)) {
	for(int qi = 0; qi < mod_cnt; qi++) {
		printf(qi == 0 ? "{" : ",");
		printf("%u", a[qi]);
	}
	printf("}");
}

bool rns_eq(int mod_cnt, const RNS(a), const RNS(b)) {
	for(int qi = 0; qi < mod_cnt; qi++) {
		if(a[qi] != b[qi])
			return false;
	}
	return true;
}

void rnspoly_from_ipoly(RnsPoly* rp, const IPoly* ip, int mod_cnt) {
	rp->deg = ip->degree;
	rp->mod_cnt = mod_cnt;
	rp->coefs = malloc(rp->deg * mod_cnt * sizeof(uint32_t));
	for(int i = 0; i < rp->deg; i++) {
		rns_from_mpz(mod_cnt, rp->coefs + i*mod_cnt, ip->coefs[i]);
	}
}

void rnspoly_to_ipoly(IPoly* ip, RnsPoly* rp) {
	assert(ip->degree == rp->deg);
	rns_to_mpz(rp->deg, ip->coefs, rp->mod_cnt, rp->coefs);
	free(rp->coefs);
}

void rnspoly_free(RnsPoly* rp) {
	free(rp->coefs);
}

void rnspoly_log(const RnsPoly* rp) {
	for(int i = 0; i < rp->deg; i++) {
		printf(i == 0 ? "[" : ", ");
		rns_log(rp->mod_cnt, rp->coefs + i*rp->mod_cnt);
	}
	printf("]\n");
}

void rns_ntt(RnsPoly* p, uint32_t tmp[RNS_BUF_SIZE]) {
	int mod_cnt = p->mod_cnt;
	int max_d = p->deg;
	int d = 1;
	int s = max_d;
	int l = ANTRAG_D / max_d;
	
	uint32_t *restrict t1 = tmp;
	uint32_t *restrict t2 = tmp + mod_cnt;
	uint32_t *restrict t3 = tmp + mod_cnt * 2;
	uint32_t *restrict t4 = tmp + mod_cnt * 3;
	
	while(d < max_d) {
		if(s % 2 == 0) {
			d *= 2;
			s /= 2;
			if(d == 2 && ANTRAG_LOG3_D) {
				uint32_t *w1 = t1, *w2 = t2;
				rns_root_pow(mod_cnt, w1, s*l, t3);
				rns_root_pow(mod_cnt, w2, ANTRAG_M - s*l, t3);
				for(int o = 0; o < s; o++) {
					uint32_t* pa = p->coefs + o*mod_cnt;
					uint32_t* pb = p->coefs + (o+s)*mod_cnt;
					
					uint32_t *bw1 = t3, *bw2 = t4;
					rns_mul(mod_cnt, bw1, pb, w1);
					rns_mul(mod_cnt, bw2, pb, w2);
					
					rns_add(mod_cnt, pb, pa, bw2);
					rns_add(mod_cnt, pa, pa, bw1);
				}
			} else {
				for(int i = 0; i*2 < d; i++) {
					uint32_t *w = t1;
					rns_root_pow(mod_cnt, w, ROOT_IDX[i*s*2*l]*s*l, t2);
					for(int o = 0; o < s; o++) {
						int k = o + 2*s*i;
						uint32_t* pa = p->coefs + k*mod_cnt;
						uint32_t* pb = p->coefs + (k+s)*mod_cnt;
						
						uint32_t *bw = t2;
						rns_mul(mod_cnt, bw, pb, w);
						
						rns_sub(mod_cnt, pb, pa, bw);
						rns_add(mod_cnt, pa, pa, bw);
					}
				}
			}
		} else {
			d *= 3;
			s /= 3;
			uint32_t* j = t1;
			rns_root_pow(mod_cnt, j, ANTRAG_M/3, t2);
			for(int i = 0; i*3 < d; i++) {
				uint32_t *w = t2;
				rns_root_pow(mod_cnt, w, ROOT_IDX[i*s*3*l]*s*l, t3);
				for(int o = 0; o < s; o++) {
					int k = o + 3*s*i;
					uint32_t* pa = p->coefs + k*mod_cnt;
					uint32_t* pb = p->coefs + (k+s)*mod_cnt;
					uint32_t* pc = p->coefs + (k+2*s)*mod_cnt;
					
					uint32_t *b = t3, *c = t4;
					rns_mul(mod_cnt, b, pb, w);
					rns_mul(mod_cnt, c, pc, w);
					rns_mul(mod_cnt, c, c, w);
					rns_cpy(mod_cnt, pb, pa);
					rns_cpy(mod_cnt, pc, pa);
					
					rns_add(mod_cnt, pa, pa, b);
					rns_add(mod_cnt, pa, pa, c);
					
					rns_mul(mod_cnt, b, b, j);
					rns_mul(mod_cnt, c, c, j);
					rns_add(mod_cnt, pb, pb, b);
					rns_add(mod_cnt, pc, pc, c);
					
					rns_mul(mod_cnt, b, b, j);
					rns_mul(mod_cnt, c, c, j);
					rns_add(mod_cnt, pb, pb, c);
					rns_add(mod_cnt, pc, pc, b);
				}
			}
		}
	}
}

void rns_intt(RnsPoly* p, uint32_t tmp[RNS_BUF_SIZE]) {
	int mod_cnt = p->mod_cnt;
	int max_d = p->deg;
	int d = max_d;
	int s = 1;
	int l = ANTRAG_D / max_d;
	
	uint32_t *restrict t1 = tmp;
	uint32_t *restrict t2 = tmp + mod_cnt;
	uint32_t *restrict t3 = tmp + mod_cnt * 2;
	uint32_t *restrict t4 = tmp + mod_cnt * 3;
	
	while(d > 1) {
		if(d % 3 == 0) {
			uint32_t *j = t1;
			rns_root_pow(mod_cnt, j, ANTRAG_M/3, t2);
			for(int i = 0; i*3 < d; i++) {
				uint32_t *inv_w = t2;
				rns_root_pow(mod_cnt, inv_w, ANTRAG_M - ROOT_IDX[i*s*3*l]*s*l, t3);
				for(int o = 0; o < s; o++) {
					int k = o + 3*s*i;
					uint32_t* pa = p->coefs + k*mod_cnt;
					uint32_t* pb = p->coefs + (k+s)*mod_cnt;
					uint32_t* pc = p->coefs + (k+2*s)*mod_cnt;
					
					uint32_t *b = t3, *c = t4;
					rns_cpy(mod_cnt, b, pb);
					rns_cpy(mod_cnt, c, pc);
					rns_cpy(mod_cnt, pb, pa);
					rns_cpy(mod_cnt, pc, pa);
					
					rns_add(mod_cnt, pa, pa, b);
					rns_add(mod_cnt, pa, pa, c);
					
					rns_mul(mod_cnt, b, b, j);
					rns_mul(mod_cnt, c, c, j);
					rns_add(mod_cnt, pb, pb, c);
					rns_add(mod_cnt, pc, pc, b);
					
					rns_mul(mod_cnt, b, b, j);
					rns_mul(mod_cnt, c, c, j);
					rns_add(mod_cnt, pb, pb, b);
					rns_add(mod_cnt, pc, pc, c);
					
					rns_mul(mod_cnt, pb, pb, inv_w);
					rns_mul(mod_cnt, pc, pc, inv_w);
					rns_mul(mod_cnt, pc, pc, inv_w);
					
					rns_mul(mod_cnt, pa, pa, RNS_INV3);
					rns_mul(mod_cnt, pb, pb, RNS_INV3);
					rns_mul(mod_cnt, pc, pc, RNS_INV3);
				}
			}
			d /= 3;
			s *= 3;
		} else {
			if(d == 2 && ANTRAG_LOG3_D) {
				uint32_t *w1 = t1;
				rns_root_pow(mod_cnt, w1, s*l, t2);
				for(int o = 0; o < s; o++) {
					uint32_t* pa = p->coefs + o*mod_cnt;
					uint32_t* pb = p->coefs + (o+s)*mod_cnt;
					
					uint32_t *c = t2;
					rns_sub(mod_cnt, c, pa, pb);
					rns_mul(mod_cnt, c, c, RNS_INV_ROOT6);
					rns_cpy(mod_cnt, pb, c);
					rns_mul(mod_cnt, c, c, w1);
					rns_sub(mod_cnt, pa, pa, c);
				}
			} else {
				for(int i = 0; i*2 < d; i++) {
					uint32_t *inv_w = t1;
					rns_root_pow(mod_cnt, inv_w, ANTRAG_M - ROOT_IDX[i*s*2*l]*s*l, t2);
					for(int o = 0; o < s; o++) {
						int k = o + 2*s*i;
						uint32_t* pa = p->coefs + k*mod_cnt;
						uint32_t* pb = p->coefs + (k+s)*mod_cnt;
						
						uint32_t *a = t2, *b = t3;
						rns_cpy(mod_cnt, a, pa);
						rns_cpy(mod_cnt, b, pb);
						
						rns_add(mod_cnt, pa, a, b);
						rns_sub(mod_cnt, pb, a, b);
						rns_mul(mod_cnt, pb, pb, inv_w);
						
						rns_mul(mod_cnt, pa, pa, RNS_INV2);
						rns_mul(mod_cnt, pb, pb, RNS_INV2);
					}
				}
			}
			d /= 2;
			s *= 2;
		}
	}
}

void rns_elem_mul(RnsPoly* p1, const RnsPoly* p2) {
	assert(p1->deg == p2->deg && p1->mod_cnt == p2->mod_cnt);
	int deg = p1->deg;
	int mod_cnt = p1->mod_cnt;
	for(int off = 0; off < deg * mod_cnt; off += mod_cnt) {
		rns_mul(mod_cnt, p1->coefs + off, p1->coefs + off, p2->coefs + off);
	}
}
