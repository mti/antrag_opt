#ifndef NTRU_RNS_H
#define NTRU_RNS_H

#include <stdint.h>

#include <gmp.h>

#include "ipoly.h"

#include "../gen/rns.h"

// large integers are represented in a Residue Number System (RNS)
// mod_cnt is the number of moduli used in the representation

// short notation for a sized array parameter representing a number in RNS form
#define RNS(id) uint32_t (id)[mod_cnt]

// size of required scratch buffer (in uint32_t)
#define RNS_BUF_SIZE (4*MAX_MODULI)

// convert between large GNU MP integers and RNS form
void rns_from_mpz(int mod_cnt, RNS(ri), const mpz_t mi);
void rns_to_mpz(int cnt, mpz_t mi[cnt], int mod_cnt, const uint32_t ri[mod_cnt*cnt]);

// operations on RNS form integers
void rns_cpy(int mod_cnt, RNS(dst), const RNS(src));
void rns_add(int mod_cnt, RNS(dst), const RNS(a), const RNS(b));
void rns_sub(int mod_cnt, RNS(dst), const RNS(a), const RNS(b));
void rns_mul(int mod_cnt, RNS(dst), const RNS(a), const RNS(b));

void rns_log(int mod_cnt, const RNS(a));
bool rns_eq(int mod_cnt, const RNS(a), const RNS(b));

// polynomial with coefficients in RNS form
typedef struct {
	int deg; // number of coefficients
	int mod_cnt; // number of moduli per coefficient
	uint32_t* coefs; // residuals for a given coefficient are contiguous
} RnsPoly;

void rnspoly_from_ipoly(RnsPoly* rp, const IPoly* ip, int mod_cnt);
void rnspoly_to_ipoly(IPoly* ip, RnsPoly* rp);
void rnspoly_free(RnsPoly* rp);

void rnspoly_log(const RnsPoly* rp);

void rns_ntt(RnsPoly* p, uint32_t tmp[RNS_BUF_SIZE]);
void rns_intt(RnsPoly* p, uint32_t tmp[RNS_BUF_SIZE]);
void rns_elem_mul(RnsPoly* p1, const RnsPoly* p2);

#endif
