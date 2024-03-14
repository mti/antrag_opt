/*
 * Falcon signature verification.
 *
 * ==========================(LICENSE BEGIN)============================
 *
 * Copyright (c) 2017-2019  Falcon Project
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * ===========================(LICENSE END)=============================
 *
 * @author   Thomas Pornin <thomas.pornin@nccgroup.com>
 */

#include "inner.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "../gen/ntt.h"

/*
 * Reduce a small signed integer modulo q. The source integer MUST
 * be between -q/2 and +q/2.
 */
static inline uint32_t mq_conv_small(int x)
{
	/*
	 * If x < 0, the cast to uint32_t will set the high bit to 1.
	 */
	uint32_t y;

	y = (uint32_t)x;
	y += ANTRAG_Q & -(y >> 31);
	return y;
}

/*
 * Addition modulo q. Operands must be in the 0..q-1 range.
 */
static inline uint32_t mq_add(uint32_t x, uint32_t y)
{
	/*
	 * We compute x + y - q. If the result is negative, then the
	 * high bit will be set, and 'd >> 31' will be equal to 1;
	 * thus '-(d >> 31)' will be an all-one pattern. Otherwise,
	 * it will be an all-zero pattern. In other words, this
	 * implements a conditional addition of q.
	 */
	uint32_t d;

	d = x + y - ANTRAG_Q;
	d += ANTRAG_Q & -(d >> 31);
	return d;
}

/*
 * Subtraction modulo q. Operands must be in the 0..q-1 range.
 */
static inline uint32_t mq_sub(uint32_t x, uint32_t y)
{
	/*
	 * As in mq_add(), we use a conditional addition to ensure the
	 * result is in the 0..q-1 range.
	 */
	uint32_t d;

	d = x - y;
	d += ANTRAG_Q & -(d >> 31);
	return d;
}

/* (Used in next function)
 * Montgomery multiplication modulo q. If we set R = 2^16 mod q, then
 * this function computes: x * y / R mod q
 * Operands must be in the 0..q-1 range.
 */
uint32_t mq_montymul(uint32_t x, uint32_t y)
{
	/*
	 * We compute x*y + k*q with a value of k chosen so that the 16
	 * low bits of the result are 0. We can then shift the value.
	 * After the shift, result may still be larger than q, but it
	 * will be lower than 2*q, so a conditional subtraction works.
	 */
	uint32_t z = x * y;
	uint32_t w = ((z * M_INV_Q) & 0xFFFF) * ANTRAG_Q;

	/*
	 * When adding z and w, the result will have its low 16 bits
	 * equal to 0. Since x, y and z are lower than q, the sum will
	 * be no more than (2^15 - 1) * q + (q - 1)^2, which will
	 * fit on 29 bits.
	 */
	z = (z + w) >> 16;

	/*
	 * After the shift, analysis shows that the value will be less
	 * than 2q. We do a subtraction then conditional addition to
	 * ensure the result is in the expected range.
	 */
	z -= ANTRAG_Q;
	z += ANTRAG_Q & -(z >> 31);
	return z;
}

/* Computes multiplicative inverse of x, by computing x^(q-2)
 * Input and output are in Montgomery representation.
 */
static inline uint32_t mq_montyinv(uint32_t x)
{
	uint32_t chain[INV_CHAIN_LEN] = { x };
	for (int i = 0; i < INV_CHAIN_LEN - 1; i++) {
		chain[i+1] = mq_montymul(chain[i], chain[INV_CHAIN[i]]);
	}
	return chain[INV_CHAIN_LEN-1];
}

static inline uint32_t mq_montydiv(uint32_t x, uint32_t y)
{
	return mq_montymul(x, mq_montyinv(y));
}

static inline uint16_t to_monty(uint16_t x)
{
	return mq_montymul(x, R2);
}

static inline uint16_t from_monty(uint16_t x)
{
	return mq_montymul(x, 1);
}

static inline fqx fqx_read(const uint16_t* v, int i)
{
	#if NTT_EXT
		return (fqx) { v[i], v[i+NTT_SIZE] };
	#else
		return v[i];
	#endif
}

static inline void fqx_write(uint16_t* v, int i, fqx x)
{
	#if NTT_EXT
		v[i] = x.c0;
		v[i+NTT_SIZE] = x.c1;
	#else
		v[i] = x;
	#endif
}

static inline fqx fqx_add(fqx a, fqx b)
{
	#if NTT_EXT
		return (fqx) { mq_add(a.c0, b.c0), mq_add(a.c1, b.c1) };
	#else
		return mq_add(a, b);
	#endif
}

static inline fqx fqx_sub(fqx a, fqx b)
{
	#if NTT_EXT
		return (fqx) { mq_sub(a.c0, b.c0), mq_sub(a.c1, b.c1) };
	#else
		return mq_sub(a, b);
	#endif
}

static inline fqx fqx_mul(fqx a, fqx b)
{
	#if NTT_EXT
		uint16_t c2 = mq_montymul(a.c1, b.c1);
		uint16_t c1 = mq_add(mq_montymul(a.c0, b.c1), mq_montymul(a.c1, b.c0));
		c1 = mq_sub(c1, mq_montymul(c2, EXT_REM1));
		uint16_t c0 = mq_montymul(a.c0, b.c0);
		c0 = mq_sub(c0, mq_montymul(c2, EXT_REM0));
		return (fqx) { c0, c1 };
	#else
		return mq_montymul(a, b);
	#endif
}


#if NTT_EXT
static inline fqx fqx_conj(fqx a)
{
	return (fqx) {
		mq_add(a.c0, mq_montymul(a.c1, EXT_CONJ0)),
		mq_montymul(a.c1, EXT_CONJ1)
	};
}
#endif

static inline fqx fqx_inv(fqx x)
{
	#if NTT_EXT
		// norm(x) = x conj(x) is in fq
		// x^-1 = conj(x) / norm(x)
		fqx conj = fqx_conj(x);
		uint16_t norm = fqx_mul(x, conj).c0;
		return (fqx) { mq_montydiv(conj.c0, norm), mq_montydiv(conj.c1, norm) };
	#else
		return mq_montyinv(x);
	#endif
}

/* Divide x by y in F_(q^k).
 * x and y are in Montgomery representation.
 */
static inline fqx fqx_div(fqx x, fqx y)
{
	return fqx_mul(x, fqx_inv(y));
}

static _Bool fqx_is_zero(fqx x)
{
	#if NTT_EXT
		return (x.c0 | x.c1) == 0;
	#else
		return x == 0;
	#endif
}

static fqx to_fqx(uint16_t x) {
	#if NTT_EXT
		return (fqx) { x, 0 };
	#else
		return x;
	#endif
}

void logic_error(const char* msg) {
	fprintf(stderr, "Logic error: %s\n", msg);
	abort();
}

static uint16_t from_fqx(fqx x) {
	#if NTT_EXT
		if(x.c1 != 0) logic_error("iNTT result not in F_q");
		return x.c0;
	#else
		return x;
	#endif
}

/*
 * Compute NTT on a ring element.
 * Result is in Montgomery representation.
 */
void Zf(NTT)(uint16_t v[ANTRAG_D])
{
	// Put in Montgomery representation
	for (int i = 0; i < ANTRAG_D; i++) {
		v[i] = to_monty(v[i]);
	}

	int d = 1;
	int s = ANTRAG_D; // ANTRAG_M/m = ANTRAG_D/d
	while (d < ANTRAG_D) {
		if (s % 2 == 0) {
			d *= 2;
			s /= 2;
		} else {
			d *= 3;
			s /= 3;
		}

		for (int o = 0; o < s; o++) { // offset
			if (d == 2) {
				fqx a = to_fqx(v[o]), b = to_fqx(v[o+s]);
				fqx w1 = MTH_ROOTS[1*s];
				fqx_write(v,o, fqx_add(a, fqx_mul(b, w1)));
				#if !NTT_EXT
					fqx w2 = MTH_ROOTS[ANTRAG_M-1*s];
					fqx_write(v,o+s, fqx_add(a, fqx_mul(b, w2)));
				#endif
				
			} else if (d % 3 == 0) {
				fqx j = MTH_ROOTS[ANTRAG_M/3];
				for (int i = 0; i*3*(NTT_EXT+1) < d; i++) {
					fqx w = MTH_ROOTS[ROOT_IDX[i*s*3]*s];
					fqx w2 = MTH_ROOTS[ROOT_IDX[i*s*3]*s * 2 % ANTRAG_M];
					int k = o + 3*s*i;
					fqx a = fqx_read(v,k);
					fqx b = fqx_mul(fqx_read(v,k+s), w);
					fqx c = fqx_mul(fqx_read(v,k+2*s), w2);
					fqx bj  = fqx_mul(b,j),  cj  = fqx_mul(c,j);
					fqx bj2 = fqx_mul(bj,j), cj2 = fqx_mul(cj,j);
					fqx_write(v,k,     fqx_add(a, fqx_add(b,   c)));
					fqx_write(v,k+s,   fqx_add(a, fqx_add(bj,  cj2)));
					fqx_write(v,k+2*s, fqx_add(a, fqx_add(bj2, cj)));
				}

			} else if (d % 2 == 0) {
				for (int i = 0; i*2*(NTT_EXT+1) < d; i++) {
					fqx w = MTH_ROOTS[ROOT_IDX[i*s*2]*s];
					int k = o + 2*s*i;
					fqx a = fqx_read(v,k);
					fqx b = fqx_mul(fqx_read(v,k+s), w);
					fqx_write(v,k,   fqx_add(a, b));
					fqx_write(v,k+s, fqx_sub(a, b));
				}

			} else {
				logic_error("D not even or multiple of 3");
			}
		}
	}
}

/*
 * Compute the inverse NTT on a ring element.
 */
void Zf(iNTT)(uint16_t v[ANTRAG_D])
{
	int d = ANTRAG_D;
	int s = 1;
	while (d > 1) {
		for (int o = 0; o < s; o++) { // offset
			if (d == 2) {
				fqx a = fqx_read(v,o);
				#if NTT_EXT
					fqx b = fqx_conj(a);
				#else
					fqx b = fqx_read(v,o+s);
				#endif
				fqx w1 = MTH_ROOTS[1*s];
				fqx c = fqx_mul(fqx_sub(a, b), INV_ROOT_DIFF);
				v[o] = from_fqx(fqx_sub(a, fqx_mul(w1, c)));
				v[o+s] = from_fqx(c);

			} else if (d % 3 == 0) {
				fqx j = MTH_ROOTS[ANTRAG_M/3];
				for (int i = 0; i*3*(NTT_EXT+1) < d; i++) {
					fqx inv_w = MTH_ROOTS[ANTRAG_M - ROOT_IDX[i*s*3]*s];
					fqx inv_w2 = MTH_ROOTS[ANTRAG_M - ROOT_IDX[i*s*3]*s * 2 % ANTRAG_M];
					int k = o + 3*s*i;
					fqx a = fqx_read(v,k);
					fqx b = fqx_read(v,k+s);
					fqx c = fqx_read(v,k+2*s);
					fqx bj  = fqx_mul(b,j),  cj  = fqx_mul(c,j);
					fqx bj2 = fqx_mul(bj,j), cj2 = fqx_mul(cj,j);
					fqx_write(v,k,     fqx_add(a, fqx_add(b, c)));
					fqx_write(v,k+s,   fqx_mul(fqx_add(a, fqx_add(bj2, cj)), inv_w));
					fqx_write(v,k+2*s, fqx_mul(fqx_add(a, fqx_add(bj, cj2)), inv_w2));
				}

			} else if (d % 2 == 0) {
				for (int i = 0; i*2*(NTT_EXT+1) < d; i++) {
					fqx inv_w = MTH_ROOTS[ANTRAG_M - ROOT_IDX[i*s*2]*s];
					int k = o + 2*s*i;
					fqx a = fqx_read(v,k);
					fqx b = fqx_read(v,k+s);
					fqx_write(v,k,   fqx_add(a, b));
					fqx_write(v,k+s, fqx_mul(fqx_sub(a, b), inv_w));
				}

			} else {
				logic_error("D not even or multiple of 3");
			}
		}

		if (d % 3 == 0) {
			d /= 3;
			s *= 3;
		} else {
			d /= 2;
			s *= 2;
		}
	}

	// Take out of Montgomery representation and remove D/2 factor
	uint16_t f = mq_montyinv(to_monty(ANTRAG_D/2));
	for (int i = 0; i < ANTRAG_D; i++) {
		v[i] = from_monty(mq_montymul(v[i], f));
	}
}

/* ===================================================================== */

void Zf(u16_poly_add)(uint16_t f[restrict ANTRAG_D], const uint16_t g[restrict ANTRAG_D]) {
	for(int i = 0; i < ANTRAG_D; i++) {
		f[i] = mq_add(f[i], g[i]);
	}
}
void Zf(u16_poly_sub)(uint16_t f[restrict ANTRAG_D], const uint16_t g[restrict ANTRAG_D]) {
	for(int i = 0; i < ANTRAG_D; i++) {
		f[i] = mq_sub(f[i], g[i]);
	}
}
void Zf(u16_poly_neg)(uint16_t f[restrict ANTRAG_D], const uint16_t g[restrict ANTRAG_D]) {
	for(int i = 0; i < ANTRAG_D; i++) {
		f[i] = mq_sub(0, g[i]);
	}
}

void Zf(u16_poly_mul)(uint16_t f[restrict ANTRAG_D], const uint16_t g_ntt[restrict ANTRAG_D]) {
	Zf(NTT)(f);
	for(int i = 0; i < NTT_SIZE; i++) {
		fqx_write(f,i, fqx_mul(fqx_read(f,i), fqx_read(g_ntt,i)));
	}
	Zf(iNTT)(f);
}

/* see inner.h */
int Zf(compute_public)(uint16_t *h, const int8_t *f, const int8_t *g, uint8_t *tmp)
{
	uint16_t* tt = (uint16_t*) tmp;
	for (size_t u = 0; u < ANTRAG_D; u ++) {
		tt[u] = (uint16_t) mq_conv_small(f[u]);
		h[u] = (uint16_t) mq_conv_small(g[u]);
	}

	Zf(NTT)(h);
	Zf(NTT)(tt);
	for (size_t u = 0; u < NTT_SIZE; u++) {
		if (fqx_is_zero(fqx_read(tt,u))) {
			return 0;
		}
		fqx_write(h,u, fqx_div(fqx_read(h,u), fqx_read(tt,u)));
	}
	Zf(iNTT)(h);

	return 1;
}
