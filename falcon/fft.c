/*
 * FFT code.
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

#include "../gen/fft.h"

/*
 * Rules for complex number macros:
 * --------------------------------
 *
 * Operand order is: destination, source1, source2...
 *
 * Each operand is a real and an imaginary part.
 *
 * All overlaps are allowed.
 */

/*
 * Addition of two complex numbers (d = a + b).
 */
#define FPC_ADD(d_re, d_im, a_re, a_im, b_re, b_im)   do { \
		fpr fpct_re, fpct_im; \
		fpct_re = fpr_add(a_re, b_re); \
		fpct_im = fpr_add(a_im, b_im); \
		(d_re) = fpct_re; \
		(d_im) = fpct_im; \
	} while (0)

/*
 * Subtraction of two complex numbers (d = a - b).
 */
#define FPC_SUB(d_re, d_im, a_re, a_im, b_re, b_im)   do { \
		fpr fpct_re, fpct_im; \
		fpct_re = fpr_sub(a_re, b_re); \
		fpct_im = fpr_sub(a_im, b_im); \
		(d_re) = fpct_re; \
		(d_im) = fpct_im; \
	} while (0)

/*
 * Multplication of two complex numbers (d = a * b).
 */
#define FPC_MUL(d_re, d_im, a_re, a_im, b_re, b_im)   do { \
		fpr fpct_a_re, fpct_a_im; \
		fpr fpct_b_re, fpct_b_im; \
		fpr fpct_d_re, fpct_d_im; \
		fpct_a_re = (a_re); \
		fpct_a_im = (a_im); \
		fpct_b_re = (b_re); \
		fpct_b_im = (b_im); \
		fpct_d_re = fpr_sub( \
			fpr_mul(fpct_a_re, fpct_b_re), \
			fpr_mul(fpct_a_im, fpct_b_im)); \
		fpct_d_im = fpr_add( \
			fpr_mul(fpct_a_re, fpct_b_im), \
			fpr_mul(fpct_a_im, fpct_b_re)); \
		(d_re) = fpct_d_re; \
		(d_im) = fpct_d_im; \
	} while (0)

/*
 * Squaring of a complex number (d = a * a).
 */
#define FPC_SQR(d_re, d_im, a_re, a_im)   do { \
		fpr fpct_a_re, fpct_a_im; \
		fpr fpct_d_re, fpct_d_im; \
		fpct_a_re = (a_re); \
		fpct_a_im = (a_im); \
		fpct_d_re = fpr_sub(fpr_sqr(fpct_a_re), fpr_sqr(fpct_a_im)); \
		fpct_d_im = fpr_double(fpr_mul(fpct_a_re, fpct_a_im)); \
		(d_re) = fpct_d_re; \
		(d_im) = fpct_d_im; \
	} while (0)

/*
 * Inversion of a complex number (d = 1 / a).
 */
#define FPC_INV(d_re, d_im, a_re, a_im)   do { \
		fpr fpct_a_re, fpct_a_im; \
		fpr fpct_d_re, fpct_d_im; \
		fpr fpct_m; \
		fpct_a_re = (a_re); \
		fpct_a_im = (a_im); \
		fpct_m = fpr_add(fpr_sqr(fpct_a_re), fpr_sqr(fpct_a_im)); \
		fpct_m = fpr_inv(fpct_m); \
		fpct_d_re = fpr_mul(fpct_a_re, fpct_m); \
		fpct_d_im = fpr_mul(fpr_neg(fpct_a_im), fpct_m); \
		(d_re) = fpct_d_re; \
		(d_im) = fpct_d_im; \
	} while (0)

/*
 * Division of complex numbers (d = a / b).
 */
#define FPC_DIV(d_re, d_im, a_re, a_im, b_re, b_im)   do { \
		fpr fpct_a_re, fpct_a_im; \
		fpr fpct_b_re, fpct_b_im; \
		fpr fpct_d_re, fpct_d_im; \
		fpr fpct_m; \
		fpct_a_re = (a_re); \
		fpct_a_im = (a_im); \
		fpct_b_re = (b_re); \
		fpct_b_im = (b_im); \
		fpct_m = fpr_add(fpr_sqr(fpct_b_re), fpr_sqr(fpct_b_im)); \
		fpct_m = fpr_inv(fpct_m); \
		fpct_b_re = fpr_mul(fpct_b_re, fpct_m); \
		fpct_b_im = fpr_mul(fpr_neg(fpct_b_im), fpct_m); \
		fpct_d_re = fpr_sub( \
			fpr_mul(fpct_a_re, fpct_b_re), \
			fpr_mul(fpct_a_im, fpct_b_im)); \
		fpct_d_im = fpr_add( \
			fpr_mul(fpct_a_re, fpct_b_im), \
			fpr_mul(fpct_a_im, fpct_b_re)); \
		(d_re) = fpct_d_re; \
		(d_im) = fpct_d_im; \
	} while (0)


// Convenience functions to simplify writing complex arithmetic expressions

static inline fpc fpc_read(fpr f[ANTRAG_D], int i) {
	return (fpc) { f[i], f[i+FFT_SIZE] };
}
static inline void fpc_write(fpr f[ANTRAG_D], int i, fpc x) {
	f[i] = x.re;
	f[i+FFT_SIZE] = x.im;
}

static inline fpc fpc_add(fpc x, fpc y) {
	fpc z;
	FPC_ADD(z.re, z.im, x.re, x.im, y.re, y.im);
	return z;
}
static inline fpc fpc_sub(fpc x, fpc y) {
	fpc z;
	FPC_SUB(z.re, z.im, x.re, x.im, y.re, y.im);
	return z;
}
static inline fpc fpc_mul(fpc x, fpc y) {
	fpc z;
	FPC_MUL(z.re, z.im, x.re, x.im, y.re, y.im);
	return z;
}
static inline fpc fpc_div(fpc x, fpc y) {
	fpc z;
	FPC_DIV(z.re, z.im, x.re, x.im, y.re, y.im);
	return z;
}


/*
 * Let w = exp(i*pi/M); w is a primitive M-th root of 1. We define the
 * values w_j = w^k for all j from 0 to D-1, where k is the (j+1)-th natural
 * prime with M: these are all the primitive M-th roots of 1, ie. the roots
 * of Phi_M in the field of complex numbers. A crucial property is that
 * w_{D-1-j} = conj(w_j) = 1/w_j for all j.
 *
 * FFT representation of a polynomial f (taken modulo X^N+1) is the
 * set of values f(w_j). Since f is real, conj(f(w_j)) = f(conj(w_j)),
 * thus f(w_{N-1-j}) = conj(f(w_j)). We thus store only half the values,
 * specifically the ones where j is even; the other half can be recomputed
 * easily when (if) needed. A consequence is that FFT representation has
 * the same size as normal representation: FFT_SIZE = D/2 complex numbers
 * use D real numbers (each complex number is the combination of a real and an
 * imaginary part).
 *
 * We use a specific ordering which makes computations easier. Let rev()
 * be the digit-reversal function over ANTRAG_LOG3_D trits followed by
 * ANTRAG_LOG2_D-1 bits. For l in 0..D/2-1, we store the real and imaginary
 * parts of f(w_j) in slots:
 *    Re(f(w_2l)) -> slot rev(l)
 *    Im(f(w_2l)) -> slot rev(l)+D/2
 */

/* see inner.h */
TARGET_AVX2
void
Zf(FFT)(fpr f[ANTRAG_D]) {
	int d = 1;
	int s = ANTRAG_D; // step size = ANTRAG_D/d = ANTRAG_M/m
	while (d < ANTRAG_D) {
		if (s % 2 == 0) {
			d *= 2;
			s /= 2;
		} else {
			d *= 3;
			s /= 3;
		}

		for (int o = 0; o < s; o++) {
			if (d == 2) {
				if (ANTRAG_LOG3_D) { // m = 6
					fpr a = f[o], b = f[o+s];
					fpc w = MTH_ROOTS[1*s];
					fpc_write(f,o, (fpc) {
						fpr_add(a, fpr_mul(w.re, b)),
						fpr_mul(w.im, b)
					});
				} else { // m = 4
					// evaluate in z=i -> no-op in this representation
				}

			} else if(d % 3 == 0) {
				fpc j = MTH_ROOTS[ANTRAG_M/3];
				for (int i = 0; i*6 < d; i++) {
					fpc w = MTH_ROOTS[ROOT_IDX[i*s*3]*s];
					fpc w2 = MTH_ROOTS[ROOT_IDX[i*s*3]*s * 2 % ANTRAG_M];
					int k = o + 3*s*i;
					fpc a = fpc_read(f, k);
					fpc b = fpc_mul(fpc_read(f, k+s), w);
					fpc c = fpc_mul(fpc_read(f, k+2*s), w2);
					fpc bj = fpc_mul(b,j), cj = fpc_mul(c,j);
					fpc bj2 = fpc_mul(bj,j), cj2 = fpc_mul(cj,j);
					fpc_write(f,k,     fpc_add(a, fpc_add(b,   c)));
					fpc_write(f,k+s,   fpc_add(a, fpc_add(bj,  cj2)));
					fpc_write(f,k+2*s, fpc_add(a, fpc_add(bj2, cj)));
				}

			} else if (d % 2 == 0) {
				for (int i = 0; i*4 < d; i++) {
					fpc w = MTH_ROOTS[ROOT_IDX[i*s*2]*s];
					int k = o + 2*s*i;
					fpc a = fpc_read(f, k);
					fpc b = fpc_mul(fpc_read(f, k+s), w);
					fpc_write(f,k,     fpc_add(a, b));
					fpc_write(f,k+s,   fpc_sub(a, b));
				}
			}
		}
	}
}

/* see inner.h */
TARGET_AVX2
void
Zf(iFFT)(fpr f[ANTRAG_D])
{
	int d = ANTRAG_D;
	int s = 1; // step size = ANTRAG_D/d = ANTRAG_M/m
	while (d > 1) {
		for (int o = 0; o < s; o++) {
			if (d == 2) {
				if (ANTRAG_LOG3_D) { // m = 6
					fpc z = fpc_read(f, o);
					fpc w = MTH_ROOTS[1*s];
					f[o+s] = fpr_div(z.im, w.im);
					f[o] = fpr_sub(z.re, fpr_mul(f[o+s], w.re));
				} else { // m = 4
					// no-op
				}

			} else if(d % 3 == 0) {
				fpc j = MTH_ROOTS[ANTRAG_M/3];
				for (int i = 0; i*6 < d; i++) {
					fpc w = MTH_ROOTS[ROOT_IDX[i*s*3]*s];
					fpc w2 = MTH_ROOTS[ROOT_IDX[i*s*3]*s * 2 % ANTRAG_M];
					int k = o + 3*s*i;
					fpc a = fpc_read(f, k);
					fpc b = fpc_read(f, k+s);
					fpc c = fpc_read(f, k+2*s);
					fpc bj = fpc_mul(b,j), cj = fpc_mul(c,j);
					fpc bj2 = fpc_mul(bj,j), cj2 = fpc_mul(cj,j);
					fpc_write(f,k,     fpc_add(a, fpc_add(b,   c)));
					fpc_write(f,k+s,   fpc_div(fpc_add(a, fpc_add(bj2, cj)),  w));
					fpc_write(f,k+2*s, fpc_div(fpc_add(a, fpc_add(bj,  cj2)), w2));
				}

			} else if (d % 2 == 0) {
				for (int i = 0; i*4 < d; i++) {
					fpc w = MTH_ROOTS[ROOT_IDX[i*s*2]*s];
					int k = o + 2*s*i;
					fpc a = fpc_read(f, k);
					fpc b = fpc_read(f, k+s);
					fpc_write(f,k,   fpc_add(a, b));
					fpc_write(f,k+s, fpc_div(fpc_sub(a, b), w));
				}
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

	for (int i = 0; i < ANTRAG_D; i++) {
		f[i] = fpr_div(f[i], FPR(FFT_SIZE));
	}
}

/* see inner.h */
TARGET_AVX2
void
Zf(poly_add)(
	fpr a[restrict ANTRAG_D], const fpr b[restrict ANTRAG_D])
{
	size_t u;
#if FALCON_AVX2 // yyyAVX2+1
	for (u = 0; u < ANTRAG_D; u += 4) {
		_mm256_storeu_pd(&a[u].v,
			_mm256_add_pd(
				_mm256_loadu_pd(&a[u].v),
				_mm256_loadu_pd(&b[u].v)));
	}
#else // yyyAVX2+0
	for (u = 0; u < ANTRAG_D; u ++) {
		a[u] = fpr_add(a[u], b[u]);
	}
#endif // yyyAVX2-
}

/* see inner.h */
TARGET_AVX2
void
Zf(poly_sub)(
	fpr a[restrict ANTRAG_D], const fpr b[restrict ANTRAG_D])
{
	size_t u;

#if FALCON_AVX2 // yyyAVX2+1
	for (u = 0; u < ANTRAG_D; u += 4) {
		_mm256_storeu_pd(&a[u].v,
			_mm256_sub_pd(
				_mm256_loadu_pd(&a[u].v),
				_mm256_loadu_pd(&b[u].v)));
	}
#else // yyyAVX2+0
	for (u = 0; u < ANTRAG_D; u ++) {
		a[u] = fpr_sub(a[u], b[u]);
	}
#endif // yyyAVX2-
}

/* see inner.h */
TARGET_AVX2
void
Zf(poly_adj_fft)(fpr a[ANTRAG_D])
{
	size_t u;
#if FALCON_AVX2 // yyyAVX2+1
	__m256d s;
	s = _mm256_set1_pd(-0.0);
	for (u = FFT_SIZE; u < ANTRAG_D; u += 4) {
		_mm256_storeu_pd(&a[u].v,
			_mm256_xor_pd(_mm256_loadu_pd(&a[u].v), s));
	}
#else // yyyAVX2+0
	for (u = FFT_SIZE; u < ANTRAG_D; u++) {
		a[u] = fpr_neg(a[u]);
	}
#endif // yyyAVX2-
}

/* see inner.h */
TARGET_AVX2
void
Zf(poly_mul_fft)(
	fpr a[restrict ANTRAG_D], const fpr b[restrict ANTRAG_D])
{
	size_t hn, u;
	hn = FFT_SIZE;
#if FALCON_AVX2 // yyyAVX2+1
	for (u = 0; u < hn; u += 4) {
		__m256d a_re, a_im, b_re, b_im, c_re, c_im;

		a_re = _mm256_loadu_pd(&a[u].v);
		a_im = _mm256_loadu_pd(&a[u + hn].v);
		b_re = _mm256_loadu_pd(&b[u].v);
		b_im = _mm256_loadu_pd(&b[u + hn].v);
		c_re = FMSUB(
			a_re, b_re, _mm256_mul_pd(a_im, b_im));
		c_im = FMADD(
			a_re, b_im, _mm256_mul_pd(a_im, b_re));
		_mm256_storeu_pd(&a[u].v, c_re);
		_mm256_storeu_pd(&a[u + hn].v, c_im);
	}
#else // yyyAVX2+0
	for (u = 0; u < hn; u++) {
		fpr a_re, a_im, b_re, b_im;

		a_re = a[u];
		a_im = a[u + hn];
		b_re = b[u];
		b_im = b[u + hn];
		FPC_MUL(a[u], a[u + hn], a_re, a_im, b_re, b_im);
	}
#endif // yyyAVX2-
}

/* see inner.h */
TARGET_AVX2
void
Zf(poly_muladj_fft)(
	fpr a[restrict ANTRAG_D], const fpr b[restrict ANTRAG_D])
{
	size_t hn, u;
	hn = FFT_SIZE;
#if FALCON_AVX2 // yyyAVX2+1
	for (u = 0; u < hn; u += 4) {
		__m256d a_re, a_im, b_re, b_im, c_re, c_im;

		a_re = _mm256_loadu_pd(&a[u].v);
		a_im = _mm256_loadu_pd(&a[u + hn].v);
		b_re = _mm256_loadu_pd(&b[u].v);
		b_im = _mm256_loadu_pd(&b[u + hn].v);
		c_re = FMADD(
			a_re, b_re, _mm256_mul_pd(a_im, b_im));
		c_im = FMSUB(
			a_im, b_re, _mm256_mul_pd(a_re, b_im));
		_mm256_storeu_pd(&a[u].v, c_re);
		_mm256_storeu_pd(&a[u + hn].v, c_im);
	}
#else // yyyAVX2+0
	for (u = 0; u < hn; u ++) {
		fpr a_re, a_im, b_re, b_im;

		a_re = a[u];
		a_im = a[u + hn];
		b_re = b[u];
		b_im = fpr_neg(b[u + hn]);
		FPC_MUL(a[u], a[u + hn], a_re, a_im, b_re, b_im);
	}
#endif // yyyAVX2-
}

/* see inner.h */
TARGET_AVX2
void
Zf(poly_mulselfadj_fft)(fpr a[ANTRAG_D])
{
	/*
	 * Since each coefficient is multiplied with its own conjugate,
	 * the result contains only real values.
	 */
	size_t hn, u;
	hn = FFT_SIZE;
#if FALCON_AVX2 // yyyAVX2+1
	__m256d zero;
	zero = _mm256_setzero_pd();
	for (u = 0; u < hn; u += 4) {
		__m256d a_re, a_im;

		a_re = _mm256_loadu_pd(&a[u].v);
		a_im = _mm256_loadu_pd(&a[u + hn].v);
		_mm256_storeu_pd(&a[u].v,
			FMADD(a_re, a_re,
				_mm256_mul_pd(a_im, a_im)));
		_mm256_storeu_pd(&a[u + hn].v, zero);
	}
#else // yyyAVX2+0
	for (u = 0; u < hn; u ++) {
		fpr a_re, a_im;

		a_re = a[u];
		a_im = a[u + hn];
		a[u] = fpr_add(fpr_sqr(a_re), fpr_sqr(a_im));
		a[u + hn] = fpr_zero;
	}
#endif // yyyAVX2-
}

/* see inner.h */
TARGET_AVX2
void
Zf(poly_div_fft)(
	fpr a[restrict ANTRAG_D], const fpr b[restrict ANTRAG_D])
{
	size_t hn, u;
	hn = FFT_SIZE;
#if FALCON_AVX2 // yyyAVX2+1
	__m256d one;
	one = _mm256_set1_pd(1.0);
	for (u = 0; u < hn; u += 4) {
		__m256d a_re, a_im, b_re, b_im, c_re, c_im, t;

		a_re = _mm256_loadu_pd(&a[u].v);
		a_im = _mm256_loadu_pd(&a[u + hn].v);
		b_re = _mm256_loadu_pd(&b[u].v);
		b_im = _mm256_loadu_pd(&b[u + hn].v);
		t = _mm256_div_pd(one,
			FMADD(b_re, b_re,
				_mm256_mul_pd(b_im, b_im)));
		b_re = _mm256_mul_pd(b_re, t);
		b_im = _mm256_mul_pd(b_im, t);
		c_re = FMADD(
			a_re, b_re, _mm256_mul_pd(a_im, b_im));
		c_im = FMSUB(
			a_im, b_re, _mm256_mul_pd(a_re, b_im));
		_mm256_storeu_pd(&a[u].v, c_re);
		_mm256_storeu_pd(&a[u + hn].v, c_im);
	}
#else // yyyAVX2+0
	for (u = 0; u < hn; u ++) {
		fpr a_re, a_im, b_re, b_im;

		a_re = a[u];
		a_im = a[u + hn];
		b_re = b[u];
		b_im = b[u + hn];
		FPC_DIV(a[u], a[u + hn], a_re, a_im, b_re, b_im);
	}
#endif // yyyAVX2-
}
