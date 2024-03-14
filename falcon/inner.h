#ifndef FALCON_INNER_H__
#define FALCON_INNER_H__

/*
 * Internal functions for Falcon. This is not the API intended to be
 * used by applications; instead, this internal API provides all the
 * primitives on which wrappers build to provide external APIs.
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

/*
 * IMPORTANT API RULES
 * -------------------
 *
 * This API has some non-trivial usage rules:
 *
 *
 *  - All public functions (i.e. the non-static ones) must be referenced
 *    with the Zf() macro (e.g. Zf(verify_raw) for the verify_raw()
 *    function). That macro adds a prefix to the name, which is
 *    configurable with the FALCON_PREFIX macro. This allows compiling
 *    the code into a specific "namespace" and potentially including
 *    several versions of this code into a single application (e.g. to
 *    have an AVX2 and a non-AVX2 variants and select the one to use at
 *    runtime based on availability of AVX2 opcodes).
 *
 *  - Functions that need temporary buffers expects them as a final
 *    tmp[] array of type uint8_t*, with a size which is documented for
 *    each function. However, most have some alignment requirements,
 *    because they will use the array to store 16-bit, 32-bit or 64-bit
 *    values (e.g. uint64_t or double). The caller must ensure proper
 *    alignment. What happens on unaligned access depends on the
 *    underlying architecture, ranging from a slight time penalty
 *    to immediate termination of the process.
 *
 *  - Some functions rely on specific rounding rules and precision for
 *    floating-point numbers. On some systems (in particular 32-bit x86
 *    with the 387 FPU), this requires setting an hardware control
 *    word. The caller MUST use set_fpu_cw() to ensure proper precision:
 *
 *      oldcw = set_fpu_cw(2);
 *      Zf(sign_dyn)(...);
 *      set_fpu_cw(oldcw);
 *
 *    On systems where the native floating-point precision is already
 *    proper, or integer-based emulation is used, the set_fpu_cw()
 *    function does nothing, so it can be called systematically.
 */

// yyyPQCLEAN+0 yyyNIST+0 yyySUPERCOP+0
#include "config.h"
// yyyPQCLEAN- yyyNIST- yyySUPERCOP-
// yyySUPERCOP+1
// yyyCONF*
// yyySUPERCOP-

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "../gen/const.h"

#if defined FALCON_AVX2 && FALCON_AVX2 // yyyAVX2+1
/*
 * This implementation uses AVX2 and optionally FMA intrinsics.
 */
#include <immintrin.h>
#ifndef FALCON_LE
#define FALCON_LE   1
#endif
#ifndef FALCON_UNALIGNED
#define FALCON_UNALIGNED   1
#endif
#if defined __GNUC__
#if defined FALCON_FMA && FALCON_FMA
#define TARGET_AVX2   __attribute__((target("avx2,fma")))
#else
#define TARGET_AVX2   __attribute__((target("avx2")))
#endif
#elif defined _MSC_VER && _MSC_VER
#pragma warning( disable : 4752 )
#endif
#if defined FALCON_FMA && FALCON_FMA
#define FMADD(a, b, c)   _mm256_fmadd_pd(a, b, c)
#define FMSUB(a, b, c)   _mm256_fmsub_pd(a, b, c)
#else
#define FMADD(a, b, c)   _mm256_add_pd(_mm256_mul_pd(a, b), c)
#define FMSUB(a, b, c)   _mm256_sub_pd(_mm256_mul_pd(a, b), c)
#endif
#endif // yyyAVX2-

// yyyNIST+0 yyyPQCLEAN+0
/*
 * On MSVC, disable warning about applying unary minus on an unsigned
 * type: this is perfectly defined standard behaviour and we do it
 * quite often.
 */
#if defined _MSC_VER && _MSC_VER
#pragma warning( disable : 4146 )
#endif

// yyySUPERCOP+0
/*
 * Enable ARM assembly on any ARMv7m platform (if it was not done before).
 */
#ifndef FALCON_ASM_CORTEXM4
#if (defined __ARM_ARCH_7EM__ && __ARM_ARCH_7EM__) \
	&& (defined __ARM_FEATURE_DSP && __ARM_FEATURE_DSP)
#define FALCON_ASM_CORTEXM4   1
#else
#define FALCON_ASM_CORTEXM4   0
#endif
#endif
// yyySUPERCOP-

#if defined __i386__ || defined _M_IX86 \
	|| defined __x86_64__ || defined _M_X64 || \
	(defined _ARCH_PWR8 && \
		(defined __LITTLE_ENDIAN || defined __LITTLE_ENDIAN__))

#ifndef FALCON_LE
#define FALCON_LE     1
#endif
#ifndef FALCON_UNALIGNED
#define FALCON_UNALIGNED   1
#endif

#elif defined FALCON_ASM_CORTEXM4 && FALCON_ASM_CORTEXM4

#ifndef FALCON_LE
#define FALCON_LE     1
#endif
#ifndef FALCON_UNALIGNED
#define FALCON_UNALIGNED   0
#endif

#elif (defined __LITTLE_ENDIAN__ && __LITTLE_ENDIAN__) \
	|| (defined __BYTE_ORDER__ && defined __ORDER_LITTLE_ENDIAN__ \
		&& __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)

#ifndef FALCON_LE
#define FALCON_LE     1
#endif
#ifndef FALCON_UNALIGNED
#define FALCON_UNALIGNED   0
#endif

#else

#ifndef FALCON_LE
#define FALCON_LE     0
#endif
#ifndef FALCON_UNALIGNED
#define FALCON_UNALIGNED   0
#endif

#endif

/*
 * We ensure that both FALCON_FPEMU and FALCON_FPNATIVE are defined,
 * with compatible values (exactly one of them must be non-zero).
 * If none is defined, then default FP implementation is 'native'
 * except on ARM Cortex M4.
 */
#if !defined FALCON_FPEMU && !defined FALCON_FPNATIVE

#if (defined __ARM_FP && ((__ARM_FP & 0x08) == 0x08)) \
	|| (!defined __ARM_FP && defined __ARM_VFPV2__)
#define FALCON_FPEMU      0
#define FALCON_FPNATIVE   1
#elif defined FALCON_ASM_CORTEXM4 && FALCON_ASM_CORTEXM4
#define FALCON_FPEMU      1
#define FALCON_FPNATIVE   0
#else
#define FALCON_FPEMU      0
#define FALCON_FPNATIVE   1
#endif

#elif defined FALCON_FPEMU && !defined FALCON_FPNATIVE

#if FALCON_FPEMU
#define FALCON_FPNATIVE   0
#else
#define FALCON_FPNATIVE   1
#endif

#elif defined FALCON_FPNATIVE && !defined FALCON_FPEMU

#if FALCON_FPNATIVE
#define FALCON_FPEMU   0
#else
#define FALCON_FPEMU   1
#endif

#endif

#if (FALCON_FPEMU && FALCON_FPNATIVE) || (!FALCON_FPEMU && !FALCON_FPNATIVE)
#error Exactly one of FALCON_FPEMU and FALCON_FPNATIVE must be selected
#endif

// yyySUPERCOP+0
/*
 * For seed generation from the operating system:
 *  - On Linux and glibc-2.25+, FreeBSD 12+ and OpenBSD, use getentropy().
 *  - On Unix-like systems, use /dev/urandom (including as a fallback
 *    for failed getentropy() calls).
 *  - On Windows, use CryptGenRandom().
 */

#ifndef FALCON_RAND_GETENTROPY
#if (defined __linux__ && defined __GLIBC__ \
	&& (__GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 25))) \
	|| (defined __FreeBSD__ && __FreeBSD__ >= 12) \
	|| defined __OpenBSD__
#define FALCON_RAND_GETENTROPY   1
#else
#define FALCON_RAND_GETENTROPY   0
#endif
#endif

#ifndef FALCON_RAND_URANDOM
#if defined _AIX \
	|| defined __ANDROID__ \
	|| defined __FreeBSD__ \
	|| defined __NetBSD__ \
	|| defined __OpenBSD__ \
	|| defined __DragonFly__ \
	|| defined __linux__ \
	|| (defined __sun && (defined __SVR4 || defined __svr4__)) \
	|| (defined __APPLE__ && defined __MACH__)
#define FALCON_RAND_URANDOM   1
#else
#define FALCON_RAND_URANDOM   0
#endif
#endif

#ifndef FALCON_RAND_WIN32
#if defined _WIN32 || defined _WIN64
#define FALCON_RAND_WIN32   1
#else
#define FALCON_RAND_WIN32   0
#endif
#endif
// yyySUPERCOP-

/*
 * For still undefined compile-time macros, define them to 0 to avoid
 * warnings with -Wundef.
 */
#ifndef FALCON_AVX2
#define FALCON_AVX2   0
#endif
#ifndef FALCON_FMA
#define FALCON_FMA   0
#endif
#ifndef FALCON_KG_CHACHA20
#define FALCON_KG_CHACHA20   0
#endif
// yyyNIST- yyyPQCLEAN-

// yyyPQCLEAN+0 yyySUPERCOP+0
/*
 * "Naming" macro used to apply a consistent prefix over all global
 * symbols.
 */
#ifndef FALCON_PREFIX
#define FALCON_PREFIX   falcon_inner
#endif
#define Zf(name)             Zf_(FALCON_PREFIX, name)
#define Zf_(prefix, name)    Zf__(prefix, name)
#define Zf__(prefix, name)   prefix ## _ ## name  
// yyyPQCLEAN- yyySUPERCOP-

// yyyAVX2+1
/*
 * We use the TARGET_AVX2 macro to tag some functions which, in some
 * configurations, may use AVX2 and FMA intrinsics; this depends on
 * the compiler. In all other cases, we just define it to emptiness
 * (i.e. it will have no effect).
 */
#ifndef TARGET_AVX2
#define TARGET_AVX2
#endif
// yyyAVX2-

/*
 * Some computations with floating-point elements, in particular
 * rounding to the nearest integer, rely on operations using _exactly_
 * the precision of IEEE-754 binary64 type (i.e. 52 bits). On 32-bit
 * x86, the 387 FPU may be used (depending on the target OS) and, in
 * that case, may use more precision bits (i.e. 64 bits, for an 80-bit
 * total type length); to prevent miscomputations, we define an explicit
 * function that modifies the precision in the FPU control word.
 *
 * set_fpu_cw() sets the precision to the provided value, and returns
 * the previously set precision; callers are supposed to restore the
 * previous precision on exit. The correct (52-bit) precision is
 * configured with the value "2". On unsupported compilers, or on
 * targets other than 32-bit x86, or when the native 'double' type is
 * not used, the set_fpu_cw() function does nothing at all.
 */
#if FALCON_FPNATIVE  // yyyFPNATIVE+1
#if defined __GNUC__ && defined __i386__
static inline unsigned
set_fpu_cw(unsigned x)
{
	unsigned short t;
	unsigned old;

	__asm__ __volatile__ ("fstcw %0" : "=m" (t) : : );
	old = (t & 0x0300u) >> 8;
	t = (unsigned short)((t & ~0x0300u) | (x << 8));
	__asm__ __volatile__ ("fldcw %0" : : "m" (t) : );
	return old;
}
#elif defined _M_IX86
static inline unsigned
set_fpu_cw(unsigned x)
{
	unsigned short t;
	unsigned old;

	__asm { fstcw t }
	old = (t & 0x0300u) >> 8;
	t = (unsigned short)((t & ~0x0300u) | (x << 8));
	__asm { fldcw t }
	return old;
}
#else
static inline unsigned
set_fpu_cw(unsigned x)
{
	return x;
}
#endif
#else  // yyyFPNATIVE+0
static inline unsigned
set_fpu_cw(unsigned x)
{
	return x;
}
#endif  // yyyFPNATIVE-

#if FALCON_FPNATIVE && !FALCON_AVX2  // yyyFPNATIVE+1 yyyAVX2+0
/*
 * If using the native 'double' type but not AVX2 code, on an x86
 * machine with SSE2 activated for maths, then we will use the
 * SSE2 intrinsics.
 */
#if defined __GNUC__ && defined __SSE2_MATH__
#include <immintrin.h>
#endif
#endif  // yyyFPNATIVE- yyyAVX2-

#if FALCON_FPNATIVE  // yyyFPNATIVE+1
/*
 * For optimal reproducibility of values, we need to disable contraction
 * of floating-point expressions; otherwise, on some architectures (e.g.
 * PowerPC), the compiler may generate fused-multiply-add opcodes that
 * may round differently than two successive separate opcodes. C99 defines
 * a standard pragma for that, but GCC-6.2.2 appears to ignore it,
 * hence the GCC-specific pragma (that Clang does not support).
 */
#if defined __clang__
#pragma STDC FP_CONTRACT OFF
#elif defined __GNUC__
#pragma GCC optimize ("fp-contract=off")
#endif
#endif  // yyyFPNATIVE-

// yyyPQCLEAN+0
/*
 * MSVC 2015 does not know the C99 keyword 'restrict'.
 */
#if defined _MSC_VER && _MSC_VER
#ifndef restrict
#define restrict   __restrict
#endif
#endif
// yyyPQCLEAN-

/* ==================================================================== */
/*
 * SHAKE256 implementation (shake.c).
 *
 * API is defined to be easily replaced with the fips202.h API defined
 * as part of PQClean.
 */

// yyyPQCLEAN+0
typedef struct {
	union {
		uint64_t A[25];
		uint8_t dbuf[200];
	} st;
	uint64_t dptr;
} inner_shake256_context;

#define inner_shake256_init      Zf(i_shake256_init)
#define inner_shake256_inject    Zf(i_shake256_inject)
#define inner_shake256_flip      Zf(i_shake256_flip)
#define inner_shake256_extract   Zf(i_shake256_extract)

void Zf(i_shake256_init)(
	inner_shake256_context *sc);
void Zf(i_shake256_inject)(
	inner_shake256_context *sc, const uint8_t *in, size_t len);
void Zf(i_shake256_flip)(
	inner_shake256_context *sc);
void Zf(i_shake256_extract)(
	inner_shake256_context *sc, uint8_t *out, size_t len);

/*
// yyyPQCLEAN+1

#include "fips202.h"

#define inner_shake256_context                shake256incctx
#define inner_shake256_init(sc)               shake256_inc_init(sc)
#define inner_shake256_inject(sc, in, len)    shake256_inc_absorb(sc, in, len)
#define inner_shake256_flip(sc)               shake256_inc_finalize(sc)
#define inner_shake256_extract(sc, out, len)  shake256_inc_squeeze(out, len, sc)

// yyyPQCLEAN+0
 */
// yyyPQCLEAN-

/* ==================================================================== */
/*
 * NTT-based functions for key generation (ntt.c).
 */

void Zf(NTT)(uint16_t v[ANTRAG_D]);
void Zf(iNTT)(uint16_t v[ANTRAG_D]);

void Zf(u16_poly_add)(uint16_t f[restrict ANTRAG_D], const uint16_t g[restrict ANTRAG_D]);
void Zf(u16_poly_sub)(uint16_t f[restrict ANTRAG_D], const uint16_t g[restrict ANTRAG_D]);
void Zf(u16_poly_neg)(uint16_t f[restrict ANTRAG_D], const uint16_t g[restrict ANTRAG_D]);
void Zf(u16_poly_mul)(uint16_t f[restrict ANTRAG_D], const uint16_t g_ntt[restrict ANTRAG_D]);

/*
 * Compute the public key h[], given the private key elements f[] and
 * g[]. This computes h = g/f mod phi mod q, where phi is the polynomial
 * modulus. This function returns 1 on success, 0 on error (an error is
 * reported if f is not invertible mod phi mod q).
 *
 * The tmp[] array must have room for at least 2*ANTRAG_D elements.
 * tmp[] must have 16-bit alignment.
 */
int Zf(compute_public)(uint16_t *h,
	const int8_t *f, const int8_t *g, uint8_t *tmp);

/* ==================================================================== */
/*
 * Implementation of floating-point real numbers (fpr.h, fpr.c).
 */

/*
 * Real numbers are implemented by an extra header file, included below.
 * This is meant to support pluggable implementations. The default
 * implementation relies on the C type 'double'.
 *
 * The included file must define the following types, functions and
 * constants:
 *
 *   fpr
 *         type for a real number
 *
 *   fpr fpr_of(int64_t i)
 *         cast an integer into a real number; source must be in the
 *         -(2^63-1)..+(2^63-1) range
 *
 *   fpr fpr_scaled(int64_t i, int sc)
 *         compute i*2^sc as a real number; source 'i' must be in the
 *         -(2^63-1)..+(2^63-1) range
 *
 *   fpr fpr_ldexp(fpr x, int e)
 *         compute x*2^e
 *
 *   int64_t fpr_rint(fpr x)
 *         round x to the nearest integer; x must be in the -(2^63-1)
 *         to +(2^63-1) range
 *
 *   int64_t fpr_trunc(fpr x)
 *         round to an integer; this rounds towards zero; value must
 *         be in the -(2^63-1) to +(2^63-1) range
 *
 *   fpr fpr_add(fpr x, fpr y)
 *         compute x + y
 *
 *   fpr fpr_sub(fpr x, fpr y)
 *         compute x - y
 *
 *   fpr fpr_neg(fpr x)
 *         compute -x
 *
 *   fpr fpr_half(fpr x)
 *         compute x/2
 *
 *   fpr fpr_double(fpr x)
 *         compute x*2
 *
 *   fpr fpr_mul(fpr x, fpr y)
 *         compute x * y
 *
 *   fpr fpr_sqr(fpr x)
 *         compute x * x
 *
 *   fpr fpr_inv(fpr x)
 *         compute 1/x
 *
 *   fpr fpr_div(fpr x, fpr y)
 *         compute x/y
 *
 *   fpr fpr_sqrt(fpr x)
 *         compute the square root of x
 *
 *   int fpr_lt(fpr x, fpr y)
 *         return 1 if x < y, 0 otherwise
 *
 *   uint64_t fpr_expm_p63(fpr x)
 *         return exp(x), assuming that 0 <= x < log(2). Returned value
 *         is scaled to 63 bits (i.e. it really returns 2^63*exp(-x),
 *         rounded to the nearest integer). Computation should have a
 *         precision of at least 45 bits.
 *
 * Constants of type 'fpr':
 *
 *   fpr fpr_q                 12289
 *   fpr fpr_inverse_of_q      1/12289
 *   fpr fpr_inv_2sqrsigma0    1/(2*(1.8205^2))
 *   fpr fpr_inv_sigma[]       1/sigma (indexed by logn, 1 to 10)
 *   fpr fpr_sigma_min[]       1/sigma_min (indexed by logn, 1 to 10)
 *   fpr fpr_log2              log(2)
 *   fpr fpr_inv_log2          1/log(2)
 *   fpr fpr_bnorm_max         16822.4121
 *   fpr fpr_zero              0
 *   fpr fpr_one               1
 *   fpr fpr_two               2
 *   fpr fpr_onehalf           0.5
 *   fpr fpr_ptwo31            2^31
 *   fpr fpr_ptwo31m1          2^31-1
 *   fpr fpr_mtwo31m1          -(2^31-1)
 *   fpr fpr_ptwo63m1          2^63-1
 *   fpr fpr_mtwo63m1          -(2^63-1)
 *   fpr fpr_ptwo63            2^63
 */
#include "fpr.h"

/* ==================================================================== */
/*
 * RNG (rng.c).
 *
 * A PRNG based on ChaCha20 is implemented; it is seeded from a SHAKE256
 * context (flipped) and is used for bulk pseudorandom generation.
 * A system-dependent seed generator is also provided.
 */

/*
 * Obtain a random seed from the system RNG.
 *
 * Returned value is 1 on success, 0 on error.
 */
int Zf(get_seed)(void *seed, size_t seed_len);

/*
 * Structure for a PRNG. This includes a large buffer so that values
 * get generated in advance. The 'state' is used to keep the current
 * PRNG algorithm state (contents depend on the selected algorithm).
 *
 * The unions with 'dummy_u64' are there to ensure proper alignment for
 * 64-bit direct access.
 */
typedef struct {
	union {
		uint8_t d[512]; /* MUST be 512, exactly */
		uint64_t dummy_u64;
	} buf;
	size_t ptr;
	union {
		uint8_t d[256];
		uint64_t dummy_u64;
	} state;
	int type;
} prng;

/*
 * Instantiate a PRNG. That PRNG will feed over the provided SHAKE256
 * context (in "flipped" state) to obtain its initial state.
 */
void Zf(prng_init)(prng *p, inner_shake256_context *src);

/*
 * Refill the PRNG buffer. This is normally invoked automatically, and
 * is declared here only so that prng_get_u64() may be inlined.
 */
void Zf(prng_refill)(prng *p);

/*
 * Get some bytes from a PRNG.
 */
void Zf(prng_get_bytes)(prng *p, void *dst, size_t len);

/*
 * Get a 64-bit random value from a PRNG.
 */
static inline uint64_t
prng_get_u64(prng *p)
{
	size_t u;

	/*
	 * If there are less than 9 bytes in the buffer, we refill it.
	 * This means that we may drop the last few bytes, but this allows
	 * for faster extraction code. Also, it means that we never leave
	 * an empty buffer.
	 */
	u = p->ptr;
	if (u >= (sizeof p->buf.d) - 9) {
		Zf(prng_refill)(p);
		u = 0;
	}
	p->ptr = u + 8;

	/*
	 * On systems that use little-endian encoding and allow
	 * unaligned accesses, we can simply read the data where it is.
	 */
#if FALCON_LE && FALCON_UNALIGNED  // yyyLEU+1
	return *(uint64_t *)(p->buf.d + u);
#else  // yyyLEU+0
	return (uint64_t)p->buf.d[u + 0]
		| ((uint64_t)p->buf.d[u + 1] << 8)
		| ((uint64_t)p->buf.d[u + 2] << 16)
		| ((uint64_t)p->buf.d[u + 3] << 24)
		| ((uint64_t)p->buf.d[u + 4] << 32)
		| ((uint64_t)p->buf.d[u + 5] << 40)
		| ((uint64_t)p->buf.d[u + 6] << 48)
		| ((uint64_t)p->buf.d[u + 7] << 56);
#endif  // yyyLEU-
}

/*
 * Get an 8-bit random value from a PRNG.
 */
static inline unsigned
prng_get_u8(prng *p)
{
	unsigned v;

	v = p->buf.d[p->ptr ++];
	if (p->ptr == sizeof p->buf.d) {
		Zf(prng_refill)(p);
	}
	return v;
}

/* ==================================================================== */
/*
 * FFT (falcon-fft.c).
 *
 * A real polynomial is represented as an array of N 'fpr' elements.
 * The FFT representation of a real polynomial contains N/2 complex
 * elements; each is stored as two real numbers, for the real and
 * imaginary parts, respectively. See falcon-fft.c for details on the
 * internal representation.
 */

/*
 * Compute FFT in-place: the source array should contain a real
 * polynomial (N coefficients); its storage area is reused to store
 * the FFT representation of that polynomial (N/2 complex numbers).
 *
 * 'logn' MUST lie between 1 and 10 (inclusive).
 */
void Zf(FFT)(fpr f[ANTRAG_D]);

/*
 * Compute the inverse FFT in-place: the source array should contain the
 * FFT representation of a real polynomial (N/2 elements); the resulting
 * real polynomial (N coefficients of type 'fpr') is written over the
 * array.
 *
 * 'logn' MUST lie between 1 and 10 (inclusive).
 */
void Zf(iFFT)(fpr f[ANTRAG_D]);

/*
 * Add polynomial b to polynomial a. a and b MUST NOT overlap. This
 * function works in both normal and FFT representations.
 */
void Zf(poly_add)(fpr a[restrict ANTRAG_D], const fpr b[restrict ANTRAG_D]);

/*
 * Subtract polynomial b from polynomial a. a and b MUST NOT overlap. This
 * function works in both normal and FFT representations.
 */
void Zf(poly_sub)(fpr a[restrict ANTRAG_D], const fpr b[restrict ANTRAG_D]);

/*
 * Compute adjoint of polynomial a. This function works only in FFT
 * representation.
 */
void Zf(poly_adj_fft)(fpr a[ANTRAG_D]);

/*
 * Multiply polynomial a with polynomial b. a and b MUST NOT overlap.
 * This function works only in FFT representation.
 */
void Zf(poly_mul_fft)(fpr a[restrict ANTRAG_D], const fpr b[restrict ANTRAG_D]);

/*
 * Multiply polynomial a with the adjoint of polynomial b. a and b MUST NOT
 * overlap. This function works only in FFT representation.
 */
void Zf(poly_muladj_fft)(fpr a[restrict ANTRAG_D], const fpr b[restrict ANTRAG_D]);

/*
 * Multiply polynomial with its own adjoint. This function works only in FFT
 * representation.
 */
void Zf(poly_mulselfadj_fft)(fpr a[ANTRAG_D]);

/*
 * Divide polynomial a by polynomial b, modulo X^N+1 (FFT representation).
 * a and b MUST NOT overlap.
 */
void Zf(poly_div_fft)(fpr a[restrict ANTRAG_D], const fpr b[restrict ANTRAG_D]);

/* ==================================================================== */
/*
 * Key pair generation.
 */

/*
 * Required sizes of the temporary buffer (in bytes).
 *
 * This size is 28*2^logn bytes, except for degrees 2 and 4 (logn = 1
 * or 2) where it is slightly greater.
 */
#define FALCON_KEYGEN_TEMP_1      136
#define FALCON_KEYGEN_TEMP_2      272
#define FALCON_KEYGEN_TEMP_3      224
#define FALCON_KEYGEN_TEMP_4      448
#define FALCON_KEYGEN_TEMP_5      896
#define FALCON_KEYGEN_TEMP_6     1792
#define FALCON_KEYGEN_TEMP_7     3584
#define FALCON_KEYGEN_TEMP_8     7168
#define FALCON_KEYGEN_TEMP_9    14336
#define FALCON_KEYGEN_TEMP_10   28672

/*
 * Generate a new key pair. Randomness is extracted from the provided
 * SHAKE256 context, which must have already been seeded and flipped.
 * The tmp[] array must have suitable size (see FALCON_KEYGEN_TEMP_*
 * macros) and be aligned for the uint32_t, uint64_t and fpr types.
 *
 * The private key elements are written in f, g, F and G, and the
 * public key is written in h. Either or both of G and h may be NULL,
 * in which case the corresponding element is not returned (they can
 * be recomputed from f, g and F).
 *
 * tmp[] must have 64-bit alignment.
 * This function uses floating-point rounding (see set_fpu_cw()).
 */
int Zf(solve_NTRU)(unsigned logn, int8_t *F, int8_t *G,
	const int8_t *f, const int8_t *g, int lim, uint32_t *tmp);

void Zf(keygen)(inner_shake256_context *rng,
	int8_t *f, int8_t *g, int8_t *F, int8_t *G, uint16_t *h,
	unsigned logn, uint8_t *tmp);

/* ==================================================================== */
/*
 * Signature generation.
 */

/*
 * Expand a private key into the B0 matrix in FFT representation and
 * the LDL tree. All the values are written in 'expanded_key', for
 * a total of (8*logn+40)*2^logn bytes.
 *
 * The tmp[] array must have room for at least 48*2^logn bytes.
 *
 * tmp[] must have 64-bit alignment.
 * This function uses floating-point rounding (see set_fpu_cw()).
 */
void Zf(expand_privkey)(fpr *restrict expanded_key,
	const int8_t *f, const int8_t *g, const int8_t *F, const int8_t *G,
	unsigned logn, uint8_t *restrict tmp);

/*
 * Compute a signature over the provided hashed message (hm); the
 * signature value is one short vector. This function uses an
 * expanded key (as generated by Zf(expand_privkey)()).
 *
 * The sig[] and hm[] buffers may overlap.
 *
 * On successful output, the start of the tmp[] buffer contains the s1
 * vector (as int16_t elements).
 *
 * The minimal size (in bytes) of tmp[] is 48*2^logn bytes.
 *
 * tmp[] must have 64-bit alignment.
 * This function uses floating-point rounding (see set_fpu_cw()).
 */
void Zf(sign_tree)(int16_t *sig, inner_shake256_context *rng,
	const fpr *restrict expanded_key,
	const uint16_t *hm, unsigned logn, uint8_t *tmp);

/*
 * Compute a signature over the provided hashed message (hm); the
 * signature value is one short vector. This function uses a raw
 * key and dynamically recompute the B0 matrix and LDL tree; this
 * saves RAM since there is no needed for an expanded key, but
 * increases the signature cost.
 *
 * The sig[] and hm[] buffers may overlap.
 *
 * On successful output, the start of the tmp[] buffer contains the s1
 * vector (as int16_t elements).
 *
 * The minimal size (in bytes) of tmp[] is 72*2^logn bytes.
 *
 * tmp[] must have 64-bit alignment.
 * This function uses floating-point rounding (see set_fpu_cw()).
 */
void Zf(sign_dyn)(int16_t *sig, inner_shake256_context *rng,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint16_t *hm, unsigned logn, uint8_t *tmp);

/*
 * Internal sampler engine. Exported for tests.
 *
 * sampler_context wraps around a source of random numbers (PRNG) and
 * the sigma_min value (nominally dependent on the degree).
 *
 * sampler() takes as parameters:
 *   ctx      pointer to the sampler_context structure
 *   mu       center for the distribution
 *   isigma   inverse of the distribution standard deviation
 * It returns an integer sampled along the Gaussian distribution centered
 * on mu and of standard deviation sigma = 1/isigma.
 *
 * gaussian0_sampler() takes as parameter a pointer to a PRNG, and
 * returns an integer sampled along a half-Gaussian with standard
 * deviation sigma0 = 1.8205 (center is 0, returned value is
 * nonnegative).
 */

typedef struct {
	prng p;
	fpr sigma_min;
} sampler_context;

TARGET_AVX2
int Zf(sampler)(void *ctx, fpr mu, fpr isigma);

TARGET_AVX2
int Zf(gaussian0_sampler)(prng *p);

/* ==================================================================== */


#endif
