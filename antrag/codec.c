#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "api.h"
#include "../gen/const.h"
#include "poly.h"

inline static size_t encode_i8_vector(size_t off, uint8_t* buf, const int8_t* v) {
	for(int i = 0; i < ANTRAG_D; i++) {
		assert(off < ANTRAG_SK_SIZE);
		buf[off++] = (uint8_t) v[i]; // conversion is 2's complement
	}
	return off;
}

uint64_t encode_double(double x) {
	assert(x == x); // check NaN
	uint64_t s, e, m;
	if(x < 0) {
		s = 1;
		x = -x;
	} else {
		s = 0;
	}
	int e2;
	double f = frexp(x, &e2);
	f *= 2;
	e2--;
	if(isinf(f) || e2 > 1023) {
		m = 0;
		e2 = 1024;
	} else if(f == 0 || e2 < -1022) {
		m = 0;
		e2 = -1023;
	} else {
		uint64_t b52 = (uint64_t) 1 << 52;
		m = round(f * b52);
		m -= b52;
	}
	e = e2 + 1023;
	return (s << 63) | (e << 52) | m;
}

inline static size_t encode_poly(size_t off, uint8_t* buf, const poly* p) {
	for(int i = 0; i < ANTRAG_D; i++) {
		double x = p->coeffs[i].v;
		uint64_t y;
		#if SIMPLE_DOUBLE_CODEC
		memcpy(&y, &x, 8);
		#else
		y = encode_double(x);
		#endif
		assert(off + 8 <= ANTRAG_SK_SIZE);
		// likely optimized to memcpy on little-endian platforms
		for(int j = 0; j < 8; j++) {
			buf[off++] = y & 0xff;
			y >>= 8;
		}
	}
	return off;
}

void encode_sk(uint8_t* buf, const secret_key* sk) {
	size_t off = 0;
	off = encode_i8_vector(off, buf, sk->f);
	off = encode_i8_vector(off, buf, sk->g);
	off = encode_i8_vector(off, buf, sk->F);
	off = encode_i8_vector(off, buf, sk->G);
	off = encode_poly(off, buf, &sk->b10);
	off = encode_poly(off, buf, &sk->b11);
	off = encode_poly(off, buf, &sk->b20);
	off = encode_poly(off, buf, &sk->b21);
	off = encode_poly(off, buf, &sk->GSO_b10);
	off = encode_poly(off, buf, &sk->GSO_b11);
	off = encode_poly(off, buf, &sk->GSO_b20);
	off = encode_poly(off, buf, &sk->GSO_b21);
	off = encode_poly(off, buf, &sk->beta10);
	off = encode_poly(off, buf, &sk->beta11);
	off = encode_poly(off, buf, &sk->beta20);
	off = encode_poly(off, buf, &sk->beta21);
	off = encode_poly(off, buf, &sk->sigma1);
	off = encode_poly(off, buf, &sk->sigma2);
}

inline static size_t decode_i8_vector(size_t off, const uint8_t* buf, int8_t* v) {
	assert(off + ANTRAG_D <= ANTRAG_SK_SIZE);
	for(int i = 0; i < ANTRAG_D; i++) {
		uint8_t b = buf[off++];
		v[i] = b < 0x80 ? b : b - 0x100; // b is promoted to int -> no overflow
	}
	return off;
}

double decode_double(uint64_t y) {
	int s = y >> 63;
	uint32_t e = (y >> 52) & 2047;
	uint64_t b52 = (uint64_t) 1 << 52;
	uint64_t m = y & (b52 - 1);
	int e2 = (int) e - 1023;
	assert(e2 != 1024 || m == 0); // not NaN
	if(e2 == -1023) {
		e2 = -1022;
	} else {
		m += b52;
	}
	double x = ldexp((double) m, e2 - 52);
	assert(isfinite(x)); // check NaN
	if(s) x = -x;
	return x;
}

inline static size_t decode_poly(size_t off, const uint8_t* buf, poly* p) {
	for(int i = 0; i < ANTRAG_D; i++) {
		uint64_t y = 0;
		assert(off + 8 <= ANTRAG_SK_SIZE);
		for(int j = 7; j >= 0; j--) {
			y <<= 8;
			y |= buf[off + j];
		}
		off += 8;
		
		double x;
		#if SIMPLE_DOUBLE_CODEC
		memcpy(&x, &y, 8);
		#else
		x = decode_double(y);
		#endif
		
		p->coeffs[i].v = x;
	}
	return off;
}

void decode_sk(secret_key* sk, const uint8_t* buf) {
	size_t off = 0;
	off = decode_i8_vector(off, buf, sk->f);
	off = decode_i8_vector(off, buf, sk->g);
	off = decode_i8_vector(off, buf, sk->F);
	off = decode_i8_vector(off, buf, sk->G);
	off = decode_poly(off, buf, &sk->b10);
	off = decode_poly(off, buf, &sk->b11);
	off = decode_poly(off, buf, &sk->b20);
	off = decode_poly(off, buf, &sk->b21);
	off = decode_poly(off, buf, &sk->GSO_b10);
	off = decode_poly(off, buf, &sk->GSO_b11);
	off = decode_poly(off, buf, &sk->GSO_b20);
	off = decode_poly(off, buf, &sk->GSO_b21);
	off = decode_poly(off, buf, &sk->beta10);
	off = decode_poly(off, buf, &sk->beta11);
	off = decode_poly(off, buf, &sk->beta20);
	off = decode_poly(off, buf, &sk->beta21);
	off = decode_poly(off, buf, &sk->sigma1);
	off = decode_poly(off, buf, &sk->sigma2);
}

typedef struct {
	uint32_t bit_buf;
	int buf_len;
} BitStream;
BitStream empty_stream = { 0, 0 };

static inline void queue_bits(BitStream* bst, uint16_t n, int bit_cnt) {
	assert(bst->buf_len + bit_cnt <= 32);
	uint32_t mask = (1 << bit_cnt) - 1;
	bst->bit_buf |= (n & mask) << bst->buf_len;
	bst->buf_len += bit_cnt;
}
static inline uint16_t dequeue_bits(BitStream* bst, int bit_cnt) {
	assert(bst->buf_len >= bit_cnt);
	uint32_t mask = (1 << bit_cnt) - 1;
	uint16_t n = bst->bit_buf & mask;
	bst->bit_buf >>= bit_cnt;
	bst->buf_len -= bit_cnt;
	return n;
}

static inline void queue_byte(BitStream* bst, const uint8_t* buf, size_t* off) {
	queue_bits(bst, buf[(*off)++], 8);
}
static inline bool has_queued_byte(BitStream* bst, bool pad) {
	if(pad && bst->buf_len < 8 && bst->buf_len > 0) {
		bst->buf_len = 8;
	}
	return bst->buf_len >= 8;
}
static inline void dequeue_byte(BitStream* bst, uint8_t* buf, size_t* off) {
	buf[(*off)++] = dequeue_bits(bst, 8);
}

void encode_pk(uint8_t* buf, const public_key* pk) {
	BitStream bst = empty_stream;
	size_t off = 0;
	for(int i = 0; i < ANTRAG_D; i++) {
		uint16_t n = pk->h_ntt.coeffs[i];
		assert(n < ANTRAG_Q);
		queue_bits(&bst, n, ANTRAG_Q_BITS);
		while(has_queued_byte(&bst, i == ANTRAG_D-1)) {
			assert(off < ANTRAG_PK_SIZE);
			dequeue_byte(&bst, buf, &off);
		}
	}
}

void decode_pk(public_key* pk, const uint8_t* buf) {
	BitStream bst = empty_stream;
	size_t off = 0;
	for(int i = 0; i < ANTRAG_D; i++) {
		while(bst.buf_len < ANTRAG_Q_BITS) {
			assert(off < ANTRAG_PK_SIZE);
			queue_byte(&bst, buf, &off);
		}
		uint16_t n = dequeue_bits(&bst, ANTRAG_Q_BITS);
		assert(n < ANTRAG_Q);
		pk->h_ntt.coeffs[i] = n;
	}
}

static inline uint16_t decenter(int16_t x) {
  uint16_t xu = x;
  return xu + (ANTRAG_Q & -(xu >> 15));
}
static inline int16_t recenter(uint16_t x) {
	uint16_t t = x - ANTRAG_Q/2;
	return x - (ANTRAG_Q & ((t >> 15) - 1));
}

int encode_sig_vec(uint8_t* buf, const uint16_t s[ANTRAG_D]) {
	int large_val = 1 << SIG_LARGE_BITS;

	BitStream bst = empty_stream;
	size_t off = 0;
	for(int i = 0; i < ANTRAG_D; i++) {
		int16_t k = recenter(s[i]);
		int sign;
		int16_t km;
		if(k < 0) {
			sign = 1;
			km = -k;
		} else {
			sign = 0;
			km = k;
		}

		if(km >= large_val) {
			return 1;
		}
		queue_bits(&bst, sign, 1); // sign bit
		queue_bits(&bst, km, SIG_SMALL_BITS); // low bits in binary
		km >>= SIG_SMALL_BITS;
		queue_bits(&bst, 1 << km, km+1); // extra bits in unary (at most 2^(LARGE_BITS-SMALL_BITS)+1)

		while(has_queued_byte(&bst, i == ANTRAG_D-1)) {
			if(off >= SIG_VECTOR_SIZE) return 1;
			dequeue_byte(&bst, buf, &off);
		}
	}

	while(off < SIG_VECTOR_SIZE) {
		buf[off++] = 0;
	}
	return 0;
}

int encode_sig(uint8_t* buf, const signature* s) {
	if(encode_sig_vec(buf, s->s1.coeffs)) return 1;
	for(int i = 0; i < SALT_BYTES; i++) {
		buf[SIG_VECTOR_SIZE + i] = s->r[i];
	}
	return 0;
}

int decode_sig_vec(const uint8_t* buf, uint16_t s[ANTRAG_D]) {
	int large_val = 1 << SIG_LARGE_BITS;

	BitStream bst = empty_stream;
	size_t off = 0;
	for(int i = 0; i < ANTRAG_D; i++) {
		while(bst.buf_len < 2 + SIG_SMALL_BITS) {
			if(off >= SIG_VECTOR_SIZE) return 1;
			queue_byte(&bst, buf, &off);
		}
		int sign = dequeue_bits(&bst, 1);
		uint16_t km = dequeue_bits(&bst, SIG_SMALL_BITS);

		while(dequeue_bits(&bst, 1) == 0) {
			km += 1 << SIG_SMALL_BITS;
			if(bst.buf_len == 0) {
				if(off >= SIG_VECTOR_SIZE) return 1;
				queue_byte(&bst, buf, &off);
			}
		}

		int16_t k = km * (sign ? -1 : 1);
		if(k < -large_val || k >= large_val) return 1;
		s[i] = decenter(k);
	}

	while(off < SIG_VECTOR_SIZE) {
		if(buf[off++] != 0) return 1;
	}
	return 0;
}

int decode_sig(signature* s, const uint8_t* buf) {
	if(decode_sig_vec(buf, s->s1.coeffs)) return 1;
	for(int i = 0; i < SALT_BYTES; i++) {
		s->r[i] = buf[SIG_VECTOR_SIZE + i];
	}
	return 0;
}
