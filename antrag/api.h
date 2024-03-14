#ifndef API_H
#define API_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include "poly.h"

int keygen_fg(secret_key *sk);
int keygen_full(secret_key *sk, public_key *pk);
void sign(size_t m_len, const uint8_t m[m_len], const secret_key* sk, signature* s);
int verify(size_t m_len, const uint8_t m[m_len], const public_key* pk, const signature* s);
void sampler(const secret_key* sk, const u16_poly* c2i, u16_poly* v1i, u16_poly* v2i);


void encode_sk(uint8_t* buf, const secret_key* sk);
void decode_sk(secret_key* sk, const uint8_t* buf);
void encode_pk(uint8_t* buf, const public_key* pk);
void decode_pk(public_key* pk, const uint8_t* buf);
/* Returns 1 if signature does not fit in buffer */
int encode_sig(uint8_t* buf, const signature* s);
/* Returns 1 if signature could not be decoded */
int decode_sig(signature* s, const uint8_t* buf);

/* Constant-time macros */
#define LSBMASK(c)      (-((c)&1))
#define CMUX(x,y,c)     (((x)&(LSBMASK(c)))^((y)&(~LSBMASK(c))))
#define CFLIP(x,c)      CMUX(x,-(x),c)
#define CZERO64(x)      ((~(x)&((x)-1))>>63)

#endif
