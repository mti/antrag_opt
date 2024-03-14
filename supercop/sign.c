#include "crypto_sign.h"
#include "randombytes.h"

#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include "antrag_api.h"

uint64_t get64() {
	uint64_t n;
	randombytes((uint8_t*) &n, 8);
	return n;
}
uint8_t get8() {
	uint8_t n;
	randombytes(&n, 1);
	return n;
}

int crypto_sign_keypair(uint8_t *pkb, uint8_t *skb) {
	public_key pk;
	secret_key sk;
	keygen_full(&sk, &pk);

	encode_pk(pkb, &pk);
	encode_sk(skb, &sk);
	return 0;
}

int crypto_sign(
	uint8_t *sm, unsigned long long *smlen,
	const uint8_t *m, unsigned long long mlen,
	const uint8_t *skb
) {
	secret_key sk;
	decode_sk(&sk, skb);

	signature s;
	sign(mlen, m, &sk, &s);

	memcpy(sm, m, mlen);
	encode_sig(sm + mlen, &s);
	*smlen = mlen + ANTRAG_SIG_SIZE;

	return 0;
}

int crypto_sign_open(
	uint8_t *m, unsigned long long *mlen,
	const uint8_t *sm, unsigned long long smlen,
	const uint8_t *pkb
) {
	public_key pk;
	decode_pk(&pk, pkb);

	size_t m_len;
	*mlen = m_len = smlen - ANTRAG_SIG_SIZE;
	memcpy(m, sm, m_len);
	signature s;
	decode_sig(&s, sm + m_len);

	bool success = verify(m_len, m, &pk, &s);

	return success ? 0 : -1;
}
