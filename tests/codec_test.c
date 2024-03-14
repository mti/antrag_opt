#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "../antrag/codec.c"
#include "../antrag/randombytes.h"

int main() {
	seed_rng();

	secret_key sk;
	public_key pk;
	keygen_full(&sk, &pk);

	uint8_t msg[MSG_BYTES] = "hey bob it's alice (not mallory)";
	
	uint8_t sk_buf[ANTRAG_SK_SIZE];
	encode_sk(sk_buf, &sk);
	secret_key sk2;
	decode_sk(&sk2, sk_buf);
	if(memcmp(&sk, &sk2, sizeof(secret_key)) != 0) {
		printf("Secret key does not round trip!\n");
		return 1;
	}
	printf("Secret key round trips.\n");
	
	uint8_t pk_buf[ANTRAG_PK_SIZE];
	encode_pk(pk_buf, &pk);
	public_key pk2;
	decode_pk(&pk2, pk_buf);
	if(memcmp(&pk, &pk2, sizeof(public_key)) != 0) {
		printf("Public key does not round trip!\n");
		return 1;
	}
	printf("Public key round trips.\n");

	signature s;
	uint8_t s_buf[ANTRAG_SIG_SIZE];
	
	do {
		printf("Generating signature...\n");
		sign(MSG_BYTES, msg, &sk, &s);
	} while(encode_sig(s_buf, &s));

	signature s2;
	if(decode_sig(&s2, s_buf)) {
		printf("Signature decoding does not work!\n");
		return 1;
	}

	if(memcmp(&s.s1, &s2.s1, sizeof(s.s1)) != 0) {
		printf("Signature vectors do not round trip!\n");
		return 1;
	}
	if(memcmp(&s.r, &s2.r, sizeof(s.r)) != 0) {
		printf("Signature salt does not round trip!\n");
		return 1;
	}
	printf("Signature round trips.\n");

	if(!verify(MSG_BYTES, msg, &pk2, &s2)) {
		printf("Signature verification fails after serialization!\n");
		return 1;
	}
	printf("Signature verification succeeds.\n");

	printf("Testing rejection rate... (SIG_VECTOR_SIZE = %d)\n", SIG_VECTOR_SIZE);
	int failures = 0;
	for(int i = 0; i < 10000; i++) {
		sign(MSG_BYTES, msg, &sk, &s);
		if(encode_sig(s_buf, &s)) failures++;
	}
	printf("Percentage of rejections: %.2f%%\n", failures / 100.0);
	
	return 0;
}