#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "../antrag/api.h"
#include "../antrag/randombytes.h"

#ifndef CST_TIME
#error "cst_time_test.c should be compiled with the CST_TIME macro"
#endif

int main() {
	seed_rng();

	secret_key sk;
	public_key pk;
	keygen_full(&sk, &pk);

	uint8_t msg[MSG_BYTES] = "hey bob it's alice (not mallory)";

	poison(&sk, sizeof(secret_key));
	poisoned_rng = true;

	signature s;
	sign(MSG_BYTES, msg, &sk, &s);
    
	return 0;
}