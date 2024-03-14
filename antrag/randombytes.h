#ifndef RNG_H
#define RNG_H
#include <stdlib.h>
#include <stdint.h>
#include "../falcon/inner.h"

/* `declassify(<condition>)` prevents a branch from triggering non-constant
 * time detection (cf. tests/cst_time_test.c). This is useful to avoid false
 * positives when using rejection sampling, where newly-generated secret data
 * only ends up being used if it passes the test.
**/

#ifdef CST_TIME
#include <valgrind/memcheck.h>

#define poison(addr, len) VALGRIND_MAKE_MEM_UNDEFINED(addr, len)
#define unpoison(addr, len) VALGRIND_MAKE_MEM_DEFINED(addr, len)

static inline int declassify(int cond) {
    unpoison(&cond, sizeof(int));
    return cond;
}

extern int poisoned_rng;
#else

static inline int declassify(int cond) {
    return cond;
}

#endif

void seed_rng();

void randombytes(uint8_t *buf, unsigned long n);
uint64_t get64();
uint8_t get8();


#endif
