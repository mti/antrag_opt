#ifndef NTRU_H
#define NTRU_H

#include <stdint.h>
#include <stdbool.h>

#include "../gen/const.h"

bool solve_ntru(const int8_t f[ANTRAG_D], const int8_t g[ANTRAG_D], int8_t F[ANTRAG_D], int8_t G[ANTRAG_D]);

#endif
