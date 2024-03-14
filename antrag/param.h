#ifndef PARAM_H
#define PARAM_H

#include "../gen/const.h"

#define MSG_BYTES 32

/* sigma^2 = r^2 * alpha^2 * q */
/* gamma^2 = slack^2 * sigma^2 * 2d */
/* slack = 1.042 */

#define ANTRAG_SLACK 1.042
#define SIGMA_SQUARE (R_SQUARE * ANTRAG_ALPHA * ANTRAG_ALPHA * ANTRAG_Q)
#define GAMMA_SQUARE (ANTRAG_SLACK * ANTRAG_SLACK * SIGMA_SQUARE * ANTRAG_D * 2)


#endif
