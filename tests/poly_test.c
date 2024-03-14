#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../antrag/poly.h"
#include "../gen/fft.h"

#define TURN 6.2831853071795865

double complex eval(poly* p, double complex z) {
    double complex zn = 1;
    double complex res = 0;
    for(int i = 0; i < ANTRAG_D; i++) {
        res += p->coeffs[i].v * zn;
        zn *= z;
    }
    return res;
}

bool close_enough(double x, double y) {
    return fabs(x - y) < 1e-4;
}

void print_vec(fpr f[ANTRAG_D]) {
	printf("[");
	for (int i = 0; i < ANTRAG_D; i++) {
		if(i != 0) printf(", ");
		printf("%.2f", f[i].v);
	}
	printf("]\n");
}

int main() {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    srand(ts.tv_sec * 1000 + ts.tv_nsec / 1000000);

    // Generate random polynomials
    poly p1, p2;
    for(int i = 0; i < ANTRAG_D; i++) {
        p1.coeffs[i].v = rand() % 100 - 50;
        p2.coeffs[i].v = rand() % 100 - 50;
    }

    poly p3 = p1;
    poly_add(&p3, &p2);
    for(int i = 0; i < ANTRAG_D; i++) {
        if(p3.coeffs[i].v != p1.coeffs[i].v + p2.coeffs[i].v) {
            printf("Polynomial addition does not work!\n");
            exit(EXIT_FAILURE);
        }
    }
    printf("Polynomial addition works.\n");

    poly v1 = p1;
    FFT(&v1);

    for(int i = 0; i < FFT_SIZE; i++) {
        int k = ROOT_IDX[i];
        double complex w = cexp(I * TURN / ANTRAG_M * k);
        double complex pw = eval(&p1, w);
        if(!close_enough(creal(pw), v1.coeffs[i].v)
            || !close_enough(cimag(pw), v1.coeffs[i + FFT_SIZE].v)) {
            printf("FFT does not work!\n");
            exit(EXIT_FAILURE);
        }
    }
    printf("FFT works.\n");

    poly p1rt = v1;
    invFFT(&p1rt);

    for(int i = 0; i < ANTRAG_D; i++) {
        if(!close_enough(p1.coeffs[i].v, p1rt.coeffs[i].v)) {
            printf("FFT does not round trip!\n");
            exit(EXIT_FAILURE);
        }
    }
    printf("FFT round trips.\n");

    p3 = p2;
    FFT(&p3);
    pointwise_mul(&p3, &v1);
    invFFT(&p3);

    // Compute product directly (quadratic algorithm)
    double p3_arr[ANTRAG_D*2-1];
    for(int i = 0; i < ANTRAG_D*2-1; i++) {
        p3_arr[i] = 0;
    }
    for(int i = 0; i < ANTRAG_D; i++) {
        for(int j = 0; j < ANTRAG_D; j++) {
            double prod = p1.coeffs[i].v * p2.coeffs[j].v;
            p3_arr[i+j] += prod;
        }
    }
    // Apply modulo
    // Phi = X^n + 1 in binary case
    // Phi = X^n - X^(n/2) + 1) in ternary case
    for(int i = ANTRAG_D*2-2; i >= ANTRAG_D; i--) {
        p3_arr[i-ANTRAG_D] -= p3_arr[i];
        if(ANTRAG_LOG3_D) {
            p3_arr[i-ANTRAG_D/2] += p3_arr[i];
        }
        p3_arr[i] = 0;
    }

    poly p3_cor;
    for(int i = 0; i < ANTRAG_D; i++) {
        p3_cor.coeffs[i].v = p3_arr[i];
    }

    for(int i = 0; i < ANTRAG_D; i++) {
        if(!close_enough(p3.coeffs[i].v, p3_cor.coeffs[i].v)) {
            printf("FFT multiplication does not work!\n");
            exit(EXIT_FAILURE);
        }
    }
    printf("FFT multiplication works.\n");

    return EXIT_SUCCESS;
}
