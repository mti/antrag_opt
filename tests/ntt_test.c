#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../falcon/ntt.c"

void print_p(uint16_t* p) {
    printf("=> [");
    for(int i = 0; i < ANTRAG_D; i++) {
        if(i != 0) printf(", ");
        printf("%d", p[i]);
    }
    printf("]\n");
}

int main() {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    srand(ts.tv_sec * 1000 + ts.tv_nsec / 1000000);

    // Generate random polynomial
    uint16_t p[ANTRAG_D];
    for(int i = 0; i < ANTRAG_D; i++) {
        p[i] = rand() % ANTRAG_Q;
    }

    // Compute NTT
    uint16_t v[ANTRAG_D];
    memcpy(v, p, sizeof(p));
    Zf(NTT)(v);

    // Compute NTT directly (quadratic algorithm)
    uint16_t v_cor[ANTRAG_D];
    for(int i = 0; i < NTT_SIZE; i++) {
        int gpi = ROOT_IDX[i];
        fqx x = to_fqx(0);
        for (int j = 0; j < ANTRAG_D; j++) {
            fqx c = to_fqx(to_monty(p[j]));
            fqx wij = MTH_ROOTS[gpi * j % ANTRAG_M];
            x = fqx_add(x, fqx_mul(wij, c));
        }
        fqx_write(v_cor,i, x);
    }

    // Compare results
    for(int i = 0; i < ANTRAG_D; i++) {
        if(v_cor[i] != v[i]) {
            printf("Direct NTT does not work!\n");
            exit(EXIT_FAILURE);
        }
    }
    printf("Direct NTT works.\n");
    
    // Inverse NTT
    uint16_t p_rt[ANTRAG_D];
    memcpy(p_rt, v, sizeof(v));
    Zf(iNTT)(p_rt);

    // Check we get back the original polynomial
    for(int i = 0; i < ANTRAG_D; i++) {
        if(p[i] != p_rt[i]) {
            printf("NTT does not round trip!\n");
            exit(EXIT_FAILURE);
        }
    }
    printf("NTT round trips.\n");
    
    // Generate second random polynomial
    uint16_t p2[ANTRAG_D];
    for(int i = 0; i < ANTRAG_D; i++) {
        p2[i] = rand() % ANTRAG_Q;
    }
    
    uint16_t v2[ANTRAG_D];
    memcpy(v2, p2, sizeof(p2));
    Zf(NTT)(v2);

    // Compute product of p1 and p2 via NTT
    uint16_t p3[ANTRAG_D];
    for(int i = 0; i < NTT_SIZE; i++) {
        fqx_write(p3,i, fqx_mul(fqx_read(v,i), fqx_read(v2,i)));
    }
    Zf(iNTT)(p3);

    // Compute product directly (quadratic algorithm)
    uint16_t p3_cor[ANTRAG_D*2-1];
    for(int i = 0; i < ANTRAG_D*2-1; i++) {
        p3_cor[i] = 0;
    }
    for(int i = 0; i < ANTRAG_D; i++) {
        for(int j = 0; j < ANTRAG_D; j++) {
            uint32_t prod = ((uint32_t) p[i] * p2[j]) % ANTRAG_Q;
            p3_cor[i+j] = (p3_cor[i+j] + prod) % ANTRAG_Q;
        }
    }
    // Apply modulo
    // Phi = X^n + 1 in binary case
    // Phi = X^n - X^(n/2) + 1) in ternary case
    for(int i = ANTRAG_D*2-2; i >= ANTRAG_D; i--) {
        if(p3_cor[i] != 0) {
            p3_cor[i-ANTRAG_D] = (p3_cor[i-ANTRAG_D] + ANTRAG_Q - p3_cor[i]) % ANTRAG_Q;
            if(ANTRAG_LOG3_D) {
                p3_cor[i-ANTRAG_D/2] = (p3_cor[i-ANTRAG_D/2] + p3_cor[i]) % ANTRAG_Q;
            }
        }
    }

    // Compare results
    for(int i = 0; i < ANTRAG_D; i++) {
        if(p3[i] != p3_cor[i]) {
            printf("NTT multiplication does not work!\n");
            exit(EXIT_FAILURE);
        }
    }
    printf("NTT multiplication works.\n");
    
    return EXIT_SUCCESS;
}
