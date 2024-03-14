#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>  
#include <math.h>  

#include "api.h"
#include "param.h"
#include "poly.h"
#include "cpucycles.h"
#include "normaldist.h"
#include "samplerZ.h"
#include "randombytes.h"

int intcmp(const void *x, const void *y)
{
  int ix = *(int*)x, iy = *(int*)y;
  return (ix>iy) - (ix<iy);
}
int uint64cmp(const void *x, const void *y)
{
  uint64_t ix = *(uint64_t*)x, iy = *(uint64_t*)y;
  return (ix>iy) - (ix<iy);
}
int doublecmp(const void *x, const void *y)
{
  double ix = *(double*)x, iy = *(double*)y;
  return (ix>iy) - (ix<iy);
}

void speed(){
  secret_key sk;
  public_key pk;
  signature s;

  uint8_t m[MSG_BYTES] = {0x46,0xb6,0xc4,0x83,0x3f,0x61,0xfa,0x3e,0xaa,0xe9,0xad,0x4a,0x68,0x8c,0xd9,0x6e,0x22,0x6d,0x93,0x3e,0xde,0xc4,0x64,0x9a,0xb2,0x18,0x45,0x2,0xad,0xf3,0xc,0x61};
  
  keygen_full(&sk, &pk);

  int iter = 10000;
  clock_t start_time = clock();
  uint64_t start = cpucycles();
  
  for(int i=0; i < iter; ++i){
    sign(MSG_BYTES, m, &sk, &s);
  }
  uint64_t stop = cpucycles();
  clock_t stop_time = clock();
  
  double delta = (double)(stop - start)/iter;
  double delta_time = (double)(stop_time - start_time)/iter*1.0e6/CLOCKS_PER_SEC;
  printf("Average number of cycles per signing: %f (%f us)\n", delta, delta_time);

  start_time = clock();
  start = cpucycles();
  for(int i=0; i < iter; ++i){
    verify(MSG_BYTES, m, &pk, &s);
  }
  stop = cpucycles();
  stop_time = clock();
  delta = (double)(stop - start)/iter;
  delta_time = (double)(stop_time - start_time)/iter*1.0e6/CLOCKS_PER_SEC;
  printf("Average number of cycles per verification: %f (%f us)\n", delta, delta_time);
}

void speed_keygen(){
  secret_key sk;
  public_key pk;

  int iter = 100;
  clock_t start_time = clock();
  uint64_t start = cpucycles();
  
  for(int i=0; i < iter; ++i){
    keygen_full(&sk, &pk);
  }
  uint64_t stop = cpucycles();
  clock_t stop_time = clock();
  
  double delta = (double)(stop - start)/iter;
  double delta_time = (double)(stop_time - start_time)/iter*1.0e3/CLOCKS_PER_SEC;
  printf("Average number of cycles per keygen: %f (%f ms)\n", delta, delta_time);
}

#define SIGNVERIF_TESTS 200
#define KEYGEN_TESTS 200
#define BENCHMARK_ITER 200

struct {
  double kcycles[BENCHMARK_ITER];
  double us[BENCHMARK_ITER];
  double kcycles_avg;
  double us_avg;

  int iter;
  clock_t start_time;
  uint64_t start;
  clock_t stop_time;
  uint64_t stop;
} bench;

void bench_start() {
  bench.kcycles_avg = 0.;
  bench.us_avg = 0.;
  bench.iter = 0;
}
// Use macros for least timing overhead
#define bench_before() { \
  bench.start_time = clock(); \
  bench.start = cpucycles(); \
}
#define bench_after() { \
  bench.stop = cpucycles(); \
  bench.stop_time = clock(); \
}
void bench_store() {
  int i = bench.iter;
  bench.kcycles[i]   = (double)(bench.stop - bench.start)/1.0e3;
  bench.us[i]       = (double)(bench.stop_time - bench.start_time)*1.0e6/CLOCKS_PER_SEC;
  bench.kcycles_avg += bench.kcycles[i];
  bench.us_avg     += bench.us[i];
  bench.iter++;
}
void bench_results(const char* bench_name, bool is_slow) {
  int n = bench.iter;
  bench.kcycles_avg /= n;
  bench.us_avg     /= n;

  qsort(bench.kcycles, n, sizeof(uint64_t), doublecmp);
  qsort(bench.us,      n, sizeof(double),   doublecmp);
  
  if(is_slow) {
    for(int i = 0; i < n; i++) {
      bench.kcycles[i] /= 1.0e3;
      bench.us[i] /= 1.0e3;
    }
    bench.kcycles_avg /= 1.0e3;
    bench.us_avg /= 1.0e3;
  }

  printf(
    is_slow
      ? "%11s Mcycles      %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n"
      : "%11s kcycles      %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f\n",
    bench_name,
    bench.kcycles[0],              bench.kcycles[n/4],
    bench.kcycles[n/2], bench.kcycles[3*n/4],
    bench.kcycles[n-1], bench.kcycles_avg);
  printf(
    is_slow
      ? "%11s speed (ms)   %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n"
      : "%11s speed (µs)   %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f\n",
    bench_name,
    bench.us[0],              bench.us[n/4],
    bench.us[n/2], bench.us[3*n/4],
    bench.us[n-1], bench.us_avg);
}

int main(){
  srand(time(0));
  seed_rng();
  printf("Hello world, signature is Antrag-%u with q=%u\n", ANTRAG_D, ANTRAG_Q);
  secret_key sk;
  public_key pk;
  signature s;

  uint8_t m[MSG_BYTES] = {0x46,0xb6,0xc4,0x83,0x3f,0x61,0xfa,0x3e,0xaa,0xe9,0xad,0x4a,0x68,0x8c,0xd9,0x6e,0x22,0x6d,0x93,0x3e,0xde,0xc4,0x64,0x9a,0xb2,0x18,0x45,0x2,0xad,0xf3,0xc,0x61};
  
  printf("\n* Generate initial key pair.\n");
  keygen_full(&sk, &pk);
  printf("  ...done.\n\n");

  printf("* Test correctness of the scheme.\n");
  int correct = 0;
  for(int i=0; i<SIGNVERIF_TESTS; i++) {
    sign(MSG_BYTES, m, &sk, &s);
    correct += verify(MSG_BYTES, m, &pk, &s);
  }
  printf("  %d/%d correct signatures. (%s).\n\n", correct, SIGNVERIF_TESTS,
    (correct == SIGNVERIF_TESTS)?"ok":"ERROR!");

  printf("* Test keygen repetitions (alpha=%.2f, radius=%.3f, tests=%d).\n\n",
    ANTRAG_ALPHA, ANTRAG_RADIUS, KEYGEN_TESTS);

  printf("                       min  lowq  med.  uppq   max  avg.\n");
  printf("--------------------------------------------------------\n");
  int trials[KEYGEN_TESTS];
  double trialsavg = 0.;
  for(int i=0; i<KEYGEN_TESTS; i++) {
    trials[i] = keygen_fg(&sk);
    trialsavg+= trials[i];
  }
  qsort(trials, KEYGEN_TESTS, sizeof(int), intcmp);
  trialsavg /= KEYGEN_TESTS;

  printf("keygen_fg repetitions %4d %5d %5d %5d %5d %5.2f\n",
    trials[0], trials[KEYGEN_TESTS/4], trials[KEYGEN_TESTS/2],
    trials[3*KEYGEN_TESTS/4], trials[KEYGEN_TESTS-1], trialsavg);

  trialsavg = 0.;
  for(int i=0; i<KEYGEN_TESTS; i++) {
    trials[i] = keygen_full(&sk, &pk);
    trialsavg+= trials[i];
  }
  qsort(trials, KEYGEN_TESTS, sizeof(int), intcmp);
  trialsavg /= KEYGEN_TESTS;
  printf("keygen repetitions   %5d %5d %5d %5d %5d %5.2f\n",
    trials[0], trials[KEYGEN_TESTS/4], trials[KEYGEN_TESTS/2],
    trials[3*KEYGEN_TESTS/4], trials[KEYGEN_TESTS-1], trialsavg);
  printf("----------------------------------------------------------\n\n");

  printf("* Benchmarking the scheme.\n\n");
  printf("                           min  lowq  med.  uppq   max  avg.\n");
  printf("------------------------------------------------------------\n");

  bench_start();
  for(int i=0; i < BENCHMARK_ITER; ++i){
    bench_before();
    keygen_full(&sk, &pk);
    bench_after();
    bench_store();
  }
  bench_results("keygen", true);

  bench_start();
  for(int i=0; i < BENCHMARK_ITER; ++i){
    bench_before();
    sign(MSG_BYTES, m, &sk, &s);
    bench_after();
    bench_store();
  }
  bench_results("sign", false);
  
  bench_start();
  for(int i=0; i < BENCHMARK_ITER; ++i){
    bench_before();
    verify(MSG_BYTES, m, &pk, &s);
    bench_after();
    bench_store();
  }
  bench_results("verif", false);

  uint8_t sk_buf[ANTRAG_SK_SIZE];
  encode_sk(sk_buf, &sk);
  uint8_t pk_buf[ANTRAG_SK_SIZE];
  encode_pk(pk_buf, &pk);
  uint8_t sig_buf[ANTRAG_SIG_SIZE];

  bench_start();
  for(int i=0; i < BENCHMARK_ITER; ++i){
    bench_before();
    decode_sk(&sk, sk_buf);
    encode_sig(sig_buf, &s);
    bench_after();
    bench_store();
  }
  bench_results("sign_codec", false);

  bench_start();
  for(int i=0; i < BENCHMARK_ITER; ++i){
    bench_before();
    decode_pk(&pk, pk_buf);
    decode_sig(&s, sig_buf);
    bench_after();
    bench_store();
  }
  bench_results("verif_codec", false);
  
  printf("-------------------------------------------------------------\n\n");

  return 0;
}
