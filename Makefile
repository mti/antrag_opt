SRCS = $(filter-out antrag/main.c, $(wildcard antrag/*.c)) $(wildcard falcon/*.c) $(wildcard ntru/*.c)
OBJS = $(SRCS:%.c=build/%.o)
MAIN_SRCS = antrag/main.c $(wildcard tests/*.c)
DEPS = $(OBJS:%.o=%.d) $(MAIN_SRCS:%.c=build/%.d)

TARGET  = ./main

CC      = gcc

CFLAGS  += -Wall -Wextra -Wshadow -Wundef
LDLIBS  = -lm -lgmp

#CFLAGS  += -fsanitize=address
#LDFLAGS += -fsanitize=address
CFLAGS  += -DNDEBUG

#CFLAGS  += -DDEBUG_NTRU

CFLAGS += -g

CFLAGS  += -march=native -O3
CFLAGS  += -flto
LDFLAGS += -flto

.PHONY: default all clean test-main test-ntt test-poly test-codec test-cst-time benchmark build-supercop gen1f gen1s gen2f gen2s gen3f gen3s gen4f gen4s gen5f gen5s antrag1f antrag1s antrag2f antrag2s antrag3f antrag3s antrag4f antrag4s antrag5f antrag5s

default: $(TARGET)
all: default

gen:
	mkdir -p gen

gen1f: scripts/gen_headers.sage | gen
	sage scripts/gen_headers.sage 1f
gen1s: scripts/gen_headers.sage | gen
	sage scripts/gen_headers.sage 1s
gen2f: scripts/gen_headers.sage | gen
	sage scripts/gen_headers.sage 2f
gen2s: scripts/gen_headers.sage | gen
	sage scripts/gen_headers.sage 2s
gen3f: scripts/gen_headers.sage | gen
	sage scripts/gen_headers.sage 3f
gen3s: scripts/gen_headers.sage | gen
	sage scripts/gen_headers.sage 3s
gen4f: scripts/gen_headers.sage | gen
	sage scripts/gen_headers.sage 4f
gen4s: scripts/gen_headers.sage | gen
	sage scripts/gen_headers.sage 4s
gen5f: scripts/gen_headers.sage | gen
	sage scripts/gen_headers.sage 5f
gen5s: scripts/gen_headers.sage | gen
	sage scripts/gen_headers.sage 5s

antrag1f: gen1f | $(TARGET)
antrag1s: gen1s | $(TARGET)
antrag2f: gen2f | $(TARGET)
antrag2s: gen2s | $(TARGET)
antrag3f: gen3f | $(TARGET)
antrag3s: gen3s | $(TARGET)
antrag4f: gen4f | $(TARGET)
antrag4s: gen4s | $(TARGET)
antrag5f: gen5f | $(TARGET)
antrag5s: gen5s | $(TARGET)

gen/const.h gen/fft.h gen/ntru.h gen/ntt.h: scripts/gen_headers.sage | gen
	sage scripts/gen_headers.sage

build:
	mkdir -p build
	mkdir -p build/antrag
	mkdir -p build/tests
	mkdir -p build/falcon
	mkdir -p build/ntru

-include $(DEPS)

build/%.o: %.c | gen/const.h build
	$(CC) $(CFLAGS) -MMD -c $< -o $@

$(TARGET): build/antrag/main.o $(OBJS)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

build/tests/ntt_test: build/tests/ntt_test.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

build/tests/poly_test: build/tests/poly_test.o build/antrag/poly.o build/falcon/fft.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

build/tests/codec_test: build/tests/codec_test.o $(filter-out build/antrag/codec.o, $(OBJS))
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

build/tests/cst_time_test: CFLAGS += -DCST_TIME
build/tests/cst_time_test: build/tests/cst_time_test.o $(OBJS)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

build/tests/rns_test: build/tests/rns_test.o build/ntru/rns.o build/ntru/ipoly.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

build/tests/benchmarks: build/tests/benchmarks.o $(OBJS)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

test-main: $(TARGET)
	$(TARGET)
test-ntt: build/tests/ntt_test
	build/tests/ntt_test
test-poly: build/tests/poly_test
	build/tests/poly_test
test-codec: build/tests/codec_test
	build/tests/codec_test
test-cst-time: build/tests/cst_time_test
	valgrind --track-origins=yes build/tests/cst_time_test
test-rns: build/tests/rns_test
	build/tests/rns_test
benchmark: build/tests/benchmarks
	build/tests/benchmarks

build-supercop: gen/const.h
	bash supercop/build.sh

clean:
	-rm -rf build gen supercop-build
	-rm -f $(TARGET)
	-rm -f tests/ntt_test

