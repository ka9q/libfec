# Makefile for FEC library
# Copyright 2014 Phil Karn, KA9Q
# The 2014 version does a lot of housecleaning. Apple killed the PowerPC,
# and 64-bit machines have taken over so all the Altivec and early SSE code
# is obsolete. We can pretty much assume 64 bits, SSE 4.2 and popcount

# May be used under the terms of the GNU Lesser General Public License (LGPL)

CFLAGS=-O3 -march=native -mtune=native -g -Wall -Wunreachable-code
LDFLAGS=-g -march=native -mtune=native
CC=gcc

OBJ=hybridtest fanotest vtest27sse vtest27port vtest224sse vtest224port simtest

LIB_PORT=libfec_port.so
LIB_SSE=libfec_sse.so

# select a library for use with test programs
LIB=libfec_sse.so
#LIB=libfec_port.so


#LIB_PORT_OBJ=sim.o viterbi27_port.o viterbi224_port.o

LIB_PORT_OBJ=viterbi27_port.o viterbi224_port.o fano.o metrics.o sim.o

LIB_SSE_OBJ=viterbi27_sse2.o viterbi224_sse2.o fano.o metrics.o sim.o


#LIB_SSE_OBJ=sim.o viterbi27_sse2.o viterbi29_port.o viterbi39_sse2.o \
#	viterbi615_sse2.o encode_rs_char.o encode_rs_int.o encode_rs_8.o \
#	decode_rs_char.o decode_rs_int.o decode_rs_8.o \
#	init_rs_char.o init_rs_int.o ccsds_tab.o \
#	encode_rs_ccsds.o decode_rs_ccsds.o ccsds_tal.o \
#	dotprod_port.o peakval_port.o sumsq_port.o


#all: $(LIB_PORT) $(LIB_SSE) vtest27sse vtest27port
all: vtest27sse vtest27port

$(LIB_PORT): $(LIB_PORT_OBJ)
	gcc -shared -Xlinker -soname=$@ -o $@ -Wl,-whole-archive $^ -Wl,-no-whole-archive -lc

$(LIB_SSE): $(LIB_SSE_OBJ)
	gcc -shared -Xlinker -soname=$@ -o $@ -Wl,-whole-archive $^ -Wl,-no-whole-archive -lc


hybridtest: hybridtest.o encode.o viterbi224_sse2.o fano.o metrics.o sim.o
	gcc $(LDFLAGS) -o $@ $^ -lm

fanotest: fanotest.o encode.o fano.o metrics.o sim.o
	gcc $(LDFLAGS) -o $@ $^ -lm

simtest: simtest.o sim.o
	gcc $(LDFLAGS) -o $@ $^ -lm

vtest224sse: vtest224.o encode.o viterbi224_sse2.o sim.o
	gcc $(LDFLAGS) -o $@ $^ -lm

vtest224port: vtest224.o encode.o viterbi224_port.o sim.o
	gcc $(LDFLAGS) -o $@ $^ -lm

vtest27sse: vtest27.o viterbi27_sse2.o sim.o
	gcc $(LDFLAGS) -o $@ $^ -lm

vtest27port: vtest27.o viterbi27_port.o sim.o
	gcc $(LDFLAGS) -o $@ $^ -lm

vtest27.o: vtest27.c fec.h

viterbi27_port.o: viterbi27_port.c fec.h

viterbi27_sse2.o: viterbi27_sse2.c fec.h

vtest224.o: vtest224.c sim.h viterbi224.h

fano.o: fano.c fano.h

viterbi224_sse2.o: viterbi224_sse2.c viterbi224.h

viterbi224_port.o: viterbi224_port.c viterbi224.h

hybridtest.o: hybridtest.c fano.h viterbi224.h

encode.o: encode.c code.h



clean:
	rm -f *.o $(LIB_PORT) $(LIB_SSE) $(OBJ) core
