// K=7 r=1/2 Viterbi decoder in portable C
// Copyright 2004-2014, Phil Karn, KA9Q
// May be used under the terms of the GNU Lesser General Public License (LGPL)
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <limits.h>
#include <stdint.h>
#include "fec.h"

// Only needs to handle 16 bits
#if 0
static inline int parity(unsigned int x){
  x ^= x >> 16;
  x ^= x >> 8;
  x ^= x >> 4;
  x ^= x >> 2;
  return (x ^ (x >> 1)) & 1;
}
#else
static inline int parity(int x){
  return __builtin_parity(x);
}
#endif

int encode_27(
   unsigned char *symbols,	// Output buffer, 2*8*nbytes
   const unsigned char *data,	// Input buffer, nbytes
   unsigned int nbytes)		// Number of bytes in data
{
  unsigned int encstate;
  int i;
  
  encstate = 0;       // Encoder shift register is bits 15-8
  while(nbytes-- != 0){
    encstate |= *data++;
    for(i=0;i<8;i++){  // Transmit 8 bits, MSB first
      encstate <<= 1;
      *symbols++ = parity(encstate & (V27POLYB << 8));
      *symbols++ = !parity(encstate & (V27POLYA << 8));
    }
  }
  if((encstate & (63 << 8)) != 0)
    return -1;             // Warn if frame wasn't tailed to 0
  return 0;
}

typedef union { uint32_t w[64]; } metric_t;
typedef union { uint32_t w[2]; uint64_t ll; } decision_t;
static union branchtab27 { uint8_t c[32]; } Branchtab27[2] __attribute__ ((aligned(16)));

// State info for instance of Viterbi decoder
struct v27 {
  metric_t metrics1;       // path metric buffer 1
  metric_t metrics2;       // path metric buffer 2
  decision_t *dp;          // Pointer to current decision
  metric_t *old_metrics,*new_metrics; // Pointers to path metrics, swapped on every bit
  decision_t *decisions;   // Beginning of decisions for block
};

// Initialize Viterbi decoder for start of new frame
int init_viterbi27(void *p,int starting_state){
  struct v27 *vp = p;
  int i;

  if(p == NULL)
    return -1;
  for(i=0;i<64;i++)
    vp->metrics1.w[i] = 63;

  vp->old_metrics = &vp->metrics1;
  vp->new_metrics = &vp->metrics2;
  vp->dp = vp->decisions;
  vp->old_metrics->w[starting_state & 63] = 0; /* Bias known start state */
  return 0;
}

// Create a new instance of a Viterbi decoder
void *create_viterbi27(int len){
  struct v27 *vp;
  int state;

  for(state=0;state < 32;state++){
    Branchtab27[0].c[state] = parity((2*state) & V27POLYB) ? 255 : 0;
    Branchtab27[1].c[state] = !parity((2*state) & V27POLYA) ? 255 : 0;
  }

  if((vp = malloc(sizeof(struct v27))) == NULL)
     return NULL;
  if((vp->decisions = malloc(len*sizeof(decision_t))) == NULL){
    free(vp);
    return NULL;
  }
  init_viterbi27(vp,0);
  return vp;
}

// Viterbi chainback
int chainback_viterbi27(
      void *p,
      unsigned char *data,    // Decoded output data
      unsigned int nbits,     // Number of data bits
      unsigned int endstate){ // Terminal encoder state
  struct v27 *vp = p;
  decision_t *d;
  unsigned char dbyte = 0;

  if(p == NULL)
    return -1;

  d = vp->decisions;
  endstate &= 63;
  while(nbits-- != 0){
    int bit;

    // Accumulate decoded data bits as they fall off the right end of endstate
    dbyte = ((endstate & 1) << 7) | (dbyte >> 1);
    if((nbits % 8) == 0)
      data[nbits/8] = dbyte;

    // Extract decision and push onto left (tail) end of encoder
    bit = (d[nbits].w[endstate/32] >> (endstate%32)) & 1;
    endstate = (bit << 5) | (endstate >> 1);
  }
  return 0;
}

// Delete instance of a Viterbi decoder
void delete_viterbi27(void *p){
  struct v27 *vp = p;

  if(vp != NULL){
    free(vp->decisions);
    free(vp);
  }
}

// C-language butterfly
#define BFLY(i) {\
    unsigned int metric,m0,m1; unsigned long long decision; \
    metric = (Branchtab27[0].c[i] ^ sym0) + (Branchtab27[1].c[i] ^ sym1);\
    m0 = vp->old_metrics->w[i] + metric;\
    m1 = vp->old_metrics->w[i+32] + (510 - metric);\
    decision = (signed int)(m0-m1) > 0;\
    vp->new_metrics->w[2*i] = decision ? m1 : m0;\
    d->ll |= decision << (2*i);\
    m0 -= (metric+metric-510);\
    m1 += (metric+metric-510);\
    decision = (signed int)(m0-m1) > 0;\
    vp->new_metrics->w[2*i+1] = decision ? m1 : m0;\
    d->ll |= decision << (2*i+1);\
}

// Update decoder with a block of demodulated symbols
// Note that nbits is the number of decoded data bits, not the number
// of symbols!
int update_viterbi27_blk(void *p,unsigned char *syms,int nbits){
  struct v27 *vp = p;
  void *tmp;
  decision_t *d;

  if(p == NULL)
    return -1;
  d = (decision_t *)vp->dp;
  while(nbits--){
    unsigned char sym0,sym1;

    d->ll = 0;
    sym0 = *syms++;
    sym1 = *syms++;
    
#if 1
    BFLY(0);
    BFLY(1);
    BFLY(2);
    BFLY(3);
    BFLY(4);
    BFLY(5);
    BFLY(6);
    BFLY(7);
    BFLY(8);
    BFLY(9);
    BFLY(10);
    BFLY(11);
    BFLY(12);
    BFLY(13);
    BFLY(14);
    BFLY(15);
    BFLY(16);
    BFLY(17);
    BFLY(18);
    BFLY(19);
    BFLY(20);
    BFLY(21);
    BFLY(22);
    BFLY(23);
    BFLY(24);
    BFLY(25);
    BFLY(26);
    BFLY(27);
    BFLY(28);
    BFLY(29);
    BFLY(30);
    BFLY(31);
#else
    for(i=0;i<32;i++)
      BFLY(i);
#endif

    d++;
    // Swap pointers to old and new metrics
    tmp = vp->old_metrics;
    vp->old_metrics = vp->new_metrics;
    vp->new_metrics = tmp;
  }    
  vp->dp = d;
  return 0;
}
