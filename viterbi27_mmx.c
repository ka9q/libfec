/* K=7 r=1/2 Viterbi decoder for MMX
 * Copyright Feb 2004, Phil Karn, KA9Q
 */
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <mmintrin.h>
#include "fec.h"

typedef union { char c[64]; __m64 v[8];} decision_t;
typedef union { unsigned char c[64]; __m64 v[8];} metric_t;

unsigned char Mettab27_1[256][32] __attribute__ ((aligned(16)));
unsigned char Mettab27_2[256][32] __attribute__ ((aligned(16)));
static int Init = 0;

/* State info for instance of Viterbi decoder
 * Don't change this without also changing references in mmxbfly27.s!
 */
struct v27 {
  metric_t metrics1; /* path metric buffer 1 */
  metric_t metrics2; /* path metric buffer 2 */
  decision_t *dp;          /* Pointer to current decision */
  metric_t *old_metrics,*new_metrics; /* Pointers to path metrics, swapped on every bit */
  decision_t *decisions;   /* Beginning of decisions for block */
};

/* Initialize Viterbi decoder for start of new frame */
int init_viterbi27_mmx(void *p,int starting_state){
  struct v27 *vp = (struct v27 *)p;
  int i;

  for(i=0;i<64;i++)
    vp->metrics1.c[i] = 63;

  vp->old_metrics = &vp->metrics1;
  vp->new_metrics = &vp->metrics2;
  vp->dp = vp->decisions;
  vp->old_metrics->c[starting_state & 63] = 0; /* Bias known start state */
  return 0;
}

/* Create a new instance of a Viterbi decoder */
void *create_viterbi27_mmx(int len){
  int state;
  struct v27 *vp;

  if(Init == 0){
    /* Initialize branch tables */
    for(state=0;state < 32;state++){
      int symbol,bit1,bit2;
      bit1 = parity((2*state) & V27POLYA);
      bit2 = parity((2*state) & V27POLYB);
      for(symbol = 0;symbol < 256;symbol++){
	Mettab27_1[symbol][state] = (bit1 ? (255-symbol):symbol) / 16;
	Mettab27_2[symbol][state] = (bit2 ? (255-symbol):symbol) / 16;
      }
    }
    Init++;
  }
  vp = (struct v27 *)malloc(sizeof(struct v27));
  vp->decisions = (decision_t *)malloc((len+6)*sizeof(decision_t));
  init_viterbi27_mmx(vp,0);
  return vp;
}

/* Viterbi chainback */
int chainback_viterbi27_mmx(
      void *p,
      unsigned char *data, /* Decoded output data */
      unsigned int nbits, /* Number of data bits */
      unsigned int endstate){ /* Terminal encoder state */

  struct v27 *vp = (struct v27 *)p;
  decision_t *d = (decision_t *)vp->decisions;

  endstate &= 63;
  d += 6; /* Look past tail */
  while(nbits-- != 0){
    int k;

    k = d[nbits].c[endstate>>2] & 1;
    data[nbits>>3] = endstate = (endstate >> 1) | (k << 7);
  }
  return 0;
}

/* Delete instance of a Viterbi decoder */
void delete_viterbi27_mmx(void *p){
  struct v27 *vp = p;

  if(vp != NULL){
    free(vp->decisions);
    free(vp);
  }
}
