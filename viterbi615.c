/* K=15 r=1/6 Viterbi decoder with optional Intel or PowerPC SIMD
 * Copyright Feb 2004, Phil Karn, KA9Q
 */
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "fec.h"

/* Create a new instance of a Viterbi decoder */
void *create_viterbi615(int len){

  find_cpu_mode();

  switch(Cpu_mode){
  case PORT:
  default:
    return create_viterbi615_port(len);
#ifdef __VEC__
  case ALTIVEC:
    return create_viterbi615_av(len);
#endif
#ifdef __i386__
  case MMX:
    return create_viterbi615_mmx(len);
  case SSE:
    return create_viterbi615_sse(len);
  case SSE2:
    return create_viterbi615_sse2(len);
#endif
  }
}

/* Initialize Viterbi decoder for start of new frame */
int init_viterbi615(void *p,int starting_state){
    switch(Cpu_mode){
    case PORT:
    default:
      return init_viterbi615_port(p,starting_state);
#ifdef __VEC__
    case ALTIVEC:
      return init_viterbi615_av(p,starting_state);
#endif
#ifdef __i386__
    case MMX:
      return init_viterbi615_mmx(p,starting_state);
    case SSE:
      return init_viterbi615_sse(p,starting_state);
    case SSE2:
      return init_viterbi615_sse2(p,starting_state);
#endif
    }
}

/* Viterbi chainback */
int chainback_viterbi615(
      void *p,
      unsigned char *data, /* Decoded output data */
      unsigned int nbits, /* Number of data bits */
      unsigned int endstate){ /* Terminal encoder state */

    switch(Cpu_mode){
    case PORT:
    default:
      return chainback_viterbi615_port(p,data,nbits,endstate);
#ifdef __VEC__
    case ALTIVEC:
      return chainback_viterbi615_av(p,data,nbits,endstate);
#endif
#ifdef __i386__
    case MMX:
      return chainback_viterbi615_mmx(p,data,nbits,endstate);
    case SSE:
      return chainback_viterbi615_sse(p,data,nbits,endstate);
    case SSE2:
      return chainback_viterbi615_sse2(p,data,nbits,endstate);
#endif
    }
}

/* Delete instance of a Viterbi decoder */
void delete_viterbi615(void *p){
    switch(Cpu_mode){
    case PORT:
    default:
      delete_viterbi615_port(p);
      break;
#ifdef __VEC__
    case ALTIVEC:
      delete_viterbi615_av(p);
      break;
#endif
#ifdef __i386__
    case MMX:
      delete_viterbi615_mmx(p);
      break;
    case SSE:
      delete_viterbi615_sse(p);
      break;
    case SSE2:
      delete_viterbi615_sse2(p);
      break;
#endif
    }
}

/* Update decoder with a block of demodulated symbols
 * Note that nbits is the number of decoded data bits, not the number
 * of symbols!
 */
int update_viterbi615_blk(void *p,unsigned char syms[],int nbits){
    switch(Cpu_mode){
    case PORT:
    default:
      return update_viterbi615_blk_port(p,syms,nbits);
#ifdef __VEC__
    case ALTIVEC:
      return update_viterbi615_blk_av(p,syms,nbits);
#endif
#ifdef __i386__
    case MMX:
      return update_viterbi615_blk_mmx(p,syms,nbits);
    case SSE:
      return update_viterbi615_blk_sse(p,syms,nbits);
    case SSE2:
      return update_viterbi615_blk_sse2(p,syms,nbits);
#endif
    }
}

