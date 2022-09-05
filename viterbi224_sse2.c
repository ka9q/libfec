// K=24 r=1/2 Viterbi decoder for ICE
// Copyright April 2014, Phil Karn, KA9Q
#include <smmintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include <assert.h>

#include "viterbi224.h"
#include "code.h"

static inline int parity(int x){
  return __builtin_parity(x);
}


typedef union { uint32_t w[1<<18]; uint16_t s[1<<19];} decision_t;
typedef union { int16_t s[1<<23]; __m128i v[1<<20];} metric_t;

static union branchtab224 { uint16_t s[1<<22]; __m128i v[1<<19];} Branchtab224[2];

// State info for instance of Viterbi decoder
struct v224 {
  metric_t metrics1; // path metric buffer 1
  metric_t metrics2; // path metric buffer 2
  void *dp;          // Pointer to current decision
  metric_t *old_metrics,*new_metrics; // Pointers to path metrics, swapped on every bit
  void *decisions;   // Beginning of decisions for block
};

// Initialize Viterbi decoder for start of new frame
int init_viterbi224(void *p,int starting_state){
  struct v224 *vp = p;
  int i;

  if(p == NULL)
    return -1;

  for(i=0;i<(1<<(K-1));i++)
    vp->metrics1.s[i] = (SHRT_MIN+5000);

  vp->old_metrics = &vp->metrics1;
  vp->new_metrics = &vp->metrics2;
  vp->dp = vp->decisions;
  vp->old_metrics->s[starting_state & ((1<<(K-1))-1)] = SHRT_MIN; // Bias known start state
  return 0;
}

// Create a new instance of a Viterbi decoder
void *create_viterbi224(int len){
  void *p;
  struct v224 *vp;
  int i;
  int state;

  // Ordinary malloc() only returns 8-byte alignment, we need 16
  i = posix_memalign(&p, sizeof(__m128i),sizeof(struct v224));
  assert(i == 0);

  vp = (struct v224 *)p;
  p = malloc(len*sizeof(decision_t));
  assert(p != NULL);
  //    free(vp);
  //    return NULL;

  vp->decisions = (decision_t *)p;

  for(state=0;state < (1<<(K-2));state++){
    Branchtab224[0].s[state] = G1FLIP ^ parity((2*state) & POLY1) ? 255 : 0;
    Branchtab224[1].s[state] = G2FLIP ^ parity((2*state) & POLY2) ? 255 : 0;
  }
  init_viterbi224(vp,0);
  return vp;
}

// Viterbi chainback
int chainback_viterbi224(
      void *p,
      unsigned char *data, /* Decoded output data */
      unsigned int nbits, /* Number of data bits */
      unsigned int endstate){ /* Terminal encoder state */
  struct v224 *vp = p;
  decision_t *d = (decision_t *)vp->decisions;

  if(d == NULL)
    return -1;

  endstate &= (1<<(K-1))-1;

#if 1
  {
  unsigned char dbyte = 0;
  while(nbits-- > 0){
    int bit;

    // Accumulate decoded data bits as they fall off the right end of endstate
    dbyte = ((endstate & 1) << 7) | (dbyte >> 1);
    if((nbits & 7) == 0)
      data[nbits>>3] = dbyte;
    bit = (d[nbits].w[endstate>>5] >> (endstate & 31)) & 1; // these constants do NOT change with K
    endstate = (bit << (K-2)) | (endstate >> 1);
  }
  }
#else
  // alternative that accumulates the data byte in the low 8 bits of endstate
  // This requires endstate to be at least K+8 bits wide
  endstate <<= 8;
  while(nbits-- > 0){
    int bit;

    bit = (d[nbits].w[endstate>>13] >> ((endstate >> 8)& 31)) & 1; // Constants do NOT change with K
    endstate = (bit << (K+6)) | (endstate >> 1);
    if((nbits & 7) == 0)
      data[nbits>>3] = endstate & 0xff;
  }
#endif
  return 0;

}

// Delete instance of a Viterbi decoder
void delete_viterbi224(void *p){
  struct v224 *vp = p;

  if(vp != NULL){
    free(vp->decisions);
    free(vp);
  }
}


// Process received symbols
int update_viterbi224_blk(void *p,const unsigned char *syms,int nbits){
  struct v224 *vp = p;
  decision_t *d = (decision_t *)vp->dp;

  while(nbits--){
    __m128i sym0v,sym1v;
    void *tmp;
    int i;

    // Splat the 0th symbol across sym0v, the 1st symbol across sym1v, etc
    sym0v = _mm_set1_epi16(syms[0]);
    sym1v = _mm_set1_epi16(syms[1]);
    syms += 2;

    // SSE2 doesn't support minimum of unsigned words but does minimum of signed words.
    // So we use signed shorts. SSE 4.1 does do minimum of unsigned words (PMINUW)
    // but there  doesn't seem to be any particular advantage to it
    for(i=0; i < 1<<(K-5); i++){
      __m128i decision0,decision1,metric,m_metric,m0,m1,m2,m3,survivor0,survivor1;

      // Form branch metrics
      // Because Branchtab takes on values 0 and 255, and the values of sym?v are offset binary in the range 0-255,
      // the XOR operations constitute conditional negation.
      // metric and m_metric (-metric) are in the range 0-510
      metric = _mm_add_epi16(_mm_xor_si128(Branchtab224[0].v[i],sym0v),_mm_xor_si128(Branchtab224[1].v[i],sym1v));
      m_metric = _mm_sub_epi16(_mm_set1_epi16(510),metric);
    
      // Add branch metrics to path metrics using saturating signed addition
      m0 = _mm_adds_epi16(vp->old_metrics->v[i],metric);
      m3 = _mm_adds_epi16(vp->old_metrics->v[(1<<(K-5))+i],metric);
      m1 = _mm_adds_epi16(vp->old_metrics->v[(1<<(K-5))+i],m_metric);
      m2 = _mm_adds_epi16(vp->old_metrics->v[i],m_metric);
    
#if 0
      _mm_prefetch((void *)&Branchtab224[0].v[i+1],3);
      _mm_prefetch((void *)&Branchtab224[1].v[i+1],3);
      _mm_prefetch((void *)&vp->old_metrics->v[i+1],0);
      _mm_prefetch((void *)&vp->old_metrics->v[(1<<(K-5))+i],0);
#endif

      // Determine winners in one of two ways:
      // Direct comparison ought to be faster by not having to wait for the minimum operation,
      // but it doesn't seem to make any difference on the machines I've tried (Core i7-3612QM; Xeon E5620),
      // perhaps because memory bandwidth is the bottleneck anyway.
      // One minor difference is that the first method breaks ties in favor of the 1-branch
      // while the second breaks ties in favor of the 0-branch. That probably doesn't matter in practice.
#if 0
      // Find winning metrics. Note: signed, there is an unsigned version on SSE 4.1
      survivor0 = _mm_min_epi16(m0,m1);
      survivor1 = _mm_min_epi16(m2,m3);
      // Compare survivors with metrics to see who won
      decision0 = _mm_cmpeq_epi16(survivor0,m1);
      decision1 = _mm_cmpeq_epi16(survivor1,m3);
#else
      // Do this before computing the surviving metrics
      decision0 = _mm_cmpgt_epi16(m0,m1);
      decision1 = _mm_cmpgt_epi16(m2,m3);
      // Find winning metrics. Note: signed, there is an unsigned version on SSE 4.1
      survivor0 = _mm_min_epi16(m0,m1);
      survivor1 = _mm_min_epi16(m2,m3);
#endif
 
      // Pack each set of decisions into 8 8-bit bytes, then interleave and compress into 16 bits
      d->s[i] = _mm_movemask_epi8(_mm_unpacklo_epi8(_mm_packs_epi16(decision0,_mm_setzero_si128()),_mm_packs_epi16(decision1,_mm_setzero_si128())));

      // Store surviving metrics
      vp->new_metrics->v[2*i] = _mm_unpacklo_epi16(survivor0,survivor1);
      vp->new_metrics->v[2*i+1] = _mm_unpackhi_epi16(survivor0,survivor1);
    }

    // Eyeball experiments show the metric spread for this code to start at 5000-6000 or so due to the large
    // starting bias for state 0. After 1 constraint length it falls to 1000-2000
#if 0
    {
      int minmetric,maxmetric;

      minmetric = 9999999;
      maxmetric = -9999999;
      for(i=0;i<(1<<(K-1));i++){
	if(minmetric > vp->new_metrics->s[i])
	  minmetric = vp->new_metrics->s[i];
	if(maxmetric < vp->new_metrics->s[i])
	  maxmetric = vp->new_metrics->s[i];
      }
      printf("metric range %d %d spread %d\n",minmetric,maxmetric,maxmetric-minmetric);
    }
#endif
    // See if we need to renormalize
    // The comparison constant is chosen empirically; a higher value causes renormalization
    // to take place less often, but it risks some other metric saturating positive
    // This number should be 32767 minus the maximum observed metric spread (see above) minus a margin
    if(vp->new_metrics->s[0] >= 25000){
      int i,adjust;
      __m128i adjustv;
      union { __m128i v; uint16_t w[8]; } t;
      
      // Find smallest metric and set adjustv to bring it down to SHRT_MIN
      adjustv = vp->new_metrics->v[0];
      for(i=1;i<(1<<(K-4));i++)
	adjustv = _mm_min_epi16(adjustv,vp->new_metrics->v[i]);

      adjustv = _mm_min_epi16(adjustv,_mm_srli_si128(adjustv,8));
      adjustv = _mm_min_epi16(adjustv,_mm_srli_si128(adjustv,4));
      adjustv = _mm_min_epi16(adjustv,_mm_srli_si128(adjustv,2));
      t.v = adjustv;
      adjust = t.w[0] - SHRT_MIN;
      adjustv = _mm_set1_epi16(adjust);

      // We cannot use a saturated subtract, because we often have to adjust by more than SHRT_MAX
      // This is okay since it can't overflow anyway
      for(i=0;i < 1<<(K-4);i++)
	vp->new_metrics->v[i] = _mm_sub_epi16(vp->new_metrics->v[i],adjustv);
#if 0
      printf("Adjust metrics by %d\n",adjust);
#endif
    }
    d++;
    // Swap pointers to old and new metrics
    tmp = vp->old_metrics;
    vp->old_metrics = vp->new_metrics;
    vp->new_metrics = tmp;
  }
  vp->dp = d;
  return 0;
}


