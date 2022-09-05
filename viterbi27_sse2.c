// K=7 r=1/2 Viterbi decoder for SSE2
// Feb 2004, Phil Karn, KA9Q
// Updated April 2014

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <xmmintrin.h>
#include <stdint.h>
#include <assert.h>
#include "fec.h"

static inline int parity(int x){
  return __builtin_parity(x);
}

int encode_27(
   unsigned char *symbols,	// Output buffer, 2*8*nbytes
   const unsigned char *data,	// Input buffer, nbytes
   unsigned int nbytes)		// Number of bytes in data
{
  unsigned int encstate;
  int i;
  
  encstate = 0;
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

typedef union { uint8_t c[64]; __m128i v[4]; } metric_t;
typedef union { uint64_t ll; uint32_t w[2]; uint8_t c[8]; uint16_t s[4]; __m64 v[1];} decision_t;
union branchtab27 { uint8_t c[32]; __m128i v[2];} Branchtab27_sse2[2];

// State info for instance of Viterbi decoder
// Don't change this without also changing references in sse2bfly27.s!
struct v27 {
  metric_t metrics;
  decision_t *dp;                     // Pointer to current decision
  decision_t *decisions;              // Beginning of decisions for block
};

// Initialize Viterbi decoder for start of new frame
int init_viterbi27(void *p,int starting_state){
  struct v27 *vp = p;
  int i;

  if(p == NULL)
    return -1;
  for(i=0;i<64;i++)
    vp->metrics.c[i] = 63;

  vp->dp = vp->decisions;
  vp->metrics.c[starting_state & 63] = 0; // Bias known start state
  return 0;
}

// Create a new instance of a Viterbi decoder
void *create_viterbi27(int len){
  void *p;
  struct v27 *vp;
  int state;

  for(state=0;state < 32;state++){
    Branchtab27_sse2[0].c[state] = parity((2*state) & V27POLYB) ? 255 : 0;
    Branchtab27_sse2[1].c[state] = !parity((2*state) & V27POLYA) ? 255 : 0;
  }

  // Ordinary malloc() only returns 8-byte alignment, we need 16
  if(posix_memalign(&p, sizeof(__m128i),sizeof(struct v27)))
    return NULL;
  vp = (struct v27 *)p;

  if((p = malloc(len*sizeof(decision_t))) == NULL){
    free(vp);
    return NULL;
  }
  vp->decisions = (decision_t *)p;
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
  //  unsigned char dbyte = 0;

  if(p == NULL)
    return -1;

  d = vp->decisions;
  //  endstate &= 63;
  endstate = (endstate & 63) << 8;
  while(nbits-- != 0){
    int bit;

    // Accumulate decoded data bits as they fall off the right end of endstate
    //    dbyte = ((endstate & 1) << 7) | (dbyte >> 1);
    //if((nbits % 8) == 0)
      //data[nbits/8] = dbyte;

    // Extract decision and push onto left (tail) end of encoder
    // We only need one bit, but reading 64 is fastest
    bit = (d[nbits].ll >> (endstate >> 8)) & 1;
    endstate = (bit << 13) | (endstate >> 1);
    if((nbits % 8) == 0)
      data[nbits/8] = endstate;

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


int update_viterbi27_blk(void *p,unsigned char *syms,int nbits){
  struct v27 *vp = p;
  decision_t *d;
  metric_t working_metrics;

  if(p == NULL)
    return 0;
  d = (decision_t *)vp->dp;

  working_metrics = vp->metrics;

  while(nbits--){
    __m128i sym0v,sym1v;
    
    // Splat the 0th symbol across sym0v, the 1st symbol across sym1v, etc
    sym0v = _mm_set1_epi8(syms[0]);
    sym1v = _mm_set1_epi8(syms[1]);
    syms += 2;

    {
      // Unwound version
      __m128i decision0,decision1,metric,m_metric,m0,m1,m2,m3,survivor0,survivor1;
      __m128i new_survivor_1;

      // First  metrics 0 & 2 are read, 0 and 1 are written
      // Form branch metrics
      metric = _mm_avg_epu8(_mm_xor_si128(Branchtab27_sse2[0].v[0],sym0v),_mm_xor_si128(Branchtab27_sse2[1].v[0],sym1v));
      // There's no packed bytes right shift in SSE2, so we use the word version and mask
      // (I'm *really* starting to like Altivec...)
      metric = _mm_srli_epi16(metric,3);
      metric = _mm_and_si128(metric,_mm_set1_epi8(31));
      m_metric = _mm_sub_epi8(_mm_set1_epi8(31),metric);
    
      // Add branch metrics to path metrics
      m0 = _mm_add_epi8(working_metrics.v[0],metric);
      m3 = _mm_add_epi8(working_metrics.v[2],metric);
      m1 = _mm_add_epi8(working_metrics.v[2],m_metric);
      m2 = _mm_add_epi8(working_metrics.v[0],m_metric);
    
      // Compare and select, using modulo arithmetic
      decision0 = _mm_cmpgt_epi8(_mm_sub_epi8(m0,m1),_mm_setzero_si128());
      decision1 = _mm_cmpgt_epi8(_mm_sub_epi8(m2,m3),_mm_setzero_si128());

      survivor0 = _mm_or_si128(_mm_and_si128(decision0,m1),_mm_andnot_si128(decision0,m0));
      survivor1 = _mm_or_si128(_mm_and_si128(decision1,m3),_mm_andnot_si128(decision1,m2));
 
      // Pack each set of decisions into 16 bits
      d->s[0] = _mm_movemask_epi8(_mm_unpacklo_epi8(decision0,decision1));
      d->s[1] = _mm_movemask_epi8(_mm_unpackhi_epi8(decision0,decision1));

      // Store surviving metrics
      working_metrics.v[0] = _mm_unpacklo_epi8(survivor0,survivor1);  // metrics 0 can get written right away
      new_survivor_1 = _mm_unpackhi_epi8(survivor0,survivor1);      // Can't write 1 yet, we need it next

      // Second
      // Second  metrics 1 & 3 are read, 2 & 3 are written
      // Form branch metrics
      metric = _mm_avg_epu8(_mm_xor_si128(Branchtab27_sse2[0].v[1],sym0v),_mm_xor_si128(Branchtab27_sse2[1].v[1],sym1v));
      // There's no packed bytes right shift in SSE2, so we use the word version and mask
      // (I'm *really* starting to like Altivec...)
      metric = _mm_srli_epi16(metric,3);
      metric = _mm_and_si128(metric,_mm_set1_epi8(31));
      m_metric = _mm_sub_epi8(_mm_set1_epi8(31),metric);
    
      // Add branch metrics to path metrics
      m0 = _mm_add_epi8(working_metrics.v[1],metric);
      m3 = _mm_add_epi8(working_metrics.v[3],metric);
      m1 = _mm_add_epi8(working_metrics.v[3],m_metric);
      m2 = _mm_add_epi8(working_metrics.v[1],m_metric);
    
      // Now safe to store surviving metrics 1 from first half
      working_metrics.v[1] = new_survivor_1;

      // Compare and select, using modulo arithmetic
      decision0 = _mm_cmpgt_epi8(_mm_sub_epi8(m0,m1),_mm_setzero_si128());
      decision1 = _mm_cmpgt_epi8(_mm_sub_epi8(m2,m3),_mm_setzero_si128());
      survivor0 = _mm_or_si128(_mm_and_si128(decision0,m1),_mm_andnot_si128(decision0,m0));
      survivor1 = _mm_or_si128(_mm_and_si128(decision1,m3),_mm_andnot_si128(decision1,m2));
 
      // Pack each set of decisions into 16 bits
      d->s[2] = _mm_movemask_epi8(_mm_unpacklo_epi8(decision0,decision1));
      d->s[3] = _mm_movemask_epi8(_mm_unpackhi_epi8(decision0,decision1));

      // Store surviving metrics
      working_metrics.v[2] = _mm_unpacklo_epi8(survivor0,survivor1);
      working_metrics.v[3] = _mm_unpackhi_epi8(survivor0,survivor1);
    }

    d++;
  }
  vp->metrics = working_metrics;
  vp->dp = d;
  return 0;
}
