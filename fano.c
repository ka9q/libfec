// Soft decision Fano sequential decoder for r=1/2 convolutional codes
// Copyright 1994, Phil Karn, KA9Q
// Updated March 2014 (!!) for r=1/2 k=24 ICE code


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "fano.h"
#include "code.h"

struct node {
  unsigned long long encstate;	// Encoder state of next node
  long gamma;		        // Cumulative metric to this node
  int metrics[4];		// Metrics indexed by all possible tx syms
  int tm[2];		        // Sorted metrics for current hypotheses
  int i;			// Current branch being tested
};

static inline int parity(unsigned long long x){
  return __builtin_parityll(x);
}

// Given an encoder state, return a rate 1/2 symbol pair.
// The POLY1 symbol goes into the next-to-LSB
// of the result and the POLY2 symbol goes into the LSB.
static inline int makesyms(unsigned long long state){
  int result;

  result = (parity(state & POLY1) << 1) ^ G1FLIP;
  result |= parity(state & POLY2) ^ G2FLIP;
  return result;
}

// Decode packet with the Fano algorithm.
// Return 0 on success, -1 on timeout
int fano(
unsigned long *metric,	// Final path metric (returned value) 
unsigned long *cycles,	// Cycle count (returned value) 
unsigned char *data,	// Decoded output data 
const unsigned char *symbols,	// Raw deinterleaved input symbols 
unsigned int nbits,	// Number of output bits, including tail
int mettab[2][256],	// Metric table, [sent sym][rx symbol] 
int delta,		// Threshold adjust parameter 
unsigned long maxcycles)// Decoding timeout in cycles per bit 
{
  struct node *nodes;		// First node 
  register struct node *np;	// Current node 
  struct node *lastnode;	// Last node 
  struct node *tail;		// First node of tail 
  long t;			// Threshold 
  long m0,m1;
  long ngamma;
  unsigned int lsym;
  unsigned long i;
  
  if((nodes = (struct node *)malloc(nbits*sizeof(struct node))) == NULL){
    fprintf(stderr,"alloc failed\n");
    return 0;
  }
  lastnode = &nodes[nbits];
  tail = &nodes[nbits-(K-1)];
  
  // Compute all possible branch metrics for each symbol pair
  // This is the only place we actually look at the raw input symbols
  for(np=nodes;np < lastnode;np++){
    np->metrics[0] = mettab[0][symbols[0]] + mettab[0][symbols[1]];
    np->metrics[1] = mettab[0][symbols[0]] + mettab[1][symbols[1]];
    np->metrics[2] = mettab[1][symbols[0]] + mettab[0][symbols[1]];
    np->metrics[3] = mettab[1][symbols[0]] + mettab[1][symbols[1]];
#if 0
    printf("k=%ld metrics %d %d %d %d\n",np-nodes,
	   np->metrics[0],np->metrics[1],np->metrics[2],np->metrics[3]);
#endif
    symbols += 2;
  }
  np = nodes;
  np->encstate = 0;
  
  // Compute and sort branch metrics from root node 
  lsym = makesyms(np->encstate);	// 0-branch (LSB is 0)

  m0 = np->metrics[lsym];
  
  // Now do the 1-branch. To save another makesyms call here and
  // inside the loop, we assume that both polynomials are odd,
  // i.e., the least significant bits are 1, providing complementary pairs of branch symbols.
  
  // This code could be sped up if a systematic code were used.
  m1 = np->metrics[3^lsym];
  if(m0 > m1){
    // 0-branch has better metric 
    np->tm[0] = m0;
    np->tm[1] = m1;
  } else {
    // 1-branch is better 
    np->tm[0] = m1;
    np->tm[1] = m0;
    np->encstate |= 1;	// Set low bit 
  }
  np->i = 0;	// Start with best branch 
  maxcycles *= nbits;
  np->gamma = t = 0;
  
  // Start the Fano decoder 
  for(i=1;i <= maxcycles;i++){

    //#define debug 1
#ifdef	debug
    fprintf(stdout,"k=%d, encoder 0x%06llx, metric=%ld, thresh=%ld, m[%d]=%d\n",
	    (int)(np-nodes),np->encstate & ((1LL<<K)-1),np->gamma,t,np->i,np->tm[np->i]);
#endif
    // Look forward 
    ngamma = np->gamma + np->tm[np->i];
    //    printf("np->gamma = %ld, ngamma = %ld\n",np->gamma,ngamma);
    if(ngamma >= t){
      // Node is acceptable 
      if(np->gamma < t + delta){
	// First time we've visited this node; tighten threshold.
	
	// This loop could be replaced with
	//   t += delta * ((ngamma - t)/delta);
	// but the multiply and divide are slower.
	while(ngamma >= t + delta)
	  t += delta;
      }
      // Move forward 
      if(++np == lastnode){
	np--;
	break;	// Done! 
      }
      np->gamma = ngamma;
      np->encstate = np[-1].encstate << 1;
      
      // Compute and sort metrics, starting with the zero branch
      lsym = makesyms(np->encstate);
      if(np >= tail){
	// The tail must be all zeroes, so don't even
	// bother computing the 1-branches there.
	np->tm[0] = np->metrics[lsym];
      } else {
	m0 = np->metrics[lsym];
	m1 = np->metrics[3^lsym];
#if 0
	printf("m0 = %ld, m1 = %ld\n",m0,m1);
#endif
	if(m0 > m1){
	  // 0-branch is better 
	  np->tm[0] = m0;
	  np->tm[1] = m1;
	} else {
	  // 1-branch is better 
	  np->tm[0] = m1;
	  np->tm[1] = m0;
	  np->encstate++;	// Set low bit 
	}
      }
      np->i = 0;	// Start with best branch 
      continue;
    }
    // Threshold violated, can't go forward 
    for(;;){
      // Look backward 
      if(np == nodes || np[-1].gamma < t){
	// Can't back up either.
	// Relax threshold and and look forward again to better branch.
	t -= delta;
	if(np->i != 0){
	  np->i = 0;
	  np->encstate ^= 1;
	}
	break;
      }
      // Back up 
      if(--np < tail && np->i != 1){
	// Search next best branch 
	np->i++;
	np->encstate ^= 1;
	break;
      } // else keep looking back 
    }
  }
  *metric =  np->gamma;	// Return final path metric 
  
  // Copy decoded data to user's buffer 
  nbits = nbits/8;	// Copy tail, which should be 0's
  np = &nodes[7]; // Start with first full byte
  while(nbits-- != 0){
    *data++ = np->encstate;
    np += 8;
  }
  
  free(nodes);
  *cycles = i;
  if(i > maxcycles)
    return -1;	// Decoder timed out 
  return 0;	// Successful completion 
}
