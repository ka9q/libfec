// User include file for libfec
// Copyright 2004-2014, Phil Karn, KA9Q
// May be used under the terms of the GNU Lesser General Public License (LGPL)

#ifndef _VITERBI_H_
#define _VITERBI_H_

static inline int parity(long long x){
  return __builtin_parityll(x);
}



void set_viterbi224_polynomial(int polys[2]);
int init_viterbi224(void *p,int starting_state);
void *create_viterbi224(int len);
int chainback_viterbi224(void *p,
      unsigned char *data, /* Decoded output data */
      unsigned int nbits, /* Number of data bits */
			 unsigned int endstate); /* Terminal encoder state */
void delete_viterbi224(void *p);
int update_viterbi224_blk(void *p,unsigned char *syms,int nbits);


#endif /* _VITERBI_H_ */



