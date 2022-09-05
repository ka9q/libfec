// User include file for r=1/2 k=24 Viterbi decoder
// Copyright 2004-2014, Phil Karn, KA9Q
// May be used under the terms of the GNU Lesser General Public License (LGPL)

#ifndef _VITERBI224_H_
#define _VITERBI224_H_

int init_viterbi224(void *p,int starting_state);
void *create_viterbi224(int len);
int chainback_viterbi224(void *p,unsigned char *data,unsigned int nbits,unsigned int endstate);
void delete_viterbi224(void *p);
int update_viterbi224_blk(void *p,const unsigned char *syms,int nbits);


#endif /* _VITERBI224_H_ */



