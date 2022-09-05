#include "code.h"


static inline int parity(unsigned long long x){
  return __builtin_parityll(x);
}


// Convolutionally encode a packet. The input data bytes are read
// high bit first and the encoded packet is written into 'symbols',
// one symbol per byte. The first symbol is generated from POLY1,
// the second from POLY2.

// Storing only one symbol per byte uses more space, but it is faster
// and easier than trying to pack them more compactly.
int encode(
   unsigned char *symbols,	// Output buffer, 2*8*nbytes
   const unsigned char *data,	// Input buffer, nbytes
   unsigned int nbytes)		// Number of bytes in data
{
  unsigned long long encstate;
  int i;
  
  encstate = 0;
  while(nbytes-- != 0){
    for(i=7;i>=0;i--){  // Transmit MSB first
      encstate = ((encstate << 1) | ((*data >> i) & 1));
      *symbols++ = G1FLIP ^ parity(encstate & POLY1);
      *symbols++ = G2FLIP ^ parity(encstate & POLY2);
    }
    data++;
  }
  if((encstate & ((1LL << K) -1)) != 0)
    return -1;  // Warn if encoder wasn't tailed back to 0
  return 0;
}
