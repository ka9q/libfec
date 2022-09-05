#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include "rs.h"

main(){
  char block[255];
  int i;
  void *rs;
  struct tms start,finish;
  double extime;
  int trials = 10000;

  for(i=0;i<223;i++)
    block[i] = 0x01;

  rs = init_rs_char(8,0x187,112,11,32,0);
  encode_rs_char(rs,block,&block[223]);


  times(&start);

  for(i=0;i<trials;i++){
#if 0
    block[0] ^= 0xff; /* Introduce an error */
    block[2] ^= 0xff; /* Introduce an error */
#endif
    decode_rs_char(rs,block,NULL,0);
  }
  times(&finish);
  extime = ((double)(finish.tms_utime-start.tms_utime))/CLK_TCK;
  printf("Execution time for %d Reed-Solomon blocks using general decoder: %.2f sec\n",trials,extime);
  printf("decoder speed: %g bits/s\n",trials*223*8/extime);


  encode_rs_8(block,&block[223],0);
  times(&start);
  for(i=0;i<trials;i++){
#if 0
    block[0] ^= 0xff; /* Introduce an error */
    block[2] ^= 0xff; /* Introduce an error */
#endif
    decode_rs_8(block,NULL,0,0);
  }
  times(&finish);
  extime = ((double)(finish.tms_utime-start.tms_utime))/CLK_TCK;
  printf("Execution time for %d Reed-Solomon blocks using CCSDS decoder: %.2f sec\n",trials,extime);
  printf("decoder speed: %g bits/s\n",trials*223*8/extime);

  exit(0);
}

