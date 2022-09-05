/* Test viterbi decoder speeds */
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <memory.h>
#include <sys/time.h>
#include <sys/resource.h>
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif
#include "viterbi224.h"
#include "code.h"
#include "sim.h"

#if HAVE_GETOPT_LONG
struct option Options[] = {
  {"frame-length",1,NULL,'l'},
  {"frame-count",1,NULL,'n'},
  {"ebn0",1,NULL,'e'},
  {"gain",1,NULL,'g'},
  {"verbose",0,NULL,'v'},
  {NULL},
};
#endif

#define RATE (0.5)
#define MAXBYTES 10000


double Gain = 24.0;
int Verbose = 0;

static inline int popcount(int x){
  return __builtin_popcount(x);
}


int main(int argc,char *argv[]){
  int i,d,tr;
  int trials = 10,errcnt,framebits=1024;
  int tot_errs=0;
  unsigned char data[MAXBYTES];
  unsigned char decoded_data[MAXBYTES];
  unsigned char xordata[MAXBYTES];
  unsigned char symbols[8*2*MAXBYTES];
  void *vp;
  extern char *optarg;
  struct rusage start,finish;
  double extime;
  double noise,esn0,ebn0;
  time_t t;
  int badframes=0;

  time(&t);
  srandom(t);
  ebn0 = -100;
#if HAVE_GETOPT_LONG
  while((d = getopt_long(argc,argv,"l:n:e:g:v",Options,NULL)) != EOF){
#else
  while((d = getopt(argc,argv,"l:n:e:g:v")) != EOF){
#endif
    switch(d){
    case 'l':
      framebits = atoi(optarg);
      break;
    case 'n':
      trials = atoi(optarg);
      break;
    case 'e':
      ebn0 = atof(optarg);
      break;
    case 'g':
      Gain = atof(optarg);
      break;
    case 'v':
      Verbose++;
      break;
    }
  }
  if(framebits > 8*MAXBYTES){
    fprintf(stderr,"Frame limited to %d bits\n",MAXBYTES*8);
    framebits = MAXBYTES*8;
  }
  if((vp = create_viterbi224(framebits)) == NULL){
    printf("create_viterbi224 failed\n");
    exit(1);
  }
  if(ebn0 != -100){
    esn0 = ebn0 + 10*log10((double)RATE); /* Es/No in dB */
    // Compute noise voltage; the sqrt(1/2) term is for BPSK
    noise = Gain * M_SQRT1_2 / pow(10.,0.05 * esn0);
    setup_channel(Gain,noise);
    
    printf("nframes = %d framesize = %d ebn0 = %.2f dB Gain = %g noise = %g\n",trials,framebits,ebn0,Gain,noise);

    for(tr=0;tr<trials;tr++){
      // Encode a frame of random data
      for(i=0;i<(framebits-K)/8;i++)
	data[i] = random() & 0xff;

      for(;i<framebits/8;i++) // leave a tail of 0's
	data[i] = 0;

      encode(symbols,data,framebits/8);

      // Add noise & scale, build histogram
      for(i=0;i<2*framebits;i++)
	symbols[i] = simulate(symbols[i]);

      // Decode it, measuring the time
      getrusage(RUSAGE_SELF,&start);
      init_viterbi224(vp,0);
      update_viterbi224_blk(vp,symbols,framebits);
      chainback_viterbi224(vp,decoded_data,framebits,0);
      getrusage(RUSAGE_SELF,&finish);
      extime = finish.ru_utime.tv_sec - start.ru_utime.tv_sec + 1e-6*(finish.ru_utime.tv_usec - start.ru_utime.tv_usec);

      // Count up errors
      errcnt = 0;
      for(i=0;i<framebits/8;i++){
	int e = popcount(xordata[i] = decoded_data[i] ^ data[i]);
	errcnt += e;
	tot_errs += e;
      }
      if(errcnt != 0)
	badframes++;
      if(Verbose > 1 && errcnt != 0){
	printf("frame %d, %d errors: ",tr,errcnt);
	for(i=0;i<framebits/8;i++){
	  printf("%02x",xordata[i]);
	}
	printf("\n");
	printf("original data: ");
	for(i=0;i<framebits/8;i++){
	  printf("%02x",data[i]);
	}
	printf("\n");
	printf("decoded data: ");
	for(i=0;i<framebits/8;i++){
	  printf("%02x",decoded_data[i]);
	}
	printf("\n");
      }
      if(Verbose)
	printf("BER %d/%d (%10.3g) FER %d/%d (%10.3g) time %.6g s (%.2lf b/s)\r",
	       tot_errs,framebits*(tr+1),tot_errs/((double)framebits*(tr+1)),
	       badframes,(tr+1),(double)badframes/(tr+1),extime,framebits/extime);
      fflush(stdout);
    }

    if(Verbose > 1)
      printf("nframes = %d framesize = %d ebn0 = %.2f dB Gain = %g\n",trials,framebits,ebn0,Gain);
    else if(Verbose == 0)
	printf("BER %d/%d (%.3g) FER %d/%d (%.3g)\n",
	       tot_errs,framebits*trials,tot_errs/((double)framebits*trials),
	       badframes,trials,(double)badframes/trials);
    else
      printf("\n");
  } else {
    /* Do time trials */
    memset(symbols,128,sizeof(symbols));
    printf("Starting time trials\n");
    getrusage(RUSAGE_SELF,&start);
    for(tr=0;tr < trials;tr++){
      /* Initialize Viterbi decoder */
      init_viterbi224(vp,0);

      /* Decode block */
      update_viterbi224_blk(vp,symbols,framebits);

      /* Do Viterbi chainback */
      chainback_viterbi224(vp,decoded_data,framebits,0);
    }
    getrusage(RUSAGE_SELF,&finish);
    extime = finish.ru_utime.tv_sec - start.ru_utime.tv_sec + 1e-6*(finish.ru_utime.tv_usec - start.ru_utime.tv_usec);
    printf("Execution time for %d %d-bit frames: %.2f sec\n",trials,
	   framebits,extime);
    printf("decoder speed: %g bits/s\n",trials*framebits/extime);
  }
  delete_viterbi224(vp);
  exit(0);
}
