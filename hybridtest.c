#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>
#include "viterbi224.h"
#include "fano.h"
#include "sim.h"

static inline int popcount(int x){
  return __builtin_popcount(x);
}

struct option Options[] = {
  {"scale",1, NULL, 'S'},
  {"delta",1, NULL, 'd'},
  {"max-cycles", 1, NULL, 'm'},
  {"frame-length",1,NULL,'l'},
  {"frame-count",1,NULL,'n'},
  {"ebn0",1,NULL,'e'},
  {"gain",1,NULL,'g'},
  {"verbose",0,NULL,'v'},
  {NULL},
};


int Nbits = 1024;                // Bits per frame, INCLUDING TAIL that must be 0's
double Signal = 30;              // Signal amplitude (scaled to 8 bit byte)
int Scale = 8;                   // Fano metric scaling factor
const double Rate = 0.5;         // Code rate = 1/2
double Ebn0 = 2.0;               // Digital signal-to-noise ratio per bit in dB
unsigned long Maxcycles = 1000;  // Maximum number of decoder moves per bit
int Trials = 1000;               // Number of frames to test
int Verbose;                     // diag diarrhea
int Zerodata;                    // Use all 0's data (no effect on linear codes)

int main(int argc,char *argv[]){
  int mettab[2][256];
  unsigned char data[Nbits/8],symbols[Nbits*2],decode_data[Nbits/8];
  unsigned char xordata[Nbits/8];
  unsigned long cycles, metric;
  int delta = 4;
  time_t t;
  int trial,i,r;
  int errcnt;
  int fano_good_frames = 0;
  int fano_failures = 0;
  int fano_frame_errors = 0;
  int fano_bit_errors = 0;
  int viterbi_attempts = 0;
  int viterbi_good_frames = 0;
  int viterbi_frame_errors = 0;
  int viterbi_bit_errors = 0;
  long totcycles = 0;
  long histogram[256];
  double noise_amp; // Actual noise amplitude, computed from Signal amplitude & Eb/N0
  void *vp;

  while((i = getopt_long(argc,argv,"d:S:l:n:e:s:m:vz",Options,NULL)) != EOF){
    switch(i){
    case 'd':
      delta = atoi(optarg);
      break;
    case 'S':
      Scale = atoi(optarg);
      break;
    case 'm':
      Maxcycles = atoi(optarg);
      break;
    case 'l':
      Nbits = atoi(optarg);
      break;
    case 'n':
      Trials = atoi(optarg);
      break;
    case 'e':
      Ebn0 = atof(optarg);
      break;
    case 's':
      Signal = atof(optarg);
      break;
    case 'v':
      Verbose++;
      break;
    case 'z':
      Zerodata++;
      break;
    default:
      printf("Usage: %s [-m maxcycles/bit] [-l bits/frame] [-n numframes] [-e Eb/No] [-s signal_amplitude] [-v] [-z]\n",
	     argv[0]);
      exit(1);
      break;
    }
  }
  if(Nbits < 64){
    printf("bits/frame must be >= 64\n");
    exit(1);
  }
  delta *= Scale;

  // Compute noise voltage. The factor of 2 accounts for BPSK seeing
  // only half the noise power, and the sqrt() converts power to voltage
  noise_amp = Signal / sqrt(2*Rate*pow(10.,Ebn0/10.));
  gen_met(mettab,Signal,noise_amp,Rate,Scale);

  // Generate channel transition matrix
  setup_channel(Signal,noise_amp);

  printf("Code rate %.2f, Nbits = %d, Maxcycles/bit %ld\n",Rate,Nbits,Maxcycles);
  printf("Eb/N0 = %.3lf dB, Signal = %lg, Noise = %lg, BER@Eb/N0 = %lg, BER@Es/N0 = %lg\n",
	 Ebn0,Signal,noise_amp,0.5*erfc(pow(10.,Ebn0/20.)),0.5*erfc(sqrt(Rate*pow(10.,Ebn0/10.))));
  
  srandom(time(&t));

  memset(histogram,0,sizeof(histogram));
  memset(data,0,sizeof(data));
  for(trial = 0; trial < Trials; trial++){

    if(!Zerodata){
      // Generate random data
      // Note last 3 bytes must be 0 to tail off the encoder
      for(i=0;i<(Nbits-64)/8;i++)      // allow room on end for max length tail
	data[i] = random() & 0xff;

      for(;i<Nbits/8;i++)
	data[i] = 0;
    }
    i =  encode(symbols,data,sizeof(data));
    assert(i == 0);
    
    // Add noise & scale, build histogram
    for(i=0;i<2*Nbits;i++){
      symbols[i] = simulate(symbols[i]);
      histogram[symbols[i]]++;
    }
    if(Verbose > 2){
      printf("Cumulative symbol histogram:\n");
      for(i=0;i<256;i++){
	printf(" %6ld",histogram[i]);
	if((i % 16) == 15)
	  putchar('\n');
      }
    }
    memset(decode_data,0,sizeof(decode_data));
    r = fano(&metric,&cycles,decode_data,symbols,Nbits,mettab,delta,Maxcycles);
    totcycles += cycles;
    if(r != 0){
      ++fano_failures;
      if(Verbose)
	printf("trial %d fano: decode failure\n",trial);
    } else { // Fano finished, look for undetected errors
      errcnt = 0;
      for(i=0;i<Nbits/8;i++){
	int e;

	e = popcount(xordata[i] = decode_data[i] ^ data[i]);
	fano_bit_errors += e;
	errcnt += e;
      }      
      if(errcnt != 0){
	fano_frame_errors++;
	fano_bit_errors += errcnt;
	if(Verbose)
	  printf("trial %d fano: metric %ld, cycles %ld, bit errors %d\n",
		 trial,metric,cycles,errcnt);

	if(Verbose > 1){
	  for(i=0;i<Nbits/8;i++)
	    printf("%02x",xordata[i]);
	  putchar('\n');
	}
      } else {
	++fano_good_frames;
	if(Verbose > 1)
	  printf("trial %d fano: metric = %ld, cycles %ld; data OK\n",trial,metric,cycles);
	continue; // Fano OK, next frame
      }
    }
    // Fano failed. Try Viterbi
    viterbi_attempts++;
    if(Verbose){
      printf("Trying Viterbi...");
      fflush(stdout);
    }
    if((vp = create_viterbi224(Nbits)) == NULL){
      printf("create_viterbi224 failed\n");
      exit(1);
    }
    init_viterbi224(vp,0);    // Initialize Viterbi decoder
    update_viterbi224_blk(vp,symbols,Nbits);    // Process symbols
    chainback_viterbi224(vp,decode_data,Nbits,0);      // Viterbi chainback
    delete_viterbi224(vp);

    // look for Viterbi errors
    errcnt = 0;
    for(i=0;i<Nbits/8;i++){
      int e = popcount(xordata[i] = data[i] ^ decode_data[i]);
      errcnt += e;
      viterbi_bit_errors += e;
    }
    if(!errcnt){
      viterbi_good_frames++;
      if(Verbose)
	  printf(" Success\n");
      continue;
    }
    // Viterbi made errors
    if(Verbose)
      printf(" %d bit errors\n",errcnt);
    viterbi_frame_errors++;
    viterbi_bit_errors += errcnt;
    if(Verbose > 1){
      for(i=0;i<Nbits/8;i++)
	printf("%02x",xordata[i]);
      putchar('\n');
    }
  }
  printf("Fano good frames: %d, decode failures %d, frame errors %d, bit errors %d cycles/bit %lf\n",
	 fano_good_frames,fano_failures,fano_frame_errors,fano_bit_errors,(double)totcycles/(trial*Nbits));
  if(viterbi_attempts != 0){
    printf("Viterbi attempts %d good frames: %d frame errors %d (%lg%%) bit errors %d (%lg%%)\n",
	   viterbi_attempts,
	   viterbi_good_frames,viterbi_frame_errors,100.*viterbi_frame_errors/viterbi_attempts,
	   viterbi_bit_errors,100.*viterbi_bit_errors/(Nbits*viterbi_attempts));
  }
  exit(0);
}
