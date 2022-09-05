#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>
#include "fano.h"
#include "sim.h"



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


int Nbits = 1024;                // Bits per frame, including tail that must be 0
double Signal = 30;              // Signal amplitude (scaled to 8 bit byte)
int Scale = 8;             // Fano metric scaling factor
const double Rate = 0.5;         // Code rate = 1/2
double Ebn0 = 2.0;               // Digital signal-to-noise ratio per bit in dB
unsigned long Maxcycles = 1000;  // Maximum number of decoder moves per bit
int Trials = 1000;               // Number of frames to test
int Verbose;                     // diag diarrhea
int Zerodata;                    // Use all 0's data (no effect on linear codes)

int main(int argc,char *argv[]){
  int mettab[2][256];
  unsigned char data[Nbits/8],symbols[Nbits*2],decode_data[Nbits/8];
  unsigned long cycles, metric;
  int delta = 4;
  time_t t;
  int trial,i,r,good=0,bad=0,undetected=0;
  long totcycles = 0,histogram[256];
  double noise_amp; // Actual noise amplitude, computed from Signal amplitude & Eb/N0

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
    
#if 0
    printf("raw data and symbols, no noise:\n");
    for(i=0;i<Nbits;i++){
      if((i % 8) == 0)
	printf("data[%d] = %02x; symbols = ",i/8,data[i/8]);
      printf("%d%d",symbols[2*i],symbols[2*i+1]);
      if((i % 8) == 7)
	putchar('\n');
    }
    putchar('\n');
#endif
    
    // Add noise & scale, build histogram
    for(i=0;i<2*Nbits;i++){
      symbols[i] = simulate(symbols[i]);
      histogram[symbols[i]]++;
    }
#if 0
    printf("data and symbols, noisy & scaled:\n");
    for(i=0;i<Nbits;i++){
      if((i % 8) == 0)
	printf("data[%d] = %02x; symbols = ",i/8,data[i/8]);
      printf(" %03d %03d",symbols[2*i],symbols[2*i+1]);
      if((i % 8) == 7)
	putchar('\n');
    }
    putchar('\n');
#endif
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
    i = memcmp(data,decode_data,sizeof(data));
    bad += (i != 0);
    undetected += (r == 0 && i != 0);
    good += (i == 0);

    if(Verbose > 1 || (Verbose && r != 0)){
      printf("trial %d fano returns %d, metric = %ld, cycles = %ld",trial,r,metric,cycles);
      if(i != 0 && (Verbose > 1 || r == 0)){
	// Error in data
	putchar(' ');
	for(i=0;i<Nbits/8;i++)
	  printf("%02x",decode_data[i] ^ data[i]);
      }
      putchar('\n');
    }
  }
  printf("trials %d avg cycles/bit %lg good %d bad %d undetected %d deletion rate %lg%%\n",
	 trial,(double)totcycles/(trial*Nbits),good,bad,undetected,100.*bad/trial);


  exit(0);
}
