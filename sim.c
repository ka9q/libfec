#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>


static int channel[2][256]; // Channel cdf lookup table, scaled to 0 - RAND_MAX
static int starts[2];   // Most likely values for 0 and 1 data bits

// Normal function integrated from -Inf to x. Range: 0-1
static double normal(double x){
  return 0.5 + 0.5*erf(x/M_SQRT2);
}

// Set up channel cdfs for memoryless AWGN BPSK with 8-bit receive samples
// The samples are offset binary, with sample 128 containing zero
void setup_channel(double signal,double noise){
  double inv_noise = 1./noise;
  int s;

  for(s=0; s < 256; s++){
    // Compute cdf at right edge of each bin
    channel[0][s] = RAND_MAX * normal((s - 128 + 0.5 + signal) * inv_noise);
    channel[1][s] = RAND_MAX * normal((s - 128 + 0.5 - signal) * inv_noise);
  }
  starts[0] = 128 - signal; // Most likely values to start search
  starts[1] = 128 + signal;
}

// Given a binary xmit symbol, simulate a soft decision receive symbol affected by noise
unsigned char simulate(int data){
  int rand,s,low,high;

  assert(channel[0][255] != 0); // ensure we're initialized
  assert(data == 0 || data == 1);

  // Generate uniform random integer and do binary search in cdf table
  rand = random();
  low = 0;
  high = 255;
  s = starts[data];

  while(high != low){
    if(rand > channel[data][s])
      low = s+1; // Must be in bin above s
    else
      high = s;  // Could be in this bin, or one below it
    s = (low + high)/2;
  }
  return s;
}


// Zignor algorithm for generating normally distributed random variables
// About 2x speed of usual rejection method
inline double DRanU(void){
  return (double)random() / RAND_MAX;
}
inline int IRanU(void){
  return random();
}


static double DRanNormalTail(double dMin, int iNegative) {
  double x, y;
  do {
    x = log(DRanU()) / dMin;
    y = log(DRanU());
  } while (-2 * y < x * x);
  return iNegative ? x - dMin : dMin - x;
}

#define ZIGNOR_C 128 // number of blocks
#define ZIGNOR_R 3.442619855899 // start of the right tail

// (R * phi(R) + Pr(X>=R)) * sqrt(2\pi)
#define ZIGNOR_V 9.91256303526217e-3
// s_adZigX holds coordinates, such that each rectangle has
// same area; s_adZigR holds s_adZigX[i + 1] / s_adZigX[i]

static double s_adZigX[ZIGNOR_C + 1], s_adZigR[ZIGNOR_C];


static void zigNorInit(int iC, double dR, double dV) {
  int i; double f;

  f = exp(-0.5 * dR * dR);
  s_adZigX[0] = dV / f; // [0] is bottom block: V / f(R)
  s_adZigX[1] = dR;
  s_adZigX[iC] = 0;
  for(i = 2; i < iC; ++i) {
      s_adZigX[i] = sqrt(-2 * log(dV / s_adZigX[i - 1] + f));
      f = exp(-0.5 * s_adZigX[i] * s_adZigX[i]);
  }
  for(i = 0; i < iC; ++i)
    s_adZigR[i] = s_adZigX[i + 1] / s_adZigX[i];
}

double DRanNormalZig(void) {
  unsigned int i;
  double x, u, f0, f1;
  static int init = 0;
  
  if(!init){
    zigNorInit(ZIGNOR_C,ZIGNOR_R,ZIGNOR_V);
    init = 1;
  }
  for (;;){
    u = 2 * DRanU() - 1;
    i = IRanU() & 0x7F;
    // first try the rectangular boxes
    if(fabs(u) < s_adZigR[i])
      return u * s_adZigX[i];
    // bottom box: sample from the tail
    if(i == 0)
      return DRanNormalTail(ZIGNOR_R, u < 0);
    // is this a sample from the wedges?
    x = u * s_adZigX[i];
    f0 = exp(-0.5 * (s_adZigX[i] * s_adZigX[i] - x * x) );
    f1 = exp(-0.5 * (s_adZigX[i+1] * s_adZigX[i+1] - x * x) );
    if(f1 + DRanU() * (f0 - f1) < 1.0)
      return x;
  }
}

// Traditional gaussian RV generator using Marsaglia polar method
double normal_rand(double mean, double std_dev){
  double fac,rsq,v1,v2;
  static double gset;
  static int iset;

  if(iset){
    // Already got one
    iset = 0;
    return mean + std_dev*gset;
  }
  // Generate two evenly distributed numbers between -1 and +1 inside the unit circle
  do {
    v1 = 2.0 * (double)random() / RAND_MAX - 1;
    v2 = 2.0 * (double)random() / RAND_MAX - 1;
    rsq = v1*v1 + v2*v2;
  } while(rsq >= 1.0 || rsq == 0.0);
  fac = sqrt(-2.0*log(rsq)/rsq);
  gset = v1*fac;
  iset++;
  return mean + std_dev*v2*fac;
}

// Scale symbol [0,1] to offset-128 and add AWGN 
unsigned char addnoise(int sym,double signal,double noise){
  int sample;
    
  assert((sym & -2) == 0); // sym must be 0 or 1
  sample = normal_rand(128. + signal*(2*sym-1),noise);
  // Clip to 0-255
  sample = sample < 0 ? 0 : sample > 255 ? 255 : sample;
  return sample;
}
