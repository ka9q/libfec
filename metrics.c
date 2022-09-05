// Generate metric tables for a soft-decision convolutional decoder
// assuming gaussian noise on a PSK channel.

// Works from "first principles" by evaluating the normal probability
// function and then computing the log-likelihood function
// for every possible received symbol value

// Copyright 1995 Phil Karn, KA9Q
// Updated March 2014 (!)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

extern int Verbose;

// Normal function integrated from -Inf to x. Range: 0-1
static double normal(double x){
  return 0.5 + 0.5*erf(x/M_SQRT2);
}

// Generate log-likelihood metrics for 8-bit soft quantized channel assuming AWGN and BPSK
void gen_met(
     int mettab[2][256], // Metric table, [sent sym][rx symbol]
     double signal,	 // Signal amplitude, units
     double noise,	 // Noise amplitude, units (absolute, no longer relative!)
     double bias,	 // Metric bias; 0 for viterbi, rate for sequential
     double scale		 // Metric scale factor */
){
  int s;
  double metrics0,metrics1,p0,p1;
  double left0,left1,right0,right1;
  double inv_noise;

  inv_noise = 1./noise;

  // Compute the channel transition probabilities, i.e., the probability of receiving each of the
  // 256 possible values when 0s and 1s were sent.
  // The bins are assumed to be centered on their nominal values:
  // Bin 0; -infinity < v < -127.5
  // Bin 1: -127.5 < v < -126.5
  // Bin 128: -0.5 < v < +0.5
  // Bin 255: +126.5 < v < +infinity
  left0 = left1 = 0.0; // area below bin 0 is zero
  for(s=0;s<256;s++){

    // Find the area below and in this bin, subtract the area below and in the previous bin,
    // leaving just the area of this bin.
    // The area above bin 255 is zero.
    right0 = (s != 255) ? normal((s - 128 + 0.5 + signal) * inv_noise) : 1.0;
    right1 = (s != 255) ? normal((s - 128 + 0.5 - signal) * inv_noise) : 1.0;

    p0 = right0 - left0;    // p0 = P(s|0), prob of receiving s given that a 0 was sent
    p1 = right1 - left1;    // p1 = P(s|1), prob of recieving s given that a 1 was sent

    left1 = right1;
    left0 = right0;

    // Compute log-likelihood ratios assuming even balance of 0's and 1's on channel
    if(p0 == p1){
      // At high SNR, extremal sample values may underflow to p0 == p1 == 0, giving
      // infinitely bad metrics for what might actually be very good samples if
      // actually encountered. Not sure what's right here, so I punt and treat both as erasures
      metrics0 = metrics1 = -bias;
    } else {
      // The smallest value from log2() is about -32, so approximate log2(0) as -33.
      // Alternatively I could represent it as -INT_MAX, the worst possible metric, but that seems excessive
      metrics0 = (p0 == 0) ? -33.0 : log2(2*p0/(p1+p0)) - bias;
      metrics1 = (p1 == 0) ? -33.0 : log2(2*p1/(p1+p0)) - bias;
      // Equivalent:
      metrics0 = (p0 == 0) ? -33.0 : 1 + log2(p0) - log2(p1+p0) - bias;
      metrics1 = (p1 == 0) ? -33.0 : 1 + log2(p1) - log2(p1+p0) - bias;
    }
    // Scale and round for table
    mettab[0][s] = lrint(metrics0 * scale);
    mettab[1][s] = lrint(metrics1 * scale);
    if(Verbose){
      printf("s=%3d P(s|0) = %-12lg P(s|1) = %-12lg P(s) = %-12lg",
	     s,p0,p1,(p1+p0)/2.);
      printf(" metrics0 = %-12lg [%4d]",metrics0,mettab[0][s]);
      printf(" metrics1 = %-12lg [%4d]\n",metrics1,mettab[1][s]);
    }
  }
}
