// User include file for libfec
// Copyright 2004-2014, Phil Karn, KA9Q
// May be used under the terms of the GNU Lesser General Public License (LGPL)

#ifndef _SIM_H_
#define _SIM_H_

// Useful utilities for simulation
double normal_rand(double mean, double std_dev);
unsigned char addnoise(int sym, double signal, double noise);
void setup_channel(double signal,double noise);
unsigned char simulate(int data);

#endif /* _SIM_H_ */



