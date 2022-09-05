#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void setup_channel(double signal, double noise);
unsigned char simulate(int data);

double signal = 40;
double esn0 = 3;     // decibels


int main(){
  double signal = 40.;
  double noise;
  int i;

  noise = signal/sqrt(0.5 * pow(10.,esn0/10.));

  setup_channel(signal,noise);

  printf("simulated 0:\n");
  for(i=0;i<1000;i++){
    printf(" %d",simulate(0));
  }
  putchar('\n');

  printf("simulated 1:\n");
  for(i=0;i<1000;i++){
    printf(" %d",simulate(1));
  }
  putchar('\n');




  exit(0);
}
