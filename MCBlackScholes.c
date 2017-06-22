/*
  Uses Monte Carlo simulation to price an option in a Black
  Scholes world.

  Author: Adam Diehl

  Date: June 7, 2017

  Compilation notes:
    $ gcc -o BlackScholes MCBlackScholes.c -O3
    $ time ./BlackScholes
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "MCBlackScholes.h"

//Define parameters
double S = 100;
double K = 102;
double r = 0.030485;
double vol = 0.3392;
double T = 0.5;
double dt = 1/252;
int N = 100000000;

//Declare function
double RNORM();

int main() {

//Seed RNG
time_t t;
srand((unsigned) time(&t));

//Calculate constants
double drift;
  drift = (r-vol*vol*0.5)*T;
double psp;
  psp = vol*sqrt(T);
double discount;
  discount = exp(-r*T);
double weight;
  weight = 1/(N + 0.000);
int progressIndicator;
  progressIndicator = floor(0.1*N);

//Calculate payoffs
double CallPayoff[N];
double PutPayoff[N];

double CallPremium = 0;
double PutPremium = 0;

int i;
int tick = -5;
for (i = 0; i < N; i = i+1) {
  double TerminalPrice;
    TerminalPrice = S*exp(drift + psp*RNORM());
  CallPayoff[i] = (TerminalPrice - K);
  PutPayoff[i] = (K - TerminalPrice);

  //Negative positions won't be exercised
  if(CallPayoff[i] < 0){
    CallPayoff[i] = 0;
  }
  if(PutPayoff[i] < 0){
    PutPayoff[i] = 0;
  }

  //Update premium
  CallPremium = CallPremium + CallPayoff[i];
  PutPremium = PutPremium + PutPayoff[i];

  //Marks program process during execution (adds ~0.6 seconds @ 100M realizations)
  if(i%progressIndicator == 0) {
    tick = tick + 10;
    printf("Status at %d percent \n", tick);
  }
}

CallPremium = CallPremium*discount*weight;
PutPremium = PutPremium*discount*weight;

printf("Call price = %f\n", CallPremium);
printf("Put price = %f\n", PutPremium);

return 0;
}
