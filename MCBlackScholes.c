/*
  Uses Monte Carlo simulation to price an option in a Black
  Scholes world.
  
  Source file is dependent on the identically named header. Parameters are hardcoded out of laziness. 

  Author: Adam Diehl

  Date: June 7, 2017
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "MCBlackScholes.h"

//Define parameters
double S = 100;
double K = 102;
double r = 0.030485;
double vol = 0.3392;
double T = 0.5;
double dt = 1/252;
int N = 1000000;

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
for (i = 0; i < N; i = i+1) {
  double TerminalPrice;
    TerminalPrice = S*exp(drift + psp*RNORM());
  CallPayoff[i] = (TerminalPrice - K);
  PutPayoff[i] = (K - TerminalPrice);

  if(CallPayoff[i] < 0){
    CallPayoff[i] = 0;
  }
  if(PutPayoff[i] < 0){
    PutPayoff[i] = 0;
  }

  CallPremium = CallPremium + CallPayoff[i];
  PutPremium = PutPremium + PutPayoff[i];
}

CallPremium = CallPremium*discount*weight;
PutPremium = PutPremium*discount*weight;

printf("Call price = %f\n", CallPremium);
printf("Put price = %f\n", PutPremium);

return 0;
}
