/*
  Uses Monte Carlo simulation to price a vanilla European option using the
  Heston model.

  Author: Adam Diehl

  Date: June 22, 2017

    $ gcc -o Heston Heston.c -O3
    $ time ./Heston
*/

/* ------------------------ Package Dependencies --------------------------- */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "Heston.h"

/* ----------------------- Parameter Declarations -------------------------- */

//Black Scholes parameters
  double S = 100; //Initial stock price
  double K = 100; //Strike price
  double r = 0.0319; //Risk-free rate
  double T = 1; //Time to maturity (in years)
  double IV = 0.010201; //Initial Volatility

//Heston parameters
  double Kappa = 6.21; //Vol reversion
  double Theta = 0.019; //Long term vol
  double Gamma = 0.019; // Vol of vol
  double Rho = -0.7; //Correlation coefficient

//Monte Carlo parameters
  double dt = 100;
  int N = 100000;

/* ------------------- Program Setup and Initializations ------------------- */

//Allocate memory for computable constants
  double CallPayoff = 0;
  double PutPayoff = 0;
  double discount = 0;
  int NCols;

//Allocate loop tickers
  int i;
  int j;

/* Declare functions */
  double RNORM();
  int FellerCondition(double Kappa, double Theta, double Gamma);

/* ----------------------- Program below this line ------------------------- */

int main() {

  //Seed rand
    srand(time(NULL));

  //Find number of necessary columns
    NCols = (int) T*dt;

  //Allocate variable dependent memory arrays
    double *dW1[N];
    double *dW2[N];
      for (i = 0; i < N; i++) {
        dW1[i] = (double *)malloc((NCols - 1) * sizeof(double));
        dW2[i] = (double *)malloc((NCols - 1) * sizeof(double));
      }

    double *PricePath[N];
    double *VolatilityPath[N];
      for(i = 0; i < N; i++) {
        PricePath[i] = (double *)malloc(NCols * sizeof(double));
        VolatilityPath[i] = (double *)malloc(NCols * sizeof(double));
      }

  //Calculate actual dt increments
    double deltaT = T/dt;
    double sqrDeltaT = sqrt(deltaT);

  //Calculate discount
    discount = exp(-r*T);

  //Status update
    printf("Generating correlated RNG...\n");

  //Generate correlated Brownian motions
    double RhoCoefSpec = sqrt(1-Rho*Rho); //Pre compute constant
    for(i = 0; i < N; i = i+1) {
      for(j = 0; j < NCols-1; j = j+1) {
        dW1[i][j] = RNORM();
        dW2[i][j] = dW1[i][j]*Rho + RhoCoefSpec*RNORM();
      }
    }

  //Status update
    printf("RNG complete.\n");

  //Check the Feller Condition
    int VolRunScheme = FellerCondition(Kappa, Theta, Gamma);
    if (VolRunScheme == 0) {
      printf("Feller condition satisfied; volatility is strictly non-negative.\n");
    }
    else {
      printf("Feller condition not satisfied; volatility process may demonstrate artificial truncation and false price convergence.\n");
    }

  //Calculate the asset and volatility paths
    //Set initial price and volatility
      for(i = 0; i < N; i = i+1) {
        PricePath[i][0] = S;
        VolatilityPath[i][0] = IV;
      }

  //Status update
    printf("Generating asset price and volatility paths...\n");

  //Iterate through time with Milstein discretization
    for(i = 0; i < N; i = i+1) {
      for(j = 1; j <= NCols-1; j = j+1) {

        //Path iteration
        VolatilityPath[i][j] = VolatilityPath[i][(j-1)] + Kappa*(Theta - VolatilityPath[i][(j-1)])*deltaT + Gamma*sqrt(VolatilityPath[i][(j-1)])*sqrDeltaT*dW1[i][(j-1)];
        PricePath[i][j] = PricePath[i][(j-1)] * exp((r - 0.5*VolatilityPath[i][(j-1)])*deltaT + sqrt(VolatilityPath[i][(j-1)]*deltaT)*dW2[i][(j-1)]);

        //Volatility truncation
        if(VolatilityPath[i][j] < 0) {
          VolatilityPath[i][j] = 0;
        }

        //Alternative discretization
        //PricePath[i][j] = PricePath[i][(j-1)] + r*PricePath[i][(j-1)]*deltaT + sqrt(VolatilityPath[i][(j-1)]*deltaT)*PricePath[i][(j-1)]*dW2[i][(j-1)] + 0.25*PricePath[i][(j-1)]*PricePath[i][(j-1)]*deltaT*(dW2[i][(j-1)]*dW2[i][(j-1)]-1);
      }
    }

  //Determine the option payoff and discount
    //Call payoff
    for(i = 0; i <= N-1; i = i+1) {
      if((PricePath[i][(NCols-1)] - K) > 0) {
        CallPayoff = CallPayoff + (PricePath[i][(NCols-1)] - K);
      }
    //Put payoff
      if(K - (PricePath[i][(NCols-1)]) > 0) {
        PutPayoff = PutPayoff + (K - PricePath[i][(NCols-1)]);
      }
    }

    CallPayoff = discount*(CallPayoff/N);
    PutPayoff = discount*(PutPayoff/N);

    printf("\n");
    printf("Call price = %f\n", CallPayoff);
    printf("Put price = %f\n", PutPayoff);

    printf("\n");
    printf("Terminating function...\n");

/* ---------------------- Memory Release and Cleanup ---------------------- */

  //Note: the arrays themselves don't need to be released, and trying to do so
    //will result in errors. The sub-array pointers do need to be
    //released, which is done below.
    for (i = 0; i < N; i++) {
      free(dW1[i]);
      free(dW2[i]);
    }

    for (i = 0; i < N; i++) {
      free(PricePath[i]);
      free(VolatilityPath[i]);
    }

    printf("Shutdown complete.\n");

  return 0;
}
