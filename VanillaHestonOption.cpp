/*

File: Prices vanilla (European) options using the Heston model using a
  parallelized Milstein discretization

Author: Adam Diehl

Date: Aug 29, 2017

Tested with GCC 7.2.0 with...
gcc-7 -o Heston VanillaHestonOption.cpp -O3 -lstdc++ -fopenmp

*/

/*---------------------------- Import Libraries ------------------------------*/

//#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include "GaussPRNG.h"

/*--------------------------- Declare Parameters -----------------------------*/

//Black Scholes parameters
  const double S = 100; //Initial stock price
  const double K = 100; //Strike price
  const double r = 0.0319; //Risk-free rate
  const double T = 1; //Time to maturity (in years)
  const double IV = 0.010201; //Initial Volatility

//Heston parameters
  const double Kappa = 6.21; //Vol reversion
  const double Theta = 0.019; //Long term vol
  const double Gamma = 0.019; // Vol of vol
  const double Rho = -0.7; //Correlation coefficient

//Monte Carlo parameters
  const double dt = 100;
  const int N = 1000000;
  const int NCols = 100; //T*dt

//Allocate memory for computable constants
  double CallPayoff = 0;
  double PutPayoff = 0;

/*---------------------------- Begin Main File -------------------------------*/

int main () {

/*----------------------------- Program Setup --------------------------------*/

  //Status update
  std::cout << "Allocating memory for arrays..." << std::endl;

  //Compute constants that depend on program parameters
  double deltaT = T/dt;
  double sqrDeltaT = sqrt(deltaT);
  double discount = exp(-r*T);
  double RhoCoefSpec = sqrt(1-Rho*Rho);

  //Allocate arrays - the number of simulations is fixed and known at compile
  std::vector<std::vector<double> > dW1(N,std::vector<double>(NCols));
  std::vector<std::vector<double> > dW2(N,std::vector<double>(NCols));
  std::vector<std::vector<double> > PricePath(N,std::vector<double>(NCols));
  std::vector<std::vector<double> > VolatilityPath(N,std::vector<double>(NCols));

  //Status update
  std::cout << "Memory allocated" << std::endl;

/*------------------------ Random Number Generation --------------------------*/

  //The built in gaussian RNG from the STL library is too slow. This uses a much
      //faster method

  //Status update
  std::cout << "Generating correlated random numbers..." << std::endl;

  //Seed random numbers
  srand (time(NULL));

  //Fill out array with correlated gaussian random numbers
  for (int i=0; i<N; ++i) {
    for (int j=0; j<NCols; ++j) {
        dW1[i][j] = RNORM();
        dW2[i][j] = dW1[i][j]*Rho + RhoCoefSpec*RNORM();
    }
  }

  //Status update
  std::cout << "Random numbers computed" << std::endl;

/*---------------------- Simulate the Path Processes -------------------------*/

  //Status update
  std::cout << "Generating price/volatility paths..." << std::endl;

  //Set initial price and volatility
  for(int i = 0; i < N; ++i) {
    PricePath[i][0] = S;
    VolatilityPath[i][0] = IV;
  }

    //Iterate through time with Milstein discretization
    #pragma omp parallel for
    for(int i = 0; i < N; ++i) {
      for(int j = 1; j < NCols; ++j) {

        //Path iteration
        VolatilityPath[i][j] = VolatilityPath[i][(j-1)] + Kappa*(Theta - VolatilityPath[i][(j-1)])*deltaT + Gamma*sqrt(VolatilityPath[i][(j-1)])*sqrDeltaT*dW1[i][(j-1)];
        PricePath[i][j] = PricePath[i][(j-1)] * exp((r - 0.5*VolatilityPath[i][(j-1)])*deltaT + sqrt(VolatilityPath[i][(j-1)]*deltaT)*dW2[i][(j-1)]);

        //Volatility truncation
        if(VolatilityPath[i][j] < 0)
          VolatilityPath[i][j] = 0;

      }
    }

    //Status update
    std::cout << "Paths computed" << std::endl;

/*---------------------- Calculate the Option Payoff -------------------------*/

  //Determine the option payoff and discount
  //Call payoff
  for(int i = 0; i < N; ++i) {
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

  std::cout << std::endl;
  std::cout << "Call price = $" << CallPayoff << std::endl;
  std::cout << "Put price = $" << PutPayoff << std::endl;

  return 0;
}
