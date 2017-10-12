/*

File: Prices vanilla (European) options using the Heston model using a
  parallelized Milstein discretization

Author: Adam Diehl

Date: oct 12, 2017

Tested with GCC 7.2.0 with...
gcc-7 -o Asian AsianEuropeanBlackScholes.cpp -O3 -lstdc++ -fopenmp

Note: try adding compiler flag: -ffast-math

*/

/*---------------------------- Import Libraries ------------------------------*/

#include <cmath>
#include <iostream>
#include <random>

/*--------------------------- Declare Parameters -----------------------------*/

//Monte Carlo parameters
  const double dt = 100;
  const int N = 1000000;
  const int NCols = 100; //T*dt

//Black Scholes parameters
  const double S = 100; //Initial stock price
  const double K = 100; //Strike price
  const double r = 0.0319; //Risk-free rate
  const double T = 1; //Time to maturity (in years)
  const double IV = 0.10201; //Initial Volatility

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
  double drift = (r-IV*IV*0.5)*T; //Risk neutral drift term
  double StockWeight = 1/(NCols + 0.00);
  double PathWeight = 1/(N + 0.00);

  //Allocate arrays - the number of simulations is fixed and known at compile
  std::vector<std::vector<double> > dW(N,std::vector<double>(NCols));
  std::vector<std::vector<double> > PricePath(N,std::vector<double>(NCols));
  std::vector<double> Average(N, 0.00);

  //Status update
  std::cout << "Memory allocated" << std::endl;

/*------------------------ Random Number Generation --------------------------*/

  //The built in gaussian RNG from the STL library is too slow. This uses a much
      //faster method

  //Status update
  std::cout << "Generating random numbers..." << std::endl;

  //Seed random numbers
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<float> RNORM(0.0, 1.0);

  //Fill out array with correlated gaussian random numbers
  for (int i=0; i<N; ++i) {
    for (int j=0; j<NCols; ++j) {
        dW[i][j] = RNORM(gen);
    }
  }

  //Status update
  std::cout << "Random numbers computed" << std::endl;

/*---------------------- Simulate the Path Processes -------------------------*/

  //Status update
  std::cout << "Generating price paths..." << std::endl;

  //Set initial price and volatility
  for(int i = 0; i < N; ++i) {
    PricePath[i][0] = S;
  }

  //Iterate through time with Milstein discretization
    #pragma omp parallel for
    for(int i = 0; i < N; ++i) {
      for(int j = 1; j < NCols; ++j) {

        //Path iteration
        PricePath[i][j] = PricePath[i][j-1] + r*PricePath[i][j-1]*deltaT + IV * PricePath[i][j-1]*sqrDeltaT*dW[i][j] + 0.5*IV*IV*deltaT*(dW[i][j]*dW[i][j]-1);

      }
    }

    //Status update
    std::cout << "Paths computed" << std::endl;

/*---------------------- Calculate the Option Payoff -------------------------*/

  //Status update
  std::cout << "Computing option prices..." << std::endl;

  //Determine the average option price
  #pragma omp parallel for
  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < NCols; ++j) {
        Average[i] += PricePath[i][j]*StockWeight;
    }
  }

  //Calculate payoffs
  for(int i = 0; i < N; ++i) {
      //Call Payoffs
      if(Average[i] - K > 0) {
        CallPayoff += (Average[i] - K);
      }
      //Put Payoffs
      if(K - Average[i] > 0) {
        PutPayoff += (K - Average[i]);
      }
  }

  CallPayoff = discount * CallPayoff * PathWeight;
  PutPayoff = discount * PutPayoff * PathWeight;

  std::cout << std::endl;
  std::cout << "Call price = $" << CallPayoff << std::endl;
  std::cout << "Put price = $" << PutPayoff << std::endl;

  return 0;
}
