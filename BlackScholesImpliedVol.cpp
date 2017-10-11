/*

File: computes Black-Scholes implied volatilities using the bisection method

Author: Adam Diehl

Date: October 10, 2017

Tested with GCC 7.2.0 with...
gcc-7 -o BSIV BlackScholesImpliedVol.cpp -O3 -lstdc++

*/

/*------------------------- Mathematical Constants ---------------------------*/

const double Pi = 3.141592653589793;

/*---------------------------- Import Libraries ------------------------------*/

#include <cmath>
#include <iostream>
#include "NormCDF.h"

/*--------------------------- Declare Parameters -----------------------------*/

//Stock price
    const double S = 100.0;
//Strike price
    const double K = 100.00;
//Risk-free rate
    const double r = 0.074;
//Time to maturity
    const double T = 45.0/252.0;
//Market option price
    const double c = 4.0457;
//Call (0) or Put (1)
    const int OptionType = 0;
//Error tolerance
    const double E = 0.000001;
//Allocate an initial IV so the function will compile
  const double ImplVol = 0.005;

/*---------------------------- Declare Functions -----------------------------*/
double BlackScholes(double S, double K, double T, double r, double ImplVol, int OptionType);

/*---------------------------- Begin Main File -------------------------------*/

int main() {

    //Function markers
    int j = 0; //Tracks convergence times

    //Binomial parameters
    double IV = 0.00; //Initial Volatility
    double LowerBound = 0.0;
    double UpperBound = 2.0;

    //Calculate an initial estimate using the Brenner-Subrahmanyam formula
    IV = c*sqrt(2*Pi) / (S*sqrt(T));

    //Expand the interval to account for moneyness
    LowerBound *= IV;
    UpperBound *= IV;

    //Utilize the bisection method to find the IV
    while(fabs(BlackScholes(S, K, T, r, IV, OptionType) - c) - E > 0) {
      if((BlackScholes(S, K, T, r, IV, OptionType) > c) && (BlackScholes(S, K, T, r, UpperBound, OptionType) > c)) {
        UpperBound = IV;
        IV = 0.5*(UpperBound + LowerBound);
      } else {
        LowerBound = IV;
        IV = 0.5*(UpperBound + LowerBound);
      }
      j += 1;
    }

    std::cout << "The implied volatility is " << IV << "." << std::endl;
    std::cout << "The algorithm converged in " << j << " steps." << std::endl;

    return 0;
}
/*---------------------------- Helper Functions ------------------------------*/
double BlackScholes(double S, double K, double T, double r, double ImplVol, int OptionType) {

  double d1 = (1/(ImplVol*sqrt(T))) * (log(S/K) + (r + 0.5*ImplVol*ImplVol)*T);
  double d2 = d1 - ImplVol*sqrt(T);

  if(OptionType == 0) { //Call
    return S*NormCDF(d1,1) - K*exp(-r*T)*NormCDF(d2,1);
  } else { //Put
    return -S*NormCDF(-d1,1) + K*exp(-r*T)*NormCDF(-d2,1);
  }
}
