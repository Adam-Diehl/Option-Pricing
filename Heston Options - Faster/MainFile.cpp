// gcc-7 -o HestonF MainFile.cpp -O3 -lstdc++ -fopenmp
/*---------------------------- Import Libraries ------------------------------*/

#include <cmath>
#include <iostream>
#include <random>
#include "HestonInnerLoop.h"

using namespace std;

/*--------------------------- Declare Parameters -----------------------------*/

//Black Scholes parameters
  const double S = 100; //Initial stock price
  const double K = 100; //Strike price
  const double r = 0.0319; //Risk-free rate
  const double T = 1; //Time to maturity (in years)
  const double IV = 0.010201; //Initial Volatility

//Heston parameters
  const double Kappa = 6.21; //Vol reversion
  const double Theta = 0.0192; //Long term vol
  const double Gamma = 0.61; // Vol of vol
  const double Rho = -0.7; //Correlation coefficient

//Monte Carlo parameters
  const int M = 100; //Time subdivisions
  const int N = 100000; //Number of paths
  const int ChunkSize = 1000; //Size of inner loop chunks

/*---------------------------- Begin Main File -------------------------------*/

int main () {

  //Print input parameters
  cout << endl;
  cout << "Generating Parameters: " << endl;
  cout << "S = $" << S << endl;
  cout << "K = $" << K << endl;
  cout << "r = " << r << endl;
  cout << "Initial Volatility = " << IV << endl;
  cout << "T = " << T << endl;
  cout << "Vol Reversion = " << Kappa << endl;
  cout << "Long Term Vol = " << Theta << endl;
  cout << "Vol of Vol = " << Gamma << endl;
  cout << "Correlation Coefficient = " << Rho << endl;
  cout << "Number of Paths = " << N << endl;
  cout << "Number of Time Increments = " << M << endl;
  cout << endl;

  //Calculate other parameters
  double discount = exp(-r*T);
  int NumChunks = N/ChunkSize;

  //Set option type
  int OptionType = 0; //0 for call, 1 for put

  //Create vars to hold intermediate values
  double OptionPayoffs[NumChunks];
  double OptionPrice;

  //Calculate option prices
  for(int i = 0; i < NumChunks; ++i) {
    OptionPayoffs[i] = HestonInnerLoop(S, K, r, T, IV, OptionType, Kappa, Theta, Gamma, Rho, ChunkSize, M);
  }
  for(int i = 0; i < NumChunks; ++i) {
    OptionPrice += OptionPayoffs[i];
  }

  OptionPrice *= discount/N;

  cout << "Option price = $" << OptionPrice << endl;

  return 0;
}
