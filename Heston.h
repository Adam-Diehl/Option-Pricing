/*
  Uses Monte Carlo simulation to price a vanilla European option using the
  Heston stochastic volatility model.

  Header file contains the random normal generator, and required functions.

  Author: Adam Diehl

  Date: June 22, 2017
*/

/* Declare formulas */
double RNORM();
int FellerCondition(double Kappa, double Theta, double Gamma);

/* Gaussian Random Normal Generator */
double RNORM() {
  static double V1, V2, S;
  static int phase = 0;
  double X;

  if(phase == 0) {
    do {
      double U1 = (double)rand() / RAND_MAX;
      double U2 = (double)rand() / RAND_MAX;

      V1 = 2 * U1 - 1;
      V2 = 2 * U2 - 1;
      S = V1 * V1 + V2 * V2;
      } while(S >= 1 || S == 0);

    X = V1 * sqrt(-2 * log(S) / S);
  } else
    X = V2 * sqrt(-2 * log(S) / S);

  phase = 1 - phase;

  return X;
}

/* Feller Condition */
int FellerCondition(double Kappa, double Theta, double Gamma) {

  //The Feller Condition says that when 2*Kappa*Theta > Gamma^2, Vol > 0.

  if(2*Kappa*Theta > Gamma*Gamma) {
    return 0; //True
  }
  else {
    return 1; //False
  }
}
