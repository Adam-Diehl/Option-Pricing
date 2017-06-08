/*
  Uses Monte Carlo simulation to price an option in a Black
  Scholes world.

  Header file contains the pseudo-random normal number generator. 

  Author: Adam Diehl

  Date: June 7, 2017
*/

//Mathematical constants
#define Pi 3.14159265358979

//NORMSINV constants
#define  C1  (-7.784894002430293e-03)
#define  C2  (-3.223964580411365e-01)
#define  C3  (-2.400758277161838e+00)
#define  C4  (-2.549732539343734e+00)
#define  C5   4.374664141464968e+00
#define  C6   2.938163982698783e+00
#define  D1   7.784695709041462e-03
#define  D2   3.224671290700398e-01
#define  D3   2.445134137142996e+00
#define  D4   3.754408661907416e+00

//Define functions
double RNORM();

//Functions
double RNORM() {
  double p;
  double e;
  double u;
  double q;
  double x;

  p = (((double)rand()+1)/((double)RAND_MAX+2));
  q = sqrt(-2*log(p));
  x = (((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
  e = 0.5 * erfc(-x/sqrt(2)) - p;
  u = e * sqrt(2*Pi) * exp(0.5*x*x);
  return (x - u/(1 + 0.5*x*u));
}
