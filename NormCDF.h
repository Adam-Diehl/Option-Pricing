//Normal CDF function
double NormCDF(double z, double sigma);

double NormCDF(double z, double sigma) {

  return 0.5 * erfc(- z / (sigma*sqrt(2)));

}
