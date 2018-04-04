//Declare function -- calculates a heston european option price
double HestonInnerLoop(double S, double K, double r, double T, double IV, int OptionType, double Kappa, double Theta, double Gamma, double Rho, int ChunkSize, int M);

double HestonInnerLoop(double S, double K, double r, double T, double IV, int OptionType, double Kappa, double Theta, double Gamma, double Rho, int ChunkSize, int M) {

  //Seed the RNG
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<float> RNORM(0.0, 1.0);

  //Allocate memory for core arrays (on the stack) and values
  double StockPaths[ChunkSize][M] = {};
  double VolatilityPaths[ChunkSize][M] = {};
  double dW1[ChunkSize][M] = {};
  double dW2[ChunkSize][M] = {};

  double Payoff = 0;

  //Calculate parameters
  double dt = T/M;
  double SqrtDt = sqrt(dt);
  double RhoCoef = sqrt(1-Rho*Rho);

  int Offset = ChunkSize/2;

  //Fill out arrays with correlated gaussian random numbers
  for (int i=0; i<Offset; ++i) {
    for (int j=0; j<M; ++j) {
        dW1[i][j] = RNORM(gen);
        dW1[i + Offset][j] = dW1[i][j]; //Antithetic paths
        dW2[i][j] = dW1[i][j]*Rho + RhoCoef*RNORM(gen);
        dW2[i + Offset][j] = dW2[i][j]; //Antithetic paths
    }
  }

  //Set initial price and volatility
  for(int i = 0; i < ChunkSize; ++i) {
    StockPaths[i][0] = S;
    VolatilityPaths[i][0] = IV;
  }

  //Iterate through time with Milstein discretization
  #pragma omp parallel for
  for(int i = 0; i < ChunkSize; ++i) {
    for(int j = 1; j < M; ++j) {

      //Path iteration
      VolatilityPaths[i][j] = VolatilityPaths[i][j-1] + Kappa * (Theta - VolatilityPaths[i][j-1]) * dt + Gamma * sqrt(VolatilityPaths[i][j-1]) * SqrtDt * dW2[i][j];
      StockPaths[i][j] = StockPaths[i][j-1] + r * StockPaths[i][j-1] * dt + sqrt(VolatilityPaths[i][j-1]) * StockPaths[i][j-1] * SqrtDt * dW1[i][j];

      //std::cout << StockPaths[i][j] << std::endl;

      //Volatility truncation
      if(VolatilityPaths[i][j] < 0)
        VolatilityPaths[i][j] = 0;

    }
  }

  //Determine the option payoffs
  if(OptionType) { //Put
    for(int i = 0; i < ChunkSize; ++i) {
      if((K - StockPaths[i][(M-1)]) > 0) {
        Payoff += (K - StockPaths[i][(M-1)]);
      }
    }
  } else { //Call
    for(int i = 0; i < ChunkSize; ++i) {
      if((StockPaths[i][(M-1)] - K) > 0) {
        Payoff += (StockPaths[i][(M-1)] - K);
      }
    }
  }

  //std::cout << exp(-r*T)*Payoff/ChunkSize << std::endl;

  return Payoff;

}
