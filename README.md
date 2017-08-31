# Option-Pricing
Option pricing algorithms in C/C++. The GaussPRNG file is used instead of the C++ standard library Gaussian random number generator due to the latter's extreme slowness. This generator is mostly accurate (should be within a cent of true prices for most cases with a high number of simulations). 
