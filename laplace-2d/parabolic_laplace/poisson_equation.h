#ifndef _POISSON
#define _POISSON

using namespace std;

#include <cmath>
#include <vector>

// Calculate electrical potential from Poisson equation, where t is time, n and p are concentration of e and h vectors
vector<double> Fourier(vector<double> T_it, double lambda, unsigned int x);

#endif // !_POISSON
