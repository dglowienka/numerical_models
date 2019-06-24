#ifndef _GAUSSIAN
#define _GAUSSIAN

#include <vector>

// Gaussian elimination function for Ax=b, where A is matrix and b, x are vectors
std::vector<double> gaussian_elimination(std::vector< std::vector<double> > A_data, std::vector<double> b_data);

#endif // !_GAUSSIAN

