#include <iostream>
#include <fstream>
#include <vector>
#include "gaussian_elimination.h"
#include "poisson_equation.h"
#include "vector.h"

int main(void)
{
	using namespace std;
	vector<vector<double>> results;
	vector<double> tmp_results;
	unsigned int x = 8;
	double t = 3;
	double dx = 2;
	unsigned int L = x / dx;
	double dt = 1;
	double lambda = 0.020875;
	
	results.resize(x, vector<double>(1)); // Dynamically increase the size of vector
	for (int t = 0; t <3; t++)
	{
		tmp_results = Fourier(get_vector(results, t), lambda, L);
		add_vector(&results, &tmp_results, t + 1);
	}
	return 0;
}