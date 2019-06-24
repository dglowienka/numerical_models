#include "poisson_equation.h"
#include "gaussian_elimination.h"

vector<double> Fourier(vector<double> T_it, double lambda, unsigned int x)
{
	double a_T, T_top, T_bottom;
	vector<vector<double>> Tit_A;
	vector<double> Tit_b;
	vector<double> tmp_T;

	Tit_A.resize(x, vector<double>(x));
	Tit_b.resize(x);

	T_top = 50;
	T_bottom = 100;

	a_T = lambda;

	//Fill the matrix A and vector b with diagonal data, A matrix from right and b vector from left part of equation
	//-lambda*T(i+1,t+1) + T(i,t+1)*(2*lambda + 1) - lambda*T(i-1,t) = T(i,t)
	int k_it = 0;
	for (unsigned int i = 1; i <= x; i++)
	{
		Tit_A[k_it][k_it] = (2*a_T + 1); // T(i,t+1)

		// Middle (no boundary around)
		if (i > 1 && i < x)
		{
			Tit_A[k_it + 1][k_it] = -lambda; // T(i-1,t+1)
			Tit_A[k_it - 1][k_it] = -lambda; // T(i+1,t+1)
			Tit_b[k_it] = T_it[i];
		}
		// Bottom
		if (i == 1)
		{
			Tit_A[k_it + 1][k_it] = -lambda; // T(i-1,t+1)
			Tit_b[k_it] = T_it[i] + lambda*T_bottom;
		}

		// Top
		if (i == x)
		{
			Tit_A[k_it - 1][k_it] = -lambda; // T(i+1,t+1)
			Tit_b[k_it] = T_it[i] + lambda*T_top;
		}
		k_it++;
	}

	// Vector x from Ax=b
	tmp_T = gaussian_elimination(Tit_A, Tit_b);

	// Put boundaries to the result
	tmp_T.insert(tmp_T.begin(), T_bottom);
	tmp_T.insert(tmp_T.end(), T_top);

	return tmp_T;
}
