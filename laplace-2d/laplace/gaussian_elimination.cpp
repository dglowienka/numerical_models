#include <vector>
#include <gsl/gsl_linalg.h>
#include "gaussian_elimination.h"

using namespace std;
vector<double> gaussian_elimination(vector<vector<double>> A_data, vector<double> b_data)
{
	int size_m = b_data.size(); // Size of b_data, as the gaussian method is only for square matrix
	
	double* tmp_A_data = new double[size_m*size_m]; // Transform 2D vector to 1D array, needed for GSL
	for (int row = 0; row < size_m; row++)
	{
		for (int col = 0; col < size_m; col++)
		{
			tmp_A_data[row * size_m + col] = A_data[row][col];
		}
	}

	double* tmp_b_data = &b_data[0]; // Vector to double* array for GSL

	gsl_matrix_view m = gsl_matrix_view_array(tmp_A_data, size_m, size_m); // Convert array to GSL library
	gsl_vector_view b = gsl_vector_view_array(tmp_b_data, size_m);
	gsl_vector *x = gsl_vector_alloc(size_m);
	
	int s;
	vector<double> tmp_vect(size_m);
	gsl_permutation * p = gsl_permutation_alloc(size_m);
	gsl_linalg_LU_decomp(&m.matrix, p, &s);
	gsl_linalg_LU_solve(&m.matrix, p, &b.vector, x);

	for (int i = 0; i < size_m; i++)
	{
		tmp_vect[i] = gsl_vector_get(x, i);
	}

	gsl_permutation_free(p);
	gsl_vector_free(x);
	delete [] tmp_A_data;
	return tmp_vect;
}