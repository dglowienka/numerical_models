#include "vector.h"

void add_vector(vector<vector<double>> * matrix_in, vector<double> * vector_in, int j)
{
	// Init matrix, in case it is not. Doesn't influence the data
	matrix_in->resize(vector_in->size());

	for (unsigned int i = 0; i < vector_in->size(); i++)
	{
		// Add j number of column
		(*matrix_in)[i].resize(j + 1);
		matrix_in->at(i).at(j) = vector_in->at(i);
	}
}

vector<double> get_vector(vector<vector<double>> matrix_in, int j)
{
	vector<double> tmp_vector;
	tmp_vector.resize(matrix_in.size());
	for (unsigned int i = 0; i< matrix_in.size(); i++)
	{
		tmp_vector[i] = matrix_in[i][j];
	}

	return tmp_vector;
}