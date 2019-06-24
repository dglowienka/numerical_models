#ifndef _VECTOR
#define _VECTOR

#include <vector>
using namespace std;

// Resize the matrix_in for new value and add the vector_in to the given j column 
void add_vector(vector<vector<double>> * matrix_in, vector<double> * vector_in, int j);

// Get one j column from 2D vector
vector<double> get_vector(vector<vector<double>> matrix_in, int j);

#endif // !_VECTOR
