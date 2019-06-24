#include <iostream>
#include <fstream>
#include <vector>
#include "gaussian_elimination.h"

int main(void)
{
	using namespace std;
	int  m, n;
	double l, r, t, b;
	vector<double> b_it;
	vector<vector<double>> A_it;
	vector<vector<double>> a;

	ofstream laplace;
	laplace.open("Laplace.txt");

	printf("\tEnter boundary conditions\n");
	printf("\tValue on left side: ");
	l = 75;
	printf("\tValue on right side: ");
	r = 50;
	printf("\tValue on top side: ");
	t = 100;
	printf("\tValue on bottom side: ");
	b = 0;

	printf("\tEnter number of steps in x and y direction: ");
	m = 8;
	printf("\tEnter number of steps in y direction: ");
	n = 10;
	
	// Init all vectors
	b_it.resize(n*m);
	A_it.resize(n*n, vector<double>(m*m));
	a.resize(n+2, vector<double>(m+2));

	//assigning boundary values begins
	for (int i = 0; i <= m+1; i++)   
	{
		a[i][0] = b;
		a[i][n+1] = t;
	}

	//assigning boundary values ends
	for (int i = 0; i <= n+1; i++)
	{
		a[0][i] = l;
		a[m+1][i] = r;
	}                         


	// Find A matrix and b vector for Ax=b linear algebric
	int k_it = 0;
	for (int j = 1; j <= n; j++)
	{
		for (int i = 1; i <= m; i++)
		{
			A_it[k_it][k_it] = -4;
			// Middle (no boundary around)
			if (i > 1 && i < m && j>1 && j < n)
			{
				A_it[k_it - n][k_it] = 1; //i,j-1
				A_it[k_it - 1][k_it] = 1; //i-1,j
				A_it[k_it + 1][k_it] = 1; //i+1,j
				A_it[k_it + n][k_it] = 1; //i,j+1
				b_it[k_it] = 0;
			}
			
			// Bottom
			if (j == 1)
			{
				
				// Bottom-left
				if (i == 1)
				{
					A_it[k_it + 1][k_it] = 1; //i+1,j
					A_it[k_it + n][k_it] = 1; //i,j+1
					b_it[k_it] = -(b + l);
				}
				// Bottom-right
				if (i == m) 
				{
					A_it[k_it - 1][k_it] = 1; //i-1,j
					A_it[k_it + n][k_it] = 1; //i,j+1
					b_it[k_it] = -(b + r);
				}
				// Bottom-middle
				if ((i > 1) && (i < m))
				{
					A_it[k_it - 1][k_it] = 1; //i-1,j
					A_it[k_it + 1][k_it] = 1; //i+1,j
					A_it[k_it + n][k_it] = 1; //i,j+1
					b_it[k_it] = -b;
				}
			}
			// Top
			if (j == n)
			{
				// Top-left
				if (i == 1) 
				{
					A_it[k_it - n][k_it] = 1; //i,j-1
					A_it[k_it + 1][k_it] = 1; //i+1,j
					b_it[k_it] = -(t + l);
				}
				// Top-right
				if (i == m) 
				{
					A_it[k_it - n][k_it] = 1; //i,j-1
					A_it[k_it - 1][k_it] = 1; //i-1,j
					b_it[k_it] = -(t + r);
				}
				// Top-middle
				if ((i > 1) && (i < m))
				{
					A_it[k_it - n][k_it] = 1; //i,j-1
					A_it[k_it - 1][k_it] = 1; //i-1,j
					A_it[k_it + 1][k_it] = 1; //i+1,j
					b_it[k_it] = -t;
				}
			}
			// Left
			if (i == 1)
			{
				// Left-bottom
				// Left-top
				// Left-middle
				if (j > 1 && j < m)
				{
					A_it[k_it - n][k_it] = 1; //i,j-1
					A_it[k_it + 1][k_it] = 1; //i+1,j
					A_it[k_it + n][k_it] = 1; //i,j+1
					b_it[k_it] = -l;
				}
			}
			//Right
			if (i == m)
			{
				// Right-top
				// Right-bottom
				// Right-middle
				if (j > 1 && j < m)
				{
					A_it[k_it - n][k_it] = 1; //i,j-1
					A_it[k_it - 1][k_it] = 1; //i-1,j
					A_it[k_it + n][k_it] = 1; //i,j+1
					b_it[k_it] = -r;
				}
			}
			k_it++;
		}
	}
	vector<double> tmp_result;

	tmp_result = gaussian_elimination(A_it, b_it);
	
	// Copy results of x to array with all the data
	int count=0;
	for (int j = 1; j <= n; j++)
	{
		for (int i = 1; i <= m; i++)
		{
			a[i][j] = tmp_result[count];
			count++;
		}
		
	}

	// Print the results
	for (int i = 0; i < a.size(); i++)
	{
		for (int j = 0; j < a[i].size(); j++)
		{
			laplace << a[i][j] << "\n";
			count++;
		}

	}
	
	laplace.close();
}