/*******************************************************************************/
/*************Program to solve Laplace's equation by finite difference method*************/
/*************************** Developed by Mahesha M G ***************************/
#include <stdio.h>
#include <iostream>
#include <fstream>
int main(void)
{
	using namespace std;
	int i, j, k, m, n, x, y;
	float a[20][20], l, r, t, b;
	ofstream laplace;
	laplace.open("Laplace.txt");

	printf("\t_______________________________________________________________\n");
	printf("\tProgram to solve Laplace's equation by finite difference method\n");
	printf("\t****************** Developed by Mahesha M G *******************\n");
	printf("\t_______________________________________________________________\n");
	printf("\tEnter boundary conditions\n");
	printf("\tValue on left side: ");
	l = 75;
	printf("\tValue on right side: ");
	r = 50;
	printf("\tValue on top side: ");
	t = 100;
	printf("\tValue on bottom side: ");
	b = 0;
	printf("\tEnter number of steps in x direction: ");
	m = 10;
	printf("\tEnter number of steps in y direction: ");
	n = 10;
	m++;
	n++; //number of mesh points is one more than number of steps
	for (i = 1; i <= m; i++)   //assigning boundary values begins
	{
		a[i][1] = b;
		a[i][n] = t;
	}
	for (i = 1; i <= n; i++)
	{
		a[1][i] = l;
		a[m][i] = r;
	}                         //assigning boundary values ends
	for (i = 2; i<m; i++)
		for (j = 2; j<n; j++)
			a[i][j] = t; //initialization of interior grid points 
	for (k = 0; k<100; k++)
	{
		for (i = 2; i<m; i++)
		{
			for (j = 2; j<n; j++)
			{
				a[i][j] = (a[i - 1][j] + a[i + 1][j] + a[i][j - 1] + a[i][j + 1]) / 4;
			}
		}
	}                     //calculation by Gauss-Seidel Method
	for (i = 1; i <= m; i++)
	{
		for (j = 1; j <= n; j++)
			laplace << a[i][j] << "\n";
	}
	laplace.close();
	printf("\nData stored\nPress any key to exit...");
}
