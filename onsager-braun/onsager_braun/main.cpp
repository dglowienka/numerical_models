#include <iostream>
#include <fstream>
#include "dissociation.h"

int main(void)
{
	using namespace std;
	double a_0 = 0.54E-9;
	double E[2000];
	double D[2000];
	double E_r = 3;
	double T = 300;
	double u_n = 2.7E-6;
	double u_p = 0;
	double k_f = 1E7;
	double f_f = 0.99;
	E[0] = 0;

	struct dissociation_params D_param1 = { 1E7, a_0, E_r, T, u_n, u_p, k_f };
	double D1 = dissociation(f_f, &D_param1);
	ofstream Onsager;
	Onsager.open("Onsager_braun.txt");

	cout << "Dissociation D: \t" << "E \t" << "\n";
	Onsager << "Dissociation D: \t" << "E \t" << "\n";
	for (int i = 0; i < 1000; i++)
	{
		E[i+1] = E[i]+2E6;
		struct dissociation_params D_param = { E[i], a_0, E_r, T, u_n, u_p, k_f };
		D[i] = dissociation(f_f, &D_param);
		cout << E[i] << "\t" << D[i] << "\n";
		Onsager << E[i] << "\t" << D[i] << "\n";
	}

	Onsager.close();
	return 0;
}