#include <iostream>
#include <fstream>
#include "dissociation.h"

int main(void)
{
	using namespace std;
	double a_0 = 0.7E-9;
	double E[2000];
	double D[2000];
	double E_r = 3;
	double T = 300;
	double u_n = 1E-4;
	double u_p = 0;
	double k_f = 1E8;
	double f_f = 0.99;
	E[0] = 0;

	//struct dissociation_params D_param1 = { a_0, E_r, T, u_n, u_p, k_f };
	//double D1 = dissociation(f_f, &D_param1);
	ofstream Onsager;
	Onsager.open("Onsager_braun.txt");

	struct dissociation_params D_param = { a_0, E_r, T, u_n, u_p, k_f };
	cout << "Dissociation D: \t" << "E \t" << "\n";
	Onsager << "Dissociation D: \t" << "E \t" << "\n";
	for (int i = 0; i < 1000; i++)
	{
		double E_scaled, kd_scaled, kd_unscaled;
		E[i + 1] = E[i] + 1E6;
		E_scaled = scaling(E[i], &D_param); // Scaling of E - electric field for input
		
		kd_scaled = dissociation(E_scaled);
		kd_unscaled = unscaling(kd_scaled, &D_param); // Unscaling of kd - dissociation rate
		D[i] = probability(kd_unscaled, k_f);

		cout << E[i] << "\t" << D[i] << "\n";
		Onsager << E[i] << "\t" << D[i] << "\n";
	}

	Onsager.close();
	return 0;
}