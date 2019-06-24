#include <iostream>
#include <fstream>
#include "dissociation.h"
#include "gsl_math.h"

// Input function for scaling of E
dissociation_scaled_params dissociation_scaling(dissociation_params * params);

// Output function for unscaling of Kd(E)
double dissociation_unscaling(double kd_scaled, void * params);

int main(void)
{
	using namespace std;
	double a_0 = 0.43E-9;
	double E[2000];
	double D[2000];
	double E_r = 3;
	double T = 300;
	double u_n = 2E-7;
	double u_p = 0;
	double k_f = 1E8;
	double f_f = 0.99;
	E[0] = 1E4;

	ofstream Onsager;
	Onsager.open("Onsager_braun.txt");

	struct dissociation_params D_param = { E[0], a_0, E_r, T, u_n, u_p, k_f, f_f };
	struct dissociation_scaled_params D_params_scaled;

	cout << "Dissociation D: \t" << "E \t" << "\n";
	Onsager << "Dissociation D: \t" << "E \t" << "\n";
	for (int i = 0; i < 1000; i++)
	{
		double D_scaled;
		E[i + 1] = E[i] + 1e6;
		D_param.E = E[i];
		D_params_scaled = dissociation_scaling(&D_param); // Scaling of E - electric field for input
		D_scaled = dissociation(&D_params_scaled);
		//D_unscaled = unscaling(kd_scaled, &D_param); // Unscaling of kd - dissociation rate
		D[i] = dissociation_unscaling(D_scaled, &D_param);// unscaling(D_scaled, &D_param);
		
		cout << E[i] << "\t" << D[i] << "\n";
		Onsager << E[i] << "\t" << D[i] << "\n";
	}

	Onsager.close();
	return 0;
}

dissociation_scaled_params dissociation_scaling(dissociation_params * params)
{
	struct dissociation_params *tmp_p = params;
	struct dissociation_scaled_params tmp_p_scaled;
	double tmp_E_scaled, tmp_k_f_scaled, u_n_scaled, u_p_scaled, tmp_alpha_scaled;
	double tmp_b, tmp_E_c, tmp_r, k_f_c;

	// Load data from unscaled structure
	double E = tmp_p->E;
	double k_f = tmp_p->k_f;
	double E_r = tmp_p->E_r;
	double T = tmp_p->T;
	double a0 = tmp_p->a_0;
	double f_f = tmp_p->f_f;
	double u_max = 1E-8; // m2 v-1 s-1

	// Scaling E
	tmp_b = (pow(C_q, 3)) / (8 * M_PI*C_E_0*E_r*pow(C_k_B, 2) * pow(T, 2));

	tmp_E_c = -1 / (2 * tmp_b);
	tmp_E_scaled = E / tmp_E_c;
	tmp_p_scaled.E_scaled = tmp_E_scaled;

	// Scaling k_f
	k_f_c = pre_exp(tmp_p)*u_max / pow(coulomb_radius(tmp_p), 3);
	tmp_k_f_scaled = k_f / k_f_c;
	tmp_p_scaled.k_f_scaled = tmp_k_f_scaled;

	// Scaling alpha
	tmp_r = coulomb_radius(tmp_p);
	tmp_alpha_scaled = -pow(tmp_r, 2) / pow(a0, 2);
	tmp_p_scaled.alpha_scaled = tmp_alpha_scaled;

	// Scaling mobility
	u_n_scaled = (tmp_p->u_n) / u_max;
	u_p_scaled = (tmp_p->u_p) / u_max;
	tmp_p_scaled.u_n_scaled = u_n_scaled;
	tmp_p_scaled.u_p_scaled = u_p_scaled;

	// Fitting factor without scaling
	tmp_p_scaled.f_f = f_f;

	return tmp_p_scaled;
}

double dissociation_unscaling(double D_scaled, void * params)
{
	struct dissociation_params *tmp_p = (struct dissociation_params *)params;
	double tmp_D_unscaled;
	// Load data from structure
	double a_0 = tmp_p->a_0;
	double E_r = tmp_p->E_r;
	double T = tmp_p->T;
	double k_f = tmp_p->k_f;
	double tmp_Nf = 4 / (sqrt(M_PI)*pow(a_0, 3));

	tmp_D_unscaled = D_scaled*tmp_Nf*pow(coulomb_radius(tmp_p), 3);
	return tmp_D_unscaled;
}