#define _USE_MATH_DEFINES // Use constants in cmath
#include <cmath>
#include <gsl/gsl_sf_bessel.h> // Bessel function
#include <gsl/gsl_integration.h> // Definite integration
#include <complex> // Complex number e.g. root of negative number
#include "dissociation.h"

double dissociation(double f_f, void * params)
{
	struct dissociation_params *p = (struct dissociation_params *)params; // pointer on structure with all the data
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(5000); // workspace for GSL function

	double tmp_D, tmp_D_int, error, tmp_Nf;
	double a_0 = p->a_0;
	double E = p->E;

	gsl_function F; // Definition of GSL function
	F.function = &dissociation_P_F; // Localization of function to integrate
	F.params = p; // Paramaters used in the function for integration

	tmp_Nf = 4 / (sqrt(M_PI)*pow(a_0, 3));
	if (E != 0)
	{
		double size;
		//gsl_integration_qags(&F, 0, 100, 0, 1e-10, 1000, w, &tmp_D_int, &error); // a to b
		gsl_integration_qagiu(&F, 0, 0, 1e-10, 5000, w, &tmp_D_int, &error); // a to infinity
		//printf("intervals       = %zu\n", w->size);
		tmp_D = f_f*tmp_Nf*tmp_D_int;
	}
	else tmp_D = 0;

	gsl_integration_workspace_free(w); // release workspace
	return tmp_D;
}

double dissociation_P_F(double a, void * params)
{
	struct dissociation_params *tmp_p = (struct dissociation_params *)params;
	// Load data from structure
	double a_0 = tmp_p->a_0;
	double E = tmp_p->E;
	double E_r = tmp_p->E_r;
	double T = tmp_p->T;
	double u_n = tmp_p->u_n;
	double u_p = tmp_p->u_p;
	double k_f = tmp_p->k_f;
	//double tmp_P_F = ftest_2(a,p);

	// Calculate P_F
	double tmp_B = bessel(E,E_r,T);
	double tmp_A = pre_exp(a,u_n,u_p,E_r);
	double tmp_E_B = coulomb(a, E_r);
	double tmp_k_d = dissociation_rate(tmp_A,tmp_E_B,tmp_B,T);
	double tmp_P_F = probability(tmp_k_d,k_f)*gauss(a,a_0);
	return tmp_P_F;
}

double probability(double k, double k_f)
{
	double tmp_p;
	tmp_p = k / (k + k_f);
	return tmp_p;
}

double dissociation_rate(double A, double E_B, double b, double T)
{
	double tmp_k_d;
	tmp_k_d = A*exp(-E_B / (C_k_B*T))*b;
	return tmp_k_d;
}

double pre_exp(double a, double u_n, double u_p, double E_r)
{
	double tmp_p_exp;
	tmp_p_exp = (u_n + u_p)*(3 * C_q) / (4 * M_PI*C_E_0*E_r*pow(a, 3));
	return tmp_p_exp;
}

double coulomb(double a, double E_r)
{
	double tmp_E_B;
	tmp_E_B = pow(C_q, 2) / (4 * M_PI*C_E_0*E_r*a);
	return tmp_E_B;
}

double bessel(double E, double E_r, double T)
{
	double tmp_B, tmp_b, tmp_sqrt_imag;
	tmp_b = (pow(C_q, 3)*E) / (8 * M_PI*C_E_0*E_r*pow(C_k_B, 2) * pow(T, 2));
	// The square root of negative number is only imaginary number
	const std::complex<double> tmp_sqrt = std::sqrt(std::complex<double>(-2 * tmp_b));
	// Take only imaginary part, because there is no real part in the result
	tmp_sqrt_imag = std::imag(tmp_sqrt);
	// Calculate Regular Modified Cylindrical Bessel Functions. It is the same as in SciLab besselj(1,x)
	tmp_B = gsl_sf_bessel_I1(2 * tmp_sqrt_imag) / tmp_sqrt_imag;
	return tmp_B;
}

double gauss(double a, double a_0)
{
	double tmp_Fa;
	tmp_Fa = pow(a, 2) * exp(-pow(a, 2) / pow(a_0, 2));
	return tmp_Fa;
}