#define _USE_MATH_DEFINES // Use constants in cmath
#include <cmath>
#include <gsl/gsl_sf_bessel.h> // Bessel function
#include <complex> // Complex number e.g. root of negative number
#include "dissociation.h"

double dissociation(double E_scaled)
{
	double tmp_D;

	if (E_scaled != 0)
	{
		tmp_D = bessel(E_scaled);
	}
	else tmp_D = 0;

	return tmp_D;
}

double dissociation_P_F(double kd_scaled, void * params)
{
	struct dissociation_params *tmp_p = (struct dissociation_params *)params;
	// Load data from structure
	double k_f = tmp_p->k_f;

	// Calculate P_F
	double tmp_P_F = probability(kd_scaled,k_f);
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

double bessel(double E_scaled)
{
	double tmp_B, tmp_sqrt_imag;
	// The square root of negative number is only imaginary number
	const std::complex<double> tmp_sqrt = std::sqrt(std::complex<double>(E_scaled));
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

double scaling(double E_unscaled, void * params)
{
	struct dissociation_params *tmp_p = (struct dissociation_params *)params;
	double tmp_E_scaled, tmp_b;
	// Load data from structure
	double E_r = tmp_p->E_r;
	double T = tmp_p->T;

	tmp_b = (pow(C_q, 3)) / (8 * M_PI*C_E_0*E_r*pow(C_k_B, 2) * pow(T, 2));
	tmp_E_scaled = -2 * E_unscaled*tmp_b;
	return tmp_E_scaled;
}

double unscaling(double kd_scaled, void * params)
{
	struct dissociation_params *tmp_p = (struct dissociation_params *)params;
	double tmp_kd_unscaled;
	// Load data from structure
	double a_0 = tmp_p->a_0;
	double E_r = tmp_p->E_r;
	double T = tmp_p->T;
	double u_n = tmp_p->u_n;
	double u_p = tmp_p->u_p;
	double k_f = tmp_p->k_f;

	// Calculate P_F
	double tmp_A = pre_exp(a_0, u_n, u_p, E_r);
	double tmp_E_B = coulomb(a_0, E_r);
	tmp_kd_unscaled = kd_scaled*tmp_A*exp(-tmp_E_B / (C_k_B*T));
	return tmp_kd_unscaled;
}
