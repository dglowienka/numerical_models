#include <gsl/gsl_sf_bessel.h> // Bessel function
#include <gsl/gsl_integration.h> // Definite integration
#include "gsl_math.h" // Constants included Pi
#include <complex> // Complex number e.g. root of negative number
#include "dissociation.h"

double dissociation(dissociation_scaled_params * params_scaled)
{
	struct dissociation_scaled_params *p = (struct dissociation_scaled_params *)params_scaled;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(5000); // workspace for GSL function
	double E_scaled = p->E_scaled;
	double f_f = p->f_f;
	double tmp_D, tmp_D_ff, error;

	gsl_function F; // Definition of GSL function
	F.function = &dissociation_P_F; // Localization of function to integrate
	F.params = p; // Paramaters used in the function for integration

	if (E_scaled != 0)
	{
		gsl_integration_qagiu(&F, 0, 0, 1e-13, 5000, w, &tmp_D, &error); // a to infinity
		tmp_D_ff = f_f*tmp_D; // multiplied by fitting factor
	}
	else tmp_D = 0;

	gsl_integration_workspace_free(w); // release workspace
	return tmp_D_ff;
}

double dissociation_P_F(double a, void * params_scaled)
{
	struct dissociation_scaled_params *p = (struct dissociation_scaled_params *)params_scaled;
	// Load data from structure
	double E_scaled = p->E_scaled;
	double k_f_scaled = p->k_f_scaled;
	double alpha_scaled = p->alpha_scaled;

	// Calculate P_F
	double tmp_k_d_scaled = dissociation_rate(a, p);
	double tmp_P_F_scaled = probability(tmp_k_d_scaled, k_f_scaled)*gauss(a, alpha_scaled);
	return tmp_P_F_scaled;
}

double probability(double k_scaled, double k_f_scaled)
{
	double tmp_p;
	tmp_p = k_scaled / (k_scaled + k_f_scaled);
	return tmp_p;
}

double dissociation_rate(double a, dissociation_scaled_params *params_scaled)
{
	struct dissociation_scaled_params *p = (struct dissociation_scaled_params *)params_scaled;
	double k_f_scaled = p->k_f_scaled;
	double E_scaled = p->E_scaled;
	double u_n_scaled = p->u_n_scaled;
	double u_p_scaled = p->u_p_scaled;
	double tmp_k_d;

	tmp_k_d = (u_n_scaled + u_p_scaled)*pow(a, -3)*exp(-pow(a, -1))*bessel(E_scaled);
	return tmp_k_d;
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

double gauss(double a, double alpha_scaled)
{
	double tmp_Fa;

	tmp_Fa = pow(a, 2) * exp(alpha_scaled*pow(a, 2));
	return tmp_Fa;
}

double pre_exp(dissociation_params * params)
{
	struct dissociation_params *tmp_p = params;
	double E_r = tmp_p->E_r;
	double tmp_p_exp;
	tmp_p_exp = (3 * C_q) / (4 * M_PI*C_E_0*E_r);
	return tmp_p_exp;
}

double coulomb_radius(dissociation_params * params)
{
	struct dissociation_params *tmp_p = params;
	double E_r = tmp_p->E_r;
	double T = tmp_p->T;
	double tmp_r;

	tmp_r = pow(C_q, 2) / (4 * M_PI*C_E_0*E_r*C_k_B*T);
	return tmp_r;
}