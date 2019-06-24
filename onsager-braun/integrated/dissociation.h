#ifndef _DISSOCIATION
#define _DISSOCIATION

#include "constants.h"

/* Structure with all the parameters needed for calculation of dissociation P(E): E, a_0, E_r, T, u_n, u_p, k_f, f_f
E - Electric field [V/m]
a_0 - Ion-pair separation distance distance [nm]
E_r - Relative permittivity [-]
T - Temperature [K]
u_n - mobility of electrons  [m2 V^(-1) s^(-1)]
u_p - mobility of holes  [m2 V^(-1) s^(-1)]
k_f - Decay rate [1/s]	
f_f - fitting factor [-] */
struct dissociation_params { double E; double a_0; double E_r; double T; double u_n; double u_p; double k_f; double f_f; };

/* Structure with all the scaled parameters needed for calculation of dissociation P'(E'): E', alpha, k_f'
E_scaled - Electric field scaled by Ec = -1/2b' [-]
k_f - decay rate scaled by k_f_c = k_dc = A' / r^3 
u_n - mobility of electrons scaled by u_max
u_p - mobility of holes  scaled by u_max
f_f - fitting factor
alpha - unitless constant for F'(a'), alpha = r^2 / a_0^2 */
struct dissociation_scaled_params { double E_scaled; double k_f_scaled; double u_n_scaled; double u_p_scaled; double f_f; double alpha_scaled; };

// SCALED - D'(E',a') is calculated with integration
double dissociation(dissociation_scaled_params * params_scaled);

// SCALED - the multiplication of P'(E',a') and F'(a') under integrate is calculated
double dissociation_P_F(double a, void * params_scaled);

// SCALED - The probability p'(F')
double probability(double k_scaled, double k_f_scaled);

// SCALED - Dissociation rate function dependent of applied electric field k_d'(E')
double dissociation_rate(double a, dissociation_scaled_params *params_scaled);

// SCALED - Argument of Bessel function J1
double bessel(double E_scaled);

// SCALED - Gaussian distribution F'(a')
double gauss(double a, double alpha_scaled);

// The preexponential factor
double pre_exp(dissociation_params * params);

// Coulomb radius
double coulomb_radius(dissociation_params * params);

#endif

