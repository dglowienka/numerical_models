#ifndef _DISSOCIATION
#define _DISSOCIATION

#include <cmath>
#define C_k_B 1.381E-23 // Boltzmann constant [J/K]
#define C_q 1.602E-19 // Elementary charge [C]
#define C_E_0 8.854E-12 // Vacuum permittivity [F/m]

/* Structure with all the parameters needed for calculation of dissociation P(E): E, a_0, E_r, T, u_n, u_p, k_f
E - Electric field V/m
a_0 - Ion-pair separation distance distance [nm]
E_r - Relative permittivity [-]
T - Temperature [K]
u_n - mobility of electrons  m2 V^(-1) s^(-1)
u_p - mobility of holes  m2 V^(-1) s^(-1)
k_f - Decay rate [1/s]	*/
struct dissociation_params { double a_0; double E_r; double T; double u_n; double u_p; double k_f; };

// D(E,s) is calculated with integration
double dissociation(double E_scaled);

// The multiplication of P(E,a) and F(a) under integrate is calculated
double dissociation_P_F(double kd_scaled, void * params);

// The probability p(F)
double probability(double k, double k_f);

// Dissociation rate function dependent of applied electric field k_d(E)
double dissociation_rate(double A, double E_B, double b, double T);

// The preexponential factor
double pre_exp(double a, double u_n, double u_p, double E_r);

// Coulomb binding energy
double coulomb(double a, double E_r);

// Argument of Bessel function
double bessel(double E_scaled);

// Gaussian distribution
double gauss(double a, double a_0);

// Input function for scaling of E
double scaling(double E_unscaled, void * params);

// Output function for unscaling of Kd(E)
double unscaling(double kd_scaled, void * params);

#endif

