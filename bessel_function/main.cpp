#define _USE_MATH_DEFINES // Use constants in cmath
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_complex.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <iomanip> //increase precision of fstream

#define C_k_B 1.381E-23 // Boltzmann constant [J/K]
#define C_q 1.602E-19 // Elementary charge [C]
#define C_E_0 8.854E-12 // Vacuum permittivity [F/m]

using namespace std;

int main(void)
{
	ofstream Bessel;
	Bessel.open("bessel1.txt");

	//double x_axis[ARRRAY_SIZE];
	//double y_axis[ARRRAY_SIZE];

	/*for (int i = 0; i < ARRRAY_SIZE; i++)
	{
		x_axis[i] = i*0.1;
		y_axis[i] = gsl_sf_bessel_J1(x_axis[i]);
		cout << x_axis[i] << "\t" << y_axis[i] << "\n";
		Bessel << x_axis[i] << "\t" << y_axis[i] << "\n";
	}*/
	double E[2000];
	E[0] = 10;
	double tmp_b[1000];
	double tmp_B[1000];
	double tmp_B_1[1000];
	double tmp_B_inside;
	double E_r = 3;
	double T = 300;
	Bessel << "E" << "\t" << "b" << "\t" << "B_GSL" << "\t" << "B_cmath" << "\n";
	Bessel << std::fixed << std::setprecision(14);
	for (int i = 0; i < 1001; i++)
	{
		E[i+1] = E[i]+1e6;
		tmp_b[i] = (pow(C_q, 3)*E[i]) / (8 * M_PI*C_E_0*E_r*pow(C_k_B, 2) * pow(T, 2));
		const std::complex<double> result = std::sqrt(std::complex<double>(-2 * tmp_b[i]));
		tmp_B_inside = std::imag(result);
		tmp_B[i] = gsl_sf_bessel_I1(2 * tmp_B_inside) / tmp_B_inside;;
		tmp_B_1[i] = _j1(2 * tmp_B_inside) / tmp_B_inside;;
		Bessel << E[i] << "\t" << tmp_b[i] << "\t" << tmp_B[i] << "\t" << tmp_B_1[i] << "\n";
		cout << E[i] << "\t" << tmp_b[i] << "\t" << tmp_B[i] << "\t" << tmp_B_1[i] << "\n";
	}
	Bessel << "end file";
	Bessel.close();
	return 0;
}
