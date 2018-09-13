/*
 * MorletGabor.cpp
 *
 *  Created on: Sep 6, 2018
 *      Author: Travis Robson
 */



// Standard Library headers
#include <vector>
#include <cmath>
#include <complex>

using namespace std;

#include "MorletGabor.h"


complex<double> Psi_FT(double f, double A, double f0, double t0, double tau, double phiBar)
{
	complex<double> j(0,1.0); // imaginary number

	double phi0 = phiBar + 2*M_PI*f0*t0; // convert phiBar to phi0

	double arg1 = pow((f - f0)*M_PI*tau, 2.0);
	double arg2 = pow((f + f0)*M_PI*tau, 2.0);
	double arg3 = 2*M_PI*f*t0 + phi0;

	//complex<double> psi = exp(2.*phi0*j);

	//complex<double> psi = exp(-arg2-j*arg3) + exp(-arg1-j*arg3+2.*phi0*j);
	complex<double> psi = exp(-arg1-j*arg3+2.*phi0*j);
//	psi *= exp(-j*arg3);
	psi *= 0.5*sqrt(M_PI)*A*tau;

	return psi;
}

complex<double> Psi(double t, double A, double f0, double t0, double tau, double phiBar)
{
	complex<double> j(0, 1.0);

	double arg1 = 2*M_PI*f0*t + phiBar;
	double arg2 = pow((t - t0)/tau, 2.0);

	complex<double> psi = A*exp(arg1*j)*exp(-arg2);

	return psi;
}

