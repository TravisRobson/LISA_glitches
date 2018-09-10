/*
 * LISA.cpp
 *
 *  Created on: Sep 6, 2018
 *      Author: Travis Robson
 */

#include "LISA.h"
#include "Constants.h"

#include <cmath>
#include <vector>

using namespace std;

namespace ls {

LISA::LISA(double T) {
	this->T = T;
}

LISA::~LISA() {

}

double LISA::get_T()
{
	return this->T;
}

double LISA::OMS_Noise(double f)
{	// single-link optical metrology noise
	double Poms = pow(1.5e-11, 2)*(1 + pow(2e-3/f, 4)); // Hz^{-1}

	return Poms;
}

double LISA::ACC_Noise(double f)
{	// single test mass acceleration noise power
	double Pacc = pow(3.0e-15, 2)*(1 + pow(0.4e-3/f, 2))*(1 + pow(f/8e-3, 4)); // m s^{-2} Hz^{-1}

	return Pacc;
}

double LISA::SnAE(double f)
{	// LISA instrument noise for the first generation TDI A, E channels
	double Poms = OMS_Noise(f);
	double Pacc = ACC_Noise(f);

	double fp = f/fstar;

	double transfer = pow(sin(fp), 2);

	double SnAE = 16/3*transfer*(Poms*(2 + cos(fp)) + 2*Pacc/pow(2*M_PI*f, 4)*(3 + 2*cos(fp) + cos(2*fp)))/pow(2*Larm, 2);

	return SnAE;
}

double LISA::SnT(double f)
{	// LISA instrument noise for the first generation TDI T channel
	double Poms = OMS_Noise(f);
	double Pacc = ACC_Noise(f);

	double fp = f/fstar;

	double transfer = pow(sin(fp), 2);

	double SnT = 16/3*transfer*(Poms*(1 - cos(fp)) + 2*Pacc/pow(2*M_PI*f, 4)*(3 - 4*cos(fp) + cos(2*fp)))/pow(2*Larm, 2);

	return SnT;
}

double LISA::SnX(double f)
{	// LISA instrument noise for the first generation TDI X, Y, Z channels
	double Poms = OMS_Noise(f);
	double Pacc = ACC_Noise(f);

	double fp = f/fstar;

	double transfer = pow(sin(fp), 2);

	double SnXYZ = 4*transfer*(4*Poms + 8*Pacc/pow(2*M_PI*f, 4)*(1 + pow(cos(fp), 2.0)))/pow(2*Larm, 2);

	return SnXYZ;
}

void LISA::SC_position_analytic(double t, vector<double> x, vector<double> y, vector<double> z)
{
	double alpha = 2*M_PI*fm*t + KAPPA;

	double beta1, beta2, beta3;
	beta1 = 0.0 + LAMBDA;
	beta2 = 2*M_PI/3 + LAMBDA;
	beta3 = 4*M_PI/3 + LAMBDA;

	double sa = sin(alpha);
	double ca = cos(alpha);

	double sb = sin(beta1);
	double cb = cos(beta1);
	x[0] = AU*(ca + ec*(sa*ca*sb - (1. + sa*sa)*cb));
	y[0] = AU*(sa + ec*(sa*ca*cb - (1. + ca*ca)*sb));
	z[0] = -sqrt(3.)*AU*ec*(ca*cb + sa*sb);

	sb = sin(beta2);
	cb = cos(beta2);
	x[0] = AU*(ca + ec*(sa*ca*sb - (1. + sa*sa)*cb));
	y[0] = AU*(sa + ec*(sa*ca*cb - (1. + ca*ca)*sb));
	z[0] = -sqrt(3.)*AU*ec*(ca*cb + sa*sb);

	sb = sin(beta3);
	cb = cos(beta3);
	x[0] = AU*(ca + ec*(sa*ca*sb - (1. + sa*sa)*cb));
	y[0] = AU*(sa + ec*(sa*ca*cb - (1. + ca*ca)*sb));
	z[0] = -sqrt(3.)*AU*ec*(ca*cb + sa*sb);
}



}
