/*
 * LISA.cpp
 *
 *  Created on: Sep 6, 2018
 *      Author: Travis Robson
 */

#include "LISA.h"

#include <cmath>

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

	double SnAE = 16/3*transfer*(Poms*(2 + cos(fp)) + 2*Pacc/pow(2*M_PI*f, 4)*(3 + 2*cos(fp) + cos(2*fp)))/pow(2*L, 2);

	return SnAE;
}

double LISA::SnT(double f)
{	// LISA instrument noise for the first generation TDI T channel
	double Poms = OMS_Noise(f);
	double Pacc = ACC_Noise(f);

	double fp = f/fstar;

	double transfer = pow(sin(fp), 2);

	double SnT = 16/3*transfer*(Poms*(1 - cos(fp)) + 2*Pacc/pow(2*M_PI*f, 4)*(3 - 4*cos(fp) + cos(2*fp)))/pow(2*L, 2);

	return SnT;
}

double LISA::SnX(double f)
{	// LISA instrument noise for the first generation TDI X, Y, Z channels

	return 1.0;
}


}
