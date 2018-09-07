/*
 * MorletGabor.h
 *
 *  Created on: Sep 6, 2018
 *      Author: travisrobson
 */



#ifndef MORLETGABOR_H_
#define MORLETGABOR_H_

using namespace std;

#include <complex>

complex<double> Psi_FT(double f, double A, double f0, double t0, double tau, double phiBar);
complex<double> Psi(double t, double A, double f0, double t0, double tau, double phiBar);


#endif /* MORLETGABOR_H_ */
