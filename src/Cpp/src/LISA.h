/*
 * LISA.h
 *
 *  Created on: Sep 6, 2018
 *      Author: Travis Robson
 */

#ifndef LISA_H_
#define LISA_H_

#include <vector>

using namespace std;

namespace ls {

#define fm 3.168753575e-8
#define fstar 0.01908538063694777
#define Clight 299792458.
#define Larm 2.5e9
#define ec 0.0048241852

#define KAPPA 0.0;
#define LAMBDA 0.0;

class LISA {
public:
	LISA(double T);
	virtual ~LISA();

	// Noise methods
	double OMS_Noise(double f);
	double ACC_Noise(double f);
	double SnAE(double f);
	double SnT(double f);
	double SnX(double f);

	void SC_position_analytic(double t, vector<double> *x, vector<double> *y, vector<double> *z);

	// getters
	double get_T();

	double T; // observation period


private:

};

}

#endif /* LISA_H_ */
