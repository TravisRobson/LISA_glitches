/*
 * LISA.h
 *
 *  Created on: Sep 6, 2018
 *      Author: Travis Robson
 */

#ifndef LISA_H_
#define LISA_H_

namespace ls {

#define fm 3.168753575e-8
#define fstar 0.01908538063694777
#define Clight 299792458.
#define L 2.5e9

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

	// getters
	double get_T();

private:
	double T; // observation period

};

}

#endif /* LISA_H_ */
