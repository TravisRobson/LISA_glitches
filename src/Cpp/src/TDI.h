/*
 * TDI.h
 *
 *  Created on: Sep 6, 2018
 *      Author: travisrobson
 */

#ifndef TDI_H_
#define TDI_H_

#include <vector>
#include <complex>
#include <list>
using namespace std;

#include "LISA.h"
using namespace ls;

namespace tdi {

class TDI {
public:
	TDI();
	virtual ~TDI();
	TDI(const TDI &src);

	TDI operator-(const TDI& a);
	TDI operator/=(double num);

	void set_N_lo(int);
	void set_N_hi(int);
	int  get_N_lo() const { return N_lo; }
	int  get_N_hi() const { return N_hi; }

	void phase_to_tdi(list<vector<complex<double>>> *phase_list, int N_lo, double T);

	vector<complex<double>> X, Y, Z;
	vector<complex<double>> A, E, T;

private:
	//double f_lo, f_hi;
	int N_lo, N_hi;

};

double nwip(TDI *a, TDI *b, LISA *lisa);

} /* namespace tdi */

#endif /* TDI_H_ */
