/*
 * Wavelet.h
 *
 *  Created on: Aug 28, 2018
 *      Author: Travis Robson
 */

#ifndef WAVELET_H_
#define WAVELET_H_

#include <string>
#include <vector>

using namespace std;

#include "TDI.h"
#include "LISA.h"
using namespace tdi;
using namespace ls;


namespace wv {

#define IDX_A    0
#define IDX_f0   1
#define IDX_t0   2
#define IDX_tau  3
#define IDX_phi0 4


class Wavelet {

public:
	Wavelet(string, vector<double>);
	virtual ~Wavelet();
	Wavelet(const Wavelet &src, int flag); // copy constructor

	void calc_TDI(LISA *lisa); // construct the TDI for a given wavelet
	void calc_OP12_TDI(double T);
	void set_snr(LISA *lisa);
	void adjust_snr(double snr_target, LISA *lisa);
	void set_Fisher(LISA *lisa);
	tuple <vector<double>,vector<vector<double>>> get_EigenBS();

	TDI tdi;
	double snr;
	vector<vector<double>> Fisher;
	int D;

private:
	string name;
	vector<double> paramsND;

};

void unpack_glitch_params(vector<double> paramsND, double *A, double *f0, double *t0, double *tau, double *phi0);



} /* namespace wv */

#endif /* WAVELET_H_ */
