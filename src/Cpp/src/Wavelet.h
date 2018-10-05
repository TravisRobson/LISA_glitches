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

#include "Constants.h"
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

#define IDX_cos_theta 5
#define IDX_phi 6
#define IDX_psi 7
#define IDX_ellip 8

#define A_scale 1.0e-20
#define f0_scale 1.0e-3
#define t0_scale WEEK
#define tau_scale HOUR


class Wavelet {

public:
	Wavelet();
	Wavelet(string, vector<double>);
	virtual ~Wavelet();
	Wavelet(const Wavelet &src, int flag); // copy constructor

	void calc_TDI(LISA *lisa); // construct the TDI for a given wavelet

	void calc_burst_TDI(LISA *lisa);

	void calc_OP12_TDI(LISA *lisa);
	void calc_OP21_TDI(LISA *lisa);
	void calc_OP13_TDI(LISA *lisa);
	void calc_OP31_TDI(LISA *lisa);
	void calc_OP23_TDI(LISA *lisa);
	void calc_OP32_TDI(LISA *lisa);

	void calc_AC12_TDI(LISA *lisa);
	void calc_AC21_TDI(LISA *lisa);
	void calc_AC13_TDI(LISA *lisa);
	void calc_AC31_TDI(LISA *lisa);
	void calc_AC23_TDI(LISA *lisa);
	void calc_AC32_TDI(LISA *lisa);

	void set_snr(LISA *lisa, int X_flag);
	void adjust_snr(double snr_target, LISA *lisa, int X_flag);
	void set_Fisher(LISA *lisa, int X_flag);
	tuple <vector<double>,vector<vector<double>>> get_EigenBS();

	void Unwrap_Phase();

	string  get_name() const { return name; }

	TDI tdi;
	double snr;
	vector<vector<double>> Fisher;
	vector<double> paramsND;
	int D;

private:
	string name;


};

void unpack_glitch_params(vector<double> paramsND, double *A, double *f0, double *t0, double *tau, double *phi0);
void unpack_burst_params(vector<double> paramsND, double *cos_theta, double *phi, double *psi, double *ellip);


} /* namespace wv */

#endif /* WAVELET_H_ */
