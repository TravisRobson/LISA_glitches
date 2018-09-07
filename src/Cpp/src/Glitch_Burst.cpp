//============================================================================
// Name        : Glitch_Burst.cpp
// Author      : Travis Robson
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <vector>
#include <cmath>

using namespace std;


#include "MorletGabor.h"
#include "Constants.h"
#include "Wavelet.h"
#include "LISA.h"
#include "TDI.h"

using namespace wv;
using namespace ls;
using namespace tdi;
using namespace mc;


int main(int argc, char **argv)
{
	auto start = chrono::system_clock::now();

	int i;

	double fi;
	complex <double> hi;

	// Find an observation period  greater than or equal to a week whose
	//		number of samples is a power of 2
	int N    = pow(2, ceil(log2(WEEK/dt)) );
	double T = N*dt;

	// Morlet-Gabor Wavelet parameters
	double A    = 1.0e-20;
	double f0   = 5.0e-3;
	double t0   = 0.5*T;
	double tau  = 5*HOUR;
	double phi0 = 1.2;

	cout << "Beginning Run" << "\n";
	cout << "T............. " << T/WEEK << "\n";

	// Store the Morlet-Gabor Wavelet Fourier Transform in a vector
//	vector<complex<double>> *h = new vector<complex<double>>(N/2);
//	for (i=0; i<N/2; i++)
//	{
//		fi = i/T;
//		h->at(i) = Psi_FT(fi, A, f0, t0, tau, phi0);
//	}
	// Clean up memory
//	delete h;

	// setup the LISA orbit
	LISA *lisa = new LISA(T);

	vector<double> paramsND = {log(A/1.0e-20), log(f0/mHz), log(t0/WEEK), log(tau/HOUR), phi0};
	Wavelet *wavelet = new Wavelet("glitch_OP12", paramsND);

	wavelet->calc_TDI(lisa);

	ofstream out_file; // Print the Glitch data out to a file
	complex<double> Ai, Ei, Ti;

	out_file.open("../Python/tests/Glitch_test.dat");
	for (i=0; i<wavelet->tdi.X.size(); i++)
	{
		fi = (i+wavelet->tdi.get_N_lo())/T;

		Ai = wavelet->tdi.A[i];
		Ei = wavelet->tdi.E[i];
		Ti = wavelet->tdi.T[i];

		out_file << scientific << setprecision(15) << fi << " " << real(Ai) << " " << imag(Ai) << " "
																<< real(Ei) << " " << imag(Ei) << " "
																<< real(Ti) << " " << imag(Ti) << "\n";
	}
	out_file.close();

	// calculate SNR and adjust it to the target SNR
	wavelet->set_snr(lisa);
	cout << "AET snr.......... " << wavelet->snr << endl;
	wavelet->adjust_snr(10.0, lisa);
	cout << "AET snr.......... " << wavelet->snr << endl;

	wavelet->set_Fisher(lisa);
	vector<double> e_vals;
	vector<vector<double>> e_vecs;
	for (i=0; i<wavelet->D; i++) vector<double> e_vecs[i];
	auto result = wavelet->get_EigenBS();
	e_vals = get<0>(result);
	e_vecs = get<1>(result);

	////////////////////////////////////////////////




	////////////////////////////////////////////////

	delete wavelet;
	delete lisa;

	// Find out how long the code took to run
	auto end = chrono::system_clock::now();
	auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
	cout << "Time elapsed..... " << (double)elapsed.count()/1000 << " seconds \n";

	return 0;
}
