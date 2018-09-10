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
#include <string>
#include <vector>
#include <cmath>
#include <list>
#include <map>

using namespace std;

#include "ProgressBar.hpp"
#include "cxxopts.hpp"

#include "MorletGabor.h"
#include "Constants.h"
#include "G_B_MCMC.h"
#include "Wavelet.h"
#include "LISA.h"
#include "TDI.h"

using namespace wv;
using namespace ls;
using namespace tdi;
using namespace mc;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


int main(int argc, char **argv)
{
	auto start = chrono::system_clock::now();

	int i, j, k;

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

	// setup the LISA orbit
	LISA *lisa = new LISA(T);

	vector<double> paramsND = {log(A/A_scale), f0/f0_scale, t0/t0_scale, log(tau/tau_scale), phi0};
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
	cout << "initial AET snr.......... " << wavelet->snr << endl;
	wavelet->adjust_snr(10.0, lisa);
	cout << "adjusted AET snr......... " << wavelet->snr << endl;

	wavelet->set_Fisher(lisa);
	vector<double> e_vals;
	vector<vector<double>> e_vecs;
	for (i=0; i<wavelet->D; i++) vector<double> e_vecs[i];
	auto result = wavelet->get_EigenBS();
	e_vals = get<0>(result);
	e_vecs = get<1>(result);

	////////////////////////////////////////////////
	Model modelX0 = Model(*wavelet);
	modelX0.set_logL(wavelet->tdi, lisa);
	cout << "modelX0 logL............. " << modelX0.logL << endl;

	// setup random number generator
	gsl_rng_env_setup();
	int seed = 1;
	const gsl_rng_type *TT = gsl_rng_default;
	gsl_rng *r = gsl_rng_alloc (TT);
	gsl_rng_set(r, seed);

	/////// Setup command line arguments
	int N_MCMC, N_BURN, N_undersample;
	int PT=1;
	cxxopts::Options options("MyProgram", "One line description of MyProgram");
	options.add_options()
		("n,N_MCMC", "Number of MCMC samples", cxxopts::value<int>(N_MCMC))
		("b,N_BURN", "Number of Burn in samples", cxxopts::value<int>(N_BURN))
		("u,N_undersample", "Number to undersample by", cxxopts::value<int>(N_undersample))
		("p,PTMCMC", "Number to undersample by", cxxopts::value<int>(PT))
		;
	auto opts = options.parse(argc, argv);

	//////////////
	int DB = 0; // detailed balance test flag
	// setup temperature ladder
	int N_Temps = 1;
	double gam = 1.2;
	double T_max;
	if (PT == 1)
	{
		double T0 = 1; // initial temperature
		double snr_eff = 2.0;
		T_max = pow(modelX0.wave.snr, 2.0)/pow(snr_eff, 2.0);
		while (T0 < T_max)
		{
			T0 *= gam;
			N_Temps += 1;
		}
	}
	cout << "Number of Temperature....... " << N_Temps << endl;
	vector<double> Temps(N_Temps);
	Temps[0] = 1.0;
	if (PT == 1)
	{
		for (i=1; i<N_Temps; i++) Temps[i] = gam*Temps[i-1];
		Temps[N_Temps-1] = T_max;
	}

	// setup e_vals_list and e_vecs_list
	vector<vector<double>> e_vals_list;
	for (i=0; i<N_Temps; i++) e_vals_list.push_back(e_vals);
	vector< vector<vector<double>> > e_vecs_list;
	for (i=0; i<N_Temps; i++) e_vecs_list.push_back(e_vecs);

	/////////// set up Prior and proposal distributions /////////
	Proposal_Fisher P_Fisher = Proposal_Fisher (0.75, r, Temps, e_vals_list, e_vecs_list);
	Prior prior = Prior (0.05, r, T);
	int hist_stride = 10;
	typedef vector<double> D1;
	typedef vector<D1> D2;
	typedef vector<D2> matrix_3D;
	matrix_3D *history = new matrix_3D(  N_MCMC/hist_stride,D2(N_Temps,D1(wavelet->D,0))  );

	Proposal_DE P_DE = Proposal_DE(0.2, r, history, hist_stride, Temps);

	struct Proposal *prop = new struct Proposal;
	prop->P_Fisher = P_Fisher;
	prop->prior = prior;
	prop->P_DE = P_DE;

	// acceptance counters
	vector<int> acc(N_Temps,0);
	vector<int> acc_cnt(N_Temps,0);
	vector<int> swap(N_Temps-1,0);
	vector<int> swap_cnt(N_Temps-1,0);
	vector<int> ID_list(N_Temps);
	for (i=0; i<N_Temps; i++) ID_list[i] = i;
	int who;

	double logH, u; // hastings ratio and decision uniform rvs
	int whoa, whob; // identifiers for PTMCMC swapping
	double Ta, Tb; // Temperatures for PRMCMC swapping

	// create models
	Model modelY = Model(modelX0);
	Model modelX = Model(modelX0);

	vector<Model*> modelX_list(N_Temps);
	for (i=0; i<N_Temps; i++)
	{
		Model *m = new Model(); // MUST INIT with empty ones... or default constructor more likely
		modelX_list.push_back(m);
	}
	if (DB==1) modelX0.logL = 1;
	for (i=0; i<N_Temps; i++)
	{
		modelX_list[i] = new Model(modelX0); // now I can copy
	}

	ofstream logL_file; // Print the Glitch data out to a file
	ofstream T1_file, THot_file;

	ProgressBar progressBar(N_MCMC + N_BURN, 70); // begin status bar
	out_file.open("../Python/tests/MCMC_Glitch_test.dat"); // open output file
	T1_file.open("../Python/tests/MCMC_Glitch_T1.dat"); // open output file
	THot_file.open("../Python/tests/MCMC_Glitch_THot.dat"); // open output file
	logL_file.open("../Python/tests/MCMC_Glitch_logL.dat"); // open output file
	for (i=-N_BURN; i<N_MCMC; i++)
	{
		++progressBar;
		for  (k=0; k<N_undersample; k++)
		{
			for(j=0; j<N_Temps; j++)
			{
				who = ID_list[j];
				acc_cnt[j]++;

				modelY.wave.paramsND = propose_params(prop, modelX_list[who]->wave.paramsND, Temps[j]);
				modelY.wave.Unwrap_Phase();

				logH = prior.logpdf(modelY.wave.paramsND);

				if (logH > -INFINITY)
				{
					logH -= prior.logpdf(modelX_list.at(who)->wave.paramsND);

					if (DB == 0)
					{
						modelY.wave.calc_TDI(lisa);
						modelY.wave.set_snr(lisa);
						modelY.set_logL(wavelet->tdi, lisa);
					}
					else modelY.logL = 1.0;

					logH += (modelY.logL - modelX_list.at(who)->logL)/Temps[j];

					u = log(gsl_ran_flat(r, 0, 1.0));
					if (u < logH)
					{
						acc[j]++;
						*modelX_list[who] = Model(modelY);
					}
				}
			}
		}
		progressBar.display();

		///////// Perform PTMCMC /////////
		for (j=0; j<N_Temps-1; j++)
		{
			swap_cnt[j]++;
			whoa = ID_list[j];
			whob = ID_list[j+1];

			Ta = Temps[j];
			Tb = Temps[j+1];

			logH  = (modelX_list.at(whoa)->logL - modelX_list.at(whob)->logL)/Tb;
			logH -= (modelX_list.at(whoa)->logL - modelX_list.at(whob)->logL)/Ta;

			u = log(gsl_ran_flat(r, 0, 1.0));
			if (u < logH)
			{
				swap[j]++;
				ID_list[j] = whob;
				ID_list[j+1] = whoa;
			}
		}
//		Temps = adapt_Temps(Temps, swap, swap_cnt, i + N_BURN);
//		prop->P_Fisher.Temps = Temps;
//		prop->P_DE.Temps = Temps;

		// update Fisher matrix proposal distribution
		if (i%20 == 0)
		{
			for (j=0; j<N_Temps; j++)
			{
				who = ID_list[j];
				modelX_list[who]->wave.set_Fisher(lisa);
				for (k=0; k<modelX_list[who]->wave.D; k++) vector<double> e_vecs[k];
				auto result = modelX_list[who]->wave.get_EigenBS();
				e_vals = get<0>(result);
				e_vecs = get<1>(result);

				prop->P_Fisher.e_vals_list[j] = e_vals;
				prop->P_Fisher.e_vecs_list[j] = e_vecs;
			}
		}


		// update the history tensor
		int loc = (i+N_BURN)%hist_stride;
		for(j=0; j<N_Temps; j++)
		{
			for(k=0; k<wavelet->D; k++)
			{
				history->at(loc).at(j).at(k) = modelX_list[ID_list[j]]->wave.paramsND[k];
			}
		}
		prop->P_DE.cnt++;

		// Store the cold chain
		who = ID_list[0]; // only store the cold chain
		out_file << scientific << setprecision(15)
							   << modelX_list.at(who)->logL << " "
							   << modelX_list.at(who)->wave.paramsND[IDX_A]    << " "
							   << modelX_list.at(who)->wave.paramsND[IDX_f0]   << " "
							   << modelX_list.at(who)->wave.paramsND[IDX_t0]   << " "
							   << modelX_list.at(who)->wave.paramsND[IDX_tau]  << " "
							   << modelX_list.at(who)->wave.paramsND[IDX_phi0] << endl;

		who = ID_list[1]; // only store the cold chain
		T1_file << scientific << setprecision(15)
							   << modelX_list.at(who)->logL << " "
							   << modelX_list.at(who)->wave.paramsND[IDX_A]    << " "
							   << modelX_list.at(who)->wave.paramsND[IDX_f0]   << " "
							   << modelX_list.at(who)->wave.paramsND[IDX_t0]   << " "
							   << modelX_list.at(who)->wave.paramsND[IDX_tau]  << " "
							   << modelX_list.at(who)->wave.paramsND[IDX_phi0] << endl;

		who = ID_list[N_Temps-1]; // only store the cold chain
		THot_file << scientific << setprecision(15)
							   << modelX_list.at(who)->logL << " "
							   << modelX_list.at(who)->wave.paramsND[IDX_A]    << " "
							   << modelX_list.at(who)->wave.paramsND[IDX_f0]   << " "
							   << modelX_list.at(who)->wave.paramsND[IDX_t0]   << " "
							   << modelX_list.at(who)->wave.paramsND[IDX_tau]  << " "
							   << modelX_list.at(who)->wave.paramsND[IDX_phi0] << endl;

		// Store the logLs
		for (j=0; j<N_Temps; j++) logL_file << scientific << setprecision(15) << modelX_list.at(ID_list[j])->logL << " ";
		logL_file << endl;

	}
	out_file.close();
	logL_file.close();
	T1_file.close();
	THot_file.close();
	progressBar.done();

	cout << endl;
	for (j=0; j<N_Temps; j++)
	{
		cout << "Acceptance rate[" << j << "]............ " << 100.*acc[j]/acc_cnt[j] << " %" << endl;
	}
	cout << endl;
	for (j=0; j<N_Temps-1; j++)
	{
		cout << "Swap Acceptance rate[" << j << "<->" << j+1 << "]....... " << 100.*acc[j]/acc_cnt[j] << " %" << endl;
	}
	cout << endl;
	for (j=0; j<N_Temps; j++) cout << "T[" << j << "]...... " << Temps[j] << endl;
	cout << endl;
	////////////////////////////////////////////////

	delete wavelet;
	delete lisa;

	for (j=0; j<N_Temps; j++) delete modelX_list[j];
	delete prop;

	// Find out how long the code took to run
	auto end = chrono::system_clock::now();
	auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
	cout << "Time elapsed..... " << (double)elapsed.count()/1000 << " seconds \n";

	return 0;
}
