//============================================================================
// Name        : Glitch_Burst.cpp
// Author      : Travis Robson
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//
//
// -ffast-math -ftree-vectorize
//
//============================================================================

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

#include "ProgressBar.hpp"
#include "cxxopts.hpp"

#include "MorletGabor.h"
#include "Constants.h"
#include "G_B_MCMC.h"
#include "Wavelet.h"
#include "LISA.h"
#include "TDI.h"
#include "IO.h"

using namespace wv;
using namespace ls;
using namespace tdi;
using namespace mc;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

int main(int argc, char **argv)
{
	auto start = chrono::system_clock::now(); // time the duration of this program

	// Find an observation period  greater than or equal to a week whose
	//		number of samples is a power of 2
	int N    = pow(2, ceil(log2(WEEK/dt)) );
	double T = N*dt;
	cout << "Beginning Run" << "\n";
	cout << "T............. " << T/WEEK << " weeks\n";

	ofstream logL_file, out_file;
	ofstream T1_file, THot_file, ID_file;

	////////////////// Setup command line arguments ///////////////////
	int N_MCMC, N_BURN, N_undersample, seed;
	int PT = 0;   // PTMCMC
	int DB = 0;  // detailed balance test flag
	int flag_adaptive_temps = 0;
	int X_flag;

	string input_file;

	cxxopts::Options options("GlitchBurstMCMC", "Calculate model posteriors via MCMC, evidence via PTMCMC");
	options.add_options()
		("n,nmcmc",      "Number of MCMC samples", cxxopts::value<int>(N_MCMC))
		("b,nburn",      "Number of Burn in samples", cxxopts::value<int>(N_BURN))
		("u,nunder",     "Number to undersample by", cxxopts::value<int>(N_undersample)->default_value("10"))
		("p,ptmcmc",     "Number to undersample by")
		("d, detbail",   "Perform detailed balance test")
		("a, adaptt",    "Adaptive temperature ladder")
		("f, inputfile", "Model Data setup file", cxxopts::value<string>(input_file))
		("s, seed",      "random number generator seed", cxxopts::value<int>(seed)->default_value("1"))
		("x, xonly",      "X channel only analysis", cxxopts::value<int>(X_flag)->default_value("0"))
		;
	auto opts = options.parse(argc, argv);

	if (opts.count("detbail")) DB = 1;
	if (opts.count("adaptt")) flag_adaptive_temps = 1;
	if (opts.count("ptmcmc")) PT = 1;

	/////////////////////////////////////////////////////////////////////

	int i, j, k;

	//complex <double> hi;


	struct Files *Files= new struct Files;
	LISA *lisa = new LISA(T); // setup the LISA orbit
	//Wavelet *wavelet = parse_config_file(input_file, T, lisa, Files);
	vector<Model*> models;
	models = parse_config_file(input_file, &T, lisa, Files, X_flag);
	Model *modelX0 = models.at(0);
	Model *model_true = models.at(1);
	cout << "T............. " << T/WEEK << " weeks\n";

//	Model modelX0 = Model(*wavelet);
//	modelX0.set_logL(wavelet->tdi, lisa);
	cout << "modelX0 logL............. " << modelX0->logL << endl;
	// setup temperature ladder
	int N_Temps = 1;
	double gam = 1.3;
	double T_max;
	if (PT == 1)
	{
		double T0 = 1; // initial temperature
		double snr_eff = 0.5;
		T_max = pow(modelX0->wave.snr, 2.0)/pow(snr_eff, 2.0);
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

	// Print temperatures (betas) to file
	logL_file.open (Files->File_logL);
	logL_file << N_Temps;
	for (j=0; j<N_Temps; j++) logL_file << scientific << setprecision(15) << " " << 1/Temps[j];
	logL_file << endl;

	modelX0->wave.set_Fisher(lisa, X_flag);
	vector<double> e_vals;
	vector<vector<double>> e_vecs;
	for (i=0; i<modelX0->wave.D; i++) vector<double> e_vecs[i];
	auto result = modelX0->wave.get_EigenBS();
	e_vals = get<0>(result);
	e_vecs = get<1>(result);
	// setup e_vals_list and e_vecs_list
	vector<vector<double>> e_vals_list;
	for (i=0; i<N_Temps; i++) e_vals_list.push_back(e_vals);
	vector< vector<vector<double>> > e_vecs_list;
	for (i=0; i<N_Temps; i++) e_vecs_list.push_back(e_vecs);

	////////////////////////////////////////////////


	// setup random number generator
	gsl_rng_env_setup();
	const gsl_rng_type *TT = gsl_rng_default;
	gsl_rng *r = gsl_rng_alloc (TT);
	gsl_rng_set(r, seed);

	////////////////// setup the infrastructure for the adaptive temperature ladder /////////////////////
	vector<double> n_up(N_Temps,0), n_down(N_Temps,0);
	vector<double> up_down_tracker(N_Temps,0);

	up_down_tracker[0] = 1; // i.e. positive, on its way up, visited cold most recently
	up_down_tracker[N_Temps-1] = -1; // negative, on its way down, visited hot most recently

	vector<double> f(N_Temps,0), f_smooth(N_Temps,0);
	for (j=0; j<N_Temps; j++) f_smooth[j] = 1 - (double)j/(N_Temps-1);
	double weight = 0.95; // how much to balance measured f by the desired f (f_smooth)
	double dir;
	int t = 1000;
	/////////////////////////////////////////////////////////////////////////////////////////////////////


	////////////////////// set up Prior and proposal distributions //////////////////////////
	double weight_prior     = 0.05;
	double weight_DE        = 0.5;
	double weight_TimeShift = 0.05;
	double weight_Target    = 0.1;
	double weight_Fisher    = 1 - weight_prior - weight_DE - weight_TimeShift - weight_Target;

	Prior prior = Prior (weight_prior, r, T);

	int hist_stride = 100;
	typedef vector<double> D1;
	typedef vector<D1> D2;
	typedef vector<D2> matrix_3D;
	matrix_3D *history = new matrix_3D(hist_stride, D2(N_Temps, D1(modelX0->wave.D,0)) );
	Proposal_DE P_DE = Proposal_DE(weight_DE, r, history, hist_stride, Temps);

	Proposal_Target P_Target = Proposal_Target(weight_Target, r, modelX0->wave.paramsND[IDX_f0],
												modelX0->wave.paramsND[IDX_tau],
												3./sqrt(modelX0->wave.Fisher[IDX_f0][IDX_f0]),
												3./sqrt(modelX0->wave.Fisher[IDX_tau][IDX_tau]));

	Proposal_TimeShift P_TimeShift = Proposal_TimeShift(weight_TimeShift, r);

	Proposal_Fisher P_Fisher = Proposal_Fisher (weight_Fisher, r, Temps, e_vals_list, e_vecs_list);

	struct Proposal *prop = new struct Proposal;
	prop->P_Fisher = P_Fisher;
	prop->prior = prior;
	prop->P_DE = P_DE;
	prop->P_Target = P_Target;
	prop->P_TimeShift = P_TimeShift;
	///////////////////////////////////////////////////////////////////////////////

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
	Model modelY = Model(*modelX0);
	Model modelX = Model(*modelX0);

	vector<Model*> modelX_list(N_Temps);
	for (i=0; i<N_Temps; i++)
	{
		Model *m = new Model(); // MUST INIT with empty ones... or default constructor more likely
		modelX_list.push_back(m);
		delete m;
	}
	if (DB==1) modelX0->logL = 1;
	for (i=0; i<N_Temps; i++)
	{
		modelX_list[i] = new Model(*modelX0); // now I can copy
	}

	double logQyx, logQxy;

	tuple<vector<double>, double, double, string> out;
	string name;


	int loc;
	ProgressBar progressBar(N_MCMC + N_BURN, 70);  // begin status bar

	out_file.open  (Files->File_cold_chain);
	T1_file.open   (Files->File_T1_chain);
	THot_file.open (Files->File_Hot_chain);
	ID_file.open   (Files->File_IDs);

	for (i=-N_BURN; i<N_MCMC; i++)
	{
		++progressBar;
		for  (k=0; k<N_undersample; k++)
		{
			for(j=0; j<N_Temps; j++)
			{
				who = ID_list[j];
				acc_cnt[j]++;

				out = propose_params(prop, modelX_list[who]->wave.paramsND, Temps[j]);
				modelY.wave.paramsND = get<0>(out);
				logQyx = get<1>(out);
				logQxy = get<2>(out);
				name = get<3>(out);

				modelY.wave.Unwrap_Phase();

				logH = prior.logpdf(modelY.wave.paramsND);

				if (logH > -INFINITY)
				{
					logH -= prior.logpdf(modelX_list.at(who)->wave.paramsND);

					if (DB == 0)
					{
						modelY.wave.calc_TDI(lisa);
						modelY.wave.set_snr(lisa, X_flag);
						//cout << Temps[j] << " before " << i << endl;
						modelY.set_logL(model_true->wave.tdi, lisa, X_flag);
						//cout << Temps[j] << " after " << i << endl << endl;
//						if (i==22)
//						{
//							cout << i << endl;
//							cout << "logL.............. " << modelY.logL << endl;
//							cout << modelY.wave.tdi.get_N_hi() - modelY.wave.tdi.get_N_lo() << endl << endl;
//							for (int h =0;h<modelY.wave.paramsND.size(); h++) cout << "p[" << h << "]  " << modelY.wave.paramsND[h] << endl;
//						}

					}
					else modelY.logL = 1.0;

					//if (modelY.logL != modelY.logL)
					logH += (modelY.logL - modelX_list.at(who)->logL)/Temps[j];

					u = log(gsl_ran_flat(r, 0, 1.0));
					if (u < logH)
					{
						acc[j]++;
						*modelX_list[who] = Model(modelY);

						if (name == prop->P_Fisher.name) prop->P_Fisher.acc++;
						else if (name == prop->prior.name) prop->prior.acc++;
						else if (name == prop->P_DE.name) prop->P_DE.acc++;
						else if (name == prop->P_Target.name) prop->P_Target.acc++;
						else if (name == prop->P_TimeShift.name) prop->P_TimeShift.acc++;
					}
				}
			}
		}
		if ((i+N_BURN)%50 == 0) progressBar.display();
		if (i==0 & PT == 1)
		{
			for (j=0; j<N_Temps-1; j++)
			{
				swap_cnt[j] = 0;
				swap[j] = 0;
			}
		}

		///////// Perform PTMCMC /////////
		if (PT == 1)
		{
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
					ID_list[j]   = whob;
					ID_list[j+1] = whoa;
				}
			}
		}
		/////////////// Update Temperature ladder //////////////////////
		for (j=0; j<N_Temps; j++)
		{
			who = ID_list[j];
			if (j == 0) up_down_tracker[who] = 1;
			else if (j==N_Temps-1) up_down_tracker[who] = -1;
		}

		for (j=0; j<N_Temps; j++)
		{
			who = ID_list[j];
			dir = up_down_tracker[who];
			if (dir == 1) n_up[j]++;
			else if (dir==-1) n_down[j]++;
		}
		if (flag_adaptive_temps==1 and i<0)
		{
//			if (i%t == 0 & i>=-N_BURN+t)
//			{
//				for (j=0; j<N_Temps; j++) f[j] = weight*(double)n_up[j]/(n_up[j] + n_down[j]) + (1-weight)*f_smooth[j];
//			    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
//			    gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, N_Temps);
//			    double f_copy[N_Temps];
//			    double Temps_copy[N_Temps];
//			    double result[N_Temps];
//			    for (j=0; j<N_Temps; j++)
//			    {	// copy to C containers...
//			    	f_copy[j] = f[N_Temps-1-j]; // reverse the order
//			    	Temps_copy[j] = Temps[N_Temps-1-j]; // reverse the order
//			    }
//			    gsl_spline_init (spline, f_copy, Temps_copy, N_Temps);
//			    for (j=1; j<N_Temps-1; j++)
//			    {
//			    	result[j] = gsl_spline_eval (spline, f_smooth[j], acc);
//			    }
//			    gsl_spline_free (spline);
//			    gsl_interp_accel_free (acc);
//
//			    for (j=1; j<N_Temps-1; j++) Temps[j] = result[j];
//
//			    if (i<0)
//			    {
//					for (j=0; j<N_Temps; j++)
//					{	// reset counters
//						n_up[j] = 0;
//						n_down[j] = 0;
//					}
//			    }
//			    //t *= 2;
//			}

			//Temps = adapt_Temps(Temps, swap, swap_cnt, i + N_BURN);
			prop->P_Fisher.Temps = Temps;
			prop->P_DE.Temps = Temps;
			//if ((i+N_BURN)%10 == 0) for (j=0; j<N_Temps-1; j++) cout << "swap rate[" << j << "]....... " << (double)swap[j]/swap_cnt[j] << endl;
//			for (j=0; j<N_Temps; j++) cout << "f[" << j << "]....... " << (double)n_up[j]/(n_up[j] + n_down[j]) << endl;
			//cout << endl;
		}

		// update Fisher matrix proposal distribution
		if (i*N_undersample%500 == 0 and i*N_undersample>N_BURN*N_undersample + 500)
		{
			for (j=0; j<N_Temps; j++)
			{
				who = ID_list[j];
				modelX_list[who]->wave.set_Fisher(lisa, X_flag);
				for (k=0; k<modelX_list[who]->wave.D; k++) vector<double> e_vecs[k];
				auto result = modelX_list[who]->wave.get_EigenBS();
				e_vals = get<0>(result);
				e_vecs = get<1>(result);

				prop->P_Fisher.e_vals_list[j] = e_vals;
				prop->P_Fisher.e_vecs_list[j] = e_vecs;
			}
		}

		// update the history tensor
		loc = (i+N_BURN)%hist_stride;
		for(j=0; j<N_Temps; j++)
		{
			for(k=0; k<modelX0->wave.D; k++) history->at(loc).at(j).at(k) = modelX_list[ID_list[j]]->wave.paramsND[k];
		}
		prop->P_DE.cnt++;

		// Store the cold chain
		if (i>=0)
		{
			who = ID_list[0]; // only store the cold chain
			out_file << scientific << setprecision(15) << modelX_list.at(who)->logL;
			for (k=0; k<modelX_list.at(who)->wave.paramsND.size(); k++)
			{
				out_file << scientific << setprecision(15) << " " << modelX_list.at(who)->wave.paramsND[k];
			}
			out_file << scientific << setprecision(15) << endl;
			out_file.flush();
		}
		if (PT == 1 and i>=0)
		{
			who = ID_list[1]; // only store the cold chain
			T1_file << scientific << setprecision(15) << modelX_list.at(who)->logL;
			for (k=0; k<modelX_list.at(who)->wave.paramsND.size(); k++)
			{
				T1_file << scientific << setprecision(15)<< " " << modelX_list.at(who)->wave.paramsND[k];
			}
			T1_file << scientific << setprecision(15) << endl;

			who = ID_list[N_Temps-1];
			THot_file << scientific << setprecision(15) << modelX_list.at(who)->logL;
			for (k=0; k<modelX_list.at(who)->wave.paramsND.size(); k++)
			{
				THot_file << scientific << setprecision(15)<< " " << modelX_list.at(who)->wave.paramsND[k];
			}
			THot_file << scientific << setprecision(15) << endl;

			// Store the logLs
			logL_file << i;
			for (j=0; j<N_Temps; j++) logL_file << scientific << setprecision(15) << " " << modelX_list.at(ID_list[j])->logL;
			logL_file << endl;
			logL_file.flush();

			for (j=0; j<N_Temps; j++) ID_file << scientific << setprecision(15) << ID_list[j] << " ";
			ID_file << endl;
		}

	}
	out_file.close();
	ID_file.close();
	logL_file.close();
	T1_file.close();
	THot_file.close();
	progressBar.done();

	///////////////////////////////// Report Summary Statistics ///////////////////////////////////////
	cout << endl;
	for (j=0; j<N_Temps; j++)
	{
		cout << "Acceptance rate[" << j << "]............ " << 100.*acc[j]/acc_cnt[j] << " %" << endl;
	}
	cout << endl;
	for (j=0; j<N_Temps-1; j++)
	{
		cout << "Swap Acceptance rate[" << j << "<->" << j+1 << "]....... " << 100.*swap[j]/swap_cnt[j] << " %" << endl;
	}
	cout << endl;
	for (j=0; j<N_Temps; j++) cout << "T[" << j << "]...... " << Temps[j] << endl;
	cout << endl;

	/////////
	ofstream Temperature_file;
	Temperature_file.open(Files->Temperature_file);
	for (j=0; j<N_Temps; j++)
	{
		Temperature_file << scientific << setprecision(15) << Temps[j] << " " << (double)n_up[j]/(n_up[j] + n_down[j]) << endl;
	}
	Temperature_file.close();

	cout << "Proposal: Fisher acceptance rate.......... " << (double)prop->P_Fisher.acc/prop->P_Fisher.acc_cnt       << endl;
	cout << "Proposal: prior acceptance rate........... " << (double)prop->prior.acc/prop->prior.acc_cnt             << endl;
	cout << "Proposal: DE acceptance rate.............. " << (double)prop->P_DE.acc/prop->P_DE.acc_cnt               << endl;
	cout << "Proposal: Target acceptance rate.......... " << (double)prop->P_Target.acc/prop->P_Target.acc_cnt       << endl;
	cout << "Proposal: TimeShift acceptance rate....... " << (double)prop->P_TimeShift.acc/prop->P_TimeShift.acc_cnt << endl;
	////////////////////////////////////////////////

	// handle leftover memory
	//delete modelX0->wave;
	delete lisa;
	for (j=0; j<N_Temps; j++)
	{
		delete modelX_list[j];
	}
	delete prop;
	gsl_rng_free (r);
	delete history;
	delete Files;
	delete modelX0;
	delete model_true;

	// Find out how long the code took to run
	auto end = chrono::system_clock::now();
	auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
	cout << "Time elapsed..... " << (double)elapsed.count()/1000 << " seconds \n";

	return 0;
}
