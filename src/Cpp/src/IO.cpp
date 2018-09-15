/*
 * IO.cpp
 *
 *  Created on: Sep 14, 2018
 *      Author: travisrobson
 */

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <iomanip>


using namespace std;

#include "IO.h"
#include "Wavelet.h"
#include "LISA.h"
#include "G_B_MCMC.h"

using namespace wv;
using namespace ls;
using namespace mc;

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/foreach.hpp>

vector<Model*> parse_config_file(string input_file, double T, LISA *lisa, struct Files *Files)
{

	ifstream config_file(input_file, ios::in);
	if (! config_file)
	{
	  cerr << "unable to open input file: " << input_file << " --bailing out! \n";
	  exit (EXIT_FAILURE);
	}
	string model_name, model_snr;
	string model_A, model_f0, model_t0, model_tau, model_phi0,
	       model_cos_theta, model_phi, model_psi, model_ellip;

	string data_name, data_snr;
	string data_A, data_f0, data_t0, data_tau, data_phi0,
	       data_cos_theta, data_phi, data_psi, data_ellip;

	string File_true_waveform;

	string line, word, tag;
	vector<string> words;

	while (config_file)
	{
		getline(config_file, line);
		cout << line << endl;

		if (line.substr(0,7) == "model g" or line.substr(0,7) == "model b")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			model_name = words[1];
		}
		if (line.substr(0,9) == "model_snr")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			model_snr = words[1];
		}
		if (line.substr(0,7) == "model_A")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			model_A = words[1];
		}
		if (line.substr(0,8) == "model_f0")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			model_f0 = words[1];
		}
		if (line.substr(0,8) == "model_t0")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			model_t0 = words[1];
		}
		if (line.substr(0,9) == "model_tau")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			model_tau = words[1];
		}
		if (line.substr(0,10) == "model_phi0")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			model_phi0 = words[1];
		}
		if (line.substr(0,15) == "model_cos_theta")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			model_cos_theta = words[1];
		}
		if (line.substr(0,9) == "model_phi")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			model_phi = words[1];
		}
		if (line.substr(0,9) == "model_psi")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			model_psi = words[1];
		}
		if (line.substr(0,11) == "model_ellip")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			model_ellip = words[1];
		}

		if (line.substr(0,6) == "data g" or line.substr(0,6) == "data b")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			data_name = words[1];
		}

		if (line.substr(0,8) == "data_snr")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			data_snr = words[1];
		}
		if (line.substr(0,6) == "data_A")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			data_A = words[1];
		}
		if (line.substr(0,7) == "data_f0")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			data_f0 = words[1];
		}
		if (line.substr(0,7) == "data_t0")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			data_t0 = words[1];
		}
		if (line.substr(0,8) == "data_tau")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			data_tau = words[1];
		}
		if (line.substr(0,9) == "data_phi0")
		{
			boost::split(words, line, boost::is_any_of(" ") );
			data_phi0 = words[1];
		}
		tag = "data_cos_theta";
		if (line.substr(0,tag.length()) == tag)
		{
			boost::split(words, line, boost::is_any_of(" ") );
			data_cos_theta = words[1];
		}
		tag = "data_phi";
		if (line.substr(0,tag.length()) == tag)
		{
			boost::split(words, line, boost::is_any_of(" ") );
			data_phi = words[1];
		}
		tag = "data_psi";
		if (line.substr(0,tag.length()) == tag)
		{
			boost::split(words, line, boost::is_any_of(" ") );
			data_psi = words[1];
		}
		tag = "data_ellip";
		if (line.substr(0,tag.length()) == tag)
		{
			boost::split(words, line, boost::is_any_of(" ") );
			data_ellip = words[1];
		}

		// extract the file names to write results to
		tag = "File_true_waveform";
		if (line.substr(0, tag.length()) == tag)
		{
			boost::split(words, line, boost::is_any_of(" ") );
			File_true_waveform = words[1];
		}
		tag = "File_Temps";
		if (line.substr(0, tag.length()) == tag)
		{
			boost::split(words, line, boost::is_any_of(" ") );
			Files->Temperature_file = words[1];
		}
	}
	cout << endl;
	config_file.close();

	double A = stod(model_A);
	double f0 = stod(model_f0);
	boost::split(words, model_t0, boost::is_any_of("*"));
	double t0 = T*stod(words[0]);
	boost::split(words, model_tau, boost::is_any_of("*"));
	double tau = HOUR*stod(words[0]);
	double phi0 = stod(model_phi0);
	vector<double> paramsND;
	if (model_name.size() == 11)
	{	// i.e. this is a glitch
		paramsND = {log(A/A_scale), f0/f0_scale, t0/t0_scale, log(tau/tau_scale), phi0};
	}
	else
	{
		double cos_theta = stod(model_cos_theta);
		double phi = stod(model_phi);
		double psi = stod(model_psi);
		double ellip = stod(model_ellip);

		paramsND = {log(A/A_scale), f0/f0_scale, t0/t0_scale, log(tau/tau_scale), phi0,
									cos_theta, phi, psi, ellip};
	}
	Wavelet *wavelet = new Wavelet(model_name, paramsND);
	wavelet->calc_TDI(lisa);
	wavelet->set_snr(lisa);
	wavelet->adjust_snr(stod(model_snr), lisa);


	Model *modelX0 = new Model(*wavelet);


	///////////////////////
	A = stod(data_A);
	f0 = stod(data_f0);
	boost::split(words, data_t0, boost::is_any_of("*"));
	t0 = T*stod(words[0]);
	boost::split(words, data_tau, boost::is_any_of("*"));
	tau = HOUR*stod(words[0]);
	phi0 = stod(data_phi0);
	vector<double> paramsND_true;
	if (model_name.size() == 11)
	{	// i.e. this is a glitch
		paramsND_true = {log(A/A_scale), f0/f0_scale, t0/t0_scale, log(tau/tau_scale), phi0};
	}
	else
	{
		double cos_theta = stod(data_cos_theta);
		double phi = stod(data_phi);
		double psi = stod(data_psi);
		double ellip = stod(data_ellip);

		paramsND_true = {log(A/A_scale), f0/f0_scale, t0/t0_scale, log(tau/tau_scale), phi0,
									cos_theta, phi, psi, ellip};
	}
	cout << "Q............. " << 2*M_PI*f0*tau << "\n";

	Wavelet *wavelet_true = new Wavelet(data_name, paramsND_true);
	wavelet_true->calc_TDI(lisa);
	wavelet_true->set_snr(lisa);
	cout << "initial AET snr.......... " << wavelet_true->snr << endl;
	wavelet_true->adjust_snr(stod(data_snr), lisa);
	cout << "adjusted AET snr......... " << wavelet_true->snr << endl;

	Model *model_true = new Model(*wavelet_true);

	model_true->set_logL(wavelet_true->tdi, lisa);
	modelX0->set_logL(wavelet_true->tdi, lisa);

	vector<Model*> models;
	models.push_back(modelX0);
	models.push_back(model_true);


	// Print the true model to file
	ofstream out_file;
	complex<double> Ai, Ei, Ti;
	double fi;
	int i;

	out_file.open(File_true_waveform);
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

	return models;
}


