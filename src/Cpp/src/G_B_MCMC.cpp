/*
 * G_B_MCMC.cpp
 *
 *  Created on: Sep 7, 2018
 *      Author: travisrobson
 */

#include <cmath>
#include <map>
#include <iostream>
#include <algorithm>

using namespace std;

#include "G_B_MCMC.h"
#include "Constants.h"
#include "Wavelet.h"
#include "LISA.h"
#include "TDI.h"

using namespace wv;
using namespace ls;
using namespace tdi;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


namespace mc {

PDF::PDF()
{
}

PDF::~PDF()
{

}

//PDF::PDF(const PDF&) {}

Model::Model(Wavelet wv)
{
	this->wave = Wavelet(wv, 1);
}

Model::Model()
{
}

Model::~Model()
{

}

Model::Model(const Model &src)
{
	this->wave = Wavelet(src.wave, 1);
	this->logL = src.logL;
}

void Model::set_logL(TDI data, LISA *lisa)
{
	double result;

	result = nwip(&data, &this->wave.tdi, lisa);

	result -= 0.5*pow(this->wave.snr, 2.0);

	this->logL = result;
}

Proposal_Fisher::Proposal_Fisher(double weight, gsl_rng *r, vector<double> Temps,  vector< vector<double> > e_vals_list, vector< vector<vector<double>> > e_vecs_list)
{
	this->weight = weight;
	this->r = r;
	this->Temps = Temps;
	this->e_vals_list = e_vals_list;
	this->e_vecs_list = e_vecs_list;
}

Proposal_Fisher::Proposal_Fisher() {}

Proposal_Fisher::Proposal_Fisher(const Proposal_Fisher& src)
{
	this->r = src.r;
	this->weight = src.weight;
	this->Temps = src.Temps;
	this->e_vals_list = src.e_vals_list;
	this->e_vecs_list = src.e_vecs_list;
}

Proposal_Fisher::~Proposal_Fisher()
{

}

double Proposal_Fisher::logpdf(vector<double> paramND)
{
	return 0.0;
}

vector<double> Proposal_Fisher::rvs(vector<double> paramsND, double T)
{
	int i, dir;
	int T_idx = -1;

	for (i=0; i<this->Temps.size(); i++)
	{	// identify which temperature we are dealing with
		if (T == this->Temps[i])
		{
			T_idx = i;
			break;
		}
	}
	if (T_idx == -1) throw invalid_argument{"Temperature T is not present in temperature ladder attached to this Fisher"};

	int D = this->e_vals_list[T_idx].size();

	vector<double> y(D);

	double u = gsl_ran_gaussian(this->r, 1.0);

	// choose a random eigenvector
	int dir_list[D];
	for (i=0; i<D; i++) dir_list[i] = i;
	double val = -1.0;
	while (val < 1.0e-10)
	{	// don't let negative eigen values get used
		gsl_ran_choose(this->r, &dir, 1, dir_list, D, sizeof(int));
		val = this->e_vals_list[T_idx][dir];
	}

	for (i=0; i<D; i++) y[i] = paramsND[i] + u/sqrt(this->e_vals_list[T_idx][dir])*this->e_vecs_list[T_idx][i][dir];

	return y;
}

Prior::Prior(double weight, gsl_rng *r, double T)
{
	this->weight = weight;
	this->r = r;
	this->T = T; // observation period
}

Prior::Prior() { }

Prior::Prior(const Prior& src)
{
	this->r = src.r;
	this->weight = src.weight;
	this->T = src.T;
}


Prior::~Prior() {}

double Prior::logpdf(vector<double> paramsND)
{
	double result = 1.0;

	result *= gsl_ran_flat_pdf(paramsND[IDX_A], A_lo, A_hi);
	result *= gsl_ran_flat_pdf(paramsND[IDX_f0], f0_lo, f0_hi);
	result *= gsl_ran_flat_pdf(paramsND[IDX_t0], t0_lo, this->T/WEEK);
	result *= gsl_ran_flat_pdf(paramsND[IDX_tau], tau_lo, tau_hi);
	result *= gsl_ran_flat_pdf(paramsND[IDX_phi0], phi0_lo, phi0_hi);

	return log(result);
}

vector<double> Prior::rvs(vector<double> paramsND, double T)
{
	vector<double> y(paramsND.size());

	double u = gsl_rng_uniform(this->r);
	y[IDX_A] = gsl_ran_flat(this->r, A_lo, A_hi);

	u = gsl_rng_uniform(this->r);
	y[IDX_f0] = gsl_ran_flat(this->r, f0_lo, f0_hi);

	u = gsl_rng_uniform(this->r);
	y[IDX_t0] = gsl_ran_flat(this->r, t0_lo, this->T/WEEK);

	u = gsl_rng_uniform(this->r);
	y[IDX_tau] = gsl_ran_flat(this->r, tau_lo, tau_hi);

	u = gsl_rng_uniform(this->r);
	y[IDX_phi0] = gsl_ran_flat(this->r, phi0_lo, phi0_hi);

	return y;
}

Proposal_DE::Proposal_DE() { }

Proposal_DE::Proposal_DE(double weight, gsl_rng *r, vector< vector<vector<double>> > *history, int hist_stride, vector<double> Temps)
{
	this->weight = weight;
	this->r = r;
	this->history = history;
	this->hist_stride = hist_stride;
	this->Temps = Temps;
	this->cnt = 0;
}

Proposal_DE::~Proposal_DE() { }



vector<double> Proposal_DE::rvs(vector<double> paramsND, double T)
{
	if (this->cnt < 2) return paramsND;

	vector<double> y(paramsND.size());

	int i;
	int T_idx = -1;

	for (i=0; i<this->Temps.size(); i++)
	{	// identify which temperature we are dealing with
		if (T == this->Temps[i])
		{
			T_idx = i;
			break;
		}
	}
	if (T_idx == -1) throw invalid_argument{"Temperature T is not present in temperature ladder attached to this Fisher"};

	int hist_size = this->history->size();

	typedef vector<double> D1;
	typedef vector<D1> D2;
	typedef vector<D2> matrix_3D;

	matrix_3D hist = *(this->history);

	int max_idx = min(this->cnt, hist_size);

	int j = (int)( (double)max_idx*gsl_rng_uniform(this->r) );
	int k = j;
	while(j == k) k = (int)( (double)max_idx*gsl_rng_uniform(this->r) );

	double alpha = 1.0;
	double beta = gsl_rng_uniform(r);

	if (beta < 0.9) alpha = gsl_ran_gaussian(this->r, 1.);
	for (i=0; i<paramsND.size(); i++)
	{
		y[i] = paramsND[i] + alpha*(hist.at(j).at(T_idx).at(i) - hist.at(k).at(T_idx).at(i));
	}

	return y;
}

double Proposal_DE::logpdf(vector<double> paramND) { return 0.0; }

vector<double> propose_params(struct Proposal *prop, vector<double> paramsND, double T)
{
	vector<double> y;

	double u = gsl_ran_flat(prop->P_Fisher.r, 0, 1); // this decides which proposal to make use of

	vector<double> which(3);
	which[0] = prop->P_Fisher.weight;
	which[1] = which[0] + prop->prior.weight;
	which[2] = which[1] + which[0] +  prop->P_DE.weight;

	if (u < which[0])
	{	// Fisher Jump
		y = prop->P_Fisher.rvs(paramsND, T);
	}
	else if (u>which[0] & u<which[1])
	{	// Prior Jump
		 y = prop->prior.rvs(paramsND, T);
	}
	else
	{  // DE Jump
		y = prop->P_DE.rvs(paramsND, T);
	}

	return y;
}

vector<double> adapt_Temps(vector<double> Temps, vector<int> swap, vector<int> swap_cnt, int t)
{
	int i;
	int N_Temps = Temps.size();

	double nu = 10;
	double t0 = 10000;
	double kappa;

	// calculate the acceptance rates
	vector<double> A(N_Temps - 1);
	for (i=0; i<N_Temps-1; i++)
	{
		if (swap_cnt[i] > 0) A[i] = (double)swap[i]/swap_cnt[i];
		else A[i] = 0.0;
	}

	vector<double> new_Temps(N_Temps);
	new_Temps[0] = 1.0;
	new_Temps[N_Temps-1] = Temps[N_Temps - 1];

	for (i=1; i<N_Temps-1; i++)
	{
		kappa = t0/(t0+t)/nu;
		new_Temps[i] = new_Temps[i-1] + (Temps[i] - Temps[i-1])*exp(kappa*(A[i-1] - A[i]));
		if (new_Temps[i]/new_Temps[i-1] < 1.1) new_Temps[i] = 1.1*new_Temps[i-1];
	}

	return new_Temps;
}

} /* namespace mcmc */
