/*
 * G_B_MCMC.cpp
 *
 *  Created on: Sep 7, 2018
 *      Author: travisrobson
 *
 *      -ffast-math -ftree-vectorize
 */

#include <cmath>
#include <map>
#include <iostream>
#include <algorithm>
#include <string>
#include <complex>


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


#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

namespace mc {
//
//PDF::PDF()
//{
//}
//
//PDF::~PDF()
//{
//
//}

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

gsl_matrix *invert_a_matrix(gsl_matrix *matrix)
{
    gsl_permutation *p = gsl_permutation_alloc(4);
    int s;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(matrix, p, &s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(4,4);
    gsl_linalg_LU_invert(matrix, p, inv);

    gsl_permutation_free(p);

    return inv;
}

void Model::set_logL(TDI data, LISA *lisa)
{
	double result;

	if (this->wave.get_name() == "burst")
	{
		vector<double> paramsND1 = {log(A_scale/A_scale), this->wave.paramsND[IDX_f0], this->wave.paramsND[IDX_t0],
								    this->wave.paramsND[IDX_tau], 0.0, this->wave.paramsND[IDX_cos_theta],
									this->wave.paramsND[IDX_phi], 0.0, 0.0};
		Wavelet *A1 = new Wavelet("burst", paramsND1);

		vector<double> paramsND2 = {log(A_scale/A_scale), this->wave.paramsND[IDX_f0], this->wave.paramsND[IDX_t0],
								    this->wave.paramsND[IDX_tau], 0.5*M_PI, this->wave.paramsND[IDX_cos_theta],
									this->wave.paramsND[IDX_phi], -0.25*M_PI, 0.0};
		Wavelet *A2 = new Wavelet("burst", paramsND2);

		vector<double> paramsND3 = {log(A_scale/A_scale), this->wave.paramsND[IDX_f0], this->wave.paramsND[IDX_t0],
								    this->wave.paramsND[IDX_tau], 0.0, this->wave.paramsND[IDX_cos_theta],
									this->wave.paramsND[IDX_phi], -0.25*M_PI, 0.0};
		Wavelet *A3 = new Wavelet("burst", paramsND3);

		vector<double> paramsND4 = {log(A_scale/A_scale), this->wave.paramsND[IDX_f0], this->wave.paramsND[IDX_t0],
								    this->wave.paramsND[IDX_tau], 0.5*M_PI, this->wave.paramsND[IDX_cos_theta],
									this->wave.paramsND[IDX_phi], 0.0, 0.0};
		Wavelet *A4 = new Wavelet("burst", paramsND4);

		A1->calc_TDI(lisa);
		A2->calc_TDI(lisa);
		A3->calc_TDI(lisa);
		A4->calc_TDI(lisa);

		vector<double> N(4);
		N[0] = nwip(&data, &A1->tdi, lisa);
		N[1] = nwip(&data, &A2->tdi, lisa);
		N[2] = nwip(&data, &A3->tdi, lisa);
		N[3] = nwip(&data, &A4->tdi, lisa);



		vector<double> v1(4);
		vector<vector<double>> M(4, v1);

		M[0][0] = nwip(&A1->tdi, &A1->tdi, lisa);
		M[0][1] = nwip(&A1->tdi, &A2->tdi, lisa);
		M[0][2] = nwip(&A1->tdi, &A3->tdi, lisa);
		M[0][3] = nwip(&A1->tdi, &A4->tdi, lisa);

		M[1][1] = nwip(&A2->tdi, &A2->tdi, lisa);
		M[1][2] = nwip(&A2->tdi, &A3->tdi, lisa);
		M[1][3] = nwip(&A2->tdi, &A4->tdi, lisa);

		M[2][2] = nwip(&A3->tdi, &A3->tdi, lisa);
		M[2][3] = nwip(&A3->tdi, &A4->tdi, lisa);

		M[3][3] = nwip(&A4->tdi, &A4->tdi, lisa);

		M[1][0] = M[0][1];
		M[2][0] = M[0][2];
		M[3][0] = M[0][3];

		M[2][1] = M[1][2];
		M[3][1] = M[1][3];

		M[3][2] = M[2][3];

//	    for (int i=0; i<4; ++i)
//	    {
//	    	cout << "N[" << i << "] " << N[i] << endl;
//	        for (int j=0; j<4; ++j)
//			{
//	        	cout << "M[" << i << "][" << j<< "] " << M[i][j] << endl;
//			}
//	        cout << endl;
//	    }

		gsl_matrix *mat = gsl_matrix_alloc(4,4);
	    for (int i=0; i<4; ++i)
	    {
	        for (int j=0; j<4; ++j) gsl_matrix_set(mat, i, j, M[i][j]);
	    }
	    gsl_matrix *inverse = invert_a_matrix(mat);

	    for (int i=0; i<4; ++i)
	    {
	        for (int j=0; j<4; ++j)
			{
	        	M[i][j] = gsl_matrix_get(inverse, i, j);
	        	//cout << M[i][j] << endl;
			}
	    }
	    gsl_matrix_free(mat);
	    gsl_matrix_free(inverse);
	    this->logL = 0.0;
	    for (int i=0; i<4; i++)
	    {
	    	for (int j=0; j<4; j++)
	    	{
	    		this->logL += M[i][j]*N[i]*N[j];
	    	}
	    }
	    this->logL *= 0.5;
	    //cout << "logL.............. " << this->logL << endl << endl;

		delete A1;
		delete A2;
		delete A3;
		delete A4;

	}
	else
	{
		result = nwip(&data, &this->wave.tdi, lisa);

		result -= 0.5*pow(this->wave.snr, 2.0);

		this->logL = result;
	}
}

Proposal_Fisher::Proposal_Fisher(double weight, gsl_rng *r, vector<double> Temps,  vector< vector<double> > e_vals_list, vector< vector<vector<double>> > e_vecs_list)
{
	this->weight = weight;
	this->r = r;
	this->Temps = Temps;
	this->e_vals_list = e_vals_list;
	this->e_vecs_list = e_vecs_list;

	this->acc = 0;
	this->acc_cnt = 0;
	this->name = "Fisher";
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
	int cnt = 0;

	while (val < 1.0e-10 or val != val)
	{	// don't let negative eigen values get used
		gsl_ran_choose(this->r, &dir, 1, dir_list, D, sizeof(int));
		val = this->e_vals_list[T_idx][dir];
		cnt++;
		if (cnt>30) return paramsND;
	}

	for (i=0; i<D; i++) y[i] = paramsND[i] + u/sqrt(this->e_vals_list[T_idx][dir])*this->e_vecs_list[T_idx][i][dir];

	return y;
}

Prior::Prior(double weight, gsl_rng *r, double T)
{
	this->weight = weight;
	this->r = r;
	this->T = T; // observation period

	this->acc = 0;
	this->acc_cnt = 0;
	this->name = "Prior";
}

Prior::Prior() { }

Prior::Prior(const Prior& src)
{
	this->r = src.r;
	this->weight = src.weight;
	this->T = src.T;

//	this->acc = src.acc;
//	this->acc_cnt = 0;
//	this->name = "Prior";
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

	if (paramsND.size()>5)
	{
		result *= gsl_ran_flat_pdf(paramsND[IDX_cos_theta], cos_theta_lo, cos_theta_hi);
		result *= gsl_ran_flat_pdf(paramsND[IDX_phi], phi_lo, phi_hi);
		result *= gsl_ran_flat_pdf(paramsND[IDX_psi], psi_lo, psi_hi);
		result *= gsl_ran_flat_pdf(paramsND[IDX_ellip], ellip_lo, ellip_hi);
	}

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

	if (paramsND.size()>5)
	{
		u = gsl_rng_uniform(this->r);
		y[IDX_cos_theta] = gsl_ran_flat(this->r, cos_theta_lo, cos_theta_hi);

		u = gsl_rng_uniform(this->r);
		y[IDX_phi]       = gsl_ran_flat(this->r, phi_lo, phi_hi);

		u = gsl_rng_uniform(this->r);
		y[IDX_psi]       = gsl_ran_flat(this->r, psi_lo, psi_hi);

		u = gsl_rng_uniform(this->r);
		y[IDX_ellip]     = gsl_ran_flat(this->r, ellip_lo, ellip_hi);
	}

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

	this->acc = 0;
	this->acc_cnt = 0;
	this->name = "DiffEv";
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

Proposal_Target::Proposal_Target() { }

Proposal_Target::Proposal_Target(double weight, gsl_rng *r, double f0, double tau, double sig_f0, double sig_tau)
{
	this->weight  = weight;
	this->r       = r;
	this->f0      = f0;
	this->tau     = tau;
	this->sig_f0  = sig_f0;
	this->sig_tau = sig_tau;

	this->weight_uni = 0.2;
	this->weight_gau = 1.0 - this->weight_uni;

	this->fact_gau = 2.0; // how many Fisher estimated sigma's wide for Gaussian?

	this->acc = 0;
	this->acc_cnt = 0;
	this->name = "Target";
}

Proposal_Target::~Proposal_Target() { }

double Proposal_Target::logpdf(vector<double> paramsND)
{
	double t1, t2;
	double result = 0.0;

	t1  = this->weight_uni*gsl_ran_flat_pdf(paramsND[IDX_f0], f0_lo, f0_hi);
	t1 += this->weight_gau*gsl_ran_gaussian_pdf(paramsND[IDX_f0] - this->f0, this->sig_f0);

	t2  = this->weight_uni*gsl_ran_flat_pdf(paramsND[IDX_tau], tau_lo, tau_hi);
	t2 += this->weight_gau*gsl_ran_gaussian_pdf(paramsND[IDX_tau] - this->tau, this->sig_tau);

	return log(t1*t2);
}

vector<double> Proposal_Target::rvs(vector<double> paramsND, double T)
{
	int D = paramsND.size();

	vector<double> y(D);

	y[IDX_A] = paramsND[IDX_A];

	double u = gsl_ran_flat(this->r, 0, 1);

	if (u < this->weight_gau)
	{
		y[IDX_f0] = paramsND[IDX_f0] + gsl_ran_gaussian(this->r, this->sig_f0);
	}
	else
	{
		y[IDX_f0] = gsl_ran_flat(this->r, f0_lo, f0_hi);
	}

	y[IDX_t0] = paramsND[IDX_t0];

	u = gsl_ran_flat(this->r, 0, 1);
	if (u < this->weight_gau)
	{
		y[IDX_tau] = paramsND[IDX_tau] + gsl_ran_gaussian(this->r, this->sig_tau);
	}
	else
	{
		y[IDX_tau] = gsl_ran_flat(this->r, tau_lo, tau_hi);
	}

	for (int i=IDX_phi0; i<D; i++) y[i] = paramsND[i];

	return y;
}

Proposal_TimeShift::Proposal_TimeShift() { }

Proposal_TimeShift::~Proposal_TimeShift() { }

Proposal_TimeShift::Proposal_TimeShift(double weight, gsl_rng *r)
{
	this->weight = weight;
	this->r = r;
	this->sig_t0 = 1*HOUR;
	this->sig_phi0 = 0.3;
	this->acc = 0;
	this->acc_cnt = 0;
	this->name = "TimeShift";
}

vector<double> Proposal_TimeShift::rvs(vector<double> paramsND, double T)
{
	int D = paramsND.size();

	vector<double> y(D);

	y[IDX_A] = paramsND[IDX_A];
	y[IDX_f0] = paramsND[IDX_f0];
	y[IDX_tau] = paramsND[IDX_tau];

	double u = gsl_ran_flat(this->r, 0, 1); // half a chance of going forward or backwards, and...
	if (u<0.25)
	{
		y[IDX_t0] = paramsND[IDX_t0] + 1./(paramsND[IDX_f0]*f0_scale)/t0_scale*gsl_ran_gaussian(this->r, this->sig_t0);
		y[IDX_phi0] = paramsND[IDX_phi0];
	}
	else if (u>0.25 & u<0.5)
	{
		y[IDX_t0] = paramsND[IDX_t0] - 1./(paramsND[IDX_f0]*f0_scale)/t0_scale*gsl_ran_gaussian(this->r, this->sig_t0);
		y[IDX_phi0] = paramsND[IDX_phi0];
	}
	else if (u>0.5 & u<0.75)
	{
		y[IDX_t0] = paramsND[IDX_t0] + 0.5/(paramsND[IDX_f0]*f0_scale)/t0_scale*gsl_ran_gaussian(this->r, this->sig_t0);
		y[IDX_phi0] = paramsND[IDX_phi0] + M_PI*gsl_ran_gaussian(this->r, this->sig_phi0);
	}
	else
	{
		y[IDX_t0] = paramsND[IDX_t0] - 0.5/(paramsND[IDX_f0]*f0_scale)/t0_scale*gsl_ran_gaussian(this->r, this->sig_t0);
		y[IDX_phi0] = paramsND[IDX_phi0] - M_PI*gsl_ran_gaussian(this->r, this->sig_phi0);
	}

	if (D>5)
	{
		for (int i=IDX_cos_theta; i<D; i++) y[i] = paramsND[i];
	}

	return y;
}


tuple<vector<double>, double, double, string> propose_params(struct Proposal *prop, vector<double> paramsND, double T)
{
	vector<double> y;
	string name;

	double u = gsl_ran_flat(prop->P_Fisher.r, 0, 1); // this decides which proposal to make use of

	vector<double> which(5);
	which[0] = prop->P_Fisher.weight;
	which[1] = which[0] + prop->prior.weight;
	which[2] = which[1] + prop->P_DE.weight;
	which[3] = which[2] + prop->P_Target.weight;
	which[4] = which[3] + prop->P_TimeShift.weight;

	double logQyx, logQxy;

	if (u < which[0])
	{	// Fisher Jump
		y = prop->P_Fisher.rvs(paramsND, T);
		logQyx = 0.0;
		logQxy = 0.0;
		name = prop->P_Fisher.name;
		prop->P_Fisher.acc_cnt++;
	}
	else if (u>which[0] & u<which[1])
	{	// Prior Jump
		 y = prop->prior.rvs(paramsND, T);
		logQyx = 0.0;
		logQxy = 0.0;
		name = prop->prior.name;
		prop->prior.acc_cnt++;
	}
	else if (u>which[1] & u<which[2])
	{  // DE Jump
		y = prop->P_DE.rvs(paramsND, T);
		logQyx = 0.0;
		logQxy = 0.0;
		name = prop->P_DE.name;
		prop->P_DE.acc_cnt++;
	}
	else if (u>which[2] & u<which[3])
	{	// Targeting Jump
		y = prop->P_Target.rvs(paramsND, T);
		logQyx = prop->P_Target.logpdf(y);
		logQxy = prop->P_Target.logpdf(paramsND);
		name = prop->P_Target.name;
		prop->P_Target.acc_cnt++;
	}
	else
	{ // TimeShift Jump
		y = prop->P_TimeShift.rvs(paramsND, T);
		logQyx = 0.0;
		logQxy = 0.0;
		name = prop->P_TimeShift.name;
		prop->P_TimeShift.acc_cnt++;
	}

	return {y, logQyx, logQxy, name};
}

vector<double> adapt_Temps(vector<double> Temps, vector<int> swap, vector<int> swap_cnt, int t)
{
	int i;
	int N_Temps = Temps.size();

	double nu = 100;
	double t0 = 1000;
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
		//if (new_Temps[i]/new_Temps[i-1] < 1.05) new_Temps[i] = 1.05*new_Temps[i-1];
	}

	return new_Temps;
}

} /* namespace mcmc */
