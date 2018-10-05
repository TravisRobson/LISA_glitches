/*
 * G_B_MCMC.h
 *
 *  Created on: Sep 7, 2018
 *      Author: travisrobson
 */

#ifndef G_B_MCMC_H_
#define G_B_MCMC_H_

#include <map>
#include <string>

#include "Wavelet.h"
#include "LISA.h"
#include "TDI.h"

using namespace wv;
using namespace tdi;
using namespace ls;

#include <gsl/gsl_rng.h>


// define limits associated with prior ranges

#define A_lo log(1.0e-23/A_scale)
#define A_hi log(1.0e-18/A_scale)

#define f0_lo 1.0e-3/f0_scale
#define f0_hi 30.0e-3/f0_scale

#define t0_lo 0.0/t0_scale

#define tau_lo log(HOUR/3600/10/tau_scale)
#define tau_hi log(4*HOUR/tau_scale)

#define phi0_lo 0.0
#define phi0_hi 2*M_PI

// burst parameters bounds
#define cos_theta_lo -1.0
#define cos_theta_hi  1.0

#define phi_lo 0.0
#define phi_hi 2*M_PI

#define psi_lo 0.0
#define psi_hi M_PI

#define ellip_lo 0.0
#define ellip_hi 1.0

namespace mc {

//class PDF
//{
//public:
//	PDF();
//	virtual ~PDF();
//	//virtual PDF(const PDF&);
//
//	gsl_rng *r;
//	double weight;
//
//	virtual vector<double> rvs(vector<double> paramsND, double T) { ; }
//	virtual double logpdf(vector<double> paramND) { ; }
//};

class Model {
public:
	Model(Wavelet a);
	Model();
	virtual ~Model();
	Model(const Model &src);

	void set_logL(TDI tdi, LISA *lisa, int X_flag);

	Wavelet wave;
	double logL;

private:

};

class Proposal_DE
{
public:
	Proposal_DE();
	Proposal_DE(double weight, gsl_rng *r, vector< vector<vector<double>> > *history, int, vector<double> Temps);
	virtual ~Proposal_DE();

	gsl_rng *r;
	double weight;
	int cnt;
	vector< vector<vector<double>> > *history;
	int hist_stride;
	vector<double> Temps;
	int acc, acc_cnt;
	string name;

	vector<double> rvs(vector<double> paramsND, double T);
	double logpdf(vector<double> paramND);
};

class Proposal_Fisher //: public PDF
{
public:
	Proposal_Fisher();
	Proposal_Fisher(double weight, gsl_rng *r, vector<double> Temps, vector< vector<double> > e_vals_list, vector< vector<vector<double>> > e_vecs_list);
	virtual ~Proposal_Fisher();
	Proposal_Fisher(const Proposal_Fisher&);

	gsl_rng *r;
	double weight;
	vector<double> Temps;
	vector< vector<double> > e_vals_list;
	vector< vector<vector<double>> > e_vecs_list;
	int acc, acc_cnt;
	string name;

	vector<double> rvs(vector<double> paramsND, double T);
	double logpdf(vector<double> paramND);
};

class Prior //: public PDF
{
public:
	Prior();
	Prior(double weight, gsl_rng *r, double T);
	virtual ~Prior();
	Prior(const Prior&);

	gsl_rng *r;
	double weight;
	double T; // observation period for t0 prior
	int acc, acc_cnt;
	string name;

	double logpdf(vector<double> paramsND);
	vector<double> rvs(vector<double> paramsND, double T);
};

class Proposal_Target
{
public:
	Proposal_Target();
	virtual ~Proposal_Target();
	Proposal_Target(double weight, gsl_rng *r, double f0, double tau, double sig_f0, double sig_tau);

	gsl_rng *r;
	double weight;
	double f0, tau, sig_f0, sig_tau;
	double weight_uni, weight_gau;
	double fact_gau;
	int acc, acc_cnt;
	string name;

	double logpdf(vector<double> paramsND);
	vector<double> rvs(vector<double> paramsND, double T);
};

class Proposal_TimeShift
{
public:
	Proposal_TimeShift();
	virtual ~Proposal_TimeShift();
	Proposal_TimeShift(double weight, gsl_rng*r);

	double weight, sig_t0, sig_phi0;
	gsl_rng *r;
	int acc, acc_cnt;
	string name;

	double logpdf(vector<double> paramsND);
	vector<double> rvs(vector<double> paramsND, double T);
};

struct Proposal
{
	Proposal_Fisher P_Fisher;
	Prior prior;
	Proposal_DE P_DE;
	Proposal_Target P_Target;
	Proposal_TimeShift P_TimeShift;
};

tuple<vector<double>, double, double, string> propose_params(struct Proposal *, vector<double> paramsND, double T);

vector<double> adapt_Temps(vector<double> Temps, vector<int> swap, vector<int> swap_cnt, int t);

} /* namespace mc */

#endif /* G_B_MCMC_H_ */
