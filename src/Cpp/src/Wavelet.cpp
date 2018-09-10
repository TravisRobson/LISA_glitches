/*
 * Wavelet.cpp
 *
 *  Created on: Aug 28, 2018
 *      Author: Travis Robson
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <list>

using namespace std;

#include "MorletGabor.h"
#include "G_B_MCMC.h"
#include "Constants.h"
#include "Wavelet.h"
#include "LISA.h"
#include "TDI.h"

using namespace tdi;
using namespace ls;
using namespace mc;


#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;

namespace wv {

Wavelet::Wavelet() {}

Wavelet::Wavelet(string name, vector<double> paramsND)
{	// Constructor

	this->name = name;
	this->paramsND = paramsND;
	this->D = paramsND.size();
}

Wavelet::~Wavelet()
{

}

Wavelet::Wavelet(const Wavelet &src, int flag)
{
   this->name = src.name;
   this->paramsND = src.paramsND;
   this->D = src.D;

   if (flag == 1)
   {	// if this flag is set copy the snr and TDI
	   this->snr = src.snr;
	   this->tdi = src.tdi;
   }
}

void Wavelet::calc_TDI(LISA *lisa)
{

	if (this->name == "glitch_OP12")
	{
		calc_OP12_TDI(lisa);
	}



	else if (this->name == "burst")
	{
		calc_burst_TDI(lisa);
	}

}

void Wavelet::calc_burst_TDI(LISA *lisa)
{	// Calculate the TDI for an optical path glitch from laser link 1->2
	this->tdi = TDI();

	double T = lisa->get_T();

	////////// extract burst parameters
	double A, f0, t0, tau, phi0;
	unpack_glitch_params(this->paramsND, &A, &f0, &t0, &tau, &phi0);

	double cos_theta, phi, psi, ellip;
	unpack_burst_params(this->paramsND, &cos_theta, &phi, &psi, &ellip);
	double sin_theta = 1. - pow(cos_theta, 2.0);
	double sin_phi = sin(phi);
	double cos_phi = cos(phi);
	double sin_2psi = sin(2*psi);
	double cos_2psi = cos(2*psi);

	////////// estimate an appropriate bandwidth
	double rho_est = 5;

	// Lowest frequency bin
	double df = 2.0/tau*pow(rho_est/5, 1.0);
	int N_lo = (int)((f0 - df)*T);
	if (N_lo < 0) N_lo = 1; // make this the lowest positive frequency (don't want to break noise curve code)
	this->tdi.set_N_lo(N_lo);

	// highest frequency bin
	int N_nyquist = (int)(0.5/dt*T);
	int N_hi = (int)((f0 + df)*T);
	if (N_hi > N_nyquist) N_hi = N_nyquist;
	this->tdi.set_N_hi(N_hi);

	vector<complex<double>> p12(N_hi - N_lo);
	vector<complex<double>> p21(N_hi - N_lo);
	vector<complex<double>> p13(N_hi - N_lo);
	vector<complex<double>> p31(N_hi - N_lo);
	vector<complex<double>> p23(N_hi - N_lo);
	vector<complex<double>> p32(N_hi - N_lo);

	vector<double> u {cos_theta*cos_phi, cos_theta*sin_phi, -sin_theta};
	vector<double> v {sin_phi, -cos_phi, 0.0};


	list<vector<complex<double>>> *phase_list = new list<vector<complex<double>>>;
	phase_list->push_back(p12);
	phase_list->push_back(p21);
	phase_list->push_back(p13);
	phase_list->push_back(p31);
	phase_list->push_back(p23);
	phase_list->push_back(p32);

	// construct the TDI channels
	this->tdi.phase_to_tdi(phase_list, N_lo, T);

	delete phase_list;

}

void Wavelet::calc_OP12_TDI(LISA *lisa)
{	// Calculate the TDI for an optical path glitch from laser link 1->2
	this->tdi = TDI();

	double T = lisa->get_T();

	// extract Morlet-Gabor parameters
	double A, f0, t0, tau, phi0;
	unpack_glitch_params(this->paramsND, &A, &f0, & t0, &tau, &phi0);

	double rho_est = sqrt(sqrt(M_PI/2)*tau/lisa->SnX(f0))*A*2*fabs(sin(f0/fstar)); //  estimate an appropriate bandwidth

	// Lowest frequency bin
	double df = 2.0/tau*pow(rho_est/5, 1.0);
	int N_lo = (int)((f0 - df)*T);
	if (N_lo < 0) N_lo = 1; // make this the lowest positive frequency (don't want to break noise curve code)
	this->tdi.set_N_lo(N_lo);

	// highest frequency bin
	int N_nyquist = (int)(0.5/dt*T);
	int N_hi = (int)((f0 + df)*T);
	if (N_hi > N_nyquist) N_hi = N_nyquist;
	this->tdi.set_N_hi(N_hi);

	// sample the wavelet
	vector<complex<double>> p12(N_hi - N_lo);
	double fi;
	for (int i=0; i<p12.size(); i++)
	{
		fi = (N_lo+i)/T;
		p12[i] = Psi_FT(fi, A, f0, t0, tau, phi0);
	}

	// Construct the other phase variables
	vector<complex<double>> p21(p12.size(), 0.0);
	vector<complex<double>> p13(p12.size(), 0.0);
	vector<complex<double>> p31(p12.size(), 0.0);
	vector<complex<double>> p23(p12.size(), 0.0);
	vector<complex<double>> p32(p12.size(), 0.0);

	list<vector<complex<double>>> *phase_list = new list<vector<complex<double>>>;
	phase_list->push_back(p12);
	phase_list->push_back(p21);
	phase_list->push_back(p13);
	phase_list->push_back(p31);
	phase_list->push_back(p23);
	phase_list->push_back(p32);

	// construct the TDI channels
	this->tdi.phase_to_tdi(phase_list, N_lo, T);

	delete phase_list;
}

void unpack_glitch_params(vector<double> paramsND, double *A, double *f0, double *t0, double *tau, double *phi0)
{
	*A    = exp(paramsND[IDX_A])*A_scale;
	*f0   = paramsND[IDX_f0]*f0_scale;
	*t0   = paramsND[IDX_t0]*t0_scale;
	*tau  = exp(paramsND[IDX_tau])*tau_scale;
	*phi0 = paramsND[IDX_phi0];
}

void unpack_burst_params(vector<double> paramsND, double *cos_theta, double *phi, double *psi, double *ellip)
{
	*cos_theta = paramsND[IDX_cos_theta];
	*phi       = paramsND[IDX_phi];
	*psi       = paramsND[IDX_psi];
	*ellip     = paramsND[IDX_ellip];
}

void Wavelet::set_snr(LISA *lisa)
{
	this->snr = sqrt(nwip(&this->tdi, &this->tdi, lisa));
}

void Wavelet::adjust_snr(double snr_target, LISA *lisa)
{
	this->paramsND[IDX_A] = log(exp(paramsND[IDX_A])*snr_target/this->snr);
	this->calc_TDI(lisa);
	this->set_snr(lisa);
}

void Wavelet::set_Fisher(LISA *lisa)
{
	int i, j;
	int D = this->paramsND.size(); // number of parameters

	double epsilon = 1.0e-7; // small number for differentiating

	this->Fisher = vector<vector<double>>(D);
	for (i=0; i<D; i++) this->Fisher[i] = vector<double>(D);

	for (i=0; i<D; i++)
	{
		Wavelet plus_RHS = Wavelet(*this, 0);
		plus_RHS.paramsND[i] += epsilon;
		plus_RHS.calc_TDI(lisa);

		Wavelet minus_RHS = Wavelet(*this, 0);
		minus_RHS.paramsND[i] -= epsilon;
		minus_RHS.calc_TDI(lisa);

		TDI tdi_diff_RHS = plus_RHS.tdi - minus_RHS.tdi;
		tdi_diff_RHS /= 2*epsilon;

		for (j=i; j<D; j++)
		{
			Wavelet plus_LHS = Wavelet(*this, 0);
			plus_LHS.paramsND[j] += epsilon;
			plus_LHS.calc_TDI(lisa);

			Wavelet minus_LHS = Wavelet(*this, 0);
			minus_LHS.paramsND[j] -= epsilon;
			minus_LHS.calc_TDI(lisa);

			TDI tdi_diff_LHS = plus_LHS.tdi - minus_LHS.tdi;

			tdi_diff_LHS /= 2*epsilon;

			Fisher[i][j] = nwip(&tdi_diff_RHS, &tdi_diff_LHS, lisa);
		}
	}

	for (i=0; i<D; i++)
	{
		for (j=i+1; j<D; j++) Fisher[j][i] = Fisher[i][j];
	}
}

tuple < vector<double>,vector<vector<double>> > Wavelet::get_EigenBS()
{	// calculate the eigenvalues and eigenvectors of a matrix
	int i, j;
	int D = this->Fisher.size();

	MatrixXd m(D, D);
	for (i=0; i<D; i++)
	{
		for (j=0; j<D; j++) m(i,j) = this->Fisher[i][j];
	}

	SelfAdjointEigenSolver<MatrixXd> es(m);

	// convert Eigen eigenvalues to a STD vector
	vector<double> e_vals(D);
	VectorXd::Map(&e_vals[0], D) = es.eigenvalues();

	// convert Eigen eigenvectors to a vector of STD vectors
	vector< vector<double> > e_vecs(D);
	for (i=0; i<D; i++) e_vecs[i] = vector<double>(D);
	for (i=0; i<D; i++) VectorXd::Map(&e_vecs[i][0], D) = es.eigenvectors().col(i);

	return {e_vals, e_vecs};
}

void Wavelet::Unwrap_Phase()
{	// unwrap the phase if it overshoots 0 to 2pi
	double phi0 = this->paramsND[IDX_phi0];
	if (phi0 > phi0_hi) phi0 = fmod(phi0,2*M_PI);
	if (phi0 < phi0_lo)
	{
		phi0 = -phi0;
		phi0 = fmod(phi0,2*M_PI);
		phi0 = 2*M_PI - phi0;
	}
	this->paramsND[IDX_phi0] = phi0;
}

} /* namespace wv */
