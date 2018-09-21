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
#include <limits>
#include <algorithm>

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
	double rho_est = sqrt(sqrt(M_PI/2)*tau/lisa->SnX(f0))*A*2*8*pow(sin(f0/fstar), 2.0); // I weirdly canceled the transfer functions...
	//cout << rho_est << endl;

	// Lowest frequency bin
	double df = 2.0/tau*pow(rho_est/5, 1.0);
	int N_lo = (int)((f0 - df)*T);
	if (N_lo <= 0) N_lo = 1; // make this the lowest positive frequency (don't want to break noise curve code)
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
	vector<double> k {-sin_theta*cos_phi, -sin_theta*sin_phi, -cos_theta};

	vector<double> v1(3);
	vector<vector<double>> e_plus(3, v1), e_cross(3, v1);
	int i, j;
	double t1, t2;
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			t1 = u[i]*u[j] - v[i]*v[j];
			t2 = u[i]*v[j] + u[j]*v[i];
			e_plus[i][j]  = cos_2psi*t1 - sin_2psi*t2;
			e_cross[i][j] = sin_2psi*t1 + cos_2psi*t2;
		}
	}

	vector<double> x(3), y(3), z(3);
	lisa->SC_position_analytic(t0, &x, &y, &z); // find positions at central time, assuming static LISA

	vector<double> r12(3), r21(3), r13(3), r31(3), r23(3), r32(3);

	r12[0] = (x[1] - x[0])/Larm;
	r12[1] = (y[1] - y[0])/Larm;
	r12[2] = (z[1] - z[0])/Larm;

	r13[0] = (x[2] - x[0])/Larm;
	r13[1] = (y[2] - y[0])/Larm;
	r13[2] = (z[2] - z[0])/Larm;

	r23[0] = (x[2] - x[1])/Larm;
	r23[1] = (y[2] - y[1])/Larm;
	r23[2] = (z[2] - z[1])/Larm;

	for (i=0; i<3; i++)
	{
		r21[i] = -r12[i];
		r31[i] = -r13[i];
		r32[i] = -r23[i];
	}

	double k_dot_r12 = 0.0;
	double k_dot_r21 = 0.0;
	double k_dot_r13 = 0.0;
	double k_dot_r31 = 0.0;
	double k_dot_r23 = 0.0;
	double k_dot_r32 = 0.0;
	for (i=0; i<3; i++)
	{
		k_dot_r12 += k[i]*r12[i];
		k_dot_r13 += k[i]*r13[i];
		k_dot_r23 += k[i]*r23[i];
	}
	k_dot_r21 = -k_dot_r12;
	k_dot_r31 = -k_dot_r13;
	k_dot_r32 = -k_dot_r23;

	double k_dot_x1 = (k[0]*x[0] + k[1]*y[0] + k[2]*z[0]);
	double k_dot_x2 = (k[0]*x[1] + k[1]*y[1] + k[2]*z[1]);
	double k_dot_x3 = (k[0]*x[2] + k[1]*y[2] + k[2]*z[2]);

	double fp12 = 0; double fp21 = 0; double fp13 = 0; double fp31 = 0; double fp23 = 0; double fp32 = 0;
	double fc12 = 0; double fc21 = 0; double fc13 = 0; double fc31 = 0; double fc23 = 0; double fc32 = 0;
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			fp12 += r12[i]*r12[j]*e_plus[i][j];
			fp13 += r13[i]*r13[j]*e_plus[i][j];
			fp23 += r23[i]*r23[j]*e_plus[i][j];

			fc12 += r12[i]*r12[j]*e_cross[i][j];
			fc13 += r13[i]*r13[j]*e_cross[i][j];
			fc23 += r23[i]*r23[j]*e_cross[i][j];
		}
	}

//	fp12 = cos(2*psi + ellip);
//	fp13 = cos(2*psi + ellip);
//	fp23 = cos(2*psi + ellip);
//
//	fc12 = cos(2*psi + ellip);
//	fc13 = cos(2*psi + ellip);
//	fc23 = cos(2*psi + ellip);

	fp21 = fp12;
	fp31 = fp13;
	fp32 = fp23;

	fc21 = fc12;
	fc31 = fc13;
	fc32 = fc23;

	fp12 /= (1 - k_dot_r12);
	fp13 /= (1 - k_dot_r13);
	fp23 /= (1 - k_dot_r23);
	fp21 /= (1 - k_dot_r21);
	fp31 /= (1 - k_dot_r31);
	fp32 /= (1 - k_dot_r32);

	fc12 /= (1 - k_dot_r12);
	fc13 /= (1 - k_dot_r13);
	fc23 /= (1 - k_dot_r23);
	fc21 /= (1 - k_dot_r21);
	fc31 /= (1 - k_dot_r31);
	fc32 /= (1 - k_dot_r32);

	complex<double> jj(0,1.0); // imaginary number

	double f;
	complex<double> t3;

	for (i=0; i<N_hi-N_lo; i++)
	{
		f = (i+N_lo)/T;

		t1 = 2*M_PI*f/Clight;//*100.;
		t3 = Psi_FT(f, A, f0, t0, tau, phi0);

		p12[i]  = (fp12 + ellip*fc12)*t3*( exp(-t1*(-Larm + k_dot_x2)*jj) - exp(-t1*k_dot_x1*jj) );
		p13[i]  = (fp13 + ellip*fc13)*t3*( exp(-t1*(-Larm + k_dot_x3)*jj) - exp(-t1*k_dot_x1*jj) );
		p23[i]  = (fp23 + ellip*fc23)*t3*( exp(-t1*(-Larm + k_dot_x3)*jj) - exp(-t1*k_dot_x2*jj) );

		p21[i]  = (fp21 + ellip*fc21)*t3*( exp(-t1*(-Larm + k_dot_x1)*jj) - exp(-t1*k_dot_x2*jj) );
		p31[i]  = (fp31 + ellip*fc31)*t3*( exp(-t1*(-Larm + k_dot_x1)*jj) - exp(-t1*k_dot_x3*jj) );
		p32[i]  = (fp32 + ellip*fc32)*t3*( exp(-t1*(-Larm + k_dot_x2)*jj) - exp(-t1*k_dot_x3*jj) );
	}

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
	if (N_lo <= 0) N_lo = 1; // make this the lowest positive frequency (don't want to break noise curve code)
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

	//////////// Use the modified Cholesky decomposition to ensure we have a positive definite matrix //////////
	vector<double> v1(D,0);
	vector<vector<double>> l(D, v1), d(D, v1), c(D, v1), A(D, v1), temp(D, v1);

	A = Fisher; // make a copy

	double mac_epsilon = numeric_limits<double>::epsilon();
	for (i=0; i<D; i++) l[i][i] = 1.0;

	// Find the largest diagonal value
	double gamma = 0.0;
	for (i=0; i<D; i++)
	{
		if (A[i][i] > gamma) gamma = A[i][i];
	}

	// find the largest off-diagonal value
	double xi = 0.0;
	for (i=0; i<D; i++)
	{
		for(j=i+1; j<D; j++)
		{
			if (fabs(A[i][j]) > xi) xi = fabs(A[i][j]);
		}
	}

	double delta = mac_epsilon*max(xi + gamma, 1.0);
	double beta = sqrt(   max( gamma, max(mac_epsilon, xi/sqrt(D*D-1.0)) )   );

	for (j=0; j<D; j++)
	{
		c[j][j] = A[j][j];
		for (int s=0; s<j; s++) c[j][j] -= d[s][s]*pow(l[j][s], 2.0);

		double dummy = max(delta, fabs(c[j][j]));

		double theta_j = 0.0;
		if (j<= D)
		{
			for (int s=j+1; s<D; s++) theta_j = max(theta_j, fabs(c[s][j]));
		}

		d[j][j] = max(dummy, pow(theta_j/beta, 2.0)); //c[j][j]; <-- plugged in for normal Cholesky decomp

		for (i=j+1; i<D; i++)
		{
			c[i][j] = A[i][j];
			for (int s=0; s<D; s++) c[i][j] -= d[s][s]*l[i][s]*l[j][s];
			l[i][j] = c[i][j]/d[j][j];
		}
	}

	for (i=0; i<D; i++)
	{
		for (j=0; j<D; j++)
		{
			for (int m=0; m<D; m++)
			{
				temp[i][j] += d[i][m]*l[j][m]; // i.e. transpose of L
			}
		}
	}
	for (i=0; i<D; i++)
	{
		for (j=0; j<D; j++)
		{
			Fisher[i][j] = 0.0;
			for (int m=0; m<D; m++)
			{
				Fisher[i][j] += l[i][m]*temp[m][j];
			}
		}
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

//	cout << "Eigenvalues.... " << es.eigenvalues() << endl << endl;

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

	// unwrap sky phi
	if (this->paramsND.size() > 5)
	{
		double phi = this->paramsND[IDX_phi];
		if (phi > phi_hi) phi = fmod(phi,2*M_PI);
		if (phi < phi_lo)
		{
			phi = -phi;
			phi = fmod(phi,2*M_PI);
			phi = 2*M_PI - phi;
		}
		this->paramsND[IDX_phi] = phi;
	}
}

} /* namespace wv */
