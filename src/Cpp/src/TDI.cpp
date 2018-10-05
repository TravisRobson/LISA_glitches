/*
 * TDI.cpp
 *
 *  Created on: Sep 6, 2018
 *      Author: travisrobson
 */

#include <list>
#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <algorithm>
using namespace std;

#include "TDI.h"
#include "LISA.h"
using namespace ls;

namespace tdi {

TDI::TDI() {

}

TDI::~TDI() {
}

TDI::TDI(const TDI &src)
{
	this->N_lo = src.get_N_lo();
	this->N_hi = src.get_N_hi();

	// deep copy the TDI data channels
	this->X = vector<complex<double>>(src.X);
	this->Y = vector<complex<double>>(src.Y);
	this->Z = vector<complex<double>>(src.Z);

	this->A = vector<complex<double>>(src.A);
	this->E = vector<complex<double>>(src.E);
	this->T = vector<complex<double>>(src.T);
}

TDI TDI::operator-(const TDI& a)
{
	// find the range of overlapping frequency bins
	int N_lo = max(a.get_N_lo(), this->get_N_lo());
	int N_hi = min(a.get_N_hi(), this->get_N_hi());

	//if (N_lo > N_hi) return 1.0; // TODO: the signals do not overlap

	int dN_lo_a    = N_lo - a.get_N_lo();
	int dN_lo_this = N_lo - this->get_N_lo();

	// construct the TDI to be returned
	TDI result;
	result.set_N_lo(N_lo);
	result.set_N_hi(N_hi);
	result.X = vector<complex<double>>(a.X);
	result.Y = vector<complex<double>>(a.Y);
	result.Z = vector<complex<double>>(a.Z);
	result.A = vector<complex<double>>(a.A);
	result.E = vector<complex<double>>(a.E);
	result.T = vector<complex<double>>(a.T);

	// form up the difference
	for (int i=0; i<(N_hi-N_lo); i++)
	{
		result.X[i] = this->X[dN_lo_this + i] - a.X[dN_lo_a + i];
		result.Y[i] = this->Y[dN_lo_this + i] - a.Y[dN_lo_a + i];
		result.Z[i] = this->Z[dN_lo_this + i] - a.Z[dN_lo_a + i];

		result.A[i] = this->A[dN_lo_this + i] - a.A[dN_lo_a + i];
		result.E[i] = this->E[dN_lo_this + i] - a.E[dN_lo_a + i];
		result.T[i] = this->T[dN_lo_this + i] - a.T[dN_lo_a + i];
	}

	return result;
}

TDI TDI::operator/=(double num)
{    //divide by num
	for (int i=0; i<this->X.size(); i++)
	{
		this->X[i] /= num;
		this->Y[i] /= num;
		this->Z[i] /= num;

		this->A[i] /= num;
		this->E[i] /= num;
		this->T[i] /= num;
	}

	return *this;
}

void TDI::set_N_lo(int N_lo)
{
	this->N_lo = N_lo;
}

void TDI::set_N_hi(int N_hi)
{
	this->N_hi = N_hi;
}

void TDI::phase_to_tdi(list<vector<complex<double>>> *phase_list, int N_lo, double T)
{	// Convert phase measurements to first generation TDI

	// extract phases
	vector<complex<double>> p12, p21, p13, p31, p23, p32;
	p12 = phase_list->front(); phase_list->pop_front();
	p21 = phase_list->front(); phase_list->pop_front();
	p13 = phase_list->front(); phase_list->pop_front();
	p31 = phase_list->front(); phase_list->pop_front();
	p23 = phase_list->front(); phase_list->pop_front();
	p32 = phase_list->front(); phase_list->pop_front();

	int N = p12.size(); // number of frequency samples

	this->X = vector<complex<double>>(N);
	this->Y = vector<complex<double>>(N);
	this->Z = vector<complex<double>>(N);

	this->A = vector<complex<double>>(N);
	this->E = vector<complex<double>>(N);
	this->T = vector<complex<double>>(N);

	complex<double> phase3, phase2, phase1;
	complex<double> j (0, 1.0); // imaginary number
	double fonfs;

	for (int i=0; i<N; i++)
	{
		fonfs  = (i+N_lo)/T/fstar;

//		phase3 = exp(-j*3.*fonfs);
//		phase2 = exp(-j*2.*fonfs);
//		phase1 = exp(-j*1.*fonfs);

		phase1 = exp(-j*fonfs);
		phase2 = phase1*phase1;
		phase3 = phase2*phase1;

		this->X[i] = (p12[i] - p13[i])*phase3 + (p21[i] - p31[i])*phase2 +
					 (p13[i] - p12[i])*phase1 + (p31[i] - p21[i]);
		this->Y[i] = (p23[i] - p21[i])*phase3 + (p32[i] - p12[i])*phase2 +
					 (p21[i] - p23[i])*phase1 + (p12[i] - p32[i]);
		this->Z[i] = (p31[i] - p32[i])*phase3 + (p13[i] - p23[i])*phase2 +
					 (p32[i] - p31[i])*phase1 + (p23[i] - p13[i]);

		this->A[i] = (2.*this->X[i] - this->Y[i] - this->Z[i])/3.;
		this->E[i] = (this->Z[i] - this->Y[i])/sqrt(3.);
		this->T[i] = (this->X[i] + this->Y[i] + this->Z[i])/3.;
	}

}

double nwip(TDI *a, TDI *b, LISA *lisa, int X_flag)
{	// Calculate the overlap between two arbitrary TDI objects

	double fi, SnAEi, SnTi, SnXi;
	double overlap = 0.0;
	double T = lisa->get_T();

	// find the range of overlapping frequency bins
	int N_lo = max(a->get_N_lo(), b->get_N_lo());
	int N_hi = min(a->get_N_hi(), b->get_N_hi());

	if (N_lo >= N_hi) return 0.0; // the signals do not overlap

	int dN_lo_a = N_lo - a->get_N_lo();
	int dN_lo_b = N_lo - b->get_N_lo();

	for (int i=0; i<(N_hi-N_lo); i++)
	{
		fi = (N_lo + i)/T;

		if (X_flag == 0)
		{
			SnAEi = lisa->SnAE(fi);
			SnTi  = lisa->SnT (fi);

			if (SnAEi != SnAEi or SnTi != SnTi)
			{
				SnAEi = INFINITY;
				SnTi  = INFINITY;
			}

			overlap += real(a->A[dN_lo_a + i]*conj(b->A[dN_lo_b + i]) + a->E[dN_lo_a + i]*conj(b->E[dN_lo_b + i]))/SnAEi;
			overlap += real(a->T[dN_lo_a + i]*conj(b->T[dN_lo_b + i]))/SnTi;
		}
		else
		{
			SnXi = lisa->SnX(fi);

			if (SnXi != SnXi) SnXi = INFINITY;

			overlap += real( a->X[dN_lo_a + i]*conj(b->X[dN_lo_b + i]) )/SnXi;
		}
	}

	return 4*overlap/T;
}

} /* namespace tdi */
