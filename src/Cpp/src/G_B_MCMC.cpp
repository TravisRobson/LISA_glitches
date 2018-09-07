/*
 * G_B_MCMC.cpp
 *
 *  Created on: Sep 7, 2018
 *      Author: travisrobson
 */

#include "G_B_MCMC.h"

namespace mc {

Model::Model(Wavelet x)
{
	this->x = x;
}

Model::~Model()
{

}

void Model::get_logL(TDI tdi, LISA *lisa)
{


}

Proposal_Fisher::Proposal_Fisher()
{

}

Proposal_Fisher::~Proposal_Fisher()
{

}

} /* namespace mcmc */
