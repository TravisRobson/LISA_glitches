/*
 * G_B_MCMC.h
 *
 *  Created on: Sep 7, 2018
 *      Author: travisrobson
 */

#ifndef G_B_MCMC_H_
#define G_B_MCMC_H_

#include "Wavelet.h"
#include "LISA.h"
#include "TDI.h"

using namespace wv;
using namespace tdi;
using namespace ls;

namespace mc {

class Model {
public:
	Model(Wavelet x);
	virtual ~Model();

	Wavelet x;

	void get_logL(TDI tdi, LISA *lisa);
private:

};

class Proposal_Fisher {
public:
	Proposal_Fisher();
	virtual ~Proposal_Fisher();
};

} /* namespace mcmc */

#endif /* G_B_MCMC_H_ */
