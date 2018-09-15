/*
 * IO.h
 *
 *  Created on: Sep 14, 2018
 *      Author: travisrobson
 */

#ifndef IO_H_
#define IO_H_

#include <string>

using namespace std;

#include "LISA.h"
#include "G_B_MCMC.h"
#include "Wavelet.h"

using namespace wv;
using namespace ls;
using namespace mc;

struct Files
{
	string Temperature_file;
};


vector<Model*> parse_config_file(string input_file, double T, LISA *lisa, struct Files *Files);




#endif /* IO_H_ */
