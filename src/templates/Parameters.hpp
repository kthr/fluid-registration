/*
 * Parameters.hpp
 *
 *  Created on: Oct 18, 2012
 *      Author: kthierbach
 */

#ifndef PARAMETERS_HPP_
#define PARAMETERS_HPP_

#include <string>

struct parameters{
	double end;
	double error;
	double alpha;
	double vortex_weight;
	double mu;
	double lambda;
	std::string boundary;
	std::string method;
	double actual_error;
	double actual_time;
};

#endif /* PARAMETERS_HPP_ */
