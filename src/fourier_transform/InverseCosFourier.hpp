/*
 * InverseCosFourier.hpp
 *
 *  Created on: Jul 12, 2012
 *      Author: kthierbach
 */

#ifndef INVERSECOSFOURIER_HPP_
#define INVERSECOSFOURIER_HPP_

#include <math.h>
#include "../../lib/fftw3.h"

void InverseCosFourier(double*out, double*in, int nx, int ny);

#endif /* INVERSECOSFOURIER_HPP_ */
