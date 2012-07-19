/*
 * InverseSinFourier.hpp
 *
 *  Created on: Jul 12, 2012
 *      Author: kthierbach
 */

#ifndef INVERSESINFOURIER_HPP_
#define INVERSESINFOURIER_HPP_

#include <math.h>
#include "../../lib/fftw3.h"

void InverseSinFourier(double*out, double*in, int nx, int ny);

#endif /* INVERSESINFOURIER_HPP_ */
