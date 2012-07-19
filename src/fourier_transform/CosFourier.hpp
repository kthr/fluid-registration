/*
 * CosFourier.hpp
 *
 *  Created on: Jul 12, 2012
 *      Author: kthierbach
 */

#ifndef COSFOURIER_HPP_
#define COSFOURIER_HPP_

#include <math.h>
#include "../../lib/fftw3.h"

void CosFourier(double *out, double*in, int nx, int ny);

#endif /* COSFOURIER_HPP_ */
