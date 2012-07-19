/*
 * SinFourier.hpp
 *
 *  Created on: Jul 12, 2012
 *      Author: kthierbach
 */

#ifndef SINFOURIER_HPP_
#define SINFOURIER_HPP_

#include <math.h>
#include "../../lib/fftw3.h"

void SinFourier(double*out, double*in, int nx, int ny);

#endif /* SINFOURIER_HPP_ */
