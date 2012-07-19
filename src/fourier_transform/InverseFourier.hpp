/*
 * InverseFourier.hpp
 *
 *  Created on: Jul 12, 2012
 *      Author: kthierbach
 */

#ifndef INVERSEFOURIER_HPP_
#define INVERSEFOURIER_HPP_

#include <math.h>
#include "../../lib/fftw3.h"

void InverseFourier(fftw_complex*out, fftw_complex*in, int nx, int ny);
void InverseFourier(double*out, fftw_complex*in, int nx, int ny);

#endif /* INVERSEFOURIER_HPP_ */
