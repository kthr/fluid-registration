/*
 * Fourier.hpp
 *
 *  Created on: Jul 12, 2012
 *      Author: kthierbach
 */

#ifndef FOURIER_HPP_
#define FOURIER_HPP_

#include <math.h>
#include "../../lib/fftw3.h"

void Fourier(fftw_complex *out, fftw_complex *in, int nx, int ny);
void Fourier(fftw_complex *out, double *in, int nx, int ny);

#endif /* FOURIER_HPP_ */
