/*
 * InverseFourier.cpp
 *
 *  Created on: Jul 12, 2012
 *      Author: kthierbach
 */

#include "InverseFourier.hpp"

void InverseFourier(fftw_complex *out, fftw_complex*in, int nx, int ny)
{
	fftw_plan plan;

	plan = fftw_plan_dft_2d(ny, nx, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
}

void InverseFourier(double *out, fftw_complex*in, int nx, int ny)
{
	fftw_plan plan;
	fftw_complex*tmp;
	double nrm;
	int i, mu = nx * ny;

	nrm = 1.0 / sqrt((double) nx * ny);
	tmp = (fftw_complex*) fftw_malloc(nx * ny * sizeof(fftw_complex));
	plan = fftw_plan_dft_2d(ny, nx, in, tmp, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
	for (i = 0; i < mu; i++)
	{
		out[i] = tmp[i][0] * nrm;
	}
	fftw_free(tmp);
}
