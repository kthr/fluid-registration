/*
 * Fourier.cpp
 *
 *  Created on: Jul 12, 2012
 *      Author: kthierbach
 */

#include "Fourier.hpp"

void Fourier(fftw_complex *out, fftw_complex *in, int nx, int ny)
{
	fftw_plan plan;

	plan = fftw_plan_dft_2d(ny, nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	if (!plan)
	{
		fftw_cleanup();
		plan = fftw_plan_dft_2d(ny, nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	}
	fftw_execute(plan);
	fftw_destroy_plan(plan);
}

void Fourier(fftw_complex *out, double *in, int nx, int ny)
{
	fftw_plan plan;
	fftw_complex*tmp;
	int i, mu = nx * ny;
	double nrm;

	nrm = 1.0 / sqrt((double) nx * ny);

	tmp = (fftw_complex*) fftw_malloc(nx * ny * sizeof(fftw_complex));
	for (i = 0; i < mu; i++)
	{
		tmp[i][0] = in[i] * nrm;
		tmp[i][1] = 0.0;
	}
	plan = fftw_plan_dft_2d(ny, nx, tmp, out, FFTW_FORWARD, FFTW_ESTIMATE);

	fftw_execute(plan);
	fftw_destroy_plan(plan);
	fftw_free(tmp);
}

