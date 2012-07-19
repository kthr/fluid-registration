/*
 * InverseCosFourier.cpp
 *
 *  Created on: Jul 12, 2012
 *      Author: kthierbach
 */

#include "InverseCosFourier.hpp"

void InverseCosFourier(double*out, double*in, int nx, int ny)
{
	fftw_plan plan;
	int i, mu = nx * ny;
	double norm;

	plan = fftw_plan_r2r_2d(ny, nx, in, out, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
	norm = 1.0 / sqrt(4.0 * nx * ny);
	for (i = 0; i < mu; i++)
		out[i] *= norm;

}
