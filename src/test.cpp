/*
 * test.cpp
 *
 *  Created on: Jul 10, 2012
 *      Author: kthierbach
 */

#include "registration/FluidCurvatureRegistration.hpp"
#include "../lib/CImg-1.5.0/CImg.h"
using namespace cimg_library;

int main()
{
	CImg<double> *tImage = new CImg<double>("./U.png");
	CImg<double> *sImage = new CImg<double>("./X.png");

	printf("tImage channels: %d dx: %d dy: %d\n", tImage->spectrum(), tImage->width(), tImage->height());
	printf("sImage channels: %d dx: %d dy: %d\n", sImage->spectrum(), sImage->width(), sImage->height());
	double 	t_end = 64.,
			dt_start = 0.001,
			alpha = 2000, //SmoothWeight
			viscosity = 1.,
			localDamping = 1.,
			vortexWeight = 0.,
			mu = 1., //LameMu
			lambda = 0.25, //LameLambda
			localError = .5,
			mismatch = 0.000015; //MismatchError
	int 	method = 7, //boundary
			returnType = 1;

	*tImage /= tImage->max();
	*sImage /= sImage->max();
	FluidCurvatureRegistration *reg = new FluidCurvatureRegistration(method, tImage->width(), tImage->height(),
			tImage->data(), (long int) tImage->size(), sImage->data(), (long int) sImage->size(), t_end, dt_start, alpha,
			viscosity, localDamping, mu, lambda, vortexWeight, localError, mismatch,returnType, NULL, 1, NULL, 0);
	reg->registerImages();
	delete tImage;
	delete sImage;
}
