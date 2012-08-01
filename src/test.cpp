/*
 * test.cpp
 *
 *  Created on: Jul 10, 2012
 *      Author: kthierbach
 */

#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include "registration/FluidCurvatureRegistration.hpp"
#include "../lib/CImg-1.5.0/CImg.h"

using namespace cimg_library;

int main(int argc, char *argv[])
{
	CImg<double> *tImage;
	CImg<double> *sImage;
	double t_end = 64., dt_start = 0.001, alpha = 2, //SmoothWeight
			viscosity = 1., localDamping = 1., vortexWeight = 0., mu = 1., //LameMu
			lambda = 0.25, //LameLambda
			localError = .5, mismatch = 0.0005; //MismatchError
	int method = 7, //boundary
			returnType = 1;
	int option;
	if(argc < 3)
	{
		fprintf(stderr,"Usage: test [OPTIONS]... templateImage sampleImage\n");
		return 1;
	}
	while ((option = getopt(argc, argv, "t:e:m:")) != EOF)
	{
		switch (option)
		{
			case 't':
				if ((t_end = atof(optarg)) == 0.)
					fprintf(stderr, "Invalid argument for option -%c.\n", optopt);
				break;
			case 'e':
				if ((mismatch = atof(optarg)) == 0.)
					fprintf(stderr, "Invalid argument for option -%c.\n", optopt);
				break;
			case 'm':
				if ((alpha = atof(optarg)) == 0.)
					fprintf(stderr, "Invalid argument for option -%c.\n", optopt);
				break;
			case '?':
				if (optopt == 't' || optopt == 'e' || optopt == 'm')
					fprintf(stderr, "Option -%c requires an argument.\n", optopt);
				else if (isprint(optopt))
					fprintf(stderr, "Unknown option `-%c'.\n", optopt);
				else
					fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
				return 1;
			default:
				abort();
		}
	}
	std::cout << "Time: " << t_end << ", MismatchError: " << mismatch << ", SmoothWeight: " << alpha << "\n";
	tImage = new CImg<double>(argv[argc-2]);
	sImage = new CImg<double>(argv[argc-1]);

	*tImage /= tImage->max();
	*sImage /= sImage->max();
	FluidCurvatureRegistration *reg = new FluidCurvatureRegistration(method, tImage->width(), tImage->height(),
			tImage->data(), (long int) tImage->size(), sImage->data(), (long int) sImage->size(), t_end, dt_start,
			alpha, viscosity, localDamping, mu, lambda, vortexWeight, localError, mismatch, returnType, NULL, 1, NULL,
			0);
	reg->registerImages();
	delete tImage;
	delete sImage;
}
