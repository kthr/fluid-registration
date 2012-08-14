/*
 * test.cpp
 *
 *  Created on: Jul 10, 2012
 *      Author: kthierbach
 */

#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include "fluidReg.hpp";
#include "registration/FluidCurvatureRegistration.hpp"
#include "../lib/CImg-1.5.0/CImg.h"

using namespace cimg_library;
using namespace std;

int main(int argc, char *argv[])
{
	CImg<double> *templateImage, *sampleImage, *referenceImage = NULL, *patternImage = NULL;
	FluidCurvatureRegistration *reg;
	double t_end = 64., dt_start = 0.001, alpha = 2, //SmoothWeight
			viscosity = 1., localDamping = 1., vortexWeight = 0., mu = 1., //LameMu
			lambda = 0.25, //LameLambda
			localError = .5, mismatch = 0.0005; //MismatchError
	int method = 7, //boundary
			returnType = 1;
	int option;
	char *flowFile = NULL, *patternFile = NULL;

	cimg::exception_mode(0);
	while ((option = getopt(argc, argv, "f:hm:p:r:s:t:")) != EOF)
	{
		switch (option)
		{
			case 'f':
				if (optarg != NULL)
				{
					flowFile = optarg;
				}
				else
				{
					fprintf(stderr, "Invalid or missing argument for option -%c.\n", optopt);
					return 1;
				}
				break;
			case 'h':
				#ifdef REVISION
					cout << "Revision: " << REVISION << "\n";
				#else
					cout << "Revision: Unknown" << "\n";
				#endif
				cout << "Usage: fluidReg [OPTIONS]... templateImage sampleImage\n";
				cout << "Example: fluidReg -t 64 -m 0.00001 -s 200 -f flow.dat -p pattern.png template.png sample.png\n\n";
				cout << "Registration parameters:\n";
				cout << "\t -f file \t the file where to save the computed flow field\n";
				cout << "\t -h \t\t prints the help message\n";
				cout << "\t -m NUM \t specifies mismatch error for the two images (default=0.0005)\n";
				cout << "\t -p file \t the file where to save the displaced sample or reference pattern\n";
				cout << "\t -r file \t the reference image\n";
				cout << "\t -s NUM \t specifies the smooth weight which controls the viscosity (default=2.)\n";
				cout << "\t -t NUM \t specifies the final time of the iteration process (default=64.)\n";
				return 0;
			case 'm':
				if ((mismatch = atof(optarg)) == 0.)
				{
					fprintf(stderr, "Invalid or missing argument for option -%c.\n", optopt);
					return 1;
				}
				break;
			case 'p':
				if (optarg != NULL)
				{
					patternFile = optarg;
				}
				else
				{
					fprintf(stderr, "Invalid or missing argument for option -%c.\n", optopt);
					return 1;
				}
				break;
			case 'r':
				try
				{
					referenceImage = new CImg<double>(optarg);
				} catch (CImgException &e)
				{
					fprintf(stderr, "Failed to open reference image, continuing without a reference image.\n");
				}
				break;
			case 's':
				if ((alpha = atof(optarg)) == 0.)
				{
					fprintf(stderr, "Invalid or missing argument for option -%c.\n", optopt);
					return 1;
				}
				break;
			case 't':
				if ((t_end = atof(optarg)) == 0.)
				{
					fprintf(stderr, "Invalid or missing argument for option -%c.\n", optopt);
					return 1;
				}
				break;
			case '?':
				if (optopt == 't' || optopt == 'm' || optopt == 's')
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
	if (argc < 3)
	{
		fprintf(stderr, "Usage: fluidReg [OPTIONS]... templateImage sampleImage\n");
		fprintf(stderr, "Try `fluidReg -h' for more information.\n");
		return 1;
	}
	try
	{
		templateImage = new CImg<double>(argv[argc - 2]);
	} catch (CImgException &e)
	{
		char message[3072];
		strlcat(message, "Failed to open template image '", 1024);
		strlcat(message, argv[argc - 2], 1024);
		strlcat(message, "'.\n", 1024);
		fprintf(stderr, "%s", message);
		return 1;
	}
	try
	{
		sampleImage = new CImg<double>(argv[argc - 1]);
	} catch (CImgException &e)
	{
		char message[3072];
		strlcat(message, "Failed to open sample image '", 1024);
		strlcat(message, argv[argc - 1], 1024);
		strlcat(message, "'.\n", 1024);
		fprintf(stderr, "%s", message);
		return 1;
	}
	std::cout << "Time: " << t_end << ", MismatchError: " << mismatch << ", SmoothWeight: " << alpha << "\n";

	/*
	 *templateImage /= templateImage->max();
	 *sampleImage /= sampleImage->max();
	 */

	if (referenceImage)
	{
		reg = new FluidCurvatureRegistration(method, templateImage->width(), templateImage->height(),
				templateImage->data(), (long int) templateImage->size(), sampleImage->data(),
				(long int) sampleImage->size(), t_end, dt_start, alpha, viscosity, localDamping, mu, lambda,
				vortexWeight, localError, mismatch, returnType, NULL, referenceImage->spectrum(),
				referenceImage->data(), (long int) referenceImage->size());
		reg->registerImages();
	}
	else
	{
		reg = new FluidCurvatureRegistration(method, templateImage->width(), templateImage->height(),
				templateImage->data(), (long int) templateImage->size(), sampleImage->data(),
				(long int) sampleImage->size(), t_end, dt_start, alpha, viscosity, localDamping, mu, lambda,
				vortexWeight, localError, mismatch, returnType, NULL, 1, NULL, 0);
		reg->registerImages();
	}

	if (referenceImage != NULL)
	{
		patternImage = new CImg<double>(reg->getReference()->bm, reg->getReference()->ny, reg->getReference()->nx);
		if (patternFile != NULL)
		{
			try
			{
				patternImage->save(patternFile);
			} catch (CImgException &e)
			{
				char message[3072];
				strlcat(message, "Failed to save reference image at '", 1024);
				strlcat(message, patternFile, 1024);
				strlcat(message, "'.\n", 1024);
				fprintf(stderr, "%s", message);
				return 1;
			}
		}
		else
		{
			patternImage->display("reference image");
		}
	}

	if (flowFile != NULL)
	{
		if (!reg->getFlowField()->save(flowFile))
		{
			char message[3072];
			strlcat(message, "Failed to save flow field at '", 1024);
			strlcat(message, flowFile, 1024);
			strlcat(message, "'.\n", 1024);
			fprintf(stderr, "%s", message);
			return 1;
		}
	}

	delete reg;
	delete patternImage;
	delete templateImage;
	delete sampleImage;
	delete referenceImage;
}
