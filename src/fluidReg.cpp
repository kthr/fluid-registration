/*
 * test.cpp
 *
 *  Created on: Jul 10, 2012
 *      Author: kthierbach
 */

#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include "fluidReg.hpp"
#include "registration/FluidCurvatureRegistration.hpp"
#include "templates/Parameters.hpp"
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
	string boundaries[] = {"Periodic","ZeroDerivative","ZeroDisplacement"};
	string methods[] = {"NavierLame","OverDampedCurvature","OverDampedDiffusion"};
	int boundary = 0;
	int method = 6;
	int option, option_index;
	char *flowFile = NULL, *patternFile = NULL;
	int verbose = 0;
	parameters param;

	struct option long_options[] =
	{
		{ "boundary-condition", required_argument, NULL, '1' },
		{ "flow-file", required_argument, NULL, '2' },
		{ "method", required_argument, NULL, '3' },
		{ "pattern-file", required_argument, NULL, '4' },
		{ "reference-image", required_argument, NULL, '5' },
		{ "verbose", no_argument, &verbose, 1 },
		{ 0, 0, 0, 0 }
	};

	cimg::exception_mode(0);

	while (1)
	{

		option = getopt_long(argc, argv, "d:e:hl:m:s:t:v:w:", long_options, &option_index);
		if (option == -1)
		{
			break;
		}
		switch (option)
		{
			case 0: //flag option
				/*option sets a flag*/
				break;
			case '1': //boundary condition
				if (optarg != NULL)
				{
					boundary = -1;
					for (int i = 0; i < 3; ++i)
					{
						if (boundaries[i].compare(optarg) == 0)
						{
							switch (i)
							{
								case 0: //Periodic
									boundary = 0;
									break;
								case 1: //ZeroDerivative
									boundary = 1;
									break;
								case 2: //ZeroDisplacement
									boundary = 2;
									break;
							}
						}
					}
					if (boundary == -1)
					{
						fprintf(stderr,
								"Invalid argument for option --%s. Possible options are 'NavierLame', 'OverDampedCurvature' and 'OverDampedDiffusion'. Continuing with default.\n",
								long_options[option_index].name);
						boundary = 0;
					}
				}
				else
				{
					fprintf(stderr, "Missing argument for option --%s.\n", long_options[option_index].name);
					return EXIT_FAILURE;
				}
				break;
			case '2': //flow file
				if (optarg != NULL)
				{
					flowFile = optarg;
				}
				else
				{
					fprintf(stderr, "Missing argument for option --%s.\n", long_options[option_index].name);
					return EXIT_FAILURE;
				}
				break;
			case '3': //method
				if (optarg != NULL)
				{
					method = -1;
					for(int i=0; i<3; ++i)
					{
						if(methods[i].compare(optarg) == 0)
						{
							switch(i)
							{
								case 0: //NavierLame
									method = 6;
									break;
								case 1: //OverDampedCurvature
									method = 7;
									break;
								case 2: //OverDampedDiffusion
									method = 10;
									break;
							}
						}
					}
					if(method == -1)
					{
						fprintf(stderr, "Invalid argument for option --%s. Possible options are 'NavierLame', 'OverDampedCurvature' and 'OverDampedDiffusion'. Continuing with default.\n", long_options[option_index].name);
						method = 6;
					}
				}
				else
				{
					fprintf(stderr, "Missing argument for option -%s.\n", long_options[option_index].name);
					return EXIT_FAILURE;
				}
				break;
			case '4': //pattern file
				if (optarg != NULL)
				{
					patternFile = optarg;
				}
				else
				{
					fprintf(stderr, "Invalid or missing argument for option --%s.\n", long_options[option_index].name);
					return EXIT_FAILURE;
				}
				break;
			case '5': //reference image
				try
				{
					referenceImage = new CImg<double>(optarg);
				}
				catch (CImgException &e)
				{
					fprintf(stderr, "Failed to open reference image, continuing without a reference image.\n");
				}
				break;
			case 'd': //local damping
				if (optarg == NULL)
				{
					fprintf(stderr, "Invalid or missing argument for option -%c.\n", option);
					return EXIT_FAILURE;
				}
				localDamping = atof(optarg);
				break;
			case 'e': //mismatch error
				if ((mismatch = atof(optarg)) == 0.)
				{
					fprintf(stderr, "Invalid or missing argument for option -%c.\n", option);
					return EXIT_FAILURE;
				}
				break;
			case 'h': //help
#ifdef REVISION
				cout << "Revision: " << REVISION << endl;
#else
				cout << "Revision: Unknown" << endl;
#endif
				cout << "Usage: fluidReg [OPTIONS]... templateImage sampleImage" << endl;
				cout
						<< "Example: fluidReg -t 64 -e 0.00001 -s 200 --flow-file flow.dat --pattern-file pattern.png template.png sample.png\n\n";
				cout << "Registration parameters:" << endl;
				//cout << "\t -d NUM \t local damping parameter (default=1.)" << endl;
				cout << "\t -e NUM \t mismatch error for the two images (default=0.0005)" << endl;
				cout << "\t -h \t\t prints the help message" << endl;
				cout << "\t -l NUM \t lambda parameter (default=.25)" << endl;
				cout << "\t -m NUM \t lame mu parameter (default=1.)" << endl;
				cout << "\t -s NUM \t smooth weight parameter which controls the viscosity (default=2.)" << endl;
				cout << "\t -t NUM \t final time of the iteration process (default=64.)" << endl;
				//cout << "\t -v NUM \t viscosity parameter (default=1.)" << endl;
				cout << "\t -w NUM \t vortex weight parameter (default=0.)" << endl;
				cout << "\t --boundary-condition CONDITION \t the boundary condition can be 'Periodic', 'ZeroDerivative' or 'ZeroDisplacement' (default=Periodic)" << endl;
				cout << "\t --flow-file FILE \t the file where to save the computed flow field" << endl;
				cout << "\t --method METHOD \t the method used for registration (default=NavierLame). Values for METHOD are:" << endl;
				cout << "\t\t NavierLame \t\t with parameters lame mu (-m) and lambda (-l)" << endl;
				cout << "\t\t OverDampedCurvature \t with parameters smooth weight (-s) and vortex weight (-w)" << endl;
				cout << "\t\t OverDampedDiffusion \t with parameters smooth weight (-s) and vortex weight (-w)" << endl;
				cout << "\t --pattern-file FILE \t the file where to save the displaced sample or reference image" << endl;
				cout << "\t --reference-image FILE  the location of the reference image" << endl;
				cout << "\t --verbose \t\t turns on verbose output (transformed sample image is shown)" << endl;
				return EXIT_SUCCESS;
			case 'l': //lambda
				if (optarg == NULL)
				{
					fprintf(stderr, "Missing argument for option -%c.\n", option);
					return EXIT_FAILURE;
				}
				lambda = atof(optarg);
				break;
			case 'm': //lame mu
				if (optarg == NULL)
				{
					fprintf(stderr, "Missing argument for option -%c.\n", option);
					return EXIT_FAILURE;
				}
				mu = atof(optarg);
				break;
			case 's': //smooth weight
				if (optarg == NULL)
				{
					fprintf(stderr, "Missing argument for option -%c.\n", option);
					return EXIT_FAILURE;
				}
				alpha = atof(optarg);
				break;
			case 't': // max time
				if (optarg == NULL)
				{
					fprintf(stderr, "Missing argument for option -%c.\n", option);
					return EXIT_FAILURE;
				}
				t_end = atof(optarg);
				break;
			case 'v': //viscosity
				if (optarg == NULL)
				{
					fprintf(stderr, "Missing argument for option -%c.\n", option);
					return EXIT_FAILURE;
				}
				viscosity = atof(optarg);
				break;
			case 'w': //vortex weight
				if (optarg == NULL)
				{
					fprintf(stderr, "Missing argument for option -%c.\n", option);
					return EXIT_FAILURE;
				}
				vortexWeight = atof(optarg);
				break;
			case '?':
				/* getopt_long already printed an error message. */
				break;
			default:
				abort();
		}
	}
	if(method != 6)
	{
		method += boundary;
	}
	if (argc < 3)
	{
		fprintf(stderr, "Usage: fluidReg [OPTIONS]... templateImage sampleImage\n");
		fprintf(stderr, "Try `fluidReg -h' for more information.\n");
		return EXIT_FAILURE;
	}
	try
	{
		templateImage = new CImg<double>(argv[argc - 2]);
	} catch (CImgException &e)
	{
		fprintf(stderr, "Failed to open template image '%s'.\n", argv[argc - 2]);
		return EXIT_FAILURE;
	}
	try
	{
		sampleImage = new CImg<double>(argv[argc - 1]);
	} catch (CImgException &e)
	{
		fprintf(stderr, "Failed to open sample image '%s'.\n", argv[argc - 1]);
		return EXIT_FAILURE;
	}
	switch(method-boundary)
	{
		case 6://NavierLame
			cout << "Method: NavierLame" << ", BoundaryCondition: " << boundaries[boundary] <<", Time: " << t_end << ", MismatchError: " << mismatch << endl;
			cout << "LameMu: " << mu << ", Lambda: " << lambda << endl;
			break;
		case 7://OverdampedCurvature
			cout << "Method: OverDampedCurvature" << ", BoundaryCondition: " << boundaries[boundary] <<", Time: " << t_end << ", MismatchError: " << mismatch << endl;
			cout << "SmoothWeight: " << alpha << ", VortexWeight: " << vortexWeight << endl;
			break;
		case 10://OverdampedDiffusion
			cout << "Method: OverDampedDiffusion" << ", BoundaryCondition: " << boundaries[boundary] <<", Time: " << t_end << ", MismatchError: " << mismatch << endl;
			cout << "SmoothWeight: " << alpha << ", VortexWeight: " << vortexWeight << endl;
			break;
	}

	*templateImage /= templateImage->max();
	*sampleImage /= sampleImage->max();

	if (referenceImage)
	{
		reg = new FluidCurvatureRegistration(method, templateImage->height(), templateImage->width(),
				templateImage->data(), (long int) templateImage->size(), sampleImage->data(),
				(long int) sampleImage->size(), t_end, dt_start, alpha, viscosity, localDamping, mu, lambda,
				vortexWeight, localError, mismatch, verbose, NULL, referenceImage->spectrum(), referenceImage->data(),
				(long int) referenceImage->size());
		reg->registerImages();
	}
	else
	{
		reg = new FluidCurvatureRegistration(method, templateImage->height(), templateImage->width(),
				templateImage->data(), (long int) templateImage->size(), sampleImage->data(),
				(long int) sampleImage->size(), t_end, dt_start, alpha, viscosity, localDamping, mu, lambda,
				vortexWeight, localError, mismatch, verbose, NULL, 1, NULL, 0);
		reg->registerImages();
	}
	//check for reference image
	if (referenceImage != NULL)
	{
		patternImage = new CImg<double>(reg->getReference()->bm, reg->getReference()->nx, reg->getReference()->ny);
	}
	else
	{
		patternImage = new CImg<double>(reg->getSample()->bm, reg->getSample()->nx, reg->getSample()->ny);
	}
	if (patternFile != NULL)
	{
		try
		{
			patternImage->save(patternFile);
		} catch (CImgException &e)
		{
			fprintf(stderr, "Failed to save reference image at '%s'.\n", patternFile);
			return 1;
		}
	}
	else
	{
		if(verbose)
			patternImage->display("result");
	}

	if (flowFile != NULL)
	{
		param.end = t_end;
		param.error = mismatch;
		param.alpha = alpha;
		param.vortex_weight = vortexWeight;
		param.mu  = mu;
		param.lambda = lambda;
		param.boundary = boundaries[boundary];
		switch(method-boundary)
		{
			case 6://NavierLame
				param.method = methods[0];
				break;
			case 7://OverdampedCurvature
				param.method = methods[1];
				break;
			case 10://OverdampedDiffusion
				param.method = methods[2];
				break;
		}
		param.actual_error = reg->getMismatchError();
		param.actual_time = reg->getMinimalTime();

		if (!reg->getFlowField()->save(flowFile, &param))
		{
			fprintf(stderr, "Failed to save flow field at '%s'.\n", flowFile);
			return 1;
		}
	}

	delete reg;
	delete patternImage;
	delete templateImage;
	delete sampleImage;
	delete referenceImage;
}
