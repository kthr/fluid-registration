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
	CImg<double> *tImage = new CImg<double>("/Users/kthierbach/Documents/data/brox/timg.png");
	CImg<double> *sImage = new CImg<double>("/Users/kthierbach/Documents/data/brox/simg.png");

	printf("tImage channels: %d dx: %d dy: %d\n", tImage->spectrum(), tImage->width(), tImage->height());
	printf("sImage channels: %d dx: %d dy: %d\n", sImage->spectrum(), sImage->width(), sImage->height());
	FluidCurvatureRegistration *reg = new FluidCurvatureRegistration(7, tImage->width(), tImage->height(),
			tImage->data(), (long int) tImage->size(), sImage->data(), (long int) sImage->size(), 64., 0.001, 1, 1.0,
			1.0, 1.0, 0.25, 0., .5, 0.00005, 1, NULL, 1, NULL, 0);
	reg->registerImages();
	delete tImage;
	delete sImage;
}
