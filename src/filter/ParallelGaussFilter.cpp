/*
 * ParallelGaussFilter.cpp
 *
 *  Created on: Jul 10, 2012
 *      Author: Konstantin Thierbach
 */

#include "ParallelGaussFilter.hpp"

template<class T>
Image<T> *parallelGaussFilter(const Image<T> *img, double*sigmaArray, int sdim, int grainSize /*= 64*/)
{
	double q, q2, q3;
	double sigma, c[4], C;
	Image<T> *outimg;

	outimg = new Image<T>(img->nx, img->ny, img->channelNo);
	sigma = sigmaArray[0];
	if (sigma < 0.5)
		q = 0.1147705018520355224609375;
	if (0.5 <= sigma && sigma <= 2.5)
		q = 3.97156 - 4.14554 * sqrt(1.0 - 0.26891 * sigma);
	else if (2.5 < sigma)
		q = 0.98711 * sigma - 0.9633;
	q2 = q * q;
	q3 = q * q2;
	c[0] = (1.57825 + (2.44413 * q) + (1.4281 * q2) + (0.422205 * q3));
	c[1] = ((2.44413 * q) + (2.85619 * q2) + (1.26661 * q3));
	c[2] = (-((1.4281 * q2) + (1.26661 * q3)));
	c[3] = ((0.422205 * q3));
	C = 1.0 - ((c[1] + c[2] + c[3]) / c[0]);
	tbb::parallel_for(tbb::blocked_range<int>(0, img->ny, grainSize), GaussFilterX < T > (outimg, img, C, c));
	sigma = sigmaArray[1 % sdim];
	if (sigma < 0.5)
		q = 0.1147705018520355224609375;
	if (0.5 <= sigma && sigma <= 2.5)
		q = 3.97156 - 4.14554 * sqrt(1.0 - 0.26891 * sigma);
	else if (2.5 < sigma)
		q = 0.98711 * sigma - 0.9633;
	q2 = q * q;
	q3 = q * q2;
	c[0] = (1.57825 + (2.44413 * q) + (1.4281 * q2) + (0.422205 * q3));
	c[1] = ((2.44413 * q) + (2.85619 * q2) + (1.26661 * q3));
	c[2] = (-((1.4281 * q2) + (1.26661 * q3)));
	c[3] = ((0.422205 * q3));
	C = 1.0 - ((c[1] + c[2] + c[3]) / c[0]);
	tbb::parallel_for(tbb::blocked_range<int>(0, img->nx, grainSize), GaussFilterY < T > (outimg, C, c));
	return outimg;
}
template<class T>
Image<T> *parallelGaussFilter(const Image<T> *img, double sigma, int grainSize /*= 64*/)
{
	double q, q2, q3;
	double c[4], C;
	Image<T> *outimg;
	outimg = new Image<T>(img->nx, img->ny, img->channelNo);

	if (sigma < 0.5)
		q = 0.1147705018520355224609375;
	if (0.5 <= sigma && sigma <= 2.5)
		q = 3.97156 - 4.14554 * sqrt(1.0 - 0.26891 * sigma);
	else if (2.5 < sigma)
		q = 0.98711 * sigma - 0.9633;
	q2 = q * q;
	q3 = q * q2;
	c[0] = (1.57825 + (2.44413 * q) + (1.4281 * q2) + (0.422205 * q3));
	c[1] = ((2.44413 * q) + (2.85619 * q2) + (1.26661 * q3));
	c[2] = (-((1.4281 * q2) + (1.26661 * q3)));
	c[3] = ((0.422205 * q3));
	C = 1.0 - ((c[1] + c[2] + c[3]) / c[0]);
	tbb::parallel_for(tbb::blocked_range<int>(0, img->ny, grainSize), GaussFilterX < T > (outimg, img, C, c));
	tbb::parallel_for(tbb::blocked_range<int>(0, img->nx, grainSize), GaussFilterY < T > (outimg, C, c));
	return outimg;
}

