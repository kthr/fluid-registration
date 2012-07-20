/*
 * ImageDifference.cpp
 *
 *  Created on: Jul 19, 2012
 *      Author: kthierbach
 */

#include "ImageDifference.hpp"

ImageDifference::ImageDifference(RKV43 *rkv43, Image<double> *templateImage, Image<double> *sampleImage,
		Image<double> *wraped, int _nx, int _ny, double _t, double _h) :
		rkv43(rkv43), templateImage(templateImage), sampleImage(sampleImage), wraped(wraped), nx(_nx), ny(_ny), t(_t), h(
				_h)
{
	u2 = new VectorArray2D(nx, ny);
}
ImageDifference::~ImageDifference()
{
	delete u2;
}
double ImageDifference::operator()(double tp)
{
	double imgDiff;
	int size = nx * ny;
	int res;
	if (tp < t - h || tp > tp)
		return MAXDOUBLE;
	res = rkv43->InterpolateRKV4(t, h, tp, u2->vx);
	if (RK_OK == res)
	{
		sampleImage->wrap(wraped, u2);
		imgDiff = 0.0;

#pragma omp parallel for reduction(+:imgDiff)
		for (int j = 0; j < ny; j++)
		{
			double rsum = 0.0;
			for (int i = 0; i < nx; i++)
			{
				double tmp = templateImage->get(i, j, 0) - wraped->get(i, j, 0);
				rsum += tmp * tmp;
			}
			imgDiff += rsum;
		}

		imgDiff = sqrt(imgDiff);
		imgDiff /= size;

	}
	else
		imgDiff = MAXDOUBLE;
	return imgDiff;
}

