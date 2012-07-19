/*
 * RetinexFilter.cpp
 *
 *  Created on: Jul 10, 2012
 *      Author: Konstantin Thierbach
 */

#include "RetinexFilter.hpp"

void retinexFilter(int mode, double sigma, int scaleNo)
{
	const double alpha = 125.;
	Image<T> *img1, *img2;
	Image<double> *out;
	double dscale, weight;
	double*scales;
	int i, ch, k;

	weight = 1.0 / scaleNo;
	out = new Image<double>(nx, ny, channelNo);

	scales = new double[scaleNo];
	if (1 == scaleNo)
		scales[0] = 0.5 * sigma;
	else if (2 == scaleNo)
	{
		scales[0] = 0.5 * sigma;
		scales[1] = sigma;
	}
	else
	{
		switch (mode)
		{
		case RETINEX_UNIFORM:
			dscale = sigma / scaleNo;
			for (i = 0; i < scaleNo; i++)
				scales[i] = 2.0 + i * dscale;
			break;
		case RETINEX_LOW:
			dscale = log(sigma - 2.0) / scaleNo;
			for (i = 0; i < scaleNo; i++)
				scales[i] = 2.0 + exp(i * dscale);
			break;
		case RETINEX_HIGH:
			dscale = log(sigma - 2.0) / scaleNo;
			for (i = 0; i < scaleNo; i++)
				scales[i] = 2.0 - exp(i * dscale);
			break;
		}
	}

	img1 = new Image<double>(nx, ny, channelNo);
	for (i = 0; i < pixelSize; i++)
		img1->bm[i] = bm[i] + 1.0;

	for (k = 0; k < scaleNo; k++)
	{
		img2 = ParallelGaussFilter<T>(img1, scales[k]);
		for (i = 0; i < pixelSize; i++)
		{
			out->bm[i] += weight * (log(img1->bm[i]) - log(img2->bm[i]));
		}
		delete img2;
	}

	delete img1;

	if (3 == channelNo)
	{
		T*imgp;
		double*outp, sum;

		for (i = 0; i < nx * ny; i++)
		{
			imgp = bm + i * channelNo;
			outp = out->bm + i * channelNo;
			sum = imgp[0] + imgp[1] + imgp[2] + 3.0;
			outp[0] *= log(1.0 + alpha * (imgp[0] + 1.0) / sum);
			outp[1] *= log(1.0 + alpha * (imgp[1] + 1.0) / sum);
			outp[2] *= log(1.0 + alpha * (imgp[2] + 1.0) / sum);
		}
	}

	{
		int maxin = max();
		for (ch = 0; ch < channelNo; ch++)
		{
			double mean, stddev;
			double scaleMin, scaleMax, range;
			T maxin = max(ch);
			out->meanAndStandarddeviation(&mean, &stddev, Max(nx / 2, ny / 2), nx / 2, ny / 2, ch);
			scaleMin = mean - 2.0 * stddev;
			scaleMax = mean + 2.0 * stddev;
			range = scaleMax - scaleMin;
			if (range < TINY)
				range = 1.0;
			range = maxin / range;
			for (i = 0; i < nx * ny; i++)
			{
				double*ptr = out->bm + channelNo * i + ch;
				double scaledValue = (*ptr - scaleMin) * range;
				scaledValue = Min(Max(scaledValue, 0.0), (double) maxin);
				*ptr = scaledValue;
			}
		}
	}

	for (i = 0; i < pixelSize; i++)
		bm[i] = (T) out->bm[i];
	delete out;
}
