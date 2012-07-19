/*
 * OptimalBinarize.cpp
 *
 *  Created on: Jul 11, 2012
 *      Author: kthierbach
 */

#include "OptimalBinarize.hpp"

Image<int>* optimalBinarize(void)
{
	Image<int> *binimg;
	T*thres, *thresNext;
	int*nObj, *nBack;
	long int*muObj, *muBack;
	int i, ch;
	double err;

	thres = new T[channelNo];
	thresNext = new T[channelNo];
	nObj = new int[channelNo];
	nBack = new int[channelNo];
	muObj = new long[channelNo];
	muBack = new long[channelNo];
	for (i = 0; i < channelNo; i++)
		thres[i] = (min(i) + max(i)) / 2;
	while (1)
	{
		for (i = 0; i < channelNo; i++)
		{
			nObj[i] = nBack[i] = 0;
			muObj[i] = muBack[i] = 0L;
		}

		for (i = 0; i < pixelSize; i++)
		{
			ch = i % channelNo;
			if (bm[i] < thres[ch])
			{
				nBack[ch]++;
				muBack[ch] += bm[i];
			}
			else
			{
				nObj[ch]++;
				muObj[ch] += bm[i];
			}
		}

		for (err = 0.0, i = 0; i < channelNo; i++)
		{

			thresNext[i] = 0;
			if (muBack[i] > 0)
				thresNext[i] += muBack[i] / (2 * nBack[i]);
			if (nObj[i] > 0)
				thresNext[i] += muObj[i] / (2 * nObj[i]);

			err = fabs((double) (thres[i] - thresNext[i]));
		}

		if (err < 0.5)
			break;
		for (i = 0; i < channelNo; i++)
			thres[i] = thresNext[i];
	}
	binimg = new Image<int>(nx, ny, channelNo);
	for (i = 0; i < pixelSize; i++)
	{
		ch = i % channelNo;
		binimg->bm[i] = (bm[i] < thresNext[ch]) ? 0 : 1;
	}
	delete[] nObj;
	delete[] nBack;
	delete[] muObj;
	delete[] muBack;
	delete[] thres;
	delete[] thresNext;
	return binimg;
}

