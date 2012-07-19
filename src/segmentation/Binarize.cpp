/*
 * Binarize.cpp
 *
 *  Created on: Jul 11, 2012
 *      Author: Konstantin Thierbach
 */

#include "Binarize.hpp"

Image<int>* binarize(double threshold) const
{
	Image<int> *binimg;
	int k;

	binimg = new Image<int>(nx, ny, channelNo);
	for (k = 0; k < pixelSize; k++)
		binimg->bm[k] = (bm[k] >= threshold ? 1 : 0);
	return binimg;
}
