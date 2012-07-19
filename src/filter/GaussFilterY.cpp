/*
 * GaussFilterY.cpp
 *
 *  Created on: Jul 11, 2012
 *      Author: Konstantin Thierbach
 */

#include "GaussFilterY.hpp"

GaussFilterY::GaussFilterY(Image<T> *_outimg, double _C, const double*_c) :
		outimg(_outimg), ny(_outimg->ny), channelNo(_outimg->channelNo), C(_C), c(_c)
{
}
void GaussFilterY::operator()(const tbb::blocked_range<int> &r) const
{
	double*w, *out;
	w = new double[ny];
	out = new double[ny];
	for (int i = r.begin(); i != r.end(); ++i)
	{
		for (int ch = 0; ch < channelNo; ch++)
		{

			w[0] = w[1] = w[2] = outimg->get(i, 0, ch);
			for (int j = 3; j < ny; j++)
				w[j] = C * outimg->get(i, j, ch) + (c[1] * w[j - 1] + c[2] * w[j - 2] + c[3] * w[j - 3]) / c[0];

			out[ny - 1] = out[ny - 2] = out[ny - 3] = w[ny - 1];
			for (int j = ny - 4; j > -1; j--)
				out[j] = C * w[j] + (c[1] * out[j + 1] + c[2] * out[j + 2] + c[3] * out[j + 3]) / c[0];

			for (int j = 0; j < ny; j++)
				outimg->setc(i, j, ch, (T) out[j]);
		}
	}
	delete[] out;
	delete[] w;
}
