/*
 * GaussFilterX.cpp
 *
 *  Created on: Jul 11, 2012
 *      Author: Konstantin Thierbach
 */

#include "GaussFilterX.hpp"

GaussFilterX::GaussFilterX(Image<T> *_outimg, const Image<T> *_img, double _C, const double*_c) :
		img(_img), outimg(_outimg), nx(_img->nx), channelNo(_img->channelNo), C(_C), c(_c)
{
}
void GaussFilterX::operator()(const tbb::blocked_range<int> &r) const
{
	double*w, *out;
	w = new double[nx];
	out = new double[nx];
	for (int j = r.begin(); j != r.end(); ++j)
	{
		for (int ch = 0; ch < channelNo; ch++)
		{

			w[0] = w[1] = w[2] = img->get(0, j, ch);
			for (int i = 3; i < nx; i++)
				w[i] = C * img->get(i, j, ch) + (c[1] * w[i - 1] + c[2] * w[i - 2] + c[3] * w[i - 3]) / c[0];

			out[nx - 1] = out[nx - 2] = out[nx - 3] = w[nx - 1];
			for (int i = nx - 4; i > -1; i--)
				out[i] = C * w[i] + (c[1] * out[i + 1] + c[2] * out[i + 2] + c[3] * out[i + 3]) / c[0];

			for (int i = 0; i < nx; i++)
			{
				outimg->setc(i, j, ch, (T) out[i]);
			}
		}
	}
	delete[] w;
	delete[] out;
}

