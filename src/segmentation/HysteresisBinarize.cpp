/*
 * HisteresisBinarize.cpp
 *
 *  Created on: Jul 11, 2012
 *      Author: kthierbach
 */

#include "HysteresisBinarize.hpp"

template <class T>
Image<int>* hysteresisBinarize(T thres1, T thres2)
{
	const int highValue = 255;
	const int medValue = 127;
	bool changed = true;
	Image<int> *res;
	res = new Image<int>(nx, ny, channelNo);
	if (thres2 > thres1)
	{
		T tmp = thres1;
		thres1 = thres2;
		thres2 = tmp;
	}
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
			for (int ch = 0; ch < channelNo; ch++)
			{
				T v = get(i, j, ch);
				if (v >= thres1)
					res->setc(i, j, ch, highValue);
				else if (v > thres2)
					res->setc(i, j, ch, medValue);
				else
					res->setc(i, j, ch, 0);
			}
	while (changed)
	{
		changed = false;

		for (int i = 1; i < nx - 1; i++)
			for (int j = 1; j < ny - 1; j++)
				for (int ch = 0; ch < channelNo; ch++)
				{
					int val = res->get(i, j, ch);
					if (highValue != val)
						continue;
					if (medValue == res->get(i + 1, j, ch))
					{
						changed = true;
						res->setc(i + 1, j, ch, highValue);
					}
					if (medValue == res->get(i - 1, j, ch))
					{
						changed = true;
						res->setc(i - 1, j, ch, highValue);
					}
					if (medValue == res->get(i, j + 1, ch))
					{
						changed = true;
						res->setc(i, j + 1, ch, highValue);
					}
					if (medValue == res->get(i, j - 1, ch))
					{
						changed = true;
						res->setc(i, j - 1, ch, highValue);
					}
					if (medValue == res->get(i + 1, j + 1, ch))
					{
						changed = true;
						res->setc(i + 1, j + 1, ch, highValue);
					}
					if (medValue == res->get(i - 1, j + 1, ch))
					{
						changed = true;
						res->setc(i - 1, j + 1, ch, highValue);
					}
					if (medValue == res->get(i - 1, j + 1, ch))
					{
						changed = true;
						res->setc(i - 1, j + 1, ch, highValue);
					}
					if (medValue == res->get(i - 1, j - 1, ch))
					{
						changed = true;
						res->setc(i - 1, j - 1, ch, highValue);
					}

				}
		if (changed)
			for (int i = nx - 2; i > 0; i--)
				for (int j = ny - 2; j > 0; j--)
					for (int ch = 0; ch < channelNo; ch++)
					{
						int val = res->get(i, j, ch);
						if (highValue != val)
							continue;
						if (medValue == res->get(i + 1, j, ch))
						{
							changed = true;
							res->setc(i + 1, j, ch, highValue);
						}
						if (medValue == res->get(i - 1, j, ch))
						{
							changed = true;
							res->setc(i - 1, j, ch, highValue);
						}
						if (medValue == res->get(i, j + 1, ch))
						{
							changed = true;
							res->setc(i, j + 1, ch, highValue);
						}
						if (medValue == res->get(i, j - 1, ch))
						{
							changed = true;
							res->setc(i, j - 1, ch, highValue);
						}
						if (medValue == res->get(i + 1, j + 1, ch))
						{
							changed = true;
							res->setc(i + 1, j + 1, ch, highValue);
						}
						if (medValue == res->get(i - 1, j + 1, ch))
						{
							changed = true;
							res->setc(i - 1, j + 1, ch, highValue);
						}
						if (medValue == res->get(i - 1, j + 1, ch))
						{
							changed = true;
							res->setc(i - 1, j + 1, ch, highValue);
						}
						if (medValue == res->get(i - 1, j - 1, ch))
						{
							changed = true;
							res->setc(i - 1, j - 1, ch, highValue);
						}
					}
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				for (int ch = 0; ch < channelNo; ch++)
				{
					int val = res->get(i, j, ch);
					if (medValue == val)
						res->setc(i, j, ch, 0);
					else if (highValue == val)
						res->setc(i, j, ch, 1);
				}
		return res;
	}

}
