/*
 * MedianFilter.cpp
 *
 *  Created on: Oct 10, 2012
 *      Author: kthierbach
 */

#include "MedianFilter.hpp"

VectorArray2D* medianFilter(VectorArray2D *vf, int neighboorhood)
{
	VectorArray2D *result;
	std::vector<double> tmpX, tmpY;
	int nx, ny;
	int x, y;

	nx = vf->nx;
	ny = vf->ny;
	result = new VectorArray2D(nx,ny);

	for(int j=0; j<ny; ++j)
	{
		for(int i=0; i<nx; ++i)
		{
			for(int l=-neighboorhood; l<=neighboorhood; ++l)
			{
				for(int k=-neighboorhood; k<=neighboorhood; ++k)
				{
					x = i+k;
					y = j+l;
					if(!((x<0 || x>=nx) && (y<0 || y>=ny)))
						tmpX.push_back(vf->getx(i+k,j+l));
				}
			}
			for (int l = -neighboorhood; l <= neighboorhood; ++l)
			{
				for (int k = -neighboorhood; k <= neighboorhood; ++k)
				{
					x = i + k;
					y = j + l;
					if (!((x < 0 || x >= nx) && (y < 0 || y >= ny)))
						tmpY.push_back(vf->gety(i + k, j + l));
				}
			}
			result->set(i,j, median(&tmpX), median(&tmpY));
			tmpX.clear();
			tmpY.clear();
		}
	}
	return result;
}

double median(std::vector<double> *values)
{
	int size;

	size = values->size();
	std::sort(values->begin(), values->end());
	if(size % 2 == 0)
		return (values->at(size/2-1) + values->at(size/2))/2;
	else
		return values->at((size-1)/2);
}
