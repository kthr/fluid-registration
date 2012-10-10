/*
 * MedianFilter.hpp
 *
 *  Created on: Oct 10, 2012
 *      Author: kthierbach
 */

#ifndef MEDIANFILTER_HPP_
#define MEDIANFILTER_HPP_

#include <vector>
#include <algorithm>

#include "../templates/VectorArray2D.hpp"

VectorArray2D* medianFilter(VectorArray2D *vf, int neighborhood);
double median(std::vector<double> *values);

#endif /* MEDIANFILTER_HPP_ */
