/*
 * ParallelGaussFilter.hpp
 *
 *  Created on: Jul 10, 2012
 *      Author: Konstantin Thierbach
 */

#ifndef PARALLELGAUSSFILTER_HPP_
#define PARALLELGAUSSFILTER_HPP_

#include "../templates/ImageTemplate.hpp"
#include "GaussFilterX.hpp"
#include "GaussFilterY.hpp"

template <class T>
Image<T> *parallelGaussFilter(const Image<T> *img, double*sigmaArray, int sdim, int grainSize = 64);
template <class T>
Image<T> *parallelGaussFilter(const Image<T> *img, double sigma, int grainSize = 64);

#endif /* PARALLELGAUSSFILTER_HPP_ */
