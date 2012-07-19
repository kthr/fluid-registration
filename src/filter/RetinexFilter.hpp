/*
 * RetinexFilter.hpp
 *
 *  Created on: Jul 10, 2012
 *      Author: Konstantin Thierbach
 */

#ifndef RETINEXFILTER_HPP_
#define RETINEXFILTER_HPP_

#include "../templates/ImageTemplate.hpp"
#include "ParallelGaussFilter.hpp"

void retinexFilter(int mode, double sigma, int scaleNo);

#endif /* RETINEXFILTER_HPP_ */
