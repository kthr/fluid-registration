/*
 * HysteresisBinarize.hpp
 *
 *  Created on: Jul 11, 2012
 *      Author: kthierbach
 */

#ifndef HYSTERESISBINARIZE_HPP_
#define HYSTERESISBINARIZE_HPP_

#include "../templates/ImageTemplate.hpp"

template <class T>
Image<int>* hysteresisBinarize(T thres1, T thres2);

#endif /* HYSTERESISBINARIZE_HPP_ */
