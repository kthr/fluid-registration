/*
 * ImageDifference.hpp
 *
 *  Created on: Jul 19, 2012
 *      Author: kthierbach
 */

#ifndef IMAGEDIFFERENCE_HPP_
#define IMAGEDIFFERENCE_HPP_

#include "../templates/ImageTemplate.hpp"
#include "../solver/RKV43.hpp"

#define MAXDOUBLE (1.7976931348623157e308)

class ImageDifference
{
	public:
		VectorArray2D*u2;
		double t, h;
		int nx, ny;
		ImageDifference(RKV43 *rkv43, Image<double> *templateImage, Image<double> *sampleImage, Image<double> *wraped, int _nx, int _ny, double _t, double _h);
		~ImageDifference();
		double operator()(double tp);
	private:
		RKV43 *rkv43;
		Image<double> *templateImage, *sampleImage, *wraped;

};
#endif /* IMAGEDIFFERENCE_HPP_ */
