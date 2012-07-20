/*
 * BracketMethod.hpp
 *
 *  Created on: Jul 20, 2012
 *      Author: kthierbach
 */

#ifndef BRACKETMETHOD_HPP_
#define BRACKETMETHOD_HPP_

#include "math.h"
#include "Utilities.hpp"
#include "ImageDifference.hpp"

#define TINY (0.1e-12)

class Bracketmethod
{
	public:
		double ax, bx, cx, fa, fb, fc;
		void bracket(const double a, const double b, ImageDifference *func);
		inline void shft2(double &a, double &b, const double c);
		inline void shft3(double &a, double &b, double &c, const double d);
		inline void mov3(double &a, double &b, double &c, const double d, const double e, const double f);
};

class Brent: public Bracketmethod
{
	public:
		double xmin, fmin;
		const double tol;
		Brent(const double toll = 3.0e-8);
		double minimize(ImageDifference *func);
};

#endif /* BRACKETMETHOD_HPP_ */
