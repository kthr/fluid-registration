/*
 * Utilities.cpp
 *
 *  Created on: Jul 20, 2012
 *      Author: kthierbach
 */

#include "Utilities.hpp"

void SplineCoefficients(double c[], double u, int i, int n)
{
	double u2;

	if (0 == i)
	{
		c[0] = 0.;
		c[1] = 1. + (-1.5 + 0.5 * u) * u;
		c[2] = (2. - 1. * u) * u;
		c[3] = (-0.5 + 0.5 * u) * u;
	}
	else if (n - 2 == i)
	{
		c[0] = 0.5 * (u - 1.) * u;
		c[1] = 1. - u * u;
		c[2] = 0.5 * (1. + u) * u;
		c[3] = 0.;
	}
	else
	{
		c[0] = u * (-0.5 + (1. - 0.5 * u) * u);
		c[1] = 1. + (u2 = u * u) * (-2.5 + 1.5 * u);
		c[2] = u * (0.5 + (2. - 1.5 * u) * u);
		c[3] = 0.5 * (u - 1.) * u2;
	}
}
