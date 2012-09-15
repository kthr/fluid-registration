/*
 * RKV43.hpp
 *
 *  Created on: Jul 11, 2012
 *      Author: Konstantin Thierbach
 */

#include "tbb/tbb.h"
//#include "../registration/FluidCurvatureRegistration.hpp"

#ifndef RKV43_HPP_
#define RKV43_HPP_

#if !defined(RK_OK)
#define RK_OK              (0)
#define RK_STEP_TOO_SMALL  (1)
#define RK_USER_STOP       (2)
#define RK_DO_NOT_ARRIVE   (3)
#define RK_NO_INTEROLATION (4)
#endif
#define savefac 0.9

static const double tiny = 2.220446049250313e-16;
static const double a1 = 3. / 8., a2 = 9. / 16., a3 = 25. / 32., a4 = 1.0;
static const double b21 = 3.0 / 8.0;
static const double b31 = 9.0 / 64., b32 = 27. / 64.;
static const double b41 = 925. / 5184., b42 = 625. / 3456., b43 = 4375. / 10368.;
static const double b51 = 63863. / 135675., b52 = -2516. / 1809, b53 = 11872. / 5427., b54 = -448.0 / 1675.;
static const double c1 = 9. / 50., c3 = 128. / 147., c4 = -1024. / 3675., c5 = 67. / 294.;
static const double d1 = -67. / 1350, d3 = 536. / 1323., d4 = -2144. / 3675., d5 = 67. / 294.;

class FluidCurvatureRegistration;
class RKV43
{
	public:
		RKV43(FluidCurvatureRegistration *fr);
		int RKVMethod43(double x0, double xe, double h, double y[], int n,
				void (FluidCurvatureRegistration::*diffy)(double, double*, double*, int),
				int (FluidCurvatureRegistration::*output)(double, double, double*, double*, int), double localerror,
				bool parallelQ = true);
		int InterpolateRKV4(double x, double h, double xp, double yp[]);
		void rkv1(double x, double h, double y[], int n,
				void (FluidCurvatureRegistration::*diffy)(double, double*, double*, int), double ynew[], double*err);
		void parallelRkv1(double x, double h, double y[], int n,
				void (FluidCurvatureRegistration::*diffy)(double, double*, double*, int), double ynew[], double*err,
				int grainSize = 8192);
		double *k1, *k2, *k3, *k4, *k5, *yy;
		int eqn_no;
		bool RKV_active;
	private:
		FluidCurvatureRegistration *fr;
};

#endif /* RKV43_HPP_ */
