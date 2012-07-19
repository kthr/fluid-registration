/*
 * RNKystroem.hpp
 *
 *  Created on: Jul 11, 2012
 *      Author: kthierbach
 */

#ifndef RKNYSTROEM_HPP_
#define RKNYSTROEM_HPP_

#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/tick_count.h"

#include "../registration/FluidCurvatureRegistration.hpp"

#define RKN_TINY               (2.220446049250313e-16)
#define RKN_OK             (0)
#define RKN_STEP_TOO_SMALL (1)
#define RKN_USER_STOP      (2)
#define RKN_DO_NOT_ARRIVE  (3)
#define RKN_OUTPUT_OK      (0)
#define RKN_OUTPUT_HALT    (1)

class FluidCurvatureRegistration;

class RKNystroemDSolve
{
	public:
		long int nFunctionCall, nStepCount, nAccept;
		RKNystroemDSolve(FluidCurvatureRegistration *fr, int dim, void (FluidCurvatureRegistration::*f)(double*, double, double*, double*, int),
				int (FluidCurvatureRegistration::*out)(double, double, double*, double*, double*, double*, int)= NULL);
		~RKNystroemDSolve(void);

		int integrate(double t0, double t1, double h, double*y, double*yp, double err, bool parallelQ);

	private:
		void (FluidCurvatureRegistration::*rhs)(double*, double, double*, double*, int);
		int (FluidCurvatureRegistration::*output)(double, double, double*, double*, double*, double*, int);
		int n;
		double *f1, *f2, *f3, *f4, *f5, *f6, *yarg, *yparg;
		double*yNew, *ypNew;
		FluidCurvatureRegistration *fr;

		double singleStep(double*yNew, double*ypNew, double t, double h, double*y, double*yp);
		double singleStepTBB(double*yNew, double*ypNew, double t, double h, double*y, double*yp, int grainSize = 2048);
};
#endif /* RKNYSTROEM_HPP_ */
