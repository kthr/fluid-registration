/*
 * RKV43.cpp
 *
 *  Created on: Jul 11, 2012
 *      Author: kthierbach
 */
#include <math.h>
#include "RKV43.hpp"

RKV43::RKV43(FluidCurvatureRegistration *fr) : fr(fr)
{
}

void RKV43::rkv1(double x, double h, double y[], int n, void (FluidCurvatureRegistration::*diffy)(double, double*, double*, int), double ynew[],
		double*err)
{
	double ycorr;
	int i;

	((fr)->*(diffy))(x, y, k1, n);
	for (i = 0; i < n; i++)
		ynew[i] = y[i] + h * b21 * k1[i];

	((fr)->*(diffy))(x + a1 * h, ynew, k2, n);
	for (i = 0; i < n; i++)
		ynew[i] = y[i] + h * (b31 * k1[i] + b32 * k2[i]);
	((fr)->*(diffy))(x + a2 * h, ynew, k3, n);
	for (i = 0; i < n; i++)
		ynew[i] = y[i] + h * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
	((fr)->*(diffy))(x + a3 * h, ynew, k4, n);
	for (i = 0; i < n; i++)
		ynew[i] = y[i] + h * (b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i]);
	((fr)->*(diffy))(x + h, ynew, k5, n);
	*err = 0.0;
	for (i = 0; i < n; i++)
	{
		ynew[i] = y[i] + h * (c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c5 * k5[i]);
		ycorr = h * (d1 * k1[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i]);
		if (fabs(ycorr) > *err)
			*err = fabs(ycorr);
	}
}
class RKVComputeArg2
{
	public:
		RKVComputeArg2(double*_ynew, double*_y, double*_k1, double _h) :
				ynew(_ynew), y(_y), k1(_k1), h(_h)
		{
		}
		void operator()(const tbb::blocked_range<int> &r) const
		{
			for (int i = r.begin(); i != r.end(); ++i)
				ynew[i] = y[i] + h * b21 * k1[i];
		}
	private:
		double* const k1;
		double* const y;
		double* const ynew;
		double h;
};
class RKVComputeArg3
{
	public:
		RKVComputeArg3(double*_ynew, double*_y, double*_k1, double*_k2, double _h) :
				ynew(_ynew), y(_y), k1(_k1), k2(_k2), h(_h)
		{
		}
		void operator()(const tbb::blocked_range<int> &r) const
		{
			for (int i = r.begin(); i != r.end(); ++i)
				ynew[i] = y[i] + h * (b31 * k1[i] + b32 * k2[i]);
		}
	private:
		double* const k1;
		double* const k2;
		double* const y;
		double* const ynew;
		double h;
};
class RKVComputeArg4
{
	public:
		RKVComputeArg4(double*_ynew, double*_y, double*_k1, double*_k2, double*_k3, double _h) :
				ynew(_ynew), y(_y), k1(_k1), k2(_k2), k3(_k3), h(_h)
		{
		}
		void operator()(const tbb::blocked_range<int> &r) const
		{
			for (int i = r.begin(); i != r.end(); ++i)
				ynew[i] = y[i] + h * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
		}
	private:
		double* const k1;
		double* const k2;
		double* const k3;
		double* const y;
		double* const ynew;
		double h;
};
class RKVComputeArg5
{
	public:
		RKVComputeArg5(double*_ynew, double*_y, double*_k1, double*_k2, double*_k3, double*_k4, double _h) :
				ynew(_ynew), y(_y), k1(_k1), k2(_k2), k3(_k3), k4(_k4), h(_h)
		{
		}
		void operator()(const tbb::blocked_range<int> &r) const
		{
			for (int i = r.begin(); i != r.end(); ++i)
				ynew[i] = y[i] + h * (b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i]);
		}
	private:
		double* const k1;
		double* const k2;
		double* const k3;
		double* const k4;
		double* const y;
		double* const ynew;
		double h;
};
class RKVComputeNext
{
	public:
		double erry;
		RKVComputeNext(double*_ynew, double*_y, double*_k1, double*_k3, double*_k4, double*_k5, double _h) :
				ynew(_ynew), y(_y), k1(_k1), k3(_k3), k4(_k4), k5(_k5), h(_h), erry(0.0)
		{
		}
		RKVComputeNext(RKVComputeNext&b, tbb::split) :
				ynew(b.ynew), y(b.y), k1(b.k1), k3(b.k3), k4(b.k4), k5(b.k5), h(b.h), erry(0.0)
		{
		}

		void operator()(const tbb::blocked_range<int> &r)
		{
			for (int i = r.begin(); i != r.end(); ++i)
			{
				double ycorr;
				ynew[i] = y[i] + h * (c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c5 * k5[i]);
				ycorr = h * (d1 * k1[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i]);
				if (fabs(ycorr) > erry)
					erry = fabs(ycorr);
			}
		}
		void join(const RKVComputeNext&b)
		{
			if (b.erry > erry)
				erry = b.erry;
		}
	private:
		double* const k1;
		double* const k3;
		double* const k4;
		double* const k5;

		double* const y;
		double*ynew;
		double h;
};

void RKV43::parallelRkv1(double x, double h, double y[], int n, void (FluidCurvatureRegistration::*diffy)(double, double*, double*, int),
		double ynew[], double*err, int grainSize /*= 8192*/)
{
	((fr)->*(diffy))(x, y, k1, n);
	tbb::parallel_for(tbb::blocked_range<int>(0, n, grainSize), RKVComputeArg2(ynew, y, k1, h));

	((fr)->*(diffy))(x + a1 * h, ynew, k2, n);
	tbb::parallel_for(tbb::blocked_range<int>(0, n, grainSize), RKVComputeArg3(ynew, y, k1, k2, h));

	((fr)->*(diffy))(x + a2 * h, ynew, k3, n);
	tbb::parallel_for(tbb::blocked_range<int>(0, n, grainSize), RKVComputeArg4(ynew, y, k1, k2, k3, h));

	((fr)->*(diffy))(x + a3 * h, ynew, k4, n);
	tbb::parallel_for(tbb::blocked_range<int>(0, n, grainSize), RKVComputeArg5(ynew, y, k1, k2, k3, k4, h));

	((fr)->*(diffy))(x + h, ynew, k5, n);
	RKVComputeNext next(ynew, y, k1, k3, k4, k5, h);
	tbb::parallel_reduce(tbb::blocked_range<int>(0, n, grainSize), next);

	*err = next.erry;

}
#undef USER_STOP_MSG
int RKV43::RKVMethod43(double x0, double xe, double h, double y[], int n, void (FluidCurvatureRegistration::*diffy)(double, double*, double*, int),
		int (FluidCurvatureRegistration::*output)(double, double, double*, double*, int), double localerror, bool parallelQ)

{
	int i;
	double x;
	double hnew, hfac;
	double ddx, error_est;
	double *ynew;
	bool succ;
	double int_length;

	if (((xe - x0) * h) < 0.0)
		return RK_DO_NOT_ARRIVE;
	RKV_active = true;
	eqn_no = n;
	localerror = fabs(localerror);
	ddx = fabs(xe - x0);
	int_length = (double) 1.0 / (xe - x0);
	k1 = new double[n];
	k2 = k5 = new double[n];
	k3 = new double[n];
	k4 = new double[n];
	yy = new double[n];
	ynew = new double[n];
	hnew = h;
	x = x0;
	if (output)
		((fr)->*(output))(x, 0.0, y, y, n);
	do
	{
		for (;;)
		{
			h = hnew;
			if (fabs(fabs(x) + 0.05 * fabs(h)) <= fabs(x) || fabs(h) <= tiny)
			{
				{
					delete[] k1;
					delete[] k2;
					delete[] k3;
					delete[] k4;
					delete[] yy;
					delete[] ynew;
				}
				RKV_active = false;
				return RK_STEP_TOO_SMALL;
			}
			if (fabs(x + h - x0) > ddx)
			{
				h = xe - x;
			}
			if (parallelQ)
				parallelRkv1(x, h, y, n, diffy, ynew, &error_est);
			else
				rkv1(x, h, y, n, diffy, ynew, &error_est);
			if (fabs(error_est) <= tiny)
				error_est = tiny;
			if (error_est > localerror)
			{
				hfac = savefac * pow(fabs(localerror / error_est), 0.333333333333);
				succ = false;
			}
			else
			{
				hfac = savefac * pow(fabs(localerror / error_est), 0.25);
				succ = true;
			}
			if (hfac < 0.1)
				hnew = 0.1 * h;
			else
				hnew = (hfac > 5.0 ? 5.0 : hfac) * h;
			if (succ)
				break;
		}
		x += h;
		for (i = 0; i < n; i++)
		{
			yy[i] = y[i];
			y[i] = ynew[i];
		}
		if (output)
			if (((fr)->*(output))(x, h, yy, y, n))
			{

				{
					delete[] k1;
					delete[] k2;
					delete[] k3;
					delete[] k4;
					delete[] yy;
					delete[] ynew;
				}

				RKV_active = false;
				return RK_USER_STOP;
			}
	} while ((fabs((x - x0) * int_length) + tiny < 1.0));
	RKV_active = false;
	{
		delete[] k1;
		delete[] k2;
		delete[] k3;
		delete[] k4;
		delete[] yy;
		delete[] ynew;
	}
	return RK_OK;
}
int RKV43::InterpolateRKV4(double x, double h, double xp, double yp[])
{
	int i;
	double u, u2;
	double bi1, bi3, bi4, bi5;

	if (!RKV_active)
	{
		return RK_NO_INTEROLATION;
	}
	x -= h;
	if (h != 0.0)
		u = (xp - x) / h;
	else
		u = -1.0;
	if (fabs(u) < tiny)
		u = 0.0;
	if (fabs(1.0 - u) < tiny)
		u = 1.0;
	if ((u > 1.0) || (u < 0.0))
	{
		return RK_NO_INTEROLATION;
	}
	u2 = u * u;
	bi1 = u * (1.0 + u * (-2.02888888888888889 + (1.77777777777777778 - 0.56888888888888889 * u) * u));
	bi3 = u2 * (7.2562358276643991 + u * (-11.0294784580498866 + 4.6439909297052154 * u));
	bi4 = u2 * (-7.5232653061224490 + (13.9319727891156463 - 6.6873469387755102 * u) * u);
	bi5 = u2 * (2.29591836734693878 + u * (-4.6802721088435374 + 2.61224489795918367 * u));
	for (i = 0; i < eqn_no; i++)
		yp[i] = yy[i] + h * (bi1 * k1[i] + bi3 * k3[i] + bi4 * k4[i] + bi5 * k5[i]);
	return RK_OK;
}
