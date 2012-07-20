/*
 * BracketMethod.cpp
 *
 *  Created on: Jul 20, 2012
 *      Author: kthierbach
 */

#include "BracketMethod.hpp"

void Bracketmethod::bracket(const double a, const double b, ImageDifference *func)
{
	const double GOLD = 1.618034;
	const double GLIMIT = 1.0; // don't grow thi interval

	ax = a;
	bx = b;
	double fu;
	fa = (*func)(ax);
	fb = (*func)(bx);
	if (fb > fa)
	{
		Swap(&ax, &bx);
		Swap(&fb, &fa);
	}
	cx = bx + GOLD * (bx - ax);
	fc = (*func)(cx);
	while (fb > fc)
	{
		double r = (bx - ax) * (fb - fc);
		double q = (bx - cx) * (fb - fa);
		double u = bx - ((bx - cx) * q - (bx - ax) * r) / (2.0 * Sign(Max(fabs(q - r), TINY), q - r));
		double ulim = bx + GLIMIT * (cx - bx);
		if ((bx - u) * (u - cx) > 0.0)
		{
			fu = (*func)(u);
			if (fu < fc)
			{
				ax = bx;
				bx = u;
				fa = fb;
				fb = fu;
				return;
			}
			else if (fu > fb)
			{
				cx = u;
				fc = fu;
				return;
			}
			u = cx + GOLD * (cx - bx);
			fu = (*func)(u);
		}
		else if ((cx - u) * (u - ulim) > 0.0)
		{
			fu = (*func)(u);
			if (fu < fc)
			{
				shft3(bx, cx, u, u + GOLD * (u - cx));
				shft3(fb, fc, fu, (*func)(u));
			}
		}
		else if ((u - ulim) * (ulim - cx) >= 0.0)
		{
			u = ulim;
			fu = (*func)(u);
		}
		else
		{
			u = cx + GOLD * (cx - bx);
			fu = (*func)(u);
		}
		shft3(ax, bx, cx, u);
		shft3(fa, fb, fc, fu);
	}
}
inline void Bracketmethod::shft2(double &a, double &b, const double c)
{
	a = b;
	b = c;
}
inline void Bracketmethod::shft3(double &a, double &b, double &c, const double d)
{
	a = b;
	b = c;
	c = d;
}
inline void Bracketmethod::mov3(double &a, double &b, double &c, const double d, const double e, const double f)
{
	a = d;
	b = e;
	c = f;
}

Brent::Brent(const double toll /*= 3.0e-8*/) :
		tol(toll)
{
}
double Brent::minimize(ImageDifference *func)
{
	const int ITMax = 100;
	const double CGOLD = 0.3819660;
	const double ZEPS = 2.220446049250313e-16 * 1.0e-3;
	double a, b, d = 0.0, etemp, fu, fv, fw, fx;
	double p, q, r, tol1, tol2, u, v, w, x, xm;
	double e = 0.0;

	a = (ax < cx ? ax : cx);
	b = (ax > cx ? ax : cx);
	x = w = v = bx;
	fw = fv = fx = (*func)(x);
	for (int iter = 0; iter < ITMax; iter++)
	{
		xm = 0.5 * (a + b);
		tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);
		if (fabs(x - xm) <= (tol2 - 0.5 * (b - a)))
		{
			fmin = fx;
			return xmin = x;
		}
		if (fabs(e) > tol1)
		{
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = 2.0 * (q - r);
			if (q > 0.0)
				p = -p;
			q = fabs(q);
			etemp = e;
			e = d;
			if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			else
			{
				d = p / q;
				u = x + d;
				if (u - a < tol2 || b - u < tol2)
					d = Sign(tol1, xm - x);
			}
		}
		else
		{
			d = CGOLD * (e = (x >= xm ? a - x : b - x));
		}
		u = (fabs(d) >= tol1 ? x + d : x + Sign(tol1, d));
		fu = (*func)(u);
		if (fu <= fx)
		{
			if (u >= x)
				a = x;
			else
				b = x;
			shft3(v, w, x, u);
			shft3(fv, fw, fx, fu);
		}
		else
		{
			if (u < x)
				a = u;
			else
				b = u;
			if (fu <= fw || w == x)
			{
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if (fu <= fv || v == x || v == w)
			{
				v = u;
				fv = fu;
			}
		}
	}
	throw("Too many iterations in brent");
}
