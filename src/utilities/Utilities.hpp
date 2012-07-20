/*
 * utilities.hpp
 *
 *  Created on: Jul 11, 2012
 *      Author: kthierbach
 */

#ifndef UTILITIES_HPP_
#define UTILITIES_HPP_

#include <stdlib.h>
#include "../templates/Point2D.hpp"

template<class T>
inline T Abs(T x)
{
	return x > 0 ? x : -x;
}

template<class T>
double computeArea(Point2D<T> *p, int n)
{
	double area = 0.0;

	for (int i = 0; i < n - 2; i++)
		area += 0.5
				* (p[1 + i].x * p[2 + i].y + p[i].x * p[1 + i].y + p[i].y * p[2 + i].x
						- (p[1 + i].y * p[2 + i].x + p[i].y * p[1 + i].x + p[i].x * p[2 + i].y));
	return area;
}
template<class T>
int Greater(const void*e1, const void*e2)
{
	T*i, *j;

	i = (T*) e1;
	j = (T*) e2;
	if (*i < *j)
		return -1;
	if (*i > *j)
		return 1;
	return 0;
}
template<class T>
T Min(T x, T y)
{
	return x < y ? x : y;
}
template<class T>
T Min(T*x, int n)
{
	T minX;

	minX = x[0];
	for (int i = 1; i < n; i++)
		if (minX > x[i])
			minX = x[i];
	return minX;
}
template<class T>
T Max(T x, T y)
{
	return x > y ? x : y;
}
template<class T>
T Max(T*x, int n)
{
	T maxX;

	maxX = x[0];
	for (int i = 1; i < n; i++)
		if (maxX < x[i])
			maxX = x[i];
	return maxX;
}
template<class T>
T minMod(T x, T y)
{
	if (x * y > 0)
		return (Sign(x) + Sign(y)) * Min(Abs(x), Abs(y)) / 2;
	return (T) 0;
}
inline int modulo(int i, int n)
{
	int k;
	k = i % n;
	while (k < 0)
		k += n;
	return k;
}
template<class T>
inline T Sign(T x)
{
	const T zero = 0;

	if (x < zero)
		return -1;
	if (x > zero)
		return 1;
	return zero;
}
template<class T>
inline T Sign(const T&a, const T&b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
void SplineCoefficients(double c[], double u, int i, int n);
template<class T>
void Swap(T*a, T*b)
{
	T tmp = *a;
	*a = *b;
	*b = tmp;
}
template<class T>
void sort(T*x, int n)
{
	qsort(x, n, sizeof(T), Greater<T>);
}
#endif /* UTILITIES_HPP_ */
