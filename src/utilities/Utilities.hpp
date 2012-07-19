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
inline T Abs(T x);
template<class T>
double computeArea(Point2D<T> *p, int n);
template<class T>
int Greater(const void*e1, const void*e2);
template<class T>
T Min(T x, T y);
template<class T>
T Min(T*x, int n);
template<class T>
T Max(T x, T y);
template<class T>
T Max(T*x, int n);
template<class T>
T minMod(T x, T y);
inline int modulo(int i, int n);
template<class T>
inline T Sign(T x);
void SplineCoefficients(double c[], double u, int i, int n);
template<class T>
void sort(T*x, int n);

#endif /* UTILITIES_HPP_ */
