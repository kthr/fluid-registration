/*
 * Point2D.hpp
 *
 *  Created on: Jul 11, 2012
 *      Author: kthierbach
 */

#ifndef POINT2D_HPP_
#define POINT2D_HPP_

#include <stdio.h>

template<class T>
class Point2D
{
	public:
		T x, y;

		Point2D() :
				x(0), y(0)
		{
		}
		Point2D(T i, T j) :
				x(i), y(j)
		{
		}
		Point2D operator+(const Point2D&p) const
		{
			Point2D q;
			q.x = x + p.x;
			q.y = y + p.y;
			return q;
		}
		Point2D& operator+=(const Point2D&v)
		{
			x += v.x;
			y += v.y;
			return *this;
		}
		Point2D operator-(const Point2D&p) const
		{
			Point2D q;
			q.x = x - p.x;
			q.y = y - p.y;
			return q;
		}
		Point2D operator-() const
		{
			Point2D p;
			p.x = -x;
			p.y = -y;
			return p;
		}
		Point2D operator-=(const Point2D&v)
		{
			x -= v.x;
			y -= v.y;
			return *this;
		}

		Point2D operator*(double f) const
		{
			Point2D q;
			q.x = f * x;
			q.y = f * y;
			return q;
		}
		Point2D& operator*=(double f)
		{
			x *= f;
			y *= f;
			return *this;
		}
		Point2D operator/(double f) const
		{
			double inv = 1.0 / f;
			Point2D q = *this;
			return q * inv;
		}
		Point2D operator/(int f) const
		{
			Point2D q = *this;
			q.x /= f;
			q.y /= f;
			return q;
		}
		Point2D& operator/=(double f)
		{
			double inv = 1.0 / f;
			x *= inv;
			y *= inv;
			return *this;
		}
		bool operator==(const Point2D&p) const
		{
			return x == p.x && y == p.y;
		}
		inline double length(void)
		{
			return sqrt((double) x * x + (double) y * y);
		}
		inline double length2(void)
		{
			return x * x + y * y;
		}
		inline Point2D perp(void) const
		{
			return Point2D<T>(y, -x);
		}
		inline Point2D unperp(void) const
		{
			return Point2D<T>(-y, x);
		}
		inline void normalize(void)
		{
			double l;
			l = length();
			if (l > 0)
			{
				l = 1.0 / l;
				x *= l;
				y *= l;
			}
		}
		inline T dot(Point2D p)
		{
			return p.x * x + p.y * y;
		}
		inline void increment(int dir)
		{
			switch (dir)
			{
				case North:
					y += 1;
					break;
				case South:
					y -= 1;
					break;
				case East:
					x += 1;
					break;
				case West:
					x -= 1;
					break;
				case SouthWest:
					x -= 1;
					y -= 1;
					break;
				case SouthEast:
					x += 1;
					y -= 1;
					break;
				case NorthWest:
					x -= 1;
					y += 1;
					break;
				case NorthEast:
					x += 1;
					y += 1;
					break;
				default:
					fprintf(stderr, "Panic .. direction is not unique %d\n", dir);
					break;
			}
		}
		inline double dPath(int dir)
		{
			const double Sqrt2 = 1.4142135623730950;
			switch (dir)
			{
				case North:
				case South:
				case East:
				case West:
					return 1.0;
				case SouthWest:
				case SouthEast:
				case NorthWest:
				case NorthEast:
					return Sqrt2;
				default:
					return 1.0;
			}
		}
	private:
		const static int SouthWest = 1, South = 2, SouthEast = 4, East = 8, NorthEast = 16, North = 32, NorthWest = 64, West = 128;
};

template<class T>
inline Point2D<T> operator*(const double f, const Point2D<T> &p)
{
	Point2D<T> q;
	q.x = f * p.x;
	q.y = f * p.y;
	return q;
}
#endif /* POINT2D_HPP_ */
