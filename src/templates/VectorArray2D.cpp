/*
 * VectorArray2D.cpp
 *
 *  Created on: Jul 10, 2012
 *      Author: Konstantin Thierbach
 */

#include "VectorArray2D.hpp"

VectorArray2D::VectorArray2D(void) :
		nx(0), ny(0), dx(1.0), dy(1.0), vx(NULL), vy(NULL)
{
}

VectorArray2D::VectorArray2D(int _nx, int _ny, double _dx /*= 1.0*/, double _dy /*= 1.0*/) :
		nx(_nx), ny(_ny), dx(_dx), dy(_dy)
{
	size_t size = nx * ny;

	vx = new double[2 * size];
	vy = vx + size;
	for (unsigned i = 0; i < size; i++)
		vx[i] = vy[i] = 0.0;
}
VectorArray2D::VectorArray2D(VectorArray2D &orig) :
		nx(orig.nx), ny(orig.ny), dx(orig.dx), dy(orig.dy)
{
	size_t size = nx * ny;
	vx = new double[2 * size];
	vy = vx + size;
	for (unsigned int i = 0; i < size; i++)
	{
		vx[i] = orig.vx[i];
		vy[i] = orig.vy[i];
	}
}
VectorArray2D::~VectorArray2D()
{
	if (vx)
		delete[] vx;
}
Vector2D VectorArray2D::get(int i, int j) const
{
	unsigned long addr;
	Vector2D v;

	addr = Address(i, j, nx, ny);
	v.x = vx[addr];
	v.y = vy[addr];
	return v;
}

double VectorArray2D::getx(int i, int j) const
{
	return vx[Address(i, j, nx, ny)];
}

double VectorArray2D::gety(int i, int j) const
{
	return vy[Address(i, j, nx, ny)];
}

void VectorArray2D::set(int i, int j, Vector2D v)
{
	unsigned long addr;

	addr = Address(i, j, nx, ny);
	vx[addr] = v.x;
	vy[addr] = v.y;
}

void VectorArray2D::set(int i, int j, Vector2D*v)
{
	unsigned long addr;

	addr = Address(i, j, nx, ny);
	vx[addr] = v->x;
	vy[addr] = v->y;
}

void VectorArray2D::set(int i, int j, double x, double y)
{
	unsigned long addr;

	addr = Address(i, j, nx, ny);
	vx[addr] = x;
	vy[addr] = y;
}

void VectorArray2D::setAll(Vector2D *v)
{
	for (int i = 0; i < nx * ny; i++)
	{
		vx[i] = v->x;
		vy[i] = v->y;
	}
}

void VectorArray2D::incAll(Vector2D*v)
{
	for (int i = 0; i < nx * ny; i++)
	{
		vx[i] += v->x;
		vy[i] += v->y;
	}
}

void VectorArray2D::setRotation(double scale, double phi, double cx, double cy)
{
	int i, j;
	double x, y, cp, sp;
	Vector2D v;

	cp = cos(phi);
	sp = sin(phi);
	for (i = 0; i < nx; i++)
	{
		x = i - cx;
		for (j = 0; j < ny; j++)
		{
			y = j - cy;
			v.x = scale * (x * cp + y * sp) + i - cx;
			v.y = scale * (-x * sp + y * cp) + j - cy;
			set(i, j, &v);
		}
	}
}

void VectorArray2D::copy(VectorArray2D*source)
{
	if (nx != source->nx || ny != source->ny)
		return;
	memcpy(vx, source->vx, sizeof(double) * nx * ny);
	memcpy(vy, source->vy, sizeof(double) * nx * ny);
}

double VectorArray2D::div(int i, int j)
{
	int ip1, im1, jp1, jm1;

	ip1 = (i + 1) % nx;
	im1 = (i - 1) % nx;
	if (im1 < 0)
		im1 = nx - 1;
	jp1 = (j + 1) % ny;
	jm1 = (j - 1) % ny;
	if (jm1 < 0)
		jm1 = ny - 1;

	return 0.5 * ((getx(ip1, j) - getx(im1, j)) / dx + (gety(i, jp1) - gety(i, jm1)) / dy);
}

double VectorArray2D::dVxdx(int i, int j)
{
	int ip1, im1, jp1, jm1;

	ip1 = (i + 1) % nx;
	im1 = (i - 1) % nx;
	if (im1 < 0)
		im1 = nx - 1;
	jp1 = (j + 1) % ny;
	jm1 = (j - 1) % ny;
	if (jm1 < 0)
		jm1 = ny - 1;

	return 0.5 * (getx(ip1, j) - getx(im1, j)) / dx;
}

double VectorArray2D::dVydx(int i, int j)
{
	int ip1, im1, jp1, jm1;

	ip1 = (i + 1) % nx;
	im1 = (i - 1) % nx;
	if (im1 < 0)
		im1 = nx - 1;
	jp1 = (j + 1) % ny;
	jm1 = (j - 1) % ny;
	if (jm1 < 0)
		jm1 = ny - 1;

	return 0.5 * (gety(ip1, j) - gety(im1, j)) / dx;
}

double VectorArray2D::dVxdy(int i, int j) const
{
	int ip1, im1, jp1, jm1;

	ip1 = (i + 1) % nx;
	im1 = (i - 1) % nx;
	if (im1 < 0)
		im1 = nx - 1;
	jp1 = (j + 1) % ny;
	jm1 = (j - 1) % ny;
	if (jm1 < 0)
		jm1 = ny - 1;

	return 0.5 * (getx(i, jp1) - getx(i, jm1)) / dy;
}

double VectorArray2D::dVydy(int i, int j) const
{
	int ip1, im1, jp1, jm1;

	ip1 = (i + 1) % nx;
	im1 = (i - 1) % nx;
	if (im1 < 0)
		im1 = nx - 1;
	jp1 = (j + 1) % ny;
	jm1 = (j - 1) % ny;
	if (jm1 < 0)
		jm1 = ny - 1;

	return 0.5 * (gety(i, jp1) - gety(i, jm1)) / dy;
}

Vector2D VectorArray2D::divComponents(int i, int j) const
{
	int ip1, im1, jp1, jm1;
	Vector2D dc;

	ip1 = (i + 1) % nx;
	im1 = (i - 1) % nx;
	if (im1 < 0)
		im1 = nx - 1;
	jp1 = (j + 1) % ny;
	jm1 = (j - 1) % ny;
	if (jm1 < 0)
		jm1 = ny - 1;

	dc.x = 0.5 * (getx(ip1, j) - getx(im1, j)) / dx;
	dc.y = 0.5 * (gety(i, jp1) - gety(i, jm1)) / dy;
	return dc;
}

double VectorArray2D::jacobian(int i, int j) const
{
	double jxx, jxy, jyx, jyy;
	int ip1, im1, jp1, jm1;

	ip1 = (i + 1) % nx;
	im1 = (i - 1) % nx;
	if (im1 < 0)
		im1 = nx - 1;
	jp1 = (j + 1) % ny;
	jm1 = (j - 1) % ny;
	if (jm1 < 0)
		jm1 = ny - 1;

	jxx = 0.5 * (getx(ip1, j) - getx(im1, j)) / dx;
	jxy = 0.5 * (getx(i, jp1) - getx(i, jm1)) / dy;
	jyx = 0.5 * (gety(ip1, j) - gety(im1, j)) / dx;
	jyy = 0.5 * (gety(i, jp1) - gety(i, jm1)) / dy;
	return (1.0 - jxx) * (1.0 - jyy) - jxy * jyx;
}

Vector2D VectorArray2D::d2x(int i, int j) const
{
	Vector2D v;
	v.x = v.y = 0.0;
	if (i >= 2 && 3 + i <= nx)
	{
		v.x = (-getx(-2 + i, j) - 30 * getx(i, j) + 16 * (getx(-1 + i, j) + getx(1 + i, j)) - getx(2 + i, j)) / 12.;
		v.y = (-gety(-2 + i, j) - 30 * gety(i, j) + 16 * (gety(-1 + i, j) + gety(1 + i, j)) - gety(2 + i, j)) / 12.;
		return v;
	}
	if (i >= 1 && 4 + i <= nx)
	{
		v.x = (11 * getx(-1 + i, j) - 20 * getx(i, j) + 6 * getx(1 + i, j) + 4 * getx(2 + i, j) - getx(3 + i, j)) / 12.;
		v.y = (11 * gety(-1 + i, j) - 20 * gety(i, j) + 6 * gety(1 + i, j) + 4 * gety(2 + i, j) - gety(3 + i, j)) / 12.;
		return v;
	}
	if (i >= 3 && 2 + i <= nx)
	{
		v.x = (-getx(-3 + i, j) + 4 * getx(-2 + i, j) + 6 * getx(-1 + i, j) - 20 * getx(i, j) + 11 * getx(1 + i, j))
				/ 12.;
		v.y = (-gety(-3 + i, j) + 4 * gety(-2 + i, j) + 6 * gety(-1 + i, j) - 20 * gety(i, j) + 11 * gety(1 + i, j))
				/ 12.;
		return v;
	}
	if (5 + i <= nx)
	{
		v.x =
				(35 * getx(i, j) - 104 * getx(1 + i, j) + 114 * getx(2 + i, j) - 56 * getx(3 + i, j)
						+ 11 * getx(4 + i, j)) / 12.;
		v.y =
				(35 * gety(i, j) - 104 * gety(1 + i, j) + 114 * gety(2 + i, j) - 56 * gety(3 + i, j)
						+ 11 * gety(4 + i, j)) / 12.;
		return v;
	}
	if (i >= 4 && 1 + i <= nx && 1 + j <= ny)
	{
		v.x = (11 * getx(-4 + i, j) - 56 * getx(-3 + i, j) + 114 * getx(-2 + i, j) - 104 * getx(-1 + i, j)
				+ 35 * getx(i, j)) / 12.;
		v.y = (11 * gety(-4 + i, j) - 56 * gety(-3 + i, j) + 114 * gety(-2 + i, j) - 104 * gety(-1 + i, j)
				+ 35 * gety(i, j)) / 12.;
		return v;
	}
	return v;
}

Vector2D VectorArray2D::d2y(int i, int j) const
{
	Vector2D v;
	v.x = v.y = 0.0;
	if (j >= 2 && 3 + j <= ny)
	{
		v.x = (-getx(i, -2 + j) - 30 * getx(i, j) + 16 * (getx(i, -1 + j) + getx(i, 1 + j)) - getx(i, 2 + j)) / 12.;
		v.y = (-gety(i, -2 + j) - 30 * gety(i, j) + 16 * (gety(i, -1 + j) + gety(i, 1 + j)) - gety(i, 2 + j)) / 12.;
		return v;
	}
	if (j >= 1 && 4 + j <= ny)
	{
		v.x = (11 * getx(i, -1 + j) - 20 * getx(i, j) + 6 * getx(i, 1 + j) + 4 * getx(i, 2 + j) - getx(i, 3 + j)) / 12.;
		v.y = (11 * gety(i, -1 + j) - 20 * gety(i, j) + 6 * gety(i, 1 + j) + 4 * gety(i, 2 + j) - gety(i, 3 + j)) / 12.;
		return v;
	}
	if (j >= 3 && 2 + j <= ny)
	{
		v.x = (-getx(i, -3 + j) + 4 * getx(i, -2 + j) + 6 * getx(i, -1 + j) - 20 * getx(i, j) + 11 * getx(i, 1 + j))
				/ 12.;
		v.y = (-gety(i, -3 + j) + 4 * gety(i, -2 + j) + 6 * gety(i, -1 + j) - 20 * gety(i, j) + 11 * gety(i, 1 + j))
				/ 12.;
		return v;
	}
	if (1 + i <= nx && 5 + j <= ny)
	{
		v.x =
				(35 * getx(i, j) - 104 * getx(i, 1 + j) + 114 * getx(i, 2 + j) - 56 * getx(i, 3 + j)
						+ 11 * getx(i, 4 + j)) / 12.;
		v.y =
				(35 * gety(i, j) - 104 * gety(i, 1 + j) + 114 * gety(i, 2 + j) - 56 * gety(i, 3 + j)
						+ 11 * gety(i, 4 + j)) / 12.;
		return v;
	}
	if (1 + i <= nx && j >= 4 && 1 + j <= ny)
	{
		v.x = (11 * getx(i, -4 + j) - 56 * getx(i, -3 + j) + 114 * getx(i, -2 + j) - 104 * getx(i, -1 + j)
				+ 35 * getx(i, j)) / 12.;
		v.y = (11 * gety(i, -4 + j) - 56 * gety(i, -3 + j) + 114 * gety(i, -2 + j) - 104 * gety(i, -1 + j)
				+ 35 * gety(i, j)) / 12.;
		return v;
	}
	return v;
}

Vector2D VectorArray2D::dxy(int i, int j) const
{
	Vector2D v;
	v.x = v.y = 0.0;
	if (i >= 2 && 3 + i <= nx && j >= 2 && 3 + j <= ny)
	{
		v.x = (getx(-2 + i, -2 + j) - getx(-2 + i, 2 + j) - 8 * (getx(-2 + i, -1 +\
 j) + getx(-1 + i, -2 + j))
				+ 64 * getx(-1 + i, -1 + j) - 64 * getx(-1 + i, 1\
 + j)
				+ 8 * (getx(-2 + i, 1 + j) + getx(-1 + i, 2 + j)) + 8 * getx(1 + i, -2\
 + j) - 64 * getx(1 + i, -1 + j)
				+ 64 * getx(1 + i, 1 + j) - 8 * getx(1 + i, 2\
 + j) - getx(2 + i, -2 + j) + 8 * getx(2 + i, -1 + j)
				- 8 * getx(2 + i, 1 +\
 j) + getx(2 + i, 2 + j)) / 144.;
		v.y = (gety(-2 + i, -2 + j) - gety(-2 + i, 2 + j) - 8 * (gety(-2 + i, -1 +\
 j) + gety(-1 + i, -2 + j))
				+ 64 * gety(-1 + i, -1 + j) - 64 * gety(-1 + i, 1\
 + j)
				+ 8 * (gety(-2 + i, 1 + j) + gety(-1 + i, 2 + j)) + 8 * gety(1 + i, -2\
 + j) - 64 * gety(1 + i, -1 + j)
				+ 64 * gety(1 + i, 1 + j) - 8 * gety(1 + i, 2\
 + j) - gety(2 + i, -2 + j) + 8 * gety(2 + i, -1 + j)
				- 8 * gety(2 + i, 1 +\
 j) + gety(2 + i, 2 + j)) / 144.;
		return v;
	}
	if (i >= 3 && 2 + i <= nx && j >= 2 && 3 + j <= ny)
	{
		v.x = (-getx(-3 + i, -2 + j) + 8 * getx(-3 + i, -1 + j) - 8 * getx(-3 + i, 1\
 + j) + getx(-3 + i, 2 + j)
				+ 6 * getx(-2 + i, -2 + j) - 48 * getx(-2 + i, -1\
 + j) + 48 * getx(-2 + i, 1 + j)
				- 6 * getx(-2 + i, 2 + j) - 18 * getx(-1 +\
 i, -2 + j) + 144 * getx(-1 + i, -1 + j)
				- 144 * getx(-1 + i, 1 + j) +\
 18 * getx(-1 + i, 2 + j) + 10 * getx(i, -2 + j) - 80 * getx(i, -1 + j)
				+\
 80 * getx(i, 1 + j) - 10 * getx(i, 2 + j)
				+ 3 * (getx(1 + i, -2 + j) -\
 8 * getx(1 + i, -1 + j) + 8 * getx(1 + i, 1 + j) - getx(1 + i, 2 +\
 j)))
				/ 144.;
		v.y = (-gety(-3 + i, -2 + j) + 8 * gety(-3 + i, -1 + j) - 8 * gety(-3 + i, 1\
 + j) + gety(-3 + i, 2 + j)
				+ 6 * gety(-2 + i, -2 + j) - 48 * gety(-2 + i, -1\
 + j) + 48 * gety(-2 + i, 1 + j)
				- 6 * gety(-2 + i, 2 + j) - 18 * gety(-1 +\
 i, -2 + j) + 144 * gety(-1 + i, -1 + j)
				- 144 * gety(-1 + i, 1 + j) +\
 18 * gety(-1 + i, 2 + j) + 10 * gety(i, -2 + j) - 80 * gety(i, -1 + j)
				+\
 80 * gety(i, 1 + j) - 10 * gety(i, 2 + j)
				+ 3 * (gety(1 + i, -2 + j) -\
 8 * gety(1 + i, -1 + j) + 8 * gety(1 + i, 1 + j) - gety(1 + i, 2 +\
 j)))
				/ 144.;
		return v;
	}
	if (i >= 2 && 3 + i <= nx && j >= 3 && 2 + j <= ny)
	{
		v.x = (-getx(-2 + i, -3 + j) + 6 * getx(-2 + i, -2 + j) - 18 * getx(-2 +\
 i, -1 + j) + 10 * getx(-2 + i, j)
				+ 3 * getx(-2 + i, 1 + j) + 8 * getx(-1 +\
 i, -3 + j) - 48 * getx(-1 + i, -2 + j)
				+ 144 * getx(-1 + i, -1 + j) -\
 80 * getx(-1 + i, j) - 24 * getx(-1 + i, 1 + j)
				- 8 * getx(1 + i, -3 + j) +\
 48 * getx(1 + i, -2 + j) - 144 * getx(1 + i, -1 + j) + 80 * getx(1 + i, j)
				+\
 24 * getx(1 + i, 1 + j) + getx(2 + i, -3 + j) - 6 * getx(2 + i, -2 + j) +\
 18 * getx(2 + i, -1 + j)
				- 10 * getx(2 + i, j) - 3 * getx(2 + i, 1 + j)) / 144.;
		v.y = (-gety(-2 + i, -3 + j) + 6 * gety(-2 + i, -2 + j) - 18 * gety(-2 +\
 i, -1 + j) + 10 * gety(-2 + i, j)
				+ 3 * gety(-2 + i, 1 + j) + 8 * gety(-1 +\
 i, -3 + j) - 48 * gety(-1 + i, -2 + j)
				+ 144 * gety(-1 + i, -1 + j) -\
 80 * gety(-1 + i, j) - 24 * gety(-1 + i, 1 + j)
				- 8 * gety(1 + i, -3 + j) +\
 48 * gety(1 + i, -2 + j) - 144 * gety(1 + i, -1 + j) + 80 * gety(1 + i, j)
				+\
 24 * gety(1 + i, 1 + j) + gety(2 + i, -3 + j) - 6 * gety(2 + i, -2 + j) +\
 18 * gety(2 + i, -1 + j)
				- 10 * gety(2 + i, j) - 3 * gety(2 + i, 1 + j)) / 144.;
		return v;
	}
	if (i >= 2 && 3 + i <= nx && j >= 1 && 4 + j <= ny)
	{
		v.x = (-3 * getx(-2 + i, -1 + j) - 10 * getx(-2 + i, j) + 18 * getx(-2 + i, 1 +\
 j) - 6 * getx(-2 + i, 2 + j)
				+ getx(-2 + i, 3 + j) + 24 * getx(-1 + i, -1 +\
 j) + 80 * getx(-1 + i, j) - 144 * getx(-1 + i, 1 + j)
				+ 48 * getx(-1 + i, 2 +\
 j) - 8 * getx(-1 + i, 3 + j) - 24 * getx(1 + i, -1 + j) - 80 * getx(1 + i, j)\

				+ 144 * getx(1 + i, 1 + j) - 48 * getx(1 + i, 2 + j) + 8 * getx(1 + i, 3 + j)\
 + 3 * getx(2 + i, -1 + j)
				+ 10 * getx(2 + i, j) - 18 * getx(2 + i, 1 + j) +\
 6 * getx(2 + i, 2 + j) - getx(2 + i, 3 + j)) / 144.;
		v.y = (-3 * gety(-2 + i, -1 + j) - 10 * gety(-2 + i, j) + 18 * gety(-2 + i, 1 +\
 j) - 6 * gety(-2 + i, 2 + j)
				+ gety(-2 + i, 3 + j) + 24 * gety(-1 + i, -1 +\
 j) + 80 * gety(-1 + i, j) - 144 * gety(-1 + i, 1 + j)
				+ 48 * gety(-1 + i, 2 +\
 j) - 8 * gety(-1 + i, 3 + j) - 24 * gety(1 + i, -1 + j) - 80 * gety(1 + i, j)\

				+ 144 * gety(1 + i, 1 + j) - 48 * gety(1 + i, 2 + j) + 8 * gety(1 + i, 3 + j)\
 + 3 * gety(2 + i, -1 + j)
				+ 10 * gety(2 + i, j) - 18 * gety(2 + i, 1 + j) +\
 6 * gety(2 + i, 2 + j) - gety(2 + i, 3 + j)) / 144.;
		return v;
	}
	if (i >= 1 && 4 + i <= nx && j >= 2 && 3 + j <= ny)
	{
		v.x = (-3 * getx(-1 + i, -2 + j) + 24 * getx(-1 + i, -1 + j) - 24 * getx(-1 +\
 i, 1 + j)
				+ 3 * getx(-1 + i, 2 + j) - 10 * getx(i, -2 + j) + 80 * getx(i, -1 +\
 j) - 80 * getx(i, 1 + j)
				+ 10 * getx(i, 2 + j) + 18 * getx(1 + i, -2 + j) -\
 144 * getx(1 + i, -1 + j)
				+ 144 * getx(1 + i, 1 + j) - 18 * getx(1 + i, 2 + j)\
 - 6 * getx(2 + i, -2 + j)
				+ 48 * getx(2 + i, -1 + j) - 48 * getx(2 + i, 1 + j)\
 + 6 * getx(2 + i, 2 + j) + getx(3 + i, -2 + j)
				- 8 * getx(3 + i, -1 + j) +\
 8 * getx(3 + i, 1 + j) - getx(3 + i, 2 + j)) / 144.;
		v.y = (-3 * gety(-1 + i, -2 + j) + 24 * gety(-1 + i, -1 + j) - 24 * gety(-1 +\
 i, 1 + j)
				+ 3 * gety(-1 + i, 2 + j) - 10 * gety(i, -2 + j) + 80 * gety(i, -1 +\
 j) - 80 * gety(i, 1 + j)
				+ 10 * gety(i, 2 + j) + 18 * gety(1 + i, -2 + j) -\
 144 * gety(1 + i, -1 + j)
				+ 144 * gety(1 + i, 1 + j) - 18 * gety(1 + i, 2 + j)\
 - 6 * gety(2 + i, -2 + j)
				+ 48 * gety(2 + i, -1 + j) - 48 * gety(2 + i, 1 + j)\
 + 6 * gety(2 + i, 2 + j) + gety(3 + i, -2 + j)
				- 8 * gety(3 + i, -1 + j) +\
 8 * gety(3 + i, 1 + j) - gety(3 + i, 2 + j)) / 144.;
		return v;
	}
	if (i >= 3 && 2 + i <= nx && j >= 3 && 2 + j <= ny)
	{
		v.x = (getx(-3 + i, -3 + j) - 6 * (getx(-3 + i, -2 + j) + getx(-2 + i, -3 +\
 j)) + 36 * getx(-2 + i, -2 + j)
				- 108 * (getx(-2 + i, -1 + j) + getx(-1 +\
 i, -2 + j)) + 324 * getx(-1 + i, -1 + j)
				- 10 * (getx(-3 + i, j) + getx(i, -3\
 + j)) + 60 * (getx(-2 + i, j) + getx(i, -2 + j))
				- 180 * (getx(-1 + i, j) +\
 getx(i, -1 + j)) + 100 * getx(i, j)
				- 3 * (getx(-3 + i, 1 + j) + getx(1 +\
 i, -3 + j))
				+ 18 * (getx(-3 + i, -1 + j) + getx(-2 + i, 1 + j) + getx(-1 +\
 i, -3 + j) + getx(1 + i, -2 + j))
				- 54 * (getx(-1 + i, 1 + j) + getx(1 +\
 i, -1 + j)) + 30 * (getx(i, 1 + j) + getx(1 + i, j))
				+ 9 * getx(1 + i, 1 +\
 j)) / 144.;
		v.y = (gety(-3 + i, -3 + j) - 6 * (gety(-3 + i, -2 + j) + gety(-2 + i, -3 +\
 j)) + 36 * gety(-2 + i, -2 + j)
				- 108 * (gety(-2 + i, -1 + j) + gety(-1 +\
 i, -2 + j)) + 324 * gety(-1 + i, -1 + j)
				- 10 * (gety(-3 + i, j) + gety(i, -3\
 + j)) + 60 * (gety(-2 + i, j) + gety(i, -2 + j))
				- 180 * (gety(-1 + i, j) +\
 gety(i, -1 + j)) + 100 * gety(i, j)
				- 3 * (gety(-3 + i, 1 + j) + gety(1 +\
 i, -3 + j))
				+ 18 * (gety(-3 + i, -1 + j) + gety(-2 + i, 1 + j) + gety(-1 +\
 i, -3 + j) + gety(1 + i, -2 + j))
				- 54 * (gety(-1 + i, 1 + j) + gety(1 +\
 i, -1 + j)) + 30 * (gety(i, 1 + j) + gety(1 + i, j))
				+ 9 * gety(1 + i, 1 +\
 j)) / 144.;
		return v;
	}
	if (i >= 3 && 2 + i <= nx && j >= 1 && 4 + j <= ny)
	{
		v.x = (-getx(-3 + i, 3 + j) - 36 * getx(-2 + i, 2 + j) + 6 * (getx(-3 + i, 2\
 + j) + getx(-2 + i, 3 + j))
				- 324 * getx(-1 + i, 1 + j) + 108 * (getx(-2 +\
 i, 1 + j) + getx(-1 + i, 2 + j)) - 100 * getx(i, j)
				+ 180 * (getx(-1 + i, j)\
 + getx(i, 1 + j)) - 60 * (getx(-2 + i, j) + getx(i, 2 + j))
				+ 10 * (getx(-3\
 + i, j) + getx(i, 3 + j)) - 9 * getx(1 + i, -1 + j)
				- 30 * (getx(i, -1 + j) +\
 getx(1 + i, j)) + 54 * (getx(-1 + i, -1 + j) + getx(1 + i, 1 + j))
				-\
 18 * (getx(-3 + i, 1 + j) + getx(-2 + i, -1 + j) + getx(-1 + i, 3 + j) +\
 getx(1 + i, 2 + j))
				+ 3 * (getx(-3 + i, -1 + j) + getx(1 + i, 3 +\
 j))) / 144.;
		v.y = (-gety(-3 + i, 3 + j) - 36 * gety(-2 + i, 2 + j) + 6 * (gety(-3 + i, 2\
 + j) + gety(-2 + i, 3 + j))
				- 324 * gety(-1 + i, 1 + j) + 108 * (gety(-2 +\
 i, 1 + j) + gety(-1 + i, 2 + j)) - 100 * gety(i, j)
				+ 180 * (gety(-1 + i, j)\
 + gety(i, 1 + j)) - 60 * (gety(-2 + i, j) + gety(i, 2 + j))
				+ 10 * (gety(-3\
 + i, j) + gety(i, 3 + j)) - 9 * gety(1 + i, -1 + j)
				- 30 * (gety(i, -1 + j) +\
 gety(1 + i, j)) + 54 * (gety(-1 + i, -1 + j) + gety(1 + i, 1 + j))
				-\
 18 * (gety(-3 + i, 1 + j) + gety(-2 + i, -1 + j) + gety(-1 + i, 3 + j) +\
 gety(1 + i, 2 + j))
				+ 3 * (gety(-3 + i, -1 + j) + gety(1 + i, 3 +\
 j))) / 144.;
		return v;
	}
	if (i >= 1 && 4 + i <= nx && j >= 3 && 2 + j <= ny)
	{
		v.x = (-9 * getx(-1 + i, 1 + j) - 100 * getx(i, j) - 30 * (getx(-1 + i, j) +\
 getx(i, 1 + j))
				- 324 * getx(1 + i, -1 + j) + 180 * (getx(i, -1 + j) +\
 getx(1 + i, j))
				+ 54 * (getx(-1 + i, -1 + j) + getx(1 + i, 1 + j)) -\
 36 * getx(2 + i, -2 + j)
				+ 108 * (getx(1 + i, -2 + j) + getx(2 + i, -1 + j))\
 - 60 * (getx(i, -2 + j) + getx(2 + i, j))
				- getx(3 + i, -3 + j) +\
 6 * (getx(2 + i, -3 + j) + getx(3 + i, -2 + j))
				- 18 * (getx(-1 + i, -2 + j)\
 + getx(1 + i, -3 + j) + getx(2 + i, 1 + j) + getx(3 + i, -1 + j))
				+\
 10 * (getx(i, -3 + j) + getx(3 + i, j)) + 3 * (getx(-1 + i, -3 + j) + getx(3\
 + i, 1 + j))) / 144.;
		v.y = (-9 * gety(-1 + i, 1 + j) - 100 * gety(i, j) - 30 * (gety(-1 + i, j) +\
 gety(i, 1 + j))
				- 324 * gety(1 + i, -1 + j) + 180 * (gety(i, -1 + j) +\
 gety(1 + i, j))
				+ 54 * (gety(-1 + i, -1 + j) + gety(1 + i, 1 + j)) -\
 36 * gety(2 + i, -2 + j)
				+ 108 * (gety(1 + i, -2 + j) + gety(2 + i, -1 + j))\
 - 60 * (gety(i, -2 + j) + gety(2 + i, j))
				- gety(3 + i, -3 + j) +\
 6 * (gety(2 + i, -3 + j) + gety(3 + i, -2 + j))
				- 18 * (gety(-1 + i, -2 + j)\
 + gety(1 + i, -3 + j) + gety(2 + i, 1 + j) + gety(3 + i, -1 + j))
				+\
 10 * (gety(i, -3 + j) + gety(3 + i, j)) + 3 * (gety(-1 + i, -3 + j) + gety(3\
 + i, 1 + j))) / 144.;
		return v;
	}
	if (i >= 1 && 4 + i <= nx && j >= 1 && 4 + j <= ny)
	{
		v.x = (9 * getx(-1 + i, -1 + j) + 30 * (getx(-1 + i, j) + getx(i, -1 + j)) +\
 100 * getx(i, j)
				- 54 * (getx(-1 + i, 1 + j) + getx(1 + i, -1 + j)) -\
 180 * (getx(i, 1 + j) + getx(1 + i, j))
				+ 324 * getx(1 + i, 1 + j) +\
 60 * (getx(i, 2 + j) + getx(2 + i, j))
				- 108 * (getx(1 + i, 2 + j) + getx(2\
 + i, 1 + j)) + 36 * getx(2 + i, 2 + j)
				- 3 * (getx(-1 + i, 3 + j) + getx(3 +\
 i, -1 + j)) - 10 * (getx(i, 3 + j) + getx(3 + i, j))
				+ 18 * (getx(-1 + i, 2 +\
 j) + getx(1 + i, 3 + j) + getx(2 + i, -1 + j) + getx(3 + i, 1 + j))
				-\
 6 * (getx(2 + i, 3 + j) + getx(3 + i, 2 + j)) + getx(3 + i, 3 + j)) / 144.;
		v.y = (9 * gety(-1 + i, -1 + j) + 30 * (gety(-1 + i, j) + gety(i, -1 + j)) +\
 100 * gety(i, j)
				- 54 * (gety(-1 + i, 1 + j) + gety(1 + i, -1 + j)) -\
 180 * (gety(i, 1 + j) + gety(1 + i, j))
				+ 324 * gety(1 + i, 1 + j) +\
 60 * (gety(i, 2 + j) + gety(2 + i, j))
				- 108 * (gety(1 + i, 2 + j) + gety(2\
 + i, 1 + j)) + 36 * gety(2 + i, 2 + j)
				- 3 * (gety(-1 + i, 3 + j) + gety(3 +\
 i, -1 + j)) - 10 * (gety(i, 3 + j) + gety(3 + i, j))
				+ 18 * (gety(-1 + i, 2 +\
 j) + gety(1 + i, 3 + j) + gety(2 + i, -1 + j) + gety(3 + i, 1 + j))
				-\
 6 * (gety(2 + i, 3 + j) + gety(3 + i, 2 + j)) + gety(3 + i, 3 + j)) / 144.;
		return v;
	}
	if (i >= 4 && 1 + i <= nx && j >= 2 && 3 + j <= ny)
	{
		v.x = (3 * getx(-4 + i, -2 + j) - 24 * getx(-4 + i, -1 + j) + 24 * getx(-4 +\
 i, 1 + j)
				- 3 * getx(-4 + i, 2 + j) - 16 * getx(-3 + i, -2 + j) +\
 128 * getx(-3 + i, -1 + j)
				- 128 * getx(-3 + i, 1 + j) + 16 * getx(-3 + i, 2 +\
 j) + 36 * getx(-2 + i, -2 + j)
				- 288 * getx(-2 + i, -1 + j) + 288 * getx(-2 +\
 i, 1 + j) - 36 * getx(-2 + i, 2 + j)
				- 48 * getx(-1 + i, -2 + j) +\
 384 * getx(-1 + i, -1 + j) - 384 * getx(-1 + i, 1 + j)
				+ 48 * getx(-1 + i, 2 +\
 j)
				+ 25 * (getx(i, -2 + j) - 8 * getx(i, -1 + j) + 8 * getx(i, 1 + j) -\
 getx(i, 2 + j))) / 144.;
		v.y = (3 * gety(-4 + i, -2 + j) - 24 * gety(-4 + i, -1 + j) + 24 * gety(-4 +\
 i, 1 + j)
				- 3 * gety(-4 + i, 2 + j) - 16 * gety(-3 + i, -2 + j) +\
 128 * gety(-3 + i, -1 + j)
				- 128 * gety(-3 + i, 1 + j) + 16 * gety(-3 + i, 2 +\
 j) + 36 * gety(-2 + i, -2 + j)
				- 288 * gety(-2 + i, -1 + j) + 288 * gety(-2 +\
 i, 1 + j) - 36 * gety(-2 + i, 2 + j)
				- 48 * gety(-1 + i, -2 + j) +\
 384 * gety(-1 + i, -1 + j) - 384 * gety(-1 + i, 1 + j)
				+ 48 * gety(-1 + i, 2 +\
 j)
				+ 25 * (gety(i, -2 + j) - 8 * gety(i, -1 + j) + 8 * gety(i, 1 + j) -\
 gety(i, 2 + j))) / 144.;
		return v;
	}
	if (i >= 2 && 3 + i <= nx && j >= 4 && 1 + j <= ny)
	{
		v.x = (3 * getx(-2 + i, -4 + j) - 16 * getx(-2 + i, -3 + j) + 36 * getx(-2 +\
 i, -2 + j)
				- 48 * getx(-2 + i, -1 + j) + 25 * getx(-2 + i, j) - 24 * getx(-1 +\
 i, -4 + j)
				+ 128 * getx(-1 + i, -3 + j) - 288 * getx(-1 + i, -2 + j) +\
 384 * getx(-1 + i, -1 + j)
				- 200 * getx(-1 + i, j) + 24 * getx(1 + i, -4 + j)\
 - 128 * getx(1 + i, -3 + j)
				+ 288 * getx(1 + i, -2 + j) - 384 * getx(1 + i, -1\
 + j) + 200 * getx(1 + i, j)
				- 3 * getx(2 + i, -4 + j) + 16 * getx(2 + i, -3 +\
 j) - 36 * getx(2 + i, -2 + j)
				+ 48 * getx(2 + i, -1 + j) - 25 * getx(2 +\
 i, j)) / 144.;
		v.y = (3 * gety(-2 + i, -4 + j) - 16 * gety(-2 + i, -3 + j) + 36 * gety(-2 +\
 i, -2 + j)
				- 48 * gety(-2 + i, -1 + j) + 25 * gety(-2 + i, j) - 24 * gety(-1 +\
 i, -4 + j)
				+ 128 * gety(-1 + i, -3 + j) - 288 * gety(-1 + i, -2 + j) +\
 384 * gety(-1 + i, -1 + j)
				- 200 * gety(-1 + i, j) + 24 * gety(1 + i, -4 + j)\
 - 128 * gety(1 + i, -3 + j)
				+ 288 * gety(1 + i, -2 + j) - 384 * gety(1 + i, -1\
 + j) + 200 * gety(1 + i, j)
				- 3 * gety(2 + i, -4 + j) + 16 * gety(2 + i, -3 +\
 j) - 36 * gety(2 + i, -2 + j)
				+ 48 * gety(2 + i, -1 + j) - 25 * gety(2 +\
 i, j)) / 144.;
		return v;
	}
	if (i >= 2 && 3 + i <= nx && 5 + j <= ny)
	{
		v.x = (-25 * getx(-2 + i, j) + 48 * getx(-2 + i, 1 + j) - 36 * getx(-2 + i, 2 +\
 j) + 16 * getx(-2 + i, 3 + j)
				- 3 * getx(-2 + i, 4 + j) + 200 * getx(-1 +\
 i, j) - 384 * getx(-1 + i, 1 + j)
				+ 288 * getx(-1 + i, 2 + j) - 128 * getx(-1\
 + i, 3 + j) + 24 * getx(-1 + i, 4 + j)
				- 200 * getx(1 + i, j) + 384 * getx(1 +\
 i, 1 + j) - 288 * getx(1 + i, 2 + j)
				+ 128 * getx(1 + i, 3 + j) - 24 * getx(1\
 + i, 4 + j) + 25 * getx(2 + i, j) - 48 * getx(2 + i, 1 + j)
				+ 36 * getx(2 +\
 i, 2 + j) - 16 * getx(2 + i, 3 + j) + 3 * getx(2 + i, 4 + j)) / 144.;
		v.y = (-25 * gety(-2 + i, j) + 48 * gety(-2 + i, 1 + j) - 36 * gety(-2 + i, 2 +\
 j) + 16 * gety(-2 + i, 3 + j)
				- 3 * gety(-2 + i, 4 + j) + 200 * gety(-1 +\
 i, j) - 384 * gety(-1 + i, 1 + j)
				+ 288 * gety(-1 + i, 2 + j) - 128 * gety(-1\
 + i, 3 + j) + 24 * gety(-1 + i, 4 + j)
				- 200 * gety(1 + i, j) + 384 * gety(1 +\
 i, 1 + j) - 288 * gety(1 + i, 2 + j)
				+ 128 * gety(1 + i, 3 + j) - 24 * gety(1\
 + i, 4 + j) + 25 * gety(2 + i, j) - 48 * gety(2 + i, 1 + j)
				+ 36 * gety(2 +\
 i, 2 + j) - 16 * gety(2 + i, 3 + j) + 3 * gety(2 + i, 4 + j)) / 144.;
		return v;
	}
	if (5 + i <= nx && j >= 2 && 3 + j <= ny)
	{
		v.x = (-25 * getx(i, -2 + j) + 200 * getx(i, -1 + j) - 200 * getx(i, 1 + j) +\
 25 * getx(i, 2 + j)
				+ 48 * getx(1 + i, -2 + j) - 384 * getx(1 + i, -1 + j) +\
 384 * getx(1 + i, 1 + j)
				- 48 * getx(1 + i, 2 + j) - 36 * getx(2 + i, -2 + j)\
 + 288 * getx(2 + i, -1 + j)
				- 288 * getx(2 + i, 1 + j) + 36 * getx(2 + i, 2 +\
 j) + 16 * getx(3 + i, -2 + j)
				- 128 * getx(3 + i, -1 + j) + 128 * getx(3 +\
 i, 1 + j) - 16 * getx(3 + i, 2 + j)
				- 3 * getx(4 + i, -2 + j) + 24 * getx(4 +\
 i, -1 + j) - 24 * getx(4 + i, 1 + j)
				+ 3 * getx(4 + i, 2 + j)) / 144.;
		v.y = (-25 * gety(i, -2 + j) + 200 * gety(i, -1 + j) - 200 * gety(i, 1 + j) +\
 25 * gety(i, 2 + j)
				+ 48 * gety(1 + i, -2 + j) - 384 * gety(1 + i, -1 + j) +\
 384 * gety(1 + i, 1 + j)
				- 48 * gety(1 + i, 2 + j) - 36 * gety(2 + i, -2 + j)\
 + 288 * gety(2 + i, -1 + j)
				- 288 * gety(2 + i, 1 + j) + 36 * gety(2 + i, 2 +\
 j) + 16 * gety(3 + i, -2 + j)
				- 128 * gety(3 + i, -1 + j) + 128 * gety(3 +\
 i, 1 + j) - 16 * gety(3 + i, 2 + j)
				- 3 * gety(4 + i, -2 + j) + 24 * gety(4 +\
 i, -1 + j) - 24 * gety(4 + i, 1 + j)
				+ 3 * gety(4 + i, 2 + j)) / 144.;
		return v;
	}
	if (i >= 4 && 1 + i <= nx && j >= 3 && 2 + j <= ny)
	{
		v.x = (-3 * getx(-4 + i, -3 + j) + 18 * getx(-4 + i, -2 + j) - 54 * getx(-4 +\
 i, -1 + j)
				+ 30 * getx(-4 + i, j) + 9 * getx(-4 + i, 1 + j) + 16 * getx(-3 +\
 i, -3 + j)
				- 96 * getx(-3 + i, -2 + j) + 288 * getx(-3 + i, -1 + j) -\
 160 * getx(-3 + i, j)
				- 48 * getx(-3 + i, 1 + j) - 36 * getx(-2 + i, -3 + j) +\
 216 * getx(-2 + i, -2 + j)
				- 648 * getx(-2 + i, -1 + j) + 360 * getx(-2 +\
 i, j) + 108 * getx(-2 + i, 1 + j)
				+ 48 * getx(-1 + i, -3 + j) - 288 * getx(-1\
 + i, -2 + j) + 864 * getx(-1 + i, -1 + j)
				- 480 * getx(-1 + i, j) -\
 144 * getx(-1 + i, 1 + j)
				- 25
						* (getx(i, -3 + j) - 6 * getx(i, -2 + j) +\
 18 * getx(i, -1 + j) - 10 * getx(i, j)
								- 3 * getx(i, 1 + j))) / 144.;
		v.y = (-3 * gety(-4 + i, -3 + j) + 18 * gety(-4 + i, -2 + j) - 54 * gety(-4 +\
 i, -1 + j)
				+ 30 * gety(-4 + i, j) + 9 * gety(-4 + i, 1 + j) + 16 * gety(-3 +\
 i, -3 + j)
				- 96 * gety(-3 + i, -2 + j) + 288 * gety(-3 + i, -1 + j) -\
 160 * gety(-3 + i, j)
				- 48 * gety(-3 + i, 1 + j) - 36 * gety(-2 + i, -3 + j) +\
 216 * gety(-2 + i, -2 + j)
				- 648 * gety(-2 + i, -1 + j) + 360 * gety(-2 +\
 i, j) + 108 * gety(-2 + i, 1 + j)
				+ 48 * gety(-1 + i, -3 + j) - 288 * gety(-1\
 + i, -2 + j) + 864 * gety(-1 + i, -1 + j)
				- 480 * gety(-1 + i, j) -\
 144 * gety(-1 + i, 1 + j)
				- 25
						* (gety(i, -3 + j) - 6 * gety(i, -2 + j) +\
 18 * gety(i, -1 + j) - 10 * gety(i, j)
								- 3 * gety(i, 1 + j))) / 144.;
		return v;
	}
	if (i >= 4 && 1 + i <= nx && j >= 1 && 4 + j <= ny)
	{
		v.x = (-9 * getx(-4 + i, -1 + j) - 30 * getx(-4 + i, j) + 54 * getx(-4 + i, 1 +\
 j) - 18 * getx(-4 + i, 2 + j)
				+ 3 * getx(-4 + i, 3 + j) + 48 * getx(-3 + i, -1\
 + j) + 160 * getx(-3 + i, j)
				- 288 * getx(-3 + i, 1 + j) + 96 * getx(-3 + i, 2\
 + j) - 16 * getx(-3 + i, 3 + j)
				- 108 * getx(-2 + i, -1 + j) - 360 * getx(-2\
 + i, j) + 648 * getx(-2 + i, 1 + j)
				- 216 * getx(-2 + i, 2 + j) + 36 * getx(-2\
 + i, 3 + j) + 144 * getx(-1 + i, -1 + j)
				+ 480 * getx(-1 + i, j) -\
 864 * getx(-1 + i, 1 + j) + 288 * getx(-1 + i, 2 + j)
				- 48 * getx(-1 + i, 3 +\
 j)
				- 25
						* (3 * getx(i, -1 + j) + 10 * getx(i, j) - 18 * getx(i, 1 + j) +\
 6 * getx(i, 2 + j)
								- getx(i, 3 + j))) / 144.;
		v.y = (-9 * gety(-4 + i, -1 + j) - 30 * gety(-4 + i, j) + 54 * gety(-4 + i, 1 +\
 j) - 18 * gety(-4 + i, 2 + j)
				+ 3 * gety(-4 + i, 3 + j) + 48 * gety(-3 + i, -1\
 + j) + 160 * gety(-3 + i, j)
				- 288 * gety(-3 + i, 1 + j) + 96 * gety(-3 + i, 2\
 + j) - 16 * gety(-3 + i, 3 + j)
				- 108 * gety(-2 + i, -1 + j) - 360 * gety(-2\
 + i, j) + 648 * gety(-2 + i, 1 + j)
				- 216 * gety(-2 + i, 2 + j) + 36 * gety(-2\
 + i, 3 + j) + 144 * gety(-1 + i, -1 + j)
				+ 480 * gety(-1 + i, j) -\
 864 * gety(-1 + i, 1 + j) + 288 * gety(-1 + i, 2 + j)
				- 48 * gety(-1 + i, 3 +\
 j)
				- 25
						* (3 * gety(i, -1 + j) + 10 * gety(i, j) - 18 * gety(i, 1 + j) +\
 6 * gety(i, 2 + j)
								- gety(i, 3 + j))) / 144.;
		return v;
	}
	if (i >= 3 && 2 + i <= nx && j >= 4 && 1 + j <= ny)
	{
		v.x = (-3 * getx(-3 + i, -4 + j) + 16 * getx(-3 + i, -3 + j) - 36 * getx(-3 +\
 i, -2 + j)
				+ 48 * getx(-3 + i, -1 + j) - 25 * getx(-3 + i, j) + 18 * getx(-2 +\
 i, -4 + j)
				- 96 * getx(-2 + i, -3 + j) + 216 * getx(-2 + i, -2 + j) -\
 288 * getx(-2 + i, -1 + j)
				+ 150 * getx(-2 + i, j) - 54 * getx(-1 + i, -4 + j)\
 + 288 * getx(-1 + i, -3 + j)
				- 648 * getx(-1 + i, -2 + j) + 864 * getx(-1 +\
 i, -1 + j) - 450 * getx(-1 + i, j)
				+ 30 * getx(i, -4 + j) - 160 * getx(i, -3 +\
 j) + 360 * getx(i, -2 + j) - 480 * getx(i, -1 + j)
				+ 250 * getx(i, j) +\
 9 * getx(1 + i, -4 + j) - 48 * getx(1 + i, -3 + j) + 108 * getx(1 + i, -2 + j)\

				- 144 * getx(1 + i, -1 + j) + 75 * getx(1 + i, j)) / 144.;
		v.y = (-3 * gety(-3 + i, -4 + j) + 16 * gety(-3 + i, -3 + j) - 36 * gety(-3 +\
 i, -2 + j)
				+ 48 * gety(-3 + i, -1 + j) - 25 * gety(-3 + i, j) + 18 * gety(-2 +\
 i, -4 + j)
				- 96 * gety(-2 + i, -3 + j) + 216 * gety(-2 + i, -2 + j) -\
 288 * gety(-2 + i, -1 + j)
				+ 150 * gety(-2 + i, j) - 54 * gety(-1 + i, -4 + j)\
 + 288 * gety(-1 + i, -3 + j)
				- 648 * gety(-1 + i, -2 + j) + 864 * gety(-1 +\
 i, -1 + j) - 450 * gety(-1 + i, j)
				+ 30 * gety(i, -4 + j) - 160 * gety(i, -3 +\
 j) + 360 * gety(i, -2 + j) - 480 * gety(i, -1 + j)
				+ 250 * gety(i, j) +\
 9 * gety(1 + i, -4 + j) - 48 * gety(1 + i, -3 + j) + 108 * gety(1 + i, -2 + j)\

				- 144 * gety(1 + i, -1 + j) + 75 * gety(1 + i, j)) / 144.;
		return v;
	}
	if (i >= 3 && 2 + i <= nx && 5 + j <= ny)
	{
		v.x = (25 * getx(-3 + i, j) - 48 * getx(-3 + i, 1 + j) + 36 * getx(-3 + i, 2 +\
 j) - 16 * getx(-3 + i, 3 + j)
				+ 3 * getx(-3 + i, 4 + j) - 150 * getx(-2 +\
 i, j) + 288 * getx(-2 + i, 1 + j)
				- 216 * getx(-2 + i, 2 + j) + 96 * getx(-2 +\
 i, 3 + j) - 18 * getx(-2 + i, 4 + j)
				+ 450 * getx(-1 + i, j) - 864 * getx(-1 +\
 i, 1 + j) + 648 * getx(-1 + i, 2 + j)
				- 288 * getx(-1 + i, 3 + j) +\
 54 * getx(-1 + i, 4 + j) - 250 * getx(i, j) + 480 * getx(i, 1 + j)
				-\
 360 * getx(i, 2 + j) + 160 * getx(i, 3 + j) - 30 * getx(i, 4 + j) - 75 * getx(1\
 + i, j)
				+ 144 * getx(1 + i, 1 + j) - 108 * getx(1 + i, 2 + j) + 48 * getx(1 +\
 i, 3 + j)
				- 9 * getx(1 + i, 4 + j)) / 144.;
		v.y = (25 * gety(-3 + i, j) - 48 * gety(-3 + i, 1 + j) + 36 * gety(-3 + i, 2 +\
 j) - 16 * gety(-3 + i, 3 + j)
				+ 3 * gety(-3 + i, 4 + j) - 150 * gety(-2 +\
 i, j) + 288 * gety(-2 + i, 1 + j)
				- 216 * gety(-2 + i, 2 + j) + 96 * gety(-2 +\
 i, 3 + j) - 18 * gety(-2 + i, 4 + j)
				+ 450 * gety(-1 + i, j) - 864 * gety(-1 +\
 i, 1 + j) + 648 * gety(-1 + i, 2 + j)
				- 288 * gety(-1 + i, 3 + j) +\
 54 * gety(-1 + i, 4 + j) - 250 * gety(i, j) + 480 * gety(i, 1 + j)
				-\
 360 * gety(i, 2 + j) + 160 * gety(i, 3 + j) - 30 * gety(i, 4 + j) - 75 * gety(1\
 + i, j)
				+ 144 * gety(1 + i, 1 + j) - 108 * gety(1 + i, 2 + j) + 48 * gety(1 +\
 i, 3 + j)
				- 9 * gety(1 + i, 4 + j)) / 144.;
		return v;
	}
	if (i >= 1 && 4 + i <= nx && j >= 4 && 1 + j <= ny)
	{
		v.x = (-9 * getx(-1 + i, -4 + j) + 48 * getx(-1 + i, -3 + j) - 108 * getx(-1 +\
 i, -2 + j)
				+ 144 * getx(-1 + i, -1 + j) - 75 * getx(-1 + i, j) -\
 30 * getx(i, -4 + j) + 160 * getx(i, -3 + j)
				- 360 * getx(i, -2 + j) +\
 480 * getx(i, -1 + j) - 250 * getx(i, j) + 54 * getx(1 + i, -4 + j)
				-\
 288 * getx(1 + i, -3 + j) + 648 * getx(1 + i, -2 + j) - 864 * getx(1 + i, -1 +\
 j)
				+ 450 * getx(1 + i, j) - 18 * getx(2 + i, -4 + j) + 96 * getx(2 + i, -3 +\
 j)
				- 216 * getx(2 + i, -2 + j) + 288 * getx(2 + i, -1 + j) - 150 * getx(2 +\
 i, j)
				+ 3 * getx(3 + i, -4 + j) - 16 * getx(3 + i, -3 + j) + 36 * getx(3 +\
 i, -2 + j)
				- 48 * getx(3 + i, -1 + j) + 25 * getx(3 + i, j)) / 144.;
		v.y = (-9 * gety(-1 + i, -4 + j) + 48 * gety(-1 + i, -3 + j) - 108 * gety(-1 +\
 i, -2 + j)
				+ 144 * gety(-1 + i, -1 + j) - 75 * gety(-1 + i, j) -\
 30 * gety(i, -4 + j) + 160 * gety(i, -3 + j)
				- 360 * gety(i, -2 + j) +\
 480 * gety(i, -1 + j) - 250 * gety(i, j) + 54 * gety(1 + i, -4 + j)
				-\
 288 * gety(1 + i, -3 + j) + 648 * gety(1 + i, -2 + j) - 864 * gety(1 + i, -1 +\
 j)
				+ 450 * gety(1 + i, j) - 18 * gety(2 + i, -4 + j) + 96 * gety(2 + i, -3 +\
 j)
				- 216 * gety(2 + i, -2 + j) + 288 * gety(2 + i, -1 + j) - 150 * gety(2 +\
 i, j)
				+ 3 * gety(3 + i, -4 + j) - 16 * gety(3 + i, -3 + j) + 36 * gety(3 +\
 i, -2 + j)
				- 48 * gety(3 + i, -1 + j) + 25 * gety(3 + i, j)) / 144.;
		return v;
	}
	if (i >= 1 && 4 + i <= nx && 5 + j <= ny)
	{
		v.x = (75 * getx(-1 + i, j) - 144 * getx(-1 + i, 1 + j) + 108 * getx(-1 + i, 2\
 + j) - 48 * getx(-1 + i, 3 + j)
				+ 9 * getx(-1 + i, 4 + j) + 250 * getx(i, j) -\
 480 * getx(i, 1 + j) + 360 * getx(i, 2 + j)
				- 160 * getx(i, 3 + j) +\
 30 * getx(i, 4 + j) - 450 * getx(1 + i, j) + 864 * getx(1 + i, 1 + j)
				-\
 648 * getx(1 + i, 2 + j) + 288 * getx(1 + i, 3 + j) - 54 * getx(1 + i, 4 + j)\
 + 150 * getx(2 + i, j)
				- 288 * getx(2 + i, 1 + j) + 216 * getx(2 + i, 2 + j) -\
 96 * getx(2 + i, 3 + j)
				+ 18 * getx(2 + i, 4 + j) - 25 * getx(3 + i, j) +\
 48 * getx(3 + i, 1 + j) - 36 * getx(3 + i, 2 + j)
				+ 16 * getx(3 + i, 3 + j) -\
 3 * getx(3 + i, 4 + j)) / 144.;
		v.y = (75 * gety(-1 + i, j) - 144 * gety(-1 + i, 1 + j) + 108 * gety(-1 + i, 2\
 + j) - 48 * gety(-1 + i, 3 + j)
				+ 9 * gety(-1 + i, 4 + j) + 250 * gety(i, j) -\
 480 * gety(i, 1 + j) + 360 * gety(i, 2 + j)
				- 160 * gety(i, 3 + j) +\
 30 * gety(i, 4 + j) - 450 * gety(1 + i, j) + 864 * gety(1 + i, 1 + j)
				-\
 648 * gety(1 + i, 2 + j) + 288 * gety(1 + i, 3 + j) - 54 * gety(1 + i, 4 + j)\
 + 150 * gety(2 + i, j)
				- 288 * gety(2 + i, 1 + j) + 216 * gety(2 + i, 2 + j) -\
 96 * gety(2 + i, 3 + j)
				+ 18 * gety(2 + i, 4 + j) - 25 * gety(3 + i, j) +\
 48 * gety(3 + i, 1 + j) - 36 * gety(3 + i, 2 + j)
				+ 16 * gety(3 + i, 3 + j) -\
 3 * gety(3 + i, 4 + j)) / 144.;
		return v;
	}
	if (5 + i <= nx && j >= 3 && 2 + j <= ny)
	{
		v.x = (25 * getx(i, -3 + j) - 150 * getx(i, -2 + j) + 450 * getx(i, -1 + j) -\
 250 * getx(i, j)
				- 75 * getx(i, 1 + j) - 48 * getx(1 + i, -3 + j) + 288 * getx(1\
 + i, -2 + j)
				- 864 * getx(1 + i, -1 + j) + 480 * getx(1 + i, j) + 144 * getx(1\
 + i, 1 + j)
				+ 36 * getx(2 + i, -3 + j) - 216 * getx(2 + i, -2 + j) +\
 648 * getx(2 + i, -1 + j)
				- 360 * getx(2 + i, j) - 108 * getx(2 + i, 1 + j) -\
 16 * getx(3 + i, -3 + j)
				+ 96 * getx(3 + i, -2 + j) - 288 * getx(3 + i, -1 +\
 j) + 160 * getx(3 + i, j)
				+ 48 * getx(3 + i, 1 + j) + 3 * getx(4 + i, -3 + j)\
 - 18 * getx(4 + i, -2 + j)
				+ 54 * getx(4 + i, -1 + j) - 30 * getx(4 + i, j) -\
 9 * getx(4 + i, 1 + j)) / 144.;
		v.y = (25 * gety(i, -3 + j) - 150 * gety(i, -2 + j) + 450 * gety(i, -1 + j) -\
 250 * gety(i, j)
				- 75 * gety(i, 1 + j) - 48 * gety(1 + i, -3 + j) + 288 * gety(1\
 + i, -2 + j)
				- 864 * gety(1 + i, -1 + j) + 480 * gety(1 + i, j) + 144 * gety(1\
 + i, 1 + j)
				+ 36 * gety(2 + i, -3 + j) - 216 * gety(2 + i, -2 + j) +\
 648 * gety(2 + i, -1 + j)
				- 360 * gety(2 + i, j) - 108 * gety(2 + i, 1 + j) -\
 16 * gety(3 + i, -3 + j)
				+ 96 * gety(3 + i, -2 + j) - 288 * gety(3 + i, -1 +\
 j) + 160 * gety(3 + i, j)
				+ 48 * gety(3 + i, 1 + j) + 3 * gety(4 + i, -3 + j)\
 - 18 * gety(4 + i, -2 + j)
				+ 54 * gety(4 + i, -1 + j) - 30 * gety(4 + i, j) -\
 9 * gety(4 + i, 1 + j)) / 144.;
		return v;
	}
	if (5 + i <= nx && j >= 1 && 4 + j <= ny)
	{
		v.x = (75 * getx(i, -1 + j) + 250 * getx(i, j) - 450 * getx(i, 1 + j) +\
 150 * getx(i, 2 + j)
				- 25 * getx(i, 3 + j) - 144 * getx(1 + i, -1 + j) -\
 480 * getx(1 + i, j) + 864 * getx(1 + i, 1 + j)
				- 288 * getx(1 + i, 2 + j) +\
 48 * getx(1 + i, 3 + j) + 108 * getx(2 + i, -1 + j)
				+ 360 * getx(2 + i, j) -\
 648 * getx(2 + i, 1 + j) + 216 * getx(2 + i, 2 + j) - 36 * getx(2 + i, 3 + j)\

				- 48 * getx(3 + i, -1 + j) - 160 * getx(3 + i, j) + 288 * getx(3 + i, 1 + j) -\
 96 * getx(3 + i, 2 + j)
				+ 16 * getx(3 + i, 3 + j) + 9 * getx(4 + i, -1 + j) +\
 30 * getx(4 + i, j) - 54 * getx(4 + i, 1 + j)
				+ 18 * getx(4 + i, 2 + j) -\
 3 * getx(4 + i, 3 + j)) / 144.;
		v.y = (75 * gety(i, -1 + j) + 250 * gety(i, j) - 450 * gety(i, 1 + j) +\
 150 * gety(i, 2 + j)
				- 25 * gety(i, 3 + j) - 144 * gety(1 + i, -1 + j) -\
 480 * gety(1 + i, j) + 864 * gety(1 + i, 1 + j)
				- 288 * gety(1 + i, 2 + j) +\
 48 * gety(1 + i, 3 + j) + 108 * gety(2 + i, -1 + j)
				+ 360 * gety(2 + i, j) -\
 648 * gety(2 + i, 1 + j) + 216 * gety(2 + i, 2 + j) - 36 * gety(2 + i, 3 + j)\

				- 48 * gety(3 + i, -1 + j) - 160 * gety(3 + i, j) + 288 * gety(3 + i, 1 + j) -\
 96 * gety(3 + i, 2 + j)
				+ 16 * gety(3 + i, 3 + j) + 9 * gety(4 + i, -1 + j) +\
 30 * gety(4 + i, j) - 54 * gety(4 + i, 1 + j)
				+ 18 * gety(4 + i, 2 + j) -\
 3 * gety(4 + i, 3 + j)) / 144.;
		return v;
	}
	if (i >= 4 && 1 + i <= nx && j >= 4 && 1 + j <= ny)
	{
		v.x = (9 * getx(-4 + i, -4 + j) + 75 * getx(-4 + i, j) - 48 * (getx(-4 + i, -3\
 + j) + getx(-3 + i, -4 + j))
				+ 256 * getx(-3 + i, -3 + j) - 400 * getx(-3 +\
 i, j)
				+ 108 * (getx(-4 + i, -2 + j) + getx(-2 + i, -4 + j))
				- 576 * (getx(-3\
 + i, -2 + j) + getx(-2 + i, -3 + j)) + 1296 * getx(-2 + i, -2 + j)
				+\
 900 * getx(-2 + i, j) - 144 * (getx(-4 + i, -1 + j) + getx(-1 + i, -4 + j))\

				+ 768 * (getx(-3 + i, -1 + j) + getx(-1 + i, -3 + j))
				- 1728 * (getx(-2 +\
 i, -1 + j) + getx(-1 + i, -2 + j)) + 2304 * getx(-1 + i, -1 + j)
				-\
 25
						* (-3 * getx(i, -4 + j) + 16 * getx(i, -3 + j) - 36 * getx(i, -2 + j)
								+\
 48 * (getx(-1 + i, j) + getx(i, -1 + j)) - 25 * getx(i, j))) / 144.;
		v.y = (9 * gety(-4 + i, -4 + j) + 75 * gety(-4 + i, j) - 48 * (gety(-4 + i, -3\
 + j) + gety(-3 + i, -4 + j))
				+ 256 * gety(-3 + i, -3 + j) - 400 * gety(-3 +\
 i, j)
				+ 108 * (gety(-4 + i, -2 + j) + gety(-2 + i, -4 + j))
				- 576 * (gety(-3\
 + i, -2 + j) + gety(-2 + i, -3 + j)) + 1296 * gety(-2 + i, -2 + j)
				+\
 900 * gety(-2 + i, j) - 144 * (gety(-4 + i, -1 + j) + gety(-1 + i, -4 + j))\

				+ 768 * (gety(-3 + i, -1 + j) + gety(-1 + i, -3 + j))
				- 1728 * (gety(-2 +\
 i, -1 + j) + gety(-1 + i, -2 + j)) + 2304 * gety(-1 + i, -1 + j)
				-\
 25
						* (-3 * gety(i, -4 + j) + 16 * gety(i, -3 + j) - 36 * gety(i, -2 + j)
								+\
 48 * (gety(-1 + i, j) + gety(i, -1 + j)) - 25 * gety(i, j))) / 144.;
		return v;
	}
	if (i >= 4 && 1 + i <= nx && 5 + j <= ny)
	{
		v.x = (-75 * getx(-4 + i, j) - 9 * getx(-4 + i, 4 + j) + 400 * getx(-3 + i, j)\
 - 256 * getx(-3 + i, 3 + j)
				+ 48 * (getx(-4 + i, 3 + j) + getx(-3 + i, 4 +\
 j)) - 900 * getx(-2 + i, j)
				- 1296 * getx(-2 + i, 2 + j) + 576 * (getx(-3 +\
 i, 2 + j) + getx(-2 + i, 3 + j))
				- 108 * (getx(-4 + i, 2 + j) + getx(-2 +\
 i, 4 + j)) + 1200 * getx(-1 + i, j)
				- 2304 * getx(-1 + i, 1 + j) +\
 1728 * (getx(-2 + i, 1 + j) + getx(-1 + i, 2 + j))
				- 768 * (getx(-3 + i, 1 +\
 j) + getx(-1 + i, 3 + j))
				+ 144 * (getx(-4 + i, 1 + j) + getx(-1 + i, 4 +\
 j))
				- 25
						* (25 * getx(i, j) - 48 * getx(i, 1 + j) + 36 * getx(i, 2 + j) -\
 16 * getx(i, 3 + j)
								+ 3 * getx(i, 4 + j))) / 144.;
		v.y = (-75 * gety(-4 + i, j) - 9 * gety(-4 + i, 4 + j) + 400 * gety(-3 + i, j)\
 - 256 * gety(-3 + i, 3 + j)
				+ 48 * (gety(-4 + i, 3 + j) + gety(-3 + i, 4 +\
 j)) - 900 * gety(-2 + i, j)
				- 1296 * gety(-2 + i, 2 + j) + 576 * (gety(-3 +\
 i, 2 + j) + gety(-2 + i, 3 + j))
				- 108 * (gety(-4 + i, 2 + j) + gety(-2 +\
 i, 4 + j)) + 1200 * gety(-1 + i, j)
				- 2304 * gety(-1 + i, 1 + j) +\
 1728 * (gety(-2 + i, 1 + j) + gety(-1 + i, 2 + j))
				- 768 * (gety(-3 + i, 1 +\
 j) + gety(-1 + i, 3 + j))
				+ 144 * (gety(-4 + i, 1 + j) + gety(-1 + i, 4 +\
 j))
				- 25
						* (25 * gety(i, j) - 48 * gety(i, 1 + j) + 36 * gety(i, 2 + j) -\
 16 * gety(i, 3 + j)
								+ 3 * gety(i, 4 + j))) / 144.;
		return v;
	}
	if (5 + i <= nx && j >= 4 && 1 + j <= ny)
	{
		v.x = (-625 * getx(i, j) - 2304 * getx(1 + i, -1 + j) + 1200 * (getx(i, -1 + j)\
 + getx(1 + i, j))
				- 1296 * getx(2 + i, -2 + j) + 1728 * (getx(1 + i, -2 + j)\
 + getx(2 + i, -1 + j))
				- 900 * (getx(i, -2 + j) + getx(2 + i, j)) -\
 256 * getx(3 + i, -3 + j)
				+ 576 * (getx(2 + i, -3 + j) + getx(3 + i, -2 +\
 j))
				- 768 * (getx(1 + i, -3 + j) + getx(3 + i, -1 + j)) + 400 * (getx(i, -3\
 + j) + getx(3 + i, j))
				- 9 * getx(4 + i, -4 + j) + 48 * (getx(3 + i, -4 + j)\
 + getx(4 + i, -3 + j))
				- 108 * (getx(2 + i, -4 + j) + getx(4 + i, -2 + j))\
 + 144 * (getx(1 + i, -4 + j) + getx(4 + i, -1 + j))
				- 75 * (getx(i, -4 + j)\
 + getx(4 + i, j))) / 144.;
		v.y = (-625 * gety(i, j) - 2304 * gety(1 + i, -1 + j) + 1200 * (gety(i, -1 + j)\
 + gety(1 + i, j))
				- 1296 * gety(2 + i, -2 + j) + 1728 * (gety(1 + i, -2 + j)\
 + gety(2 + i, -1 + j))
				- 900 * (gety(i, -2 + j) + gety(2 + i, j)) -\
 256 * gety(3 + i, -3 + j)
				+ 576 * (gety(2 + i, -3 + j) + gety(3 + i, -2 +\
 j))
				- 768 * (gety(1 + i, -3 + j) + gety(3 + i, -1 + j)) + 400 * (gety(i, -3\
 + j) + gety(3 + i, j))
				- 9 * gety(4 + i, -4 + j) + 48 * (gety(3 + i, -4 + j)\
 + gety(4 + i, -3 + j))
				- 108 * (gety(2 + i, -4 + j) + gety(4 + i, -2 + j))\
 + 144 * (gety(1 + i, -4 + j) + gety(4 + i, -1 + j))
				- 75 * (gety(i, -4 + j)\
 + gety(4 + i, j))) / 144.;
		return v;
	}
	if (5 + i <= nx && 5 + j <= ny)
	{
		v.x = (625 * getx(i, j) - 1200 * (getx(i, 1 + j) + getx(1 + i, j)) +\
 2304 * getx(1 + i, 1 + j)
				+ 900 * (getx(i, 2 + j) + getx(2 + i, j)) -\
 1728 * (getx(1 + i, 2 + j) + getx(2 + i, 1 + j))
				+ 1296 * getx(2 + i, 2 + j)\
 - 400 * (getx(i, 3 + j) + getx(3 + i, j))
				+ 768 * (getx(1 + i, 3 + j) +\
 getx(3 + i, 1 + j)) - 576 * (getx(2 + i, 3 + j) + getx(3 + i, 2 + j))
				+\
 256 * getx(3 + i, 3 + j) + 75 * (getx(i, 4 + j) + getx(4 + i, j))
				-\
 144 * (getx(1 + i, 4 + j) + getx(4 + i, 1 + j)) + 108 * (getx(2 + i, 4 + j)\
 + getx(4 + i, 2 + j))
				- 48 * (getx(3 + i, 4 + j) + getx(4 + i, 3 + j)) +\
 9 * getx(4 + i, 4 + j)) / 144.;
		v.y = (625 * gety(i, j) - 1200 * (gety(i, 1 + j) + gety(1 + i, j)) +\
 2304 * gety(1 + i, 1 + j)
				+ 900 * (gety(i, 2 + j) + gety(2 + i, j)) -\
 1728 * (gety(1 + i, 2 + j) + gety(2 + i, 1 + j))
				+ 1296 * gety(2 + i, 2 + j)\
 - 400 * (gety(i, 3 + j) + gety(3 + i, j))
				+ 768 * (gety(1 + i, 3 + j) +\
 gety(3 + i, 1 + j)) - 576 * (gety(2 + i, 3 + j) + gety(3 + i, 2 + j))
				+\
 256 * gety(3 + i, 3 + j) + 75 * (gety(i, 4 + j) + gety(4 + i, j))
				-\
 144 * (gety(1 + i, 4 + j) + gety(4 + i, 1 + j)) + 108 * (gety(2 + i, 4 + j)\
 + gety(4 + i, 2 + j))
				- 48 * (gety(3 + i, 4 + j) + gety(4 + i, 3 + j)) +\
 9 * gety(4 + i, 4 + j)) / 144.;
		return v;
	}
	return v;
}

Vector2D VectorArray2D::laplace(int i, int j)
{
	Vector2D dx2, dy2, lap;
	dx2 = d2x(i, j);
	dy2 = d2y(i, j);
	lap.x = dx2.x + dy2.x;
	lap.y = dx2.y + dy2.y;
	return lap;
}

bool VectorArray2D::save(const char*fname)
{
	FILE * f;
	int magic = VectorArrayMagic;
	f = fopen(fname, "wb");
	if (!f)
		return false;
	fwrite(&magic, sizeof(int), 1, f);
	fwrite(&dx, sizeof(double), 1, f);
	fwrite(&dy, sizeof(double), 1, f);
	fwrite(&nx, sizeof(int), 1, f);
	fwrite(&ny, sizeof(int), 1, f);
	fwrite(vx, sizeof(double), 2 * nx * ny, f);
	fclose(f);
	return true;
}

bool VectorArray2D::save(const char* fname, struct parameters *param)
{
	FILE * f;
	int magic = VectorArrayMagic2;
	f = fopen(fname, "wb");
	if (!f)
		return false;
	fwrite(&magic, sizeof(int), 1, f);
	fwrite(&dx, sizeof(double), 1, f);
	fwrite(&dy, sizeof(double), 1, f);
	fwrite(&nx, sizeof(int), 1, f);
	fwrite(&ny, sizeof(int), 1, f);
	fwrite(&param->end, sizeof(double), 1, f);
	fwrite(&param->error, sizeof(double), 1, f);
	fwrite(&param->alpha, sizeof(double), 1, f);
	fwrite(&param->vortex_weight, sizeof(double), 1, f);
	fwrite(&param->mu, sizeof(double), 1, f);
	fwrite(&param->lambda, sizeof(double), 1, f);
	fwrite(param->boundary.c_str(), sizeof(char), param->boundary.size()+1, f);
	fwrite(param->method.c_str(), sizeof(char), param->method.size()+1, f);
	fwrite(&param->actual_error, sizeof(double), 1, f);
	fwrite(&param->actual_time, sizeof(double), 1, f);
	fwrite(vx, sizeof(double), 2 * nx * ny, f);
	fclose(f);
	return true;
}

bool VectorArray2D::load(const char*fname)
{
	FILE * f;
	int magic;
	f = fopen(fname, "rb");
	if (!f)
		return false;
	fread(&magic, sizeof(int), 1, f);
	if (magic == VectorArrayMagic)
	{
		fread(&dx, sizeof(double), 1, f);
		fread(&dy, sizeof(double), 1, f);
		fread(&nx, sizeof(int), 1, f);
		fread(&ny, sizeof(int), 1, f);
		if (vx)
			delete vx;
		vx = new double[2 * nx * ny];
		fread(vx, sizeof(double), 2 * nx * ny, f);
		fclose(f);
		return true;

	}
	if(magic == VectorArrayMagic2)
	{
		fclose(f);
		return true;
	}
	fclose(f);
	return false;
}
