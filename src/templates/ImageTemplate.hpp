/*
 * imageTemplate.h
 *
 *  Created on: Jul 10, 2012
 *      Author: Konstantin Thierbach
 */

#ifndef IMAGETEMPLATE_H_
#define IMAGETEMPLATE_H_

#include <stdio.h>
#include <stdlib.h>

#include "tbb/tbb.h"

#include "../utilities/Utilities.hpp"
#include "Point2D.hpp"
#include "VectorArray2D.hpp"

#define Addr(i,j,c,nx,ny,channelNo) ((c)+(channelNo)*((i)+(nx)*(j)))
#define OddQ(i)  ((i&1)==1)
#define PixelAddr(i,j,nx,ny,channelNo) ((channelNo)*((i)+(nx)*(j)))

#define EMPTY_PIXEL 0
#define FILLED_PIXEL 15
#define MAX_CHANNELS 4
#define RETINEX_HIGH 2
#define RETINEX_LOW 1
#define RETINEX_UNIFORM 0
#define TINY (0.1e-12)

template<class T> class Image
{
public:
	int nx, ny, channelNo;
	int pixelSize, byteSize;
	T *bm;

	/*where to put this code?*/
	const static int SouthWest = 1, South = 2, SouthEast = 4, East = 8, NorthEast = 16, North = 32, NorthWest = 64, West = 128;

	Image() :
			nx(0), ny(0), channelNo(0), pixelSize(0), byteSize(0), bm(NULL)
	{
	}
	Image(T* data, int nX, int nY, int ch = 1) :
			nx(nX), ny(nY), channelNo(ch)
	{
		pixelSize = nx * ny * channelNo;
		byteSize = nx * ny * channelNo * sizeof(T);
		bm = new T[pixelSize];
		memcpy(bm, data, byteSize);
	}
	Image(int nX, int nY, int ch = 1) :
			nx(nX), ny(nY), channelNo(ch)
	{

		pixelSize = nx * ny * channelNo;
		byteSize = nx * ny * channelNo * sizeof(T);
		bm = new T[pixelSize];
		if (bm)
			memset(bm, '\0', byteSize);
		else
			fprintf(stderr, " panik can't allocate memory for the image ...\n");
	}
	Image(Image const &src) :
			nx(src.nx), ny(src.ny), channelNo(src.channelNo), pixelSize(src.pixelSize), byteSize(src.byteSize)
	{
		bm = new T[pixelSize];
		memcpy(bm, src.bm, byteSize);
	}
	~Image(void)
	{
		if (bm)
			delete[] bm;
	}

	inline bool samePixelQ(T *pixel1, T *pixel2, int channels)
	{
		bool same;

		same = pixel1[0] == pixel2[0];
		for (int k = 1; k < channels && same; k++)
			same = same && pixel1[k] == pixel2[k];
		return same;
	}
	inline bool sameRowQ(int y, T *tpixel, int channels)
	{
		bool same;

		same = samePixelQ(get(0, y), tpixel, channels);
		for (int i = 1; i < nx && same; i++)
			same = same && samePixelQ(get(i, y), tpixel, channels);
		return same;
	}
	inline void deleteRow(int y)
	{
		int addr1, addr2;

		if (ny - 1 == y)
		{
			(ny)--;
			return;
		}
		addr1 = PixelAddr(0, y, nx, ny, channelNo);
		addr2 = addr1 + channelNo * nx;
		memmove(bm + addr1, bm + addr2, sizeof(T) * channelNo * nx * (ny - y - 1));
		(ny)--;
	}
	inline void clipImageTopBottom(void)
	{
		T zero[MAX_CHANNELS];

		for (int k = 0; k < channelNo; k++)
			zero[k] = 0;
		while (sameRowQ(ny - 1, zero, channelNo) && ny > 0)
			deleteRow(ny - 1);
		while (sameRowQ(0, zero, channelNo) && ny > 0)
			deleteRow(0);
	}
	inline Image<T>* downsample(void)
	{
		const double tiny = 0.1e-14;
		const double kernel[3][3] =
		{
		{ 0.0625, 0.1250, 0.0625 },
		{ 0.1250, 0.2500, 0.1250 },
		{ 0.0625, 0.1250, 0.0625 } };
		const int kwidth = 3;
		const int kheight = 3;
		const int koriginX = 1;
		const int koriginY = 1;

		Image<T> *down = NULL;

		down = new Image<T>(nx / 2, ny / 2, channelNo);
		for (int x = 0; x < down->nx; x++)
			for (int y = 0; y < down->ny; y++)
				for (int ch = 0; ch < down->channelNo; ch++)
				{
					double val, weight, w;
					val = weight = 0.0;

					for (int i = 0; i < kwidth; i++)
						for (int j = 0; j < kheight; j++)
						{
							int nnx = 2 * x + i - koriginX;
							int nny = 2 * y + j - koriginY;
							if (nnx >= nx || nny >= ny)
								continue;
							if (nnx < 0 || nny < 0)
								continue;
							w = kernel[i][j];
							val += w * get(nnx, nny, ch);
							weight += w;
						}

					if (weight < tiny)
						continue;
					down->setc(x, y, ch, val / weight);
				}
		return down;
	}
	inline Image<T>* upsample(void)
	{
		const double tiny = 0.1e-14;
		const double kernel[3][3] =
		{
		{ 0.0625, 0.1250, 0.0625 },
		{ 0.1250, 0.2500, 0.1250 },
		{ 0.0625, 0.1250, 0.0625 } };
		const int kwidth = 3;
		const int kheight = 3;
		const int koriginX = 1;
		const int koriginY = 1;

		Image<T> *up = NULL;

		up = new Image<T>(2 * nx, 2 * ny, channelNo);
		for (int x = 0; x < up->nx; x++)
			for (int y = 0; y < up->ny; y++)
				for (int ch = 0; ch < up->channelNo; ch++)
				{
					double val, weight, w;
					val = weight = 0.0;

					for (int i = 0; i < kwidth; i++)
						for (int j = 0; j < kheight; j++)
						{
							int nnx = x / 2 + i - koriginX;
							int nny = y / 2 + j - koriginY;
							if (nnx >= nx || nny >= ny)
								continue;
							if (nnx < 0 || nny < 0)
								continue;
							w = kernel[i][j];
							val += w * get(nnx, nny, ch);
							weight += w;
						}

					if (weight < tiny)
						continue;
					up->setc(x, y, ch, val / weight);
				}
		return up;
	}
	inline Image<T> *gauss(double sigma)
	{
		int i, j, ch;
		double q, q2, q3;
		double c[4], C;
		double*w, *out;
		Image<T> *outimg;

		outimg = new Image<T>(nx, ny, channelNo);

		if (sigma < 0.5)
			q = 0.1147705018520355224609375;
		if (0.5 <= sigma && sigma <= 2.5)
			q = 3.97156 - 4.14554 * sqrt(1.0 - 0.26891 * sigma);
		else if (2.5 < sigma)
			q = 0.98711 * sigma - 0.9633;
		q2 = q * q;
		q3 = q * q2;
		c[0] = (1.57825 + (2.44413 * q) + (1.4281 * q2) + (0.422205 * q3));
		c[1] = ((2.44413 * q) + (2.85619 * q2) + (1.26661 * q3));
		c[2] = (-((1.4281 * q2) + (1.26661 * q3)));
		c[3] = ((0.422205 * q3));
		C = 1.0 - ((c[1] + c[2] + c[3]) / c[0]);

		w = new double[nx];
		out = new double[nx];
		for (ch = 0; ch < channelNo; ch++)
		{

			for (j = 0; j < ny; j++)
			{

				w[0] = w[1] = w[2] = get(0, j, ch);
				for (i = 3; i < nx; i++)
					w[i] = C * get(i, j, ch) + (c[1] * w[i - 1] + c[2] * w[i - 2] + c[3] * w[i - 3]) / c[0];

				out[nx - 1] = out[nx - 2] = out[nx - 3] = w[nx - 1];
				for (i = nx - 4; i > -1; i--)
					out[i] = C * w[i] + (c[1] * out[i + 1] + c[2] * out[i + 2] + c[3] * out[i + 3]) / c[0];

				for (i = 0; i < nx; i++)
				{
					outimg->setc(i, j, ch, (T) out[i]);
				}
			}
		}
		delete[] w;
		delete[] out;

		w = new double[ny];
		out = new double[ny];
		for (ch = 0; ch < channelNo; ch++)
		{

			for (i = 0; i < nx; i++)
			{

				w[0] = w[1] = w[2] = outimg->get(i, 0, ch);
				for (j = 3; j < ny; j++)
					w[j] = C * outimg->get(i, j, ch) + (c[1] * w[j - 1] + c[2] * w[j - 2] + c[3] * w[j - 3]) / c[0];

				out[ny - 1] = out[ny - 2] = out[ny - 3] = w[ny - 1];
				for (j = ny - 4; j > -1; j--)
					out[j] = C * w[j] + (c[1] * out[j + 1] + c[2] * out[j + 2] + c[3] * out[j + 3]) / c[0];

				for (j = 0; j < ny; j++)
					outimg->setc(i, j, ch, (T) out[j]);
			}
		}
		delete[] w;
		delete[] out;

		return outimg;
	}
	inline T interpolate(double x, double y, int c) const
	{
		int l, m, i, j;
		const double dx = 1.0;
		const double dy = 1.0;
		double u[4], v[4];
		double r[4] =
		{ 0., 0., 0., 0. };
		double val = 0.0;

		l = (int) floor(x);
		m = (int) floor(y);
		x -= l;
		y -= m;
#if 0
		if(l<0)l= 0;
		else if(l>=nx-2)l= nx-2;
		if(m<0)m= 0;
		else if(m>=ny-2)m= ny-2;
#else
		while (l < 0)
			l += nx;
		while (l >= nx - 1)
			l -= nx;
		while (m < 0)
			m += ny;
		while (m >= ny - 1)
			m -= ny;
#endif
#define G(a,b) (get(l+a,m+b,c))
#define GET(a,b) (get(a,b,c))

		if (l > 0 && l < nx - 1 && m > 0 && m < ny - 1)
		{
			int mu, lambda;

			SplineCoefficients(u, x, l, nx);
			SplineCoefficients(v, y, m, ny);

			for (j = 0; j < 4; j++)
			{
				r[j] = 0.0;
				mu = m + j - 1;
				if (mu < 0 || mu >= ny)
					continue;
				for (i = 0; i < 4; i++)
				{
					lambda = l + i - 1;
					if (lambda < 0 || lambda >= nx)
						continue;
					r[j] += u[i] * (GET(lambda,mu));
				}
				val += v[j] * r[j];
			}
		}
		else if (l > 0 && l < nx - 1 && m > 0 && m < ny - 1)
		{
			val = ((1.0 - y) * ((1.0 - x) * G(0,0) + x * G(1,0)) + y * ((1.0 - x) * G(0,1) + x * G(1,1)));
		}

#undef G
#undef GET
		return (T) val;
	}
	inline T bilinear(double x, double y, int c) const
	{
		int l, m;
		double val = 0.0;

		l = (int) floor(x);
		m = (int) floor(y);
		x -= l;
		y -= m;
#if 0
		if(l<0)l= 0;
		else if(l>=nx-2)l= nx-2;
		if(m<0)m= 0;
		else if(m>=ny-2)m= ny-2;
#else
		while (l < 0)
			l += nx;
		while (l >= nx - 1)
			l -= nx;
		while (m < 0)
			m += ny;
		while (m >= ny - 1)
			m -= ny;
#endif
		val = ((1.0 - x) * get(l, m, c) + x * get(l + 1, m, c)) * (1.0 - y)
				+ ((1.0 - x) * get(l, m + 1, c) + x * get(l + 1, m + 1, c)) * y;
		return (T) val;
	}
	void bilinear(T*val, double x, double y) const
	{
		int l, m;

		l = (int) floor(x);
		m = (int) floor(y);
		x -= l;
		y -= m;
#if 0
		if(l<0)l= 0;
		else if(l>=nx-2)l= nx-2;
		if(m<0)m= 0;
		else if(m>=ny-2)m= ny-2;
#else
		while (l < 0)
			l += nx;
		while (l >= nx)
			l -= nx;
		while (m < 0)
			m += ny;
		while (m >= ny)
			m -= ny;
#endif
		for (int c = 0; c < channelNo; c++)
		{
			val[c] = ((1.0 - x) * get(l, m, c) + x * get(l + 1, m, c)) * (1.0 - y)
					+ ((1.0 - x) * get(l, m + 1, c) + x * get(l + 1, m + 1, c)) * y;
		}
	}
	inline void wrap(Image<T> *w, VectorArray2D*u)
	{

#pragma omp parallel for
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				Vector2D delta = u->get(i, j);
				double xp = i - delta.x;
				double yp = j - delta.y;
				if (xp < 0 || xp > nx - 1)
					continue;
				if (yp < 0 || yp > ny - 1)
					continue;
				for (int ch = 0; ch < channelNo; ch++)
				{
					T val = interpolate(i - delta.x, j - delta.y, ch);
#pragma omp critical
					w->setc(i, j, ch, val);
				}
			}
	}
	inline void wrap(VectorArray2D*u)
	{
		Image<T> *w;

		w = new Image<T>(nx, ny, channelNo);
		this->wrap(w, u);
		memcpy(this->bm, w->bm, byteSize);
		delete w;
	}
	void rgbToLuv(double maxrgb)
	{
		int dataNo = nx * ny;
		if (3 == channelNo || 4 == channelNo)
		{
			for (int what = 0; what < dataNo; what++)
			{
				double rgb[3];
				for (int ch = 0; ch < 3; ch++)
					rgb[ch] = bm[what * channelNo + ch];
				RGBToLuvColor(rgb, rgb, maxrgb);
				for (int ch = 0; ch < 3; ch++)
					bm[what * channelNo + ch] = (T) rgb[ch];
			}
		}
	}
	inline void RGBToLuvColor(double*luv, const double*rgb, const double maxValue)
	{
		double xyz[3];

		RGBToXYZColor(xyz, rgb, maxValue);
		XYZToLuvColor(luv, xyz);
	}
	inline void RGBToLuvColor(int*_luv, const int*_rgb, const double maxValue)
	{
		double xyz[3];
		double rgb[3], luv[3];

		rgb[0] = _rgb[0];
		rgb[1] = _rgb[1];
		rgb[2] = _rgb[2];
		RGBToXYZColor(xyz, rgb, maxValue);
		XYZToLuvColor(luv, xyz);
		_luv[0] = luv[0];
		_luv[1] = luv[1];
		_luv[2] = luv[2];
	}
	void luvToRGB(double maxrgb)
	{
		int dataNo = nx * ny;
		if (3 == channelNo || 4 == channelNo)
		{
			for (int what = 0; what < dataNo; what++)
				LuvToRGBColor(bm + what * channelNo, bm + what * channelNo, maxrgb);
		}
	}
	inline void LuvToRGBColor(double*rgb, const double*luv, const double maxValue)
	{
		double xyz[3];

		LuvToXYZColor(xyz, luv);
		XYZToRGBColor(rgb, xyz, maxValue);
	}
	inline void LuvToRGBColor(int*_rgb, const int*_luv, const double maxValue)
	{
		double xyz[3];
		double luv[3], rgb[3];

		luv[0] = _luv[0];
		luv[1] = _luv[1];
		luv[2] = _luv[2];
		LuvToXYZColor(xyz, luv);
		XYZToRGBColor(rgb, xyz, maxValue);
		_rgb[0] = rgb[0];
		_rgb[1] = rgb[1];
		_rgb[2] = rgb[2];
	}
	inline void RGBToXYZColor(double*xyz, const double*rgb, double maxValue = 255.0)
	{

		double varRGB[3];
		for (int i = 0; i < 3; i++)
		{
			varRGB[i] = rgb[i] / maxValue;
			varRGB[i] =
					varRGB[i] > 0.04045 ?
							pow((varRGB[i] + 0.055) / 1.055, 2.4) :
							varRGB[i] / 12.92;
			varRGB[i] *= 100.0;
		}

		xyz[0] = varRGB[0] * 0.4124 + varRGB[1] * 0.3576 + varRGB[2] * 0.1805;
		xyz[1] = varRGB[0] * 0.2126 + varRGB[1] * 0.7152 + varRGB[2] * 0.0722;
		xyz[2] = varRGB[0] * 0.0193 + varRGB[1] * 0.1192 + varRGB[2] * 0.9505;
	}
	inline void XYZToRGBColor(double*rgb, const double*xyz, double maxValue = 255.0)
	{
		double varXYZ[3], varRGB[3];

		for (int i = 0; i < 3; i++)
			varXYZ[i] = xyz[i] / 100.0;

		varRGB[0] = varXYZ[0] * 3.2406 + varXYZ[1] * -1.5372 + varXYZ[2] * -0.4986;
		varRGB[1] = varXYZ[0] * -0.9689 + varXYZ[1] * 1.8758 + varXYZ[2] * 0.0415;
		varRGB[2] = varXYZ[0] * 0.0557 + varXYZ[1] * -0.2040 + varXYZ[2] * 1.0570;

		for (int i = 0; i < 3; i++)
		{
			rgb[i] =
					varRGB[i] > 0.0031308 ?
							1.055 * pow(varRGB[i], (1 / 2.4)) - 0.055 :
							12.92 * varRGB[i];
			rgb[i] *= maxValue;
		}

	}
	inline void XYZToLuvColor(double*luv, const double*xyz)
	{
		const double refXYZ[] =
		{ 95.047, 100.000, 108.883 };

		const double refUV[2] =
		{ 0.19784, 0.468336 };

		double varUVY[3], tmp;

		tmp = xyz[0] + 15. * xyz[1] + 3. * xyz[2];
		varUVY[0] = tmp > 0.0 ? 4.0 * xyz[0] / tmp : 0.0;
		varUVY[1] = tmp > 0.0 ? 9.0 * xyz[1] / tmp : 0.0;
		varUVY[2] = 0.01 * xyz[1];

		varUVY[2] =
				varUVY[2] > 0.008856 ?
						varUVY[2] = pow(varUVY[2], 1. / 3.) :
						7.787 * varUVY[2] + 16. / 116.;

		luv[0] = 116. * varUVY[2] - 16.0;
		luv[1] = 13. * luv[0] * (varUVY[0] - refUV[0]);
		luv[2] = 13. * luv[0] * (varUVY[1] - refUV[1]);
	}
	inline void LuvToXYZColor(double*xyz, const double*luv)
	{
		const double refXYZ[] =
		{ 95.047, 100.000, 108.883 };

		const double refUV[2] =
		{ 0.19784, 0.468336 };
		double varY, varU, varV;

		varY = (luv[0] + 16.0) / 116.0;
		varY = varY * varY * varY > 0.008856 ?
				varY * varY * varY : (varY - 16.0 / 116.0) / 7.787;
		if (luv[0] > 0.0)
		{
			varU = luv[1] / (13.0 * luv[0]) + refUV[0];
			varV = luv[2] / (13.0 * luv[0]) + refUV[1];
		}
		else
		{
			varU = refUV[0];
			varV = refUV[1];
		}

		xyz[1] = varY * 100.0;
		xyz[0] = -9.0 * xyz[1] * varU / ((varU - 4) * varV - varU * varV);
		xyz[2] = (9.0 * xyz[1] - 15.0 * varV * xyz[1] - varV * xyz[0])
				/ (3.0 * varV);
	}
	inline void setToBlack(void)
	{
		memset(bm, '\0', byteSize);
	}
	inline long address(int i, int j, int ch = 0) const
	{
		return Addr(i, j, ch, nx, ny, channelNo);
	}
	inline long address2d(int i, int j) const
	{
		return (i) + (nx) * (j);
	}
	inline void set(int i, int j, T*pixel)
	{
		T*img;

		img = bm + PixelAddr(i, j, nx, ny, channelNo);
		for (int c = 0; c < channelNo; c++)
			*(img++) = *(pixel++);
	}
	inline void setc(int i, int j, int c, T cval)
	{
		bm[Addr(i, j, c, nx, ny, channelNo)] = cval;
	}
	inline void incc(int i, int j, int c, T cval)
	{
		bm[Addr(i, j, c, nx, ny, channelNo)] += cval;
	}
	inline void zero(void)
	{
		memset(bm, '\0', nx * ny * channelNo * sizeof(T));
	}
	inline T* get(int i, int j) const
	{
		return bm + PixelAddr(i, j, nx, ny, channelNo);
	}
	inline T get(int i, int j, int c) const
	{
		return bm[Addr(i, j, c, nx, ny, channelNo)];
	}
	inline int blackQ(int i, int j) const
	{
		T*p;
		int c, isblack;

		p = bm + Addr(i, j, 0, nx, ny, channelNo);
		isblack = abs(*p) < TINY;
		for (c = 1; isblack && c < channelNo; c++)
			isblack = isblack && abs(*(++p)) < TINY;
		return isblack;
	}
	inline int neighborMap(int i, int j, int c = 0) const
	{
		int mapCode = 0;

		if (j > 0)
		{
			if (i > 0 && get(i - 1, j - 1, c))
				mapCode += SouthWest;
			if (get(i, j - 1, c))
				mapCode += South;
			if (i < nx - 1 && get(i + 1, j - 1, c))
				mapCode += SouthEast;
		}
		if (i > 0 && get(i - 1, j, c))
			mapCode += West;
		if (i < nx - 1 && get(i + 1, j, c))
			mapCode += East;
		if (j < ny - 1)
		{
			if (i > 0 && get(i - 1, j + 1, c))
				mapCode += NorthWest;
			if (get(i, j + 1, c))
				mapCode += North;
			if (i < nx - 1 && get(i + 1, j + 1, c))
				mapCode += NorthEast;
		}
		return mapCode;
	}
	inline void copy(const Image<T> *source)
	{
		memcpy(bm, source->bm, nx * ny * channelNo * sizeof(T));
	}
	inline void transpose(void)
	{
		T *newbm;
		int i;
		newbm = new T[pixelSize];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				T*pixel = get(i, j);
				T*newpixel = newbm + (channelNo * (j + ny * i));
				for (int k = 0; k < channelNo; k++)
					*(newpixel++) = *(pixel++);
			}
		delete[] bm;
		i = nx;
		nx = ny;
		ny = i;
		bm = newbm;
	}
	inline T max(int channel = 0) const
	{
		T mb;
		T*img;

		img = bm + channel;
		mb = *img;
		for (int i = 1; i < nx * ny; i++)
		{
			img += channelNo;
			if (mb < *img)
				mb = *img;
		}
		return mb;
	}
	inline T min(int channel) const
	{
		T mb;
		T*img;

		img = bm + channel;
		mb = *img;
		for (int i = 1; i < nx * ny; i++)
		{
			img += channelNo;
			if (mb > *img)
				mb = *img;
		}
		return mb;
	}
	inline void clipRange(int channel, T minc, T maxc)
	{
		T* img = bm + channel;
		for (int i = 0; i < nx * ny; i++)
		{
			if (*img < minc)
				*img = minc;
			if (*img > maxc)
				*img = maxc;
			img += channelNo;
		}
	}
	void invert(T maxval = 1)
	{
		int k;

		for (k = 0; k < pixelSize; k++)
		bm[k] = maxval - bm[k];
	}
	inline void normalizeRange(T maxvalue)
	{
		T minBm, maxBm;
		T*img;
		int ch, i;
		double scale;

		for (ch = 0; ch < channelNo; ch++)
		{
			minBm = min(ch);
			maxBm = max(ch);
			img = bm + ch;
			if (minBm == maxBm)
			continue;
			scale = maxvalue / (maxBm - minBm);
			for (i = 0; i < nx * ny; i++)
			{
				*img = scale * (*img - minBm);
				img += channelNo;
			}
		}
	}
	inline T dx(int i, int j, int ch = 0) const
	{
		if (0 == i)
			return (-3 * get(i, j, ch) + 4 * get(1 + i, j, ch) - get(2 + i, j, ch)) / 2;
		if (nx - 1 == i)
			return (get(-2 + i, j, ch) - 4 * get(-1 + i, j, ch) + 3 * get(i, j, ch)) / 2;
		return (get(i + 1, j, ch) - get(i - 1, j, ch)) / 2;
	}
	inline T dy(int i, int j, int ch = 0) const
	{
		if (0 == j)
			return (-3 * get(i, j, ch) + 4 * get(i, 1 + j, ch) - get(i, 2 + j, ch)) / 2;
		if (ny - 1 == j)
			return (get(i, -2 + j, ch) - 4 * get(i, -1 + j, ch) + 3 * get(i, j, ch)) / 2;
		return (get(i, j + 1, ch) - get(i, j - 1, ch)) / 2;
	}
	inline T d2x(int i, int j, int ch = 0) const
	{
		if (0 == i)
			return 2 * get(i, j, ch) - 5 * get(1 + i, j, ch) + 4 * get(2 + i, j, ch) - get(3 + i, j, ch);
		if (nx - 1 == i)
			return -get(-3 + i, j, ch) + 4 * get(-2 + i, j, ch) - 5 * get(-1 + i, j, ch) + 2 * get(i, j, ch);
		return get(i - 1, j, ch) - 2 * get(i, j, ch) + get(i + 1, j, ch);
	}
	inline T d2y(int i, int j, int ch = 0) const
	{
		if (0 == j)
			return 2 * get(i, j, ch) - 5 * get(i, 1 + j, ch) + 4 * get(i, 2 + j, ch) - get(i, 3 + j, ch);
		if (ny - 1 == j)
			return -get(i, -3 + j, ch) + 4 * get(i, -2 + j, ch) - 5 * get(i, -1 + j, ch) + 2 * get(i, j, ch);
		return get(i, j - 1, ch) - 2 * get(i, j, ch) + get(i, j + 1, ch);
	}
	inline T dxy(int i, int j, int ch = 0) const
	{
		if (0 == i)
		{
			if (0 == j)
				return (9 * get(i, j, ch) - 12 * (get(i, 1 + j, ch) + get(1 + i, j, ch)) + 16 * get(1 + i, 1 + j, ch)
						+ 3 * (get(i, 2 + j, ch) + get(2 + i, j, ch))
						- 4 * (get(1 + i, 2 + j, ch) + get(2 + i, 1 + j, ch)) + get(2 + i, 2 + j, ch)) / 4;
			if (j == ny - 1)
				return (-9 * get(i, j, ch) - 16 * get(1 + i, -1 + j, ch) + 12 * (get(i, -1 + j, ch) + get(1 + i, j, ch))
						- get(2 + i, -2 + j, ch) + 4 * (get(1 + i, -2 + j, ch) + get(2 + i, -1 + j, ch))
						- 3 * (get(i, -2 + j, ch) + get(2 + i, j, ch))) / 4;
			return (3 * (get(i, -1 + j, ch) - get(i, 1 + j, ch)) + 4 * (get(1 + i, 1 + j, ch) - get(1 + i, -1 + j, ch))
					+ get(2 + i, -1 + j, ch) - get(2 + i, 1 + j, ch)) / 4;

		}
		if (nx - 1 == i)
		{
			if (0 == j)
				return (-get(-2 + i, 2 + j, ch) - 16 * get(-1 + i, 1 + j, ch)
						+ 4 * (get(-2 + i, 1 + j, ch) + get(-1 + i, 2 + j, ch)) - 9 * get(i, j, ch)
						+ 12 * (get(-1 + i, j, ch) + get(i, 1 + j, ch)) - 3 * (get(-2 + i, j, ch) + get(i, 2 + j, ch)))
						/ 4;
			if (ny - 1 == j)
				return (get(-2 + i, -2 + j, ch) - 4 * (get(-2 + i, -1 + j, ch) + get(-1 + i, -2 + j, ch))
						+ 16 * get(-1 + i, -1 + j, ch) + 3 * (get(-2 + i, j, ch) + get(i, -2 + j, ch))
						- 12 * (get(-1 + i, j, ch) + get(i, -1 + j, ch)) + 9 * get(i, j, ch)) / 4;
			return (-get(-2 + i, -1 + j, ch) + get(-2 + i, 1 + j, ch)
					+ 4 * (get(-1 + i, -1 + j, ch) - get(-1 + i, 1 + j, ch))
					+ 3 * (get(i, 1 + j, ch) - get(i, -1 + j, ch))) / 4;
		}
		if (0 == j)
			return (3 * (get(-1 + i, j, ch) - get(1 + i, j, ch)) + 4 * (get(1 + i, 1 + j, ch) - get(-1 + i, 1 + j, ch))
					+ get(-1 + i, 2 + j, ch) - get(1 + i, 2 + j, ch)) / 4;
		if (ny - 1 == j)
			return (-get(-1 + i, -2 + j, ch) + get(1 + i, -2 + j, ch)
					+ 4 * (get(-1 + i, -1 + j, ch) - get(1 + i, -1 + j, ch))
					+ 3 * (get(1 + i, j, ch) - get(-1 + i, j, ch))) / 4;
		return (get(i + 1, j + 1, ch) - get(i + 1, j - 1, ch) - get(i - 1, j + 1, ch) + get(i - 1, j - 1, ch)) / 4;
	}
	inline T laplace(int i, int j, int ch) const
	{
		return d2x(i, j, ch) + d2y(i, j, ch);
	}
	inline T periodicLaplace(int i, int j, int ch) const
	{
		int im2, im1, ip1, ip2, jm2, jm1, jp1, jp2;
		T l, v;
		im2 = modulo(i - 2, nx);
		im1 = modulo(i - 1, nx);
		ip1 = modulo(i + 1, nx);
		ip2 = modulo(i + 2, nx);

		jm2 = modulo(j - 2, ny);
		jm1 = modulo(j - 1, ny);
		jp1 = modulo(j + 1, ny);
		jp2 = modulo(j + 2, ny);
		v = get(i, j, ch);
		l = 16 * (get(ip1, j, ch) + get(im1, j, ch)) - 30 * v - get(ip2, j, ch) - get(im2, j, ch);
		l += 16 * (get(i, jp1, ch) + get(i, jm1, ch)) - 30 * v - get(i, jp2, ch) - get(i, jm2, ch);
		return l / 12;
	}
	inline T isotropeLaplace(int i, int j, int ch) const
	{
		T lapl;

		if (i < 1 || j < i || i > nx - 2 || j > ny - 2)
			return laplace(i, j, ch);
		lapl = 4 * (get(i - 1, j, ch) + get(i + 1, j, ch) + get(i, j - 1, ch) + get(i, j + 1, ch))
				+ get(i - 1, j - 1, ch) + get(i + 1, j + 1, ch) + get(i + 1, j - 1, ch) + get(i - 1, j + 1, ch)
				- 20 * get(i, j, ch);
		return lapl / 6;
	}
	inline T DeltaMinusX(int i, int j, int ch = 0) const
	{
		if (i > 0)
			return get(i, j, ch) - get(i - 1, j, ch);
		return get(i + 1, j, ch) - get(i, j, ch);
	}
	inline T DeltaMinusY(int i, int j, int ch = 0) const
	{
		if (j > 0)
			return get(i, j, ch) - get(i, j - 1, ch);
		return get(i, j + 1, ch) - get(i, j, ch);
	}
	inline T DeltaPlusX(int i, int j, int ch = 0) const
	{
		if (i < nx - 1)
			return get(i + 1, j, ch) - get(i, j, ch);
		return get(i, j, ch) - get(i - 1, j, ch);
	}
	inline T DeltaPlusY(int i, int j, int ch = 0) const
	{
		if (j < ny - 1)
			return get(i, j + 1, ch) - get(i, j, ch);
		return get(i, j, ch) - get(i, j - 1, ch);
	}
	inline T enoDx(int i, int j, int ch) const
	{
		return minMod(DeltaPlusX(i, j, ch), DeltaMinusX(i, j, ch));
	}
	inline T enoDy(int i, int j, int ch) const
	{
		return minMod(DeltaPlusY(i, j, ch), DeltaMinusY(i, j, ch));
	}
	inline T enoAbsGrad(int i, int j, int ch) const
	{
		T ddx, ddy;

		ddx = enoDx(i, j, ch);
		ddy = enoDy(i, j, ch);
		return sqrt(ddx * ddx + ddy * ddy);
	}
	inline void hessianEigenvalues(double*l1, double*l2, int i, int j, int ch) const
	{
		T d2Idx2, d2Idxdy, d2Idy2;
		double dis, ll1, ll2;

		d2Idx2 = d2x(i, j, ch);
		d2Idy2 = d2y(i, j, ch);
		d2Idxdy = dxy(i, j, ch);
		dis = d2Idx2 - d2Idy2;
		dis = sqrt(dis * dis + 4.0 * d2Idxdy * d2Idxdy);
		ll1 = 0.5 * (d2Idx2 + d2Idy2 - dis);
		ll2 = 0.5 * (d2Idx2 + d2Idy2 + dis);
		if (ll1 > ll2)
		{
			*l1 = ll1;
			*l2 = ll2;
		}
		else
		{
			*l1 = ll2;
			*l2 = ll1;
		}
	}
	double stoppingPeronaMalik(double nu, double sqrK, int i, int j, int ch = 0) const
	{
		T gradX, gradY;
		double absGrad2;

		gradX = dx(i, j, ch);
		gradY = dy(i, j, ch);
		absGrad2 = gradX * gradX + gradY * gradY;
		return nu / (1.0 + absGrad2 / sqrK);
	}
	inline bool boundaryPixelQ(int i, int j, int c, T threshold) const
	{
		bool bq = get(i, j, c) >= threshold;
		bool ntest = false;

		if (!bq)
			return false;
		if (i == 0 || i == nx - 1)
			return true;
		if (j == 0 || j == ny - 1)
			return true;
		ntest |= get(i - 1, j, c) < threshold;
		ntest |= get(i + 1, j, c) < threshold;
		ntest |= get(i, j - 1, c) < threshold;
		ntest |= get(i, j + 1, c) < threshold;
		ntest |= get(i - 1, j - 1, c) < threshold;
		ntest |= get(i + 1, j + 1, c) < threshold;
		ntest |= get(i - 1, j + 1, c) < threshold;
		ntest |= get(i + 1, j - 1, c) < threshold;
		return ntest;
	}
	inline void meanAndStandarddeviation(double*mean, double*stddev, int w, int i, int j, int c) const
	{
		double mu = 0.0, sigma = 0.0, tmp;
		T v;
		int ii, jj, k = 0;

		for (ii = -w; ii <= w; ii++)
		{
			int px = i + ii;
			if (px < 0 || px >= nx)
				continue;
			for (jj = -w; jj <= w; jj++)
			{
				int py = j + jj;
				if (py < 0 || py >= ny)
					continue;
				mu += get(px, py, c);
				k++;
			}
		}
		mu /= k;
		*mean = mu;
		k = 0;
		for (ii = -w; ii <= w; ii++)
		{
			int px = i + ii;
			if (px < 0 || px >= nx)
				continue;
			for (jj = -w; jj <= w; jj++)
			{
				int py = j + jj;
				if (py < 0 || py >= ny)
					continue;
				v = get(px, py, c);
				tmp = v - *mean;
				sigma += tmp * tmp;
				k++;
			}
		}
		sigma /= (k - 1);
		*stddev = sqrt(sigma);
	}
	double skewness(int w, int i, int j, int c) const
	{
		double m1, m2, m3, mu2, mu3;
		int l, k;
		T*elem;

		{
			int ii, jj, px, py;

			ii = 2 * w + 1;
			elem = (T*) calloc(ii * ii, sizeof(T));
			k = 0;
			for (ii = -w; ii <= w; ii++)
			{
				px = i + ii;
				if (px < 0 || px >= nx)
					continue;
				for (jj = -w; jj <= w; jj++)
				{
					py = j + jj;
					if (py < 0 || py >= ny)
						continue;
					elem[k++] = get(px, py, c);
				}
			}
		}

		m1 = m2 = m3 = 0.0;
		for (l = 0; l < k; l++)
			m1 += elem[l];
		m1 /= k;
		for (l = 0; l < k; l++)
		{
			double tmp = (elem[l] - m1);
			m2 += tmp * tmp;
			m3 += tmp * tmp * tmp;
		}
		m2 /= k;
		m3 /= k;
		mu2 = k * m2 / (k - 1);
		mu3 = k * k * m2 / ((k - 1) * (k - 2));
		free(elem);
		return (fabs(mu2) < TINY) ? mu3 : mu3 / (mu2 * sqrt(mu2));
	}
	double excess(int w, int i, int j, int c) const
	{
		double m1, m2, m4, mu2, mu4;
		int l, k;
		T*elem;

		{
			int ii, jj, px, py;

			ii = 2 * w + 1;
			elem = (T*) calloc(ii * ii, sizeof(T));
			k = 0;
			for (ii = -w; ii <= w; ii++)
			{
				px = i + ii;
				if (px < 0 || px >= nx)
					continue;
				for (jj = -w; jj <= w; jj++)
				{
					py = j + jj;
					if (py < 0 || py >= ny)
						continue;
					elem[k++] = get(px, py, c);
				}
			}
		}

		m1 = m2 = m4 = 0.0;
		for (l = 0; l < k; l++)
			m1 += elem[l];
		m1 /= k;
		for (l = 0; l < k; l++)
		{
			double tmp = (elem[l] - m1);
			double tmp2 = tmp * tmp;

			m2 += tmp2;
			m4 += tmp2 * tmp2;
		}
		m2 /= k;
		m4 /= k;
		mu2 = k * m2 / (k - 1);
		mu4 = (m4 * k * k * k - 3 * (k - 1) * (2 * k - 3) * mu2 * mu2) / ((k - 1) * (3 + (k - 3) * k));
		free(elem);
		return fabs(mu2) < TINY ? mu4 : mu4 / (mu2 * mu2);
	}
	inline double entropy(int i, int j, int c, int w, double minVal, double maxVal) const
	{
		T*elem;
		int*histo = NULL;
		double delta, sum, q, entrop = 0.0;
		int l, k, bin, binNo;

		{
			int ii, jj, px, py;

			ii = 2 * w + 1;
			elem = (T*) calloc(ii * ii, sizeof(T));
			k = 0;
			for (ii = -w; ii <= w; ii++)
			{
				px = i + ii;
				if (px < 0 || px >= nx)
					continue;
				for (jj = -w; jj <= w; jj++)
				{
					py = j + jj;
					if (py < 0 || py >= ny)
						continue;
					elem[k++] = get(px, py, c);
				}
			}
		}

		binNo = k / 2;
		if (binNo >= 2)
		{
			delta = (maxVal - minVal) / binNo;
			histo = new int[binNo];
			for (bin = 0; bin < binNo; bin++)
				histo[bin] = 0;
			for (l = 0; l < k; l++)
			{
				q = elem[l] - minVal;
				for (bin = 0; bin < binNo; bin++)
				{
					if (q <= (bin + 1) * delta)
					{
						histo[bin] += 1;
						break;
					}
				}
			}
			for (sum = 0.0, bin = 0; bin < binNo; bin++)
				sum += histo[bin];
			sum = 1.0 / sum;
			for (bin = 0; bin < binNo; bin++)
			{
				if (!histo[bin])
					continue;
				q = histo[bin] * sum;
				entrop -= q * log(q);
			}
		}
		free(histo);
		free(elem);
		return entrop;
	}
	void minMedianMax(T*minArea, T*medianArea, T*maxArea, int w, int i, int j, int c) const
	{
		int k, m;
		T*elem;

		{
			int ii, jj, px, py;

			ii = 2 * w + 1;
			elem = (T*) calloc(ii * ii, sizeof(T));
			k = 0;
			for (ii = -w; ii <= w; ii++)
			{
				px = i + ii;
				if (px < 0 || px >= nx)
					continue;
				for (jj = -w; jj <= w; jj++)
				{
					py = j + jj;
					if (py < 0 || py >= ny)
						continue;
					elem[k++] = get(px, py, c);
				}
			}
		}

		sort<T>(elem, k);
		*minArea = elem[0];
		*maxArea = elem[k - 1];
		m = OddQ(k) ? (k >> 2) + 1 : k >> 2;
		*medianArea = elem[m];
		free(elem);
	}
	double* compactness(T threshold)
	{
		Point2D<double> p[4] =
		{ Point2D<double>(0, 0), Point2D<double>(1, 0), Point2D<double>(1, 1), Point2D<double>(0, 1) };

		const T zero = 0;
		double*compact;
		double tArea, tBoundary;

		compact = new double[channelNo];
		for (int ch = 0; ch < channelNo; ch++)
		{
			tArea = tBoundary = 0.0;

#pragma omp parallel for reduction(+:tArea,tBoundary)
			for (int y = 0; y < ny - 1; y++)
			{
				Point2D<double> c[4];
				Point2D<double> poly[6];
				double area = 0.0;
				double boundary = 0.0;
				for (int x = 0; x < nx - 1; x++)
				{
					T f00, f10, f11, f01;
					int code;
					f00 = get(x, y, ch) - threshold;
					f10 = get(x + 1, y, ch) - threshold;
					f11 = get(x + 1, y + 1, ch) - threshold;
					f01 = get(x, y + 1, ch) - threshold;
					code = (f00 >= zero ? 1 : 0) + (f10 >= zero ? 2 : 0) + (f11 >= zero ? 4 : 0)
							+ (f01 >= zero ? 8 : 0);
					if (EMPTY_PIXEL == code)
						continue;
					if (FILLED_PIXEL == code)
					{
						area += 1.0;
						if (0 == x || nx - 2 == x)
							boundary += 1.0;
						if (0 == y || ny - 2 == y)
							boundary += 1.0;
						continue;
					}

					if (Sign(f00) != Sign(f10))
					{
						c[0].x = f00 / (f00 - f10);
						c[0].y = 0.0;
					}
					if (Sign(f10) != Sign(f11))
					{
						c[1].x = 1.0;
						c[1].y = f10 / (f10 - f11);
					}
					if (Sign(f11) != Sign(f01))
					{
						c[2].x = f01 / (f01 - f11);
						c[2].y = 1.0;
					}
					if (Sign(f01) != Sign(f00))
					{
						c[3].x = 0.0;
						c[3].y = f00 / (f00 - f01);
					}

					switch (code)
					{
					case 1:
						poly[0] = p[0];
						poly[1] = c[0];
						poly[2] = c[3];
						area += computeArea(poly, 3);
						boundary += (c[0] - c[3]).length();
						if (0 == x)
							boundary += (-c[0] + p[0]).length();
						if (0 == y)
							boundary += (c[3] - p[0]).length();
						break;
					case 2:
						poly[0] = c[0];
						poly[1] = p[1];
						poly[2] = c[1];
						area += computeArea(poly, 3);
						boundary += (-c[0] + c[1]).length();
						if (0 == x)
							boundary += (c[0] - p[1]).length();
						if (ny - 2 == y)
							boundary += (-c[1] + p[1]).length();
						break;
					case 3:
						poly[0] = p[0];
						poly[1] = p[1];
						poly[2] = c[1];
						poly[3] = c[3];
						area += computeArea(poly, 4);
						boundary += (c[1] - c[3]).length();
						if (0 == x)
							boundary += (p[0] - p[1]).length();
						if (0 == y)
							boundary += (c[3] - p[0]).length();
						if (ny - 2 == y)
							boundary += (-c[1] + p[1]).length();
						break;
					case 4:
						poly[0] = c[1];
						poly[1] = p[2];
						poly[2] = c[2];
						area += computeArea(poly, 3);
						boundary += (-c[1] + c[2]).length();
						if (nx - 2 == x)
							boundary += (-c[2] + p[2]).length();
						if (ny - 2 == y)
							boundary += (c[1] - p[2]).length();
						break;
					case 5:
						poly[0] = p[0];
						poly[1] = c[0];
						poly[2] = c[1];
						poly[3] = p[2];
						poly[4] = c[2];
						poly[5] = c[3];
						area += computeArea(poly, 6);
						boundary += (c[0] - c[1]).length() + (c[2] - c[3]).length();
						if (0 == x)
							boundary += (-c[0] + p[0]).length();
						if (nx - 2 == x)
							boundary += (-c[2] + p[2]).length();
						if (0 == y)
							boundary += (c[3] - p[0]).length();
						if (ny - 2 == y)
							boundary += (c[1] - p[2]).length();
						break;
					case 6:
						poly[0] = c[0];
						poly[1] = p[1];
						poly[2] = p[2];
						poly[3] = c[2];
						area += computeArea(poly, 4);
						boundary += (-c[0] + c[2]).length();
						if (0 == x)
							boundary += (c[0] - p[1]).length();
						if (nx - 2 == x)
							boundary += (-c[2] + p[2]).length();
						if (ny - 2 == y)
							boundary += (p[1] - p[2]).length();
						break;
					case 7:
						poly[0] = p[0];
						poly[1] = p[1];
						poly[2] = p[2];
						poly[3] = c[2];
						poly[4] = c[3];
						area += computeArea(poly, 5);
						boundary += (c[2] - c[3]).length();
						if (0 == x)
							boundary += (p[0] - p[1]).length();
						if (nx - 2 == x)
							boundary += (-c[2] + p[2]).length();
						if (0 == y)
							boundary += (c[3] - p[0]).length();
						if (ny - 2 == y)
							boundary += (p[1] - p[2]).length();
						break;
					case 8:
						poly[0] = c[3];
						poly[1] = c[2];
						poly[2] = p[3];
						area += computeArea(poly, 3);
						boundary += (-c[2] + c[3]).length();
						if (nx - 2 == x)
							boundary += (c[2] - p[3]).length();
						if (0 == y)
							boundary += (-c[3] + p[3]).length();
						break;
					case 9:
						poly[0] = p[0];
						poly[1] = c[0];
						poly[2] = c[2];
						poly[3] = p[3];
						area += computeArea(poly, 4);
						boundary += (c[0] - c[2]).length();
						if (0 == x)
							boundary += (-c[0] + p[0]).length();
						if (nx - 2 == x)
							boundary += (c[2] - p[3]).length();
						if (0 == y)
							boundary += (-p[0] + p[3]).length();
						break;
					case 10:
						poly[0] = c[0];
						poly[1] = p[1];
						poly[2] = c[1];
						poly[3] = c[2];
						poly[4] = p[3];
						poly[5] = c[3];
						area += computeArea(poly, 6);
						boundary += (c[1] - c[2]).length() + (-c[0] + c[3]).length();
						if (0 == x)
							boundary += (c[0] - p[1]).length();
						if (nx - 2 == x)
							boundary += (c[2] - p[3]).length();
						if (0 == y)
							boundary += (-c[3] + p[3]).length();
						if (ny - 2 == y)
							boundary += (-c[1] + p[1]).length();
						break;
					case 11:
						poly[0] = p[0];
						poly[1] = p[1];
						poly[2] = c[1];
						poly[3] = c[2];
						poly[4] = p[3];
						area += computeArea(poly, 5);
						boundary += (c[1] - c[2]).length();
						if (0 == x)
							boundary += (p[0] - p[1]).length();
						if (nx - 2 == x)
							boundary += (c[2] - p[3]).length();
						if (0 == y)
							boundary += (-p[0] + p[3]).length();
						if (ny - 2 == y)
							boundary += (-c[1] + p[1]).length();
						break;
					case 12:
						poly[0] = c[3];
						poly[1] = c[1];
						poly[2] = p[2];
						poly[3] = p[3];
						area += computeArea(poly, 4);
						boundary += (-c[1] + c[3]).length();
						if (nx - 2 == x)
							boundary += (p[2] - p[3]).length();
						if (0 == y)
							boundary += (-c[3] + p[3]).length();
						if (ny - 2 == y)
							boundary += (c[1] - p[2]).length();
						break;
					case 13:
						poly[0] = p[0];
						poly[1] = c[0];
						poly[2] = c[1];
						poly[3] = p[2];
						poly[4] = p[3];
						area += computeArea(poly, 5);
						boundary += (c[0] - c[1]).length();
						if (0 == x)
							boundary += (-c[0] + p[0]).length();
						if (nx - 2 == x)
							boundary += (p[2] - p[3]).length();
						if (0 == y)
							boundary += (-p[0] + p[3]).length();
						if (ny - 2 == y)
							boundary += (c[1] - p[2]).length();
						break;
					case 14:
						poly[0] = c[0];
						poly[1] = p[1];
						poly[2] = p[2];
						poly[3] = p[3];
						poly[4] = c[3];
						area += computeArea(poly, 5);
						boundary += (-c[0] + c[3]).length();
						if (0 == x)
							boundary += (c[0] - p[1]).length();
						if (nx - 2 == x)
							boundary += (p[2] - p[3]).length();
						if (0 == y)
							boundary += (-c[3] + p[3]).length();
						if (ny - 2 == y)
							boundary += (p[1] - p[2]).length();
						break;
					}

				}
				tArea += area;
				tBoundary += boundary;
			}

			if (tBoundary > 0.0 && tArea > 0.0)
				compact[ch] = tBoundary * tBoundary / tArea;
			else
				compact[ch] = 0.0;
		}
		return compact;
	}
};

#endif
