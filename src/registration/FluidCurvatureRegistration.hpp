/*
 * fluidCurvatureRegistration.h
 *
 *  Created on: Jul 10, 2012
 *      Author: Konstantin Thierbach
 */

#ifndef FLUIDCURVATUREREGISTRATION_H_
#define FLUIDCURVATUREREGISTRATION_H_

#include <iostream>
#include <stdlib.h>

#include "../filter/MedianFilter.hpp"
#include "../fourier_transform/Fourier.hpp"
#include "../fourier_transform/InverseFourier.hpp"
#include "../fourier_transform/CosFourier.hpp"
#include "../fourier_transform/InverseCosFourier.hpp"
#include "../fourier_transform/SinFourier.hpp"
#include "../fourier_transform/InverseSinFourier.hpp"
#include "../templates/ImageTemplate.hpp"
#include "../templates/VectorArray2D.hpp"
#include "../templates/Vector2D.hpp"
#include "../solver/RKNystroem.hpp"
#include "../solver/RKV43.hpp"
#include "../utilities/ImageDifference.hpp"
#include "../utilities/BracketMethod.hpp"
#include "../../lib/fftw3.h"
#include "../../lib/CImg-1.5.0/CImg.h"

#define MAXDOUBLE (1.7976931348623157e308)
#define Re(c) (c)[0]
#define Im(c) (c)[1]
#define addr(i,j) ((i)+(nx)*(j))

class RKNystroemDSolve;

class FluidCurvatureRegistration
{
	public:
		FluidCurvatureRegistration(int boundary, int ny, int nx, double*tdata, long int, double*sdata, long int,
				double tmax, double dtStart, double alpha, double bharmfric, double fric, double lameMu,
				double lameLambda, double vortexw, double error, double matchError, int returnType,
				const char*vfieldfile, int refchannels, double *refpat, long rlen);
		~FluidCurvatureRegistration();

		void registerImages();
		VectorArray2D* getFlowField() const;
		Image<double>* getReference() const;
		Image<double>* getSample() const;
		double getMinimalTime() const;
		double getMismatchError() const;

	private:
		Image<double> *templateImage, *sampleImage, *__wraped, *ref;
		VectorArray2D *u, *du, *_bestU;
		RKNystroemDSolve *rknystroem;
		RKV43 *rk43;
		cimg_library::CImgDisplay display;

		int boundary, ny, nx, returnType, refchannels;
		long rlen;
		double tmax, dtStart, alpha, bharmfric, fric, lameMu, lameLambda, vortexw, error, matchError;
		double BiharmonicFriction, Friction, InverseAlpha, LameLambda, LameMu, VortexWeight;
		double _lowestError;
		static const double Pi = M_PI;
		static const double TwoPi = 2.0 * M_PI;
		double ImageMatchGoal, DampingRatio, ImageMatchError, MinimumTime;
		fftw_complex *_fx, *_fy, *_fxhat, *_fyhat;

		void (FluidCurvatureRegistration::*BiharmonicSolve)(VectorArray2D*, VectorArray2D*);
		void (FluidCurvatureRegistration::*LameSolve)(VectorArray2D*, VectorArray2D*);
		void (FluidCurvatureRegistration::*OverdampedLameSolve)(VectorArray2D*, VectorArray2D*, double);
		void (FluidCurvatureRegistration::*fluid)(double, double*, double*, int);

		void PeriodicBiharmonicSolve(VectorArray2D*solution, VectorArray2D*ff);
		void ZeroDerivativeBiharmonicSolve(VectorArray2D*solution, VectorArray2D*ff);
		void ZeroBiharmonicSolve(VectorArray2D*solution, VectorArray2D*ff);
		void PeriodicHarmonicSolve(VectorArray2D*solution, VectorArray2D*ff);
		void OverdampedZeroDerivativeHarmonicSolve(VectorArray2D*solution, VectorArray2D*ff, double ratio);
		void ZeroDerivativeHarmonicSolve(VectorArray2D*solution, VectorArray2D*ff);
		void ZeroHarmonicSolve(VectorArray2D*solution, VectorArray2D*ff);
		void PeriodicLameSolve(VectorArray2D*solution, VectorArray2D*ff);
		void OverdampedPeriodicBiharmonicSolve(VectorArray2D*solution, VectorArray2D*ff, double ratio);
		void OverdampedZeroDerivativeBiharmonicSolve(VectorArray2D*solution, VectorArray2D*ff, double ratio);
		void OverdampedZeroBiharmonicSolve(VectorArray2D*solution, VectorArray2D*ff, double ratio);
		void OverdampedPeriodicHarmonicSolve(VectorArray2D*solution, VectorArray2D*ff, double ratio);
		void OverdampedZeroHarmonicSolve(VectorArray2D*solution, VectorArray2D*ff, double ratio);

		void FluidLame(double, double*y, double*dy, int n);
		void FluidOverDamped(double, double*y, double*dy, int n);

		void FluidBiharmonic(double*d2y, double, double*y, double*dy, int n);

		int PrintFluidProgress(double t, double h, double*u, double*up, double*uNext, double*upNext, int n);
		int PrintFluidProgress1(double t, double h, double*u, double*uNext, int n);

};

#endif /* FLUIDCURVATUREREGISTRATION_H_ */
