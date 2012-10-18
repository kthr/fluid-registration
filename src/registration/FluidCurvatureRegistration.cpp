/*
 * fluidCurvatureRegistration.cpp
 *
 *  Created on: Jul 10, 2012
 *      Author: Konstantin Thierbach
 */

#include "FluidCurvatureRegistration.hpp"

FluidCurvatureRegistration::FluidCurvatureRegistration(int boundary, int ny, int nx, double*tdata, long int,
		double*sdata, long int, double tmax, double dtStart, double alpha, double bharmfric, double fric, double lameMu,
		double lameLambda, double vortexw, double error, double matchError, int returnType, const char *vfieldfile,
		int refchannels, double *refpat, long rlen) :
		boundary(boundary), ny(ny), nx(nx), tmax(tmax), dtStart(dtStart), alpha(alpha), bharmfric(bharmfric), fric(
				fric), lameMu(lameMu), lameLambda(lameLambda), vortexw(vortexw), error(error), matchError(matchError), returnType(
				returnType), refchannels(refchannels), rlen(rlen)
{
	templateImage = new Image<double>(tdata, nx, ny, 1);
	sampleImage = new Image<double>(sdata, nx, ny, 1);
	__wraped = new Image<double>(nx, ny, 1);
	if (refpat != NULL)
	{
		ref = new Image<double>(refpat, nx, ny, refchannels);
	}
	else
	{
		ref = NULL;
	}
	_fx = (fftw_complex*) fftw_malloc(nx * ny * sizeof(fftw_complex));
	_fy = (fftw_complex*) fftw_malloc(nx * ny * sizeof(fftw_complex));
	_fxhat = (fftw_complex*) fftw_malloc(nx * ny * sizeof(fftw_complex));
	_fyhat = (fftw_complex*) fftw_malloc(nx * ny * sizeof(fftw_complex));
	u = new VectorArray2D(nx, ny);
	du = new VectorArray2D(nx, ny);
	_bestU = new VectorArray2D(nx, ny);

	ImageMatchGoal = 0.0;
	DampingRatio = 1.0;
	ImageMatchError = MAXDOUBLE;
	MinimumTime = 0.0;
	_lowestError = MAXDOUBLE;
}
FluidCurvatureRegistration::~FluidCurvatureRegistration()
{
	delete templateImage;
	delete sampleImage;
	delete __wraped;
	delete ref;
	fftw_free(_fx);
	fftw_free(_fy);
	fftw_free(_fxhat);
	fftw_free(_fyhat);
	delete u;
	delete du;
	delete _bestU;
	delete rknystroem;
	delete rk43;

}
void FluidCurvatureRegistration::registerImages()
{

	int dim;
	double sMin, sMax;
	double*rMin = NULL, *rMax = NULL;

	sMin = sampleImage->min(0);
	sMax = sampleImage->max(0);

	if (rlen > 0)
	{
		rMin = new double[refchannels];
		rMax = new double[refchannels];
		for (int ch = 0; ch < refchannels; ch++)
		{
			rMin[ch] = ref->min(ch);
			rMax[ch] = ref->max(ch);
		}
	}
	InverseAlpha = 1.0 / alpha;
	ImageMatchGoal = matchError;
	BiharmonicFriction = bharmfric;
	Friction = fric;
	VortexWeight = vortexw;
	LameMu = lameMu;
	LameLambda = lameLambda;
	DampingRatio = alpha;
	dim = 2 * nx * ny;

	switch (boundary)
	{
		case 0:
			BiharmonicSolve = &FluidCurvatureRegistration::PeriodicBiharmonicSolve;
			break;
		case 1:
			BiharmonicSolve = &FluidCurvatureRegistration::ZeroDerivativeBiharmonicSolve;
			break;
		case 2:
			BiharmonicSolve = &FluidCurvatureRegistration::ZeroBiharmonicSolve;
			break;
		case 3:
			BiharmonicSolve = &FluidCurvatureRegistration::PeriodicHarmonicSolve;
			break;
		case 4:
			BiharmonicSolve = &FluidCurvatureRegistration::ZeroDerivativeHarmonicSolve;
			break;
		case 5:
			BiharmonicSolve = &FluidCurvatureRegistration::ZeroHarmonicSolve;
			break;
		case 6:
			fluid = &FluidCurvatureRegistration::FluidLame;
			LameSolve = &FluidCurvatureRegistration::PeriodicLameSolve;
			break;
		case 7:
			fluid = &FluidCurvatureRegistration::FluidOverDamped;
			OverdampedLameSolve = &FluidCurvatureRegistration::OverdampedPeriodicBiharmonicSolve;
			break;
		case 8:
			fluid = &FluidCurvatureRegistration::FluidOverDamped;
			OverdampedLameSolve = &FluidCurvatureRegistration::OverdampedZeroDerivativeBiharmonicSolve;
			break;
		case 9:
			fluid = &FluidCurvatureRegistration::FluidOverDamped;
			OverdampedLameSolve = &FluidCurvatureRegistration::OverdampedZeroBiharmonicSolve;
			break;
		case 10:
			fluid = &FluidCurvatureRegistration::FluidOverDamped;
			OverdampedLameSolve = &FluidCurvatureRegistration::OverdampedPeriodicHarmonicSolve;
			break;
		case 11:
			fluid = &FluidCurvatureRegistration::FluidOverDamped;
			OverdampedLameSolve = &FluidCurvatureRegistration::OverdampedZeroDerivativeHarmonicSolve;
			break;
		case 12:
			fluid = &FluidCurvatureRegistration::FluidOverDamped;
			OverdampedLameSolve = &FluidCurvatureRegistration::OverdampedZeroHarmonicSolve;
			break;
		default:
			fluid = &FluidCurvatureRegistration::FluidOverDamped;
			OverdampedLameSolve = &FluidCurvatureRegistration::OverdampedPeriodicBiharmonicSolve;
			boundary = 7;
			break;
	}

	if (boundary < 6)
	{
		//rknystroem = new RKNystroemDSolve(this,dim, &FluidCurvatureRegistration::FluidBiharmonic, PrintFluidProgress);
		rknystroem = new RKNystroemDSolve(this, dim, &FluidCurvatureRegistration::FluidBiharmonic, NULL);
		rknystroem->integrate(0.0, tmax, dtStart, u->vx, du->vx, error, true);
	}
	else
	{
		rk43 = new RKV43(this);
		//RKVMethod43(0.0, tmax, dtStart, u->vx, dim, fluid, PrintFluidProgress1, error, true);
		rk43->RKVMethod43(0.0, tmax, dtStart, u->vx, dim, fluid, &FluidCurvatureRegistration::PrintFluidProgress1,
				error, true);
		//rk43->RKVMethod43(0.0, tmax, dtStart, u->vx, dim, fluid, NULL, error, true);
		(this->*fluid)(tmax, u->vx, du->vx, dim);
	}
	fftw_cleanup();
	if (_lowestError < ImageMatchError)
	{
		for (int i = 0; i < nx * ny; i++)
		{
			u->vx[i] = _bestU->vx[i];
			u->vy[i] = _bestU->vy[i];
		}
		ImageMatchError = _lowestError;
	}
	if (0.0 == MinimumTime)
		MinimumTime = tmax;

	sampleImage->wrap(u);
	sampleImage->clipRange(0, sMin, sMax);
	if (ref)
	{
		ref->wrap(u);
		for (int ch = 0; ch < ref->channelNo; ch++)
			ref->clipRange(ch, rMin[ch], rMax[ch]);
	}

	std::cout << "MismatchError->" << ImageMatchError << ", MinimumTime->" << MinimumTime << std::endl;

	delete[] rMin;
	delete[] rMax;
}
void FluidCurvatureRegistration::PeriodicBiharmonicSolve(VectorArray2D*solution, VectorArray2D*ff)
{
	int nx, ny, mu;
	double dx2, dy2, dx4, dy4, nrm;
	double meanX = 0.0, meanY = 0.0;

	nx = ff->nx;
	ny = ff->ny;
	dx2 = ff->dx * ff->dx;
	dx4 = dx2 * dx2;
	dy2 = ff->dy * ff->dy;
	dy4 = dy2 * dy2;
	mu = nx * ny;

	for (int i = 0; i < mu; i++)
	{
		meanX += ff->vx[i];
		meanY += ff->vy[i];
	}
	meanX /= mu;
	meanY /= mu;

	for (int i = 0; i < mu; i++)
	{
		Re (_fx[i]) = ff->vx[i] - meanX;
		Re (_fy[i]) = ff->vy[i] - meanY;
		Im(_fx[i]) = Im(_fy[i]) = 0.0;
	}
	Fourier(_fxhat, _fx, nx, ny);
	Fourier(_fyhat, _fy, nx, ny);
#pragma omp parallel for
	for (int l = 0; l < ny; l++)
	{
		double lambda = TwoPi * l / ny;
		double cl1 = cos(lambda);
		double cl2 = cos(2.0 * lambda);
		double cl3 = cos(3.0 * lambda);
		for (int k = 0; k < nx; k++)
		{
			double kappa = TwoPi * k / nx;
			double ck1 = cos(kappa);
			double ck2 = cos(2.0 * kappa);
			double ck3 = cos(3.0 * kappa);

			double l2 = (-6 * (-28 + 39 * cl1 - 12 * cl2 + cl3) * dx4
					+ (-15 + 16 * ck1 - ck2) * (-15 + 16 * cl1 - cl2) * dx2 * dy2
					- 6 * (-28 + 39 * ck1 - 12 * ck2 + ck3) * dy4) / (18. * dx4 * dy4);

			if (fabs(l2) > TINY)
			{
				l2 = 1.0 / l2;
				int pos = addr(k,l);
				Re(_fxhat[pos]) *= l2;
				Im(_fxhat[pos]) *= l2;
				Re(_fyhat[pos]) *= l2;
				Im(_fyhat[pos]) *= l2;
			}
			else
			{

				Re(_fxhat[addr(k,l)]) /= TINY;
				Im(_fxhat[addr(k,l)]) /= TINY;
				Re(_fyhat[addr(k,l)]) /= TINY;
				Im(_fyhat[addr(k,l)]) /= TINY;
			}
		}
	}
	nrm = 1.0 / (nx * ny);
	InverseFourier(_fx, _fxhat, nx, ny);
	InverseFourier(_fy, _fyhat, nx, ny);

	for (int i = 0; i < mu; i++)
	{
		solution->vx[i] = Re(_fx[i]) * nrm + meanX;
		solution->vy[i] = Re(_fy[i]) * nrm + meanY;
	}

}
void FluidCurvatureRegistration::ZeroDerivativeBiharmonicSolve(VectorArray2D*solution, VectorArray2D*ff)
{
	double*fxhat, *fyhat, *fx, *fy, omegak, omegal;
	int nx, ny, mu;
	double meanX, meanY;
	double dx, dy, dx2, dy2, dx4, dy4;

	nx = ff->nx;
	ny = ff->ny;
	dx = ff->dx;
	dy = ff->dy;
	omegak = Pi / nx;
	omegal = Pi / ny;
	dx2 = dx * dx;
	dx4 = dx2 * dx2;
	dy2 = dy * dy;
	dy4 = dy2 * dy2;
	mu = nx * ny;
	fx = new double[mu];
	fy = new double[mu];

	fxhat = new double[mu];
	fyhat = new double[mu];
	meanX = meanY = 0.0;
	for (int i = 0; i < mu; i++)
	{
		fx[i] = ff->vx[i];
		fy[i] = ff->vy[i];
		meanX += fx[i];
		meanY += fy[i];
	}
	meanX /= mu;
	meanY /= mu;
	for (int i = 0; i < mu; i++)
	{
		fx[i] -= meanX;
		fy[i] -= meanY;
	}
	CosFourier(fxhat, fx, nx, ny);
	CosFourier(fyhat, fy, nx, ny);
#pragma omp parallel for
	for (int l = 0; l < ny; l++)
	{
		double cl1 = cos((0.5 + l) * omegal);
		double cl2 = cos((1 + 2 * l) * omegal);
		double cl3 = cos((1.5 + 3 * l) * omegal);
		for (int k = 0; k < nx; k++)
		{
			double ck1 = cos((0.5 + k) * omegak);
			double ck2 = cos((1 + 2 * k) * omegak);
			double ck3 = cos((1.5 + 3 * k) * omegak);

			double l2 = (-6 * (-28 + 39 * cl1 - 12 * cl2 + cl3) * dx4
					+ (-15 + 16 * ck1 - ck2) * (-15 + 16 * cl1 - cl2) * dx2 * dy2
					- 6 * (-28 + 39 * ck1 - 12 * ck2 + ck3) * dy4) / (18. * dx4 * dy4);
			if (fabs(l2) > TINY)
			{
				l2 = 1.0 / l2;
				int pos = addr(k,l);
				fxhat[pos] *= l2;
				fyhat[pos] *= l2;
			}
			else
			{
				fxhat[addr(k,l)] /= TINY;
				fyhat[addr(k,l)] /= TINY;
			}
		}
	}

	InverseCosFourier(solution->vx, fxhat, nx, ny);
	InverseCosFourier(solution->vy, fyhat, nx, ny);
	for (int i = 0; i < mu; i++)
	{
		solution->vx[i] += meanX;
		solution->vy[i] += meanY;
	}
	delete[] fyhat;
	delete[] fxhat;
	delete[] fy;
	delete[] fx;
}
void FluidCurvatureRegistration::ZeroBiharmonicSolve(VectorArray2D*solution, VectorArray2D*ff)
{
	double*fxhat, *fyhat, *fx, *fy, omegak, omegal;
	int nx, ny, mu;
	double meanX, meanY;
	double dx, dy, dx2, dy2, dx4, dy4;

	nx = ff->nx;
	ny = ff->ny;
	dx = ff->dx;
	dy = ff->dy;
	omegak = Pi / nx;
	omegal = Pi / ny;
	dx2 = dx * dx;
	dx4 = dx2 * dx2;
	dy2 = dy * dy;
	dy4 = dy2 * dy2;
	mu = nx * ny;
	fx = new double[mu];
	fy = new double[mu];

	fxhat = new double[mu];
	fyhat = new double[mu];
	meanX = meanY = 0.0;
	for (int i = 0; i < mu; i++)
	{
		fx[i] = ff->vx[i];
		fy[i] = ff->vy[i];
		meanX += fx[i];
		meanY += fy[i];
	}
	meanX /= mu;
	meanY /= mu;
	for (int i = 0; i < mu; i++)
	{
		fx[i] -= meanX;
		fy[i] -= meanY;
	}
	SinFourier(fxhat, fx, nx, ny);
	SinFourier(fyhat, fy, nx, ny);
#pragma omp parallel for
	for (int l = 0; l < ny; l++)
	{
		double cl1 = cos((0.5 + l) * omegal);
		double cl2 = cos((1 + 2 * l) * omegal);
		double cl3 = cos((1.5 + 3 * l) * omegal);
		for (int k = 0; k < nx; k++)
		{
			double ck1 = cos((0.5 + k) * omegak);
			double ck2 = cos((1 + 2 * k) * omegak);
			double ck3 = cos((1.5 + 3 * k) * omegak);

			double l2 = (-6 * (-28 + 39 * cl1 - 12 * cl2 + cl3) * dx4
					+ (-15 + 16 * ck1 - ck2) * (-15 + 16 * cl1 - cl2) * dx2 * dy2
					- 6 * (-28 + 39 * ck1 - 12 * ck2 + ck3) * dy4) / (18. * dx4 * dy4);

			if (fabs(l2) > TINY)
			{
				l2 = 1.0 / l2;
				int pos = addr(k,l);
				fxhat[pos] *= l2;
				fyhat[pos] *= l2;
			}
			else
			{
				fxhat[addr(k,l)] /= TINY;
				fyhat[addr(k,l)] /= TINY;
			}
		}
	}

	InverseSinFourier(solution->vx, fxhat, nx, ny);
	InverseSinFourier(solution->vy, fyhat, nx, ny);
	for (int i = 0; i < mu; i++)
	{
		solution->vx[i] += meanX;
		solution->vy[i] += meanY;
	}
	delete[] fyhat;
	delete[] fxhat;
	delete[] fy;
	delete[] fx;
}
void FluidCurvatureRegistration::PeriodicHarmonicSolve(VectorArray2D*solution, VectorArray2D*ff)
{
	int nx, ny, mu;
	double dx2, dy2, nrm;
	double meanX = 0.0, meanY = 0.0;

	nx = ff->nx;
	ny = ff->ny;
	dx2 = ff->dx * ff->dx;
	dy2 = ff->dy * ff->dy;
	mu = nx * ny;

	for (int i = 0; i < mu; i++)
	{
		meanX += ff->vx[i];
		meanY += ff->vy[i];
	}
	meanX /= mu;
	meanY /= mu;

	for (int i = 0; i < mu; i++)
	{
		Re(_fx[i]) = ff->vx[i] - meanX;
		Re(_fy[i]) = ff->vy[i] - meanY;
		Im(_fx[i]) = Im(_fy[i]) = 0.0;
	}
	Fourier(_fxhat, _fx, nx, ny);
	Fourier(_fyhat, _fy, nx, ny);
#pragma omp parallel for
	for (int l = 0; l < ny; l++)
	{
		double lambda = TwoPi * l / ny;
		double cl1 = cos(lambda);
		double cl2 = cos(2.0 * lambda);
		for (int k = 0; k < nx; k++)
		{
			double kappa = TwoPi * k / nx;
			double ck1 = cos(kappa);
			double ck2 = cos(2.0 * kappa);

			double l2 = ((15 - 16 * cl1 + cl2) * dx2 + (15 - 16 * ck1 + ck2) * dy2) / (6. * dx2 * dy2);
			if (fabs(l2) > TINY)
			{
				l2 = 1.0 / l2;
				int pos = addr(k,l);
				Re(_fxhat[pos]) *= l2;
				Im(_fxhat[pos]) *= l2;
				Re(_fyhat[pos]) *= l2;
				Im(_fyhat[pos]) *= l2;
			}
			else
			{

				Re(_fxhat[addr(k,l)]) /= TINY;
				Im(_fxhat[addr(k,l)]) /= TINY;
				Re(_fyhat[addr(k,l)]) /= TINY;
				Im(_fyhat[addr(k,l)]) /= TINY;
			}
		}
	}
	nrm = 1.0 / (nx * ny);
	InverseFourier(_fx, _fxhat, nx, ny);
	InverseFourier(_fy, _fyhat, nx, ny);

	for (int i = 0; i < mu; i++)
	{
		solution->vx[i] = Re(_fx[i]) * nrm + meanX;
		solution->vy[i] = Re(_fy[i]) * nrm + meanY;
	}

}
void FluidCurvatureRegistration::OverdampedZeroDerivativeHarmonicSolve(VectorArray2D*solution, VectorArray2D*ff,
		double ratio)
{
	double*fxhat, *fyhat, *fx, *fy, omegak, omegal;
	int nx, ny, mu;
	double meanX, meanY;
	double dx, dy, dx2, dy2;

	nx = ff->nx;
	ny = ff->ny;
	dx = ff->dx;
	dy = ff->dy;
	omegak = Pi / nx;
	omegal = Pi / ny;
	dx2 = dx * dx;
	dy2 = dy * dy;
	mu = nx * ny;
	fx = new double[mu];
	fy = new double[mu];

	fxhat = new double[mu];
	fyhat = new double[mu];
	meanX = meanY = 0.0;
	for (int i = 0; i < mu; i++)
	{
		fx[i] = ff->vx[i];
		fy[i] = ff->vy[i];
		meanX += fx[i];
		meanY += fy[i];
	}
	meanX /= mu;
	meanY /= mu;
	for (int i = 0; i < mu; i++)
	{
		fx[i] -= meanX;
		fy[i] -= meanY;
	}
	CosFourier(fxhat, fx, nx, ny);
	CosFourier(fyhat, fy, nx, ny);
#pragma omp parallel for
	for (int l = 0; l < ny; l++)
	{
		double cl1 = cos((0.5 + l) * omegal);
		double cl2 = cos((1 + 2 * l) * omegal);
		for (int k = 0; k < nx; k++)
		{
			double ck1 = cos((0.5 + k) * omegak);
			double ck2 = cos((1 + 2 * k) * omegak);

			double l2 = 1.0 - ratio * (-((15 - 16 * ck1 + ck2) / dx2) - (15 - 16 * cl1 + cl2) / dy2) / 6.;
			if (fabs(l2) > TINY)
			{
				l2 = 1.0 / l2;
				int pos = addr(k,l);
				fxhat[pos] *= l2;
				fyhat[pos] *= l2;
			}
			else
			{
				fxhat[addr(k,l)] /= TINY;
				fyhat[addr(k,l)] /= TINY;
			}
		}
	}

	InverseCosFourier(solution->vx, fxhat, nx, ny);
	InverseCosFourier(solution->vy, fyhat, nx, ny);
	for (int i = 0; i < mu; i++)
	{
		solution->vx[i] += meanX;
		solution->vy[i] += meanY;
	}
	delete[] fyhat;
	delete[] fxhat;
	delete[] fy;
	delete[] fx;
}
void FluidCurvatureRegistration::ZeroDerivativeHarmonicSolve(VectorArray2D*solution, VectorArray2D*ff)
{
	double*fxhat, *fyhat, *fx, *fy, omegak, omegal;
	int nx, ny, mu;
	double meanX, meanY;
	double dx, dy, dx2, dy2;

	nx = ff->nx;
	ny = ff->ny;
	dx = ff->dx;
	dy = ff->dy;
	omegak = Pi / nx;
	omegal = Pi / ny;
	dx2 = dx * dx;
	dy2 = dy * dy;
	mu = nx * ny;
	fx = new double[mu];
	fy = new double[mu];

	fxhat = new double[mu];
	fyhat = new double[mu];
	meanX = meanY = 0.0;
	for (int i = 0; i < mu; i++)
	{
		fx[i] = ff->vx[i];
		fy[i] = ff->vy[i];
		meanX += fx[i];
		meanY += fy[i];
	}
	meanX /= mu;
	meanY /= mu;
	for (int i = 0; i < mu; i++)
	{
		fx[i] -= meanX;
		fy[i] -= meanY;
	}
	CosFourier(fxhat, fx, nx, ny);
	CosFourier(fyhat, fy, nx, ny);
#pragma omp parallel for
	for (int l = 0; l < ny; l++)
	{
		double cl1 = cos((0.5 + l) * omegal);
		double cl2 = cos((1 + 2 * l) * omegal);
		for (int k = 0; k < nx; k++)
		{
			double ck1 = cos((0.5 + k) * omegak);
			double ck2 = cos((1 + 2 * k) * omegak);

			double l2 = -(-((15 - 16 * ck1 + ck2) / dx2) - (15 - 16 * cl1 + cl2) / dy2) / 6.;
			if (fabs(l2) > TINY)
			{
				l2 = 1.0 / l2;
				int pos = addr(k,l);
				fxhat[pos] *= l2;
				fyhat[pos] *= l2;
			}
			else
			{
				fxhat[addr(k,l)] /= TINY;
				fyhat[addr(k,l)] /= TINY;
			}
		}
	}

	InverseCosFourier(solution->vx, fxhat, nx, ny);
	InverseCosFourier(solution->vy, fyhat, nx, ny);
	for (int i = 0; i < mu; i++)
	{
		solution->vx[i] += meanX;
		solution->vy[i] += meanY;
	}
	delete[] fyhat;
	delete[] fxhat;
	delete[] fy;
	delete[] fx;
}
void FluidCurvatureRegistration::ZeroHarmonicSolve(VectorArray2D*solution, VectorArray2D*ff)
{
	double*fxhat, *fyhat, *fx, *fy, omegak, omegal;
	int nx, ny, mu;
	double meanX, meanY;
	double dx, dy, dx2, dy2;

	nx = ff->nx;
	ny = ff->ny;
	dx = ff->dx;
	dy = ff->dy;
	omegak = Pi / nx;
	omegal = Pi / ny;
	dx2 = dx * dx;
	dy2 = dy * dy;
	mu = nx * ny;
	fx = new double[mu];
	fy = new double[mu];

	fxhat = new double[mu];
	fyhat = new double[mu];
	meanX = meanY = 0.0;
	for (int i = 0; i < mu; i++)
	{
		fx[i] = ff->vx[i];
		fy[i] = ff->vy[i];
		meanX += fx[i];
		meanY += fy[i];
	}
	meanX /= mu;
	meanY /= mu;
	for (int i = 0; i < mu; i++)
	{
		fx[i] -= meanX;
		fy[i] -= meanY;
	}
	SinFourier(fxhat, fx, nx, ny);
	SinFourier(fyhat, fy, nx, ny);
#pragma omp parallel for
	for (int l = 0; l < ny; l++)
	{
		double cl1 = cos((0.5 + l) * omegal);
		double cl2 = cos((1 + 2 * l) * omegal);

		for (int k = 0; k < nx; k++)
		{
			double ck1 = cos((0.5 + k) * omegak);
			double ck2 = cos((1 + 2 * k) * omegak);
			double l2 = -(-((15 - 16 * ck1 + ck2) / dx2) - (15 - 16 * cl1 + cl2) / dy2) / 6.;

			if (fabs(l2) > TINY)
			{
				l2 = 1.0 / l2;
				int pos = addr(k,l);
				fxhat[pos] *= l2;
				fyhat[pos] *= l2;
			}
			else
			{
				fxhat[addr(k,l)] /= TINY;
				fyhat[addr(k,l)] /= TINY;
			}
		}
	}

	InverseSinFourier(solution->vx, fxhat, nx, ny);
	InverseSinFourier(solution->vy, fyhat, nx, ny);
	for (int i = 0; i < mu; i++)
	{
		solution->vx[i] += meanX;
		solution->vy[i] += meanY;
	}
	delete[] fyhat;
	delete[] fxhat;
	delete[] fy;
	delete[] fx;
}
void FluidCurvatureRegistration::PeriodicLameSolve(VectorArray2D*solution, VectorArray2D*ff)
{
	int nx, ny, mu;
	double dx2, dy2, nrm;
	double meanX = 0.0, meanY = 0.0;

	nx = ff->nx;
	ny = ff->ny;
	dx2 = ff->dx * ff->dx;
	dy2 = ff->dy * ff->dy;
	mu = nx * ny;

	for (int i = 0; i < mu; i++)
	{
		meanX += ff->vx[i];
		meanY += ff->vy[i];
	}
	meanX /= mu;
	meanY /= mu;

	for (int i = 0; i < mu; i++)
	{
		Re(_fx[i]) = ff->vx[i] - meanX;
		Re(_fy[i]) = ff->vy[i] - meanY;
		Im(_fx[i]) = Im(_fy[i]) = 0.0;
	}
	Fourier(_fxhat, _fx, nx, ny);
	Fourier(_fyhat, _fy, nx, ny);

	double mu2 = LameMu * LameMu;
	double l2 = LameLambda * LameLambda;
	double lm = LameMu * LameLambda;
	double lpm = LameMu + LameLambda;
	double lpm2 = lpm * lpm;
#pragma omp parallel for
	for (int l = 0; l < ny; l++)
	{
		double omegal = TwoPi * l / ny;
		double sl = sin(omegal), sl2 = sin(2 * omegal), cl = cos(omegal), cl2 = cos(2 * omegal), cl3 = cos(3 * omegal),
				cl4 = cos(4 * omegal);
		double mat[2][2];
		for (int k = 0; k < nx; k++)
		{
			double omegak = TwoPi * k / nx;
			double sk = sin(omegak), sk2 = sin(2 * omegak), ck = cos(omegak), ck2 = cos(2 * omegak), ck3 = cos(
					3 * omegak), ck4 = cos(4 * omegak);

			double den = ((16 * ck * (-2095 + 2288 * cl - 208 * cl2) - 16 * ck2 * (-395 + 208 * cl + 247 * cl2)
					+ 5 * (5635 - 208 * ck3 + 13 * ck4 - 6704 * cl + 1264 * cl2 - 208 * cl3 + 13 * cl4)) * l2
					+ 2
							* (111479 - 2192 * ck3 + 101 * ck4 - 103792 * cl + 16 * ck * (-6487 + 4592 * cl - 352 * cl2)
									+ 19856 * cl2 - 16 * ck2 * (-1241 + 352 * cl + 238 * cl2) - 2192 * cl3 + 101 * cl4)
							* lm
					+ (16 * (ck + 4 * ck2) * (16 * cl3 - cl4) + (16 * ck3 - ck4) * (16 * (cl + 4 * cl2 - cl3) + cl4))
							* lpm2
					+ (361391 - 5648 * ck3 + 209 * ck4 - 314608 * cl + 16 * ck * (-19663 + 11504 * cl - 784 * cl2)
							+ 60464 * cl2 - 16 * ck2 * (-3779 + 784 * cl + 211 * cl2) - 5648 * cl3 + 209 * cl4) * mu2)
					/ 4.;

			mat[0][0] = 216
					* (16 * ck * LameMu - ck2 * LameMu + (16 * cl - cl2) * (LameLambda + 2 * LameMu)
							- 15 * (LameLambda + 3 * LameMu));
			mat[0][1] = mat[1][0] = 36 * lpm * (-8 * sk + sk2) * (-8 * sl + sl2);
			mat[1][1] = 216
					* (16 * cl * LameMu - cl2 * LameMu + (16 * ck - ck2) * (LameLambda + 2 * LameMu)
							- 15 * (LameLambda + 3 * LameMu));
			if (fabs(den) > TINY)
			{
				den = 1.0 / den;
				int pos = addr(k,l);
				double reF[2];
				double imF[2];

				reF[0] = Re(_fxhat[pos]);
				reF[1] = Re(_fyhat[pos]);
				imF[0] = Im(_fxhat[pos]);
				imF[1] = Im(_fyhat[pos]);
				Re(_fxhat[pos]) = den * (mat[0][0] * reF[0] + mat[0][1] * reF[1]);
				Re(_fyhat[pos]) = den * (mat[1][0] * reF[0] + mat[1][1] * reF[1]);
				Im(_fxhat[pos]) = den * (mat[0][0] * imF[0] + mat[0][1] * imF[1]);
				Im(_fyhat[pos]) = den * (mat[1][0] * imF[0] + mat[1][1] * imF[1]);
			}
			else
			{
				Re(_fxhat[addr(k,l)]) = 0.0;
				Im(_fxhat[addr(k,l)]) = 0.0;
				Re(_fyhat[addr(k,l)]) = 0.0;
				Im(_fyhat[addr(k,l)]) = 0.0;
			}
		}
	}
	nrm = 1.0 / (nx * ny);
	InverseFourier(_fx, _fxhat, nx, ny);
	InverseFourier(_fy, _fyhat, nx, ny);

	for (int i = 0; i < mu; i++)
	{
		solution->vx[i] = Re(_fx[i]) * nrm + meanX;
		solution->vy[i] = Re(_fy[i]) * nrm + meanY;
	}

}
void FluidCurvatureRegistration::OverdampedPeriodicBiharmonicSolve(VectorArray2D*solution, VectorArray2D*ff,
		double ratio)
{
	int nx, ny, mu;
	double dx2, dy2, dx4, dy4, nrm;
	double meanX = 0.0, meanY = 0.0;

	nx = ff->nx;
	ny = ff->ny;
	dx2 = ff->dx * ff->dx;
	dx4 = dx2 * dx2;
	dy2 = ff->dy * ff->dy;
	dy4 = dy2 * dy2;
	mu = nx * ny;

	for (int i = 0; i < mu; i++)
	{
		meanX += ff->vx[i];
		meanY += ff->vy[i];
	}
	meanX /= mu;
	meanY /= mu;

	for (int i = 0; i < mu; i++)
	{
		Re(_fx[i]) = ff->vx[i] - meanX;
		Re(_fy[i]) = ff->vy[i] - meanY;
		Im(_fx[i]) = Im(_fy[i]) = 0.0;
	}
	Fourier(_fxhat, _fx, nx, ny);
	Fourier(_fyhat, _fy, nx, ny);
#pragma omp parallel for
	for (int l = 0; l < ny; l++)
	{
		double lambda = TwoPi * l / ny;
		double cl1 = cos(lambda);
		double cl2 = cos(2.0 * lambda);
		double cl3 = cos(3.0 * lambda);
		for (int k = 0; k < nx; k++)
		{
			double kappa = TwoPi * k / nx;
			double ck1 = cos(kappa);
			double ck2 = cos(2.0 * kappa);
			double ck3 = cos(3.0 * kappa);

			double l2 = 1.0
					+ ratio
							* (-6 * (-28 + 39 * cl1 - 12 * cl2 + cl3) * dx4
									+ (-15 + 16 * ck1 - ck2) * (-15 + 16 * cl1 - cl2) * dx2 * dy2
									- 6 * (-28 + 39 * ck1 - 12 * ck2 + ck3) * dy4) / (18. * dx4 * dy4);

			if (fabs(l2) > TINY)
			{
				l2 = 1.0 / l2;
				int pos = addr(k,l);
				Re(_fxhat[pos]) *= l2;
				Im(_fxhat[pos]) *= l2;
				Re(_fyhat[pos]) *= l2;
				Im(_fyhat[pos]) *= l2;
			}
			else
			{

				Re(_fxhat[addr(k,l)]) /= TINY;
				Im(_fxhat[addr(k,l)]) /= TINY;
				Re(_fyhat[addr(k,l)]) /= TINY;
				Im(_fyhat[addr(k,l)]) /= TINY;
			}
		}
	}
	nrm = 1.0 / (nx * ny);
	InverseFourier(_fx, _fxhat, nx, ny);
	InverseFourier(_fy, _fyhat, nx, ny);

	for (int i = 0; i < mu; i++)
	{
		solution->vx[i] = Re(_fx[i]) * nrm + meanX;
		solution->vy[i] = Re(_fy[i]) * nrm + meanY;
	}

}
void FluidCurvatureRegistration::OverdampedZeroDerivativeBiharmonicSolve(VectorArray2D*solution, VectorArray2D*ff,
		double ratio)
{
	double*fxhat, *fyhat, *fx, *fy, omegak, omegal;
	int nx, ny, mu;
	double meanX, meanY;
	double dx, dy, dx2, dy2, dx4, dy4;

	nx = ff->nx;
	ny = ff->ny;
	dx = ff->dx;
	dy = ff->dy;
	omegak = Pi / nx;
	omegal = Pi / ny;
	dx2 = dx * dx;
	dx4 = dx2 * dx2;
	dy2 = dy * dy;
	dy4 = dy2 * dy2;
	mu = nx * ny;
	fx = new double[mu];
	fy = new double[mu];

	fxhat = new double[mu];
	fyhat = new double[mu];
	meanX = meanY = 0.0;
	for (int i = 0; i < mu; i++)
	{
		fx[i] = ff->vx[i];
		fy[i] = ff->vy[i];
		meanX += fx[i];
		meanY += fy[i];
	}
	meanX /= mu;
	meanY /= mu;
	for (int i = 0; i < mu; i++)
	{
		fx[i] -= meanX;
		fy[i] -= meanY;
	}
	CosFourier(fxhat, fx, nx, ny);
	CosFourier(fyhat, fy, nx, ny);
#pragma omp parallel for
	for (int l = 0; l < ny; l++)
	{
		double cl1 = cos((0.5 + l) * omegal);
		double cl2 = cos((1 + 2 * l) * omegal);
		double cl3 = cos((1.5 + 3 * l) * omegal);
		for (int k = 0; k < nx; k++)
		{
			double ck1 = cos((0.5 + k) * omegak);
			double ck2 = cos((1 + 2 * k) * omegak);
			double ck3 = cos((1.5 + 3 * k) * omegak);

			double l2 = 1.0
					+ ratio
							* (-6 * (-28 + 39 * cl1 - 12 * cl2 + cl3) * dx4
									+ (-15 + 16 * ck1 - ck2) * (-15 + 16 * cl1 - cl2) * dx2 * dy2
									- 6 * (-28 + 39 * ck1 - 12 * ck2 + ck3) * dy4) / (18. * dx4 * dy4);
			if (fabs(l2) > TINY)
			{
				l2 = 1.0 / l2;
				int pos = addr(k,l);
				fxhat[pos] *= l2;
				fyhat[pos] *= l2;
			}
			else
			{
				fxhat[addr(k,l)] /= TINY;
				fyhat[addr(k,l)] /= TINY;
			}
		}
	}

	InverseCosFourier(solution->vx, fxhat, nx, ny);
	InverseCosFourier(solution->vy, fyhat, nx, ny);
	for (int i = 0; i < mu; i++)
	{
		solution->vx[i] += meanX;
		solution->vy[i] += meanY;
	}
	delete[] fyhat;
	delete[] fxhat;
	delete[] fy;
	delete[] fx;
}
void FluidCurvatureRegistration::OverdampedZeroBiharmonicSolve(VectorArray2D*solution, VectorArray2D*ff, double ratio)
{
	double *fxhat, *fyhat, *fx, *fy, omegak, omegal;
	int nx, ny, mu;
	double meanX, meanY;
	double dx, dy, dx2, dy2, dx4, dy4;

	nx = ff->nx;
	ny = ff->ny;
	dx = ff->dx;
	dy = ff->dy;
	omegak = Pi / nx;
	omegal = Pi / ny;
	dx2 = dx * dx;
	dx4 = dx2 * dx2;
	dy2 = dy * dy;
	dy4 = dy2 * dy2;
	mu = nx * ny;
	fx = new double[mu];
	fy = new double[mu];

	fxhat = new double[mu];
	fyhat = new double[mu];
	meanX = meanY = 0.0;
	for (int i = 0; i < mu; i++)
	{
		fx[i] = ff->vx[i];
		fy[i] = ff->vy[i];
		meanX += fx[i];
		meanY += fy[i];
	}
	meanX /= mu;
	meanY /= mu;
	for (int i = 0; i < mu; i++)
	{
		fx[i] -= meanX;
		fy[i] -= meanY;
	}
	SinFourier(fxhat, fx, nx, ny);
	SinFourier(fyhat, fy, nx, ny);
#pragma omp parallel for
	for (int l = 0; l < ny; l++)
	{
		double cl1 = cos((0.5 + l) * omegal);
		double cl2 = cos((1 + 2 * l) * omegal);
		double cl3 = cos((1.5 + 3 * l) * omegal);
		for (int k = 0; k < nx; k++)
		{
			double ck1 = cos((0.5 + k) * omegak);
			double ck2 = cos((1 + 2 * k) * omegak);
			double ck3 = cos((1.5 + 3 * k) * omegak);

			double l2 = 1.0
					+ ratio
							* (-6 * (-28 + 39 * cl1 - 12 * cl2 + cl3) * dx4
									+ (-15 + 16 * ck1 - ck2) * (-15 + 16 * cl1 - cl2) * dx2 * dy2
									- 6 * (-28 + 39 * ck1 - 12 * ck2 + ck3) * dy4) / (18. * dx4 * dy4);

			if (fabs(l2) > TINY)
			{
				l2 = 1.0 / l2;
				int pos = addr(k,l);
				fxhat[pos] *= l2;
				fyhat[pos] *= l2;
			}
			else
			{
				fxhat[addr(k,l)] /= TINY;
				fyhat[addr(k,l)] /= TINY;
			}
		}
	}

	InverseSinFourier(solution->vx, fxhat, nx, ny);
	InverseSinFourier(solution->vy, fyhat, nx, ny);
	for (int i = 0; i < mu; i++)
	{
		solution->vx[i] += meanX;
		solution->vy[i] += meanY;
	}
	delete[] fyhat;
	delete[] fxhat;
	delete[] fy;
	delete[] fx;
}
void FluidCurvatureRegistration::OverdampedPeriodicHarmonicSolve(VectorArray2D*solution, VectorArray2D*ff, double ratio)
{
	int nx, ny, mu;
	double dx2, dy2, nrm;
	double meanX = 0.0, meanY = 0.0;

	nx = ff->nx;
	ny = ff->ny;
	dx2 = ff->dx * ff->dx;
	dy2 = ff->dy * ff->dy;
	mu = nx * ny;

	for (int i = 0; i < mu; i++)
	{
		meanX += ff->vx[i];
		meanY += ff->vy[i];
	}
	meanX /= mu;
	meanY /= mu;

	for (int i = 0; i < mu; i++)
	{
		Re(_fx[i]) = ff->vx[i] - meanX;
		Re(_fy[i]) = ff->vy[i] - meanY;
		Im(_fx[i]) = Im(_fy[i]) = 0.0;
	}
	Fourier(_fxhat, _fx, nx, ny);
	Fourier(_fyhat, _fy, nx, ny);
#pragma omp parallel for
	for (int l = 0; l < ny; l++)
	{
		double lambda = TwoPi * l / ny;
		double cl1 = cos(lambda);
		double cl2 = cos(2.0 * lambda);
		for (int k = 0; k < nx; k++)
		{
			double kappa = TwoPi * k / nx;
			double ck1 = cos(kappa);
			double ck2 = cos(2.0 * kappa);

			double l2 = 1.0 + ratio * ((15 - 16 * cl1 + cl2) * dx2 + (15 - 16 * ck1 + ck2) * dy2) / (6. * dx2 * dy2);
			if (fabs(l2) > TINY)
			{
				l2 = 1.0 / l2;
				int pos = addr(k,l);
				Re(_fxhat[pos]) *= l2;
				Im(_fxhat[pos]) *= l2;
				Re(_fyhat[pos]) *= l2;
				Im(_fyhat[pos]) *= l2;
			}
			else
			{

				Re(_fxhat[addr(k,l)]) /= TINY;
				Im(_fxhat[addr(k,l)]) /= TINY;
				Re(_fyhat[addr(k,l)]) /= TINY;
				Im(_fyhat[addr(k,l)]) /= TINY;
			}
		}
	}
	nrm = 1.0 / (nx * ny);
	InverseFourier(_fx, _fxhat, nx, ny);
	InverseFourier(_fy, _fyhat, nx, ny);

	for (int i = 0; i < mu; i++)
	{
		solution->vx[i] = Re(_fx[i]) * nrm + meanX;
		solution->vy[i] = Re(_fy[i]) * nrm + meanY;
	}

}
void FluidCurvatureRegistration::OverdampedZeroHarmonicSolve(VectorArray2D*solution, VectorArray2D*ff, double ratio)
{
	double*fxhat, *fyhat, *fx, *fy, omegak, omegal;
	int nx, ny, mu;
	double meanX, meanY;
	double dx, dy, dx2, dy2;

	nx = ff->nx;
	ny = ff->ny;
	dx = ff->dx;
	dy = ff->dy;
	omegak = Pi / nx;
	omegal = Pi / ny;
	dx2 = dx * dx;
	dy2 = dy * dy;
	mu = nx * ny;
	fx = new double[mu];
	fy = new double[mu];

	fxhat = new double[mu];
	fyhat = new double[mu];
	meanX = meanY = 0.0;
	for (int i = 0; i < mu; i++)
	{
		fx[i] = ff->vx[i];
		fy[i] = ff->vy[i];
		meanX += fx[i];
		meanY += fy[i];
	}
	meanX /= mu;
	meanY /= mu;
	for (int i = 0; i < mu; i++)
	{
		fx[i] -= meanX;
		fy[i] -= meanY;
	}
	SinFourier(fxhat, fx, nx, ny);
	SinFourier(fyhat, fy, nx, ny);
#pragma omp parallel for
	for (int l = 0; l < ny; l++)
	{
		double cl1 = cos((0.5 + l) * omegal);
		double cl2 = cos((1 + 2 * l) * omegal);
		for (int k = 0; k < nx; k++)
		{
			double ck1 = cos((0.5 + k) * omegak);
			double ck2 = cos((1 + 2 * k) * omegak);

			double l2 = 1.0 - ratio * (-((15 - 16 * ck1 + ck2) / dx2) - (15 - 16 * cl1 + cl2) / dy2) / 6.;

			if (fabs(l2) > TINY)
			{
				l2 = 1.0 / l2;
				int pos = addr(k,l);
				fxhat[pos] *= l2;
				fyhat[pos] *= l2;
			}
			else
			{
				fxhat[addr(k,l)] /= TINY;
				fyhat[addr(k,l)] /= TINY;
			}
		}
	}

	InverseSinFourier(solution->vx, fxhat, nx, ny);
	InverseSinFourier(solution->vy, fyhat, nx, ny);
	for (int i = 0; i < mu; i++)
	{
		solution->vx[i] += meanX;
		solution->vy[i] += meanY;
	}
	delete[] fyhat;
	delete[] fxhat;
	delete[] fy;
	delete[] fx;
}

void FluidCurvatureRegistration::FluidLame(double, double*y, double*dy, int n)
{
	int nx, ny;
	VectorArray2D*u, *f;
	Image<double> *wraped;
	int size;

	nx = templateImage->nx;
	ny = templateImage->ny;
	size = nx * ny;
	u = (VectorArray2D*) malloc(sizeof(VectorArray2D));
	u->nx = nx;
	u->ny = ny;
	u->dx = 1.0;
	u->dy = 1.0;
	u->vx = y;
	u->vy = y + size;

	f = (VectorArray2D*) malloc(sizeof(VectorArray2D));
	f->nx = nx;
	f->ny = ny;
	f->dx = 1.0;
	f->dy = 1.0;
	f->vx = dy;
	f->vy = dy + size;

	wraped = new Image<double>(nx, ny, 1);
	sampleImage->wrap(wraped, u);
#pragma omp parallel for
	for (int j = 0; j < ny; j++)
		for (int i = 0; i < nx; i++)
		{
			double imagediff = templateImage->get(i, j, 0) - wraped->get(i, j, 0);
			Vector2D gradf;

			gradf.x = imagediff * wraped->dx(i, j, 0);
			gradf.y = imagediff * wraped->dy(i, j, 0);

			f->set(i, j, gradf);
		}
	(this->*LameSolve)(f, f);
	delete wraped;
	free(u);

	free(f);
}
void FluidCurvatureRegistration::FluidOverDamped(double, double*y, double*dy, int n)
{
	int nx, ny;
	VectorArray2D *u, *f;

	int size;

	nx = templateImage->nx;
	ny = templateImage->ny;
	size = nx * ny;
	u = (VectorArray2D*) malloc(sizeof(VectorArray2D));
	u->nx = nx;
	u->ny = ny;
	u->dx = 1.0;
	u->dy = 1.0;
	u->vx = y;
	u->vy = y + size;

	f = (VectorArray2D*) malloc(sizeof(VectorArray2D));
	f->nx = nx;
	f->ny = ny;
	f->dx = 1.0;
	f->dy = 1.0;
	f->vx = dy;
	f->vy = dy + size;

	sampleImage->wrap(__wraped, u);
#pragma omp parallel for
	for (int j = 0; j < ny; j++)
		for (int i = 0; i < nx; i++)
		{
			double imagediff = templateImage->get(i, j, 0) - __wraped->get(i, j, 0);
			double tmp = 64.0;
			Vector2D gradf;

			gradf.x = -tmp * imagediff * __wraped->dx(i, j, 0);
			gradf.y = -tmp * imagediff * __wraped->dy(i, j, 0);

			f->set(i, j, gradf);
		}
	if (VortexWeight > 0.0)
	{
#pragma omp parallel for
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
			{
				Vector2D rhs, d2uxx, d2uyy, d2uxy;

				rhs = f->get(i, j);

				d2uxx = u->d2x(i, j);
				d2uyy = u->d2y(i, j);
				d2uxy = u->dxy(i, j);

				rhs.x -= VortexWeight * (d2uxy.y - d2uyy.x);
				rhs.y -= VortexWeight * (d2uxy.x - d2uxx.y);
				f->set(i, j, rhs);
			}

	}
	(this->*OverdampedLameSolve)(f, f, DampingRatio);

	free(u);
	free(f);
}
void FluidCurvatureRegistration::FluidBiharmonic(double*d2y, double, double*y, double*dy, int n)
{
	int nx, ny;
	VectorArray2D*u, *du, *f;

	int size;

	nx = templateImage->nx;
	ny = templateImage->ny;
	size = nx * ny;

	u = new VectorArray2D(nx, ny, 1., 1.);
	u->vx = y;
	u->vy = y + size;

	du = new VectorArray2D(nx, ny, 1., 1.);
	du->vx = dy;
	du->vy = dy + size;

	f = new VectorArray2D(nx, ny, 1., 1.);
	f->vx = d2y;
	f->vy = d2y + size;

	sampleImage->wrap(__wraped, u);
#pragma omp parallel for
	for (int j = 0; j < ny; j++)
		for (int i = 0; i < nx; i++)
		{
			double imagediff = templateImage->get(i, j, 0) - __wraped->get(i, j, 0);
			double tmp = InverseAlpha * imagediff;
			Vector2D gradf;

			gradf.x = -tmp * __wraped->dx(i, j, 0) - BiharmonicFriction * du->getx(i, j);
			gradf.y = -tmp * __wraped->dy(i, j, 0) - BiharmonicFriction * du->gety(i, j);
			f->set(i, j, gradf);
		}
	(this->*BiharmonicSolve)(f, f);

#pragma omp parallel for
	for (int j = 0; j < ny; j++)
		for (int i = 0; i < nx; i++)
		{
			Vector2D force;
			double vx, vy;
			vx = du->getx(i, j);
			vy = du->gety(i, j);

			force = f->get(i, j);
			force.x -= Friction * vx;
			force.y -= Friction * vy;
			f->set(i, j, force);
		}

	delete u;
	delete du;
	delete f;
}
int FluidCurvatureRegistration::PrintFluidProgress(double t, double h, double*u, double*up, double*uNext, double*upNext,
		int n)
{
	static int count = 0L;
	int size, nx, ny;
	VectorArray2D *u1, *u2;

	double maxChange, imgDiff;

	if ((count++) % 16)
		return RKN_OUTPUT_OK;
	nx = templateImage->nx;
	ny = templateImage->ny;
	size = nx * ny;
	u1 = new VectorArray2D();
	u1->nx = nx;
	u1->ny = ny;
	u1->dx = u1->dy = 1.0;
	u1->vx = u;
	u1->vy = u + size;
	u2 = new VectorArray2D();
	u2->nx = nx;
	u2->ny = ny;
	u2->dx = u2->dy = 1.0;
	u2->vx = uNext;
	u2->vy = uNext + size;
	maxChange = 0.0;
	for (int i = 0; i < n; i++)
	{
		double diff;
		diff = fabs(uNext[i] - u[i]);
		if (diff > maxChange)
			maxChange = diff;
	}

	sampleImage->wrap(__wraped, u2);
	imgDiff = 0.0;

#pragma omp parallel for reduction(+:imgDiff)
	for (int j = 0; j < ny; j++)
	{
		double rsum = 0.0;
		for (int i = 0; i < nx; i++)
		{
			double tmp = templateImage->get(i, j, 0) - __wraped->get(i, j, 0);
			rsum += tmp * tmp;
		}
		imgDiff += rsum;
	}

	imgDiff = sqrt(imgDiff);
	imgDiff /= size;

#if (defined(__WIN32__)||defined(WIN32)) && !defined(__CONSOLE__)
	{
		char title[512];
		sprintf(title,"t:%2.4g (%2.4g)",t,imgDiff);
		SetWindowText(
				MLIconWindow,
				title
		);
	}
#elif defined(__CONSOLE__)
	printf("max du=%lg   [%lg] (%lg) %lg\n",maxChange,t,h,imgDiff);
#endif

	u1->vx = u1->vy = NULL;
	u2->vx = u2->vy = NULL;
	delete u1;
	delete u2;
	if (ImageMatchGoal > imgDiff)
		return RKN_OUTPUT_HALT;
	return RKN_OUTPUT_OK;
}

int FluidCurvatureRegistration::PrintFluidProgress1(double t, double h, double*u, double*uNext, int n)
{

	static int count = 0L;
	int size, nx, ny;
	VectorArray2D*u1, *u2;

	double maxChange, imgDiff;

	if (0.0 == h)
	{
		_lowestError = ImageMatchError = MAXDOUBLE;
		MinimumTime = 0.0;
		return RK_OK;
	}

	nx = templateImage->nx;
	ny = templateImage->ny;
	size = nx * ny;
	u1 = new VectorArray2D();
	u1->nx = nx;
	u1->ny = ny;
	u1->dx = u1->dy = 1.0;
	u1->vx = u;
	u1->vy = u + size;
	u2 = new VectorArray2D();
	u2->nx = nx;
	u2->ny = ny;
	u2->dx = u2->dy = 1.0;
	u2->vx = uNext;
	u2->vy = uNext + size;
	maxChange = 0.0;
	for (int i = 0; i < n; i++)
	{
		double diff;
		diff = fabs(uNext[i] - u[i]);
		if (diff > maxChange)
			maxChange = diff;
	}

	sampleImage->wrap(__wraped, u2);
	imgDiff = 0.0;

#pragma omp parallel for reduction(+:imgDiff)
	for (int j = 0; j < ny; j++)
	{
		double rsum = 0.0;
		for (int i = 0; i < nx; i++)
		{
			double tmp = templateImage->get(i, j, 0) - __wraped->get(i, j, 0);
			rsum += tmp * tmp;
		}
		imgDiff += rsum;
	}

	imgDiff = sqrt(imgDiff);
	imgDiff /= size;

	if(returnType)
	{
		const cimg_library::CImg<double> img(__wraped->bm, __wraped->nx, __wraped->ny);
		img.display(display.wait(100));
	}

	u1->vx = u1->vy = NULL;
	u2->vx = u2->vy = NULL;
	delete u1;
	delete u2;
	if (ImageMatchGoal > imgDiff)
	{
		return RK_USER_STOP;
	}
	if (imgDiff > ImageMatchError && ImageMatchError < 1.1 * _lowestError)
	{
		double xmin = t - h, fmin = ImageMatchError;

		ImageDifference *idiff = new ImageDifference(rk43, templateImage, sampleImage, __wraped, nx, ny, t, h);
		Brent *search1d = new Brent();
		search1d->bracket(t - h, t, idiff);

		xmin = search1d->minimize(idiff);

		if (_lowestError > search1d->fmin)
		{
			_lowestError = fmin = search1d->fmin;
			MinimumTime = xmin;
			rk43->InterpolateRKV4(t, h, xmin, _bestU->vx);
		}
		delete idiff;
		delete search1d;
	}

	ImageMatchError = imgDiff;

	if ((count++) % 4)
	{

#if (defined(__WIN32__)||defined(WIN32)) && !defined(__CONSOLE__)
		{
			char title[512];
			sprintf(title,"t:%2.4g (%2.4g)",t,imgDiff);
			SetWindowText(
					MLIconWindow,
					title
			);
		}
#elif defined(__CONSOLE__)
		printf("max du=%lg   [%lg] (%lg) %lg\n",maxChange,t,h,imgDiff);
#endif

	}
	return RK_OK;
}
VectorArray2D* FluidCurvatureRegistration::getFlowField() const
{
	return u;
}
Image<double>* FluidCurvatureRegistration::getReference() const
{
	return ref;
}
Image<double>* FluidCurvatureRegistration::getSample() const
{
	return sampleImage;
}
double FluidCurvatureRegistration::getMinimalTime() const
{
	return MinimumTime;
}
double FluidCurvatureRegistration::getMismatchError() const
{
	return ImageMatchError;
}
