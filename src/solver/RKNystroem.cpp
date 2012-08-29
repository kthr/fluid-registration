/*
 * RKNystroem.cpp
 *
 *  Created on: Jul 11, 2012
 *      Author: kthierbach
 */

#include "RKNystroem.hpp"

class RKNComputeArg2
{
	public:

		RKNComputeArg2(double*_yarg, double*_yparg, double*_f1, double*_y, double*_yp, double _h) :
				yarg(_yarg), yparg(_yparg), f1(_f1), y(_y), yp(_yp), h(_h)
		{
			h2 = _h * _h;
		}

		void operator()(const tbb::blocked_range<int> &r) const
		{

			static const double c[7] =
			{ 0.0, 8.0 / 39.0, 4.0 / 13.0, 5.0 / 6.0, 43.0 / 47.0, 1., 1. };
			static const double b[7] =
			{ 4817.0 / 51600.0, 0.0, 388869.0 / 1216880.0, 3276.0 / 23575.0, -1142053.0 / 22015140.0, 0.0, 0.0 };
			static const double bp[7] =
			{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53. / 240., 0. };
			static const double deltaB[7] =
			{ 8151. / 2633750., 0., -1377519. / 186334750., 586872. / 28879375., -36011118. / 2247378875., 0, 0 };
			static const double deltaBp[7] =
			{ 8151 / 2633750, 0, -5969249. / 559004250., 3521232. / 28879375., -846261273. / 4494757750., 4187.
					/ 36750., -1. / 25. };
			static const double aMatrix[7][7] =
			{
			{ 0, 0, 0, 0, 0, 0, 0 },
			{ 32. / 1521., 0, 0, 0, 0, 0, 0 },
			{ 4. / 169., 4 / 169, 0, 0, 0, 0, 0 },
			{ 175. / 5184., 0, 1625. / 5184., 0, 0, 0, 0 },
					{ -342497279. / 5618900760., 6827067. / 46824173., 35048741. / 102161832., -2201514. / 234120865.,
							0, 0, 0 },
					{ -7079 / 52152., 767. / 2173., 14027. / 52152., 30. / 2173., 0, 0, 0 },
					{ 4817 / 51600., 0, 388869. / 1216880., 3276. / 23575., -1142053. / 22015140., 0, 0 } };

			static const double apMatrix[7][7] =
					{
					{ 0, 0, 0, 0, 0, 0, 0 },
					{ 8. / 39., 0, 0, 0, 0, 0, 0 },
					{ 1. / 13., 3. / 13., 0, 0, 0, 0, 0 },
					{ 7385. / 6912., -9425. / 2304., 13325. / 3456., 0, 0, 0, 0 },
					{ 223324757. / 91364240., -174255393. / 18272848., 382840094. / 46824173., -39627252. / 234120865.,
							0, 0, 0 },
							{ 108475. / 36464., -9633. / 848., 7624604. / 806183., 8100. / 49979., -4568212.
									/ 19446707., 0, 0 },
							{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53.
									/ 240., 0 } };

			for (int i = r.begin(); i != r.end(); ++i)
			{
				yarg[i] = y[i] + c[1] * h * yp[i] + h2 * aMatrix[1][0] * f1[i];
				yparg[i] = yp[i] + h * apMatrix[1][0] * f1[i];
			}
		}

	private:
		double*yarg;
		double*yparg;
		double* const f1;
		double* const y;
		double* const yp;
		double h, h2;
};
class RKNComputeArg3
{
	public:

		RKNComputeArg3(double*_yarg, double*_yparg, double*_f1, double*_f2, double*_y, double*_yp, double _h) :
				yarg(_yarg), yparg(_yparg), f1(_f1), f2(_f2), y(_y), yp(_yp), h(_h)
		{
			h2 = _h * _h;
		}

		void operator()(const tbb::blocked_range<int> &r) const
		{

			static const double c[7] =
			{ 0.0, 8.0 / 39.0, 4.0 / 13.0, 5.0 / 6.0, 43.0 / 47.0, 1., 1. };
			static const double b[7] =
			{ 4817.0 / 51600.0, 0.0, 388869.0 / 1216880.0, 3276.0 / 23575.0, -1142053.0 / 22015140.0, 0.0, 0.0 };
			static const double bp[7] =
			{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53. / 240., 0. };
			static const double deltaB[7] =
			{ 8151. / 2633750., 0., -1377519. / 186334750., 586872. / 28879375., -36011118. / 2247378875., 0, 0 };
			static const double deltaBp[7] =
			{ 8151 / 2633750, 0, -5969249. / 559004250., 3521232. / 28879375., -846261273. / 4494757750., 4187.
					/ 36750., -1. / 25. };
			static const double aMatrix[7][7] =
			{
			{ 0, 0, 0, 0, 0, 0, 0 },
			{ 32. / 1521., 0, 0, 0, 0, 0, 0 },
			{ 4. / 169., 4 / 169, 0, 0, 0, 0, 0 },
			{ 175. / 5184., 0, 1625. / 5184., 0, 0, 0, 0 },
					{ -342497279. / 5618900760., 6827067. / 46824173., 35048741. / 102161832., -2201514. / 234120865.,
							0, 0, 0 },
					{ -7079 / 52152., 767. / 2173., 14027. / 52152., 30. / 2173., 0, 0, 0 },
					{ 4817 / 51600., 0, 388869. / 1216880., 3276. / 23575., -1142053. / 22015140., 0, 0 } };

			static const double apMatrix[7][7] =
					{
					{ 0, 0, 0, 0, 0, 0, 0 },
					{ 8. / 39., 0, 0, 0, 0, 0, 0 },
					{ 1. / 13., 3. / 13., 0, 0, 0, 0, 0 },
					{ 7385. / 6912., -9425. / 2304., 13325. / 3456., 0, 0, 0, 0 },
					{ 223324757. / 91364240., -174255393. / 18272848., 382840094. / 46824173., -39627252. / 234120865.,
							0, 0, 0 },
							{ 108475. / 36464., -9633. / 848., 7624604. / 806183., 8100. / 49979., -4568212.
									/ 19446707., 0, 0 },
							{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53.
									/ 240., 0 } };

			for (int i = r.begin(); i != r.end(); ++i)
			{
				yarg[i] = y[i] + c[2] * h * yp[i] + h2 * (aMatrix[2][0] * f1[i] + aMatrix[2][1] * f2[i]);
				yparg[i] = yp[i] + h * (apMatrix[2][0] * f1[i] + apMatrix[2][1] * f2[i]);
			}
		}

	private:
		double*yarg;
		double*yparg;
		double* const f1;
		double* const f2;
		double* const y;
		double* const yp;
		double h, h2;
};
class RKNComputeArg4
{
	public:

		RKNComputeArg4(double*_yarg, double*_yparg, double*_f1, double*_f2, double*_f3, double*_y, double*_yp,
				double _h) :
				yarg(_yarg), yparg(_yparg), f1(_f1), f2(_f2), f3(_f3), y(_y), yp(_yp), h(_h)
		{
			h2 = _h * _h;
		}

		void operator()(const tbb::blocked_range<int> &r) const
		{

			static const double c[7] =
			{ 0.0, 8.0 / 39.0, 4.0 / 13.0, 5.0 / 6.0, 43.0 / 47.0, 1., 1. };
			static const double b[7] =
			{ 4817.0 / 51600.0, 0.0, 388869.0 / 1216880.0, 3276.0 / 23575.0, -1142053.0 / 22015140.0, 0.0, 0.0 };
			static const double bp[7] =
			{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53. / 240., 0. };
			static const double deltaB[7] =
			{ 8151. / 2633750., 0., -1377519. / 186334750., 586872. / 28879375., -36011118. / 2247378875., 0, 0 };
			static const double deltaBp[7] =
			{ 8151 / 2633750, 0, -5969249. / 559004250., 3521232. / 28879375., -846261273. / 4494757750., 4187.
					/ 36750., -1. / 25. };
			static const double aMatrix[7][7] =
			{
			{ 0, 0, 0, 0, 0, 0, 0 },
			{ 32. / 1521., 0, 0, 0, 0, 0, 0 },
			{ 4. / 169., 4 / 169, 0, 0, 0, 0, 0 },
			{ 175. / 5184., 0, 1625. / 5184., 0, 0, 0, 0 },
					{ -342497279. / 5618900760., 6827067. / 46824173., 35048741. / 102161832., -2201514. / 234120865.,
							0, 0, 0 },
					{ -7079 / 52152., 767. / 2173., 14027. / 52152., 30. / 2173., 0, 0, 0 },
					{ 4817 / 51600., 0, 388869. / 1216880., 3276. / 23575., -1142053. / 22015140., 0, 0 } };

			static const double apMatrix[7][7] =
					{
					{ 0, 0, 0, 0, 0, 0, 0 },
					{ 8. / 39., 0, 0, 0, 0, 0, 0 },
					{ 1. / 13., 3. / 13., 0, 0, 0, 0, 0 },
					{ 7385. / 6912., -9425. / 2304., 13325. / 3456., 0, 0, 0, 0 },
					{ 223324757. / 91364240., -174255393. / 18272848., 382840094. / 46824173., -39627252. / 234120865.,
							0, 0, 0 },
							{ 108475. / 36464., -9633. / 848., 7624604. / 806183., 8100. / 49979., -4568212.
									/ 19446707., 0, 0 },
							{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53.
									/ 240., 0 } };

			for (int i = r.begin(); i != r.end(); ++i)
			{
				yarg[i] = y[i] + c[3] * h * yp[i]
						+ h2 * (aMatrix[3][0] * f1[i] + aMatrix[3][1] * f2[i] + aMatrix[3][2] * f3[i]);
				yparg[i] = yp[i] + h * (apMatrix[3][0] * f1[i] + apMatrix[3][1] * f2[i] + apMatrix[3][2] * f3[i]);
			}
		}

	private:
		double*yarg;
		double*yparg;
		double* const f1;
		double* const f2;
		double* const f3;
		double* const y;
		double* const yp;
		double h, h2;
};
class RKNComputeArg5
{
	public:

		RKNComputeArg5(double*_yarg, double*_yparg, double*_f1, double*_f2, double*_f3, double*_f4, double*_y,
				double*_yp, double _h) :
				yarg(_yarg), yparg(_yparg), f1(_f1), f2(_f2), f3(_f3), f4(_f4), y(_y), yp(_yp), h(_h)
		{
			h2 = _h * _h;
		}

		void operator()(const tbb::blocked_range<int> &r) const
		{
			static const double c[7] =
			{ 0.0, 8.0 / 39.0, 4.0 / 13.0, 5.0 / 6.0, 43.0 / 47.0, 1., 1. };
			static const double b[7] =
			{ 4817.0 / 51600.0, 0.0, 388869.0 / 1216880.0, 3276.0 / 23575.0, -1142053.0 / 22015140.0, 0.0, 0.0 };
			static const double bp[7] =
			{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53. / 240., 0. };
			static const double deltaB[7] =
			{ 8151. / 2633750., 0., -1377519. / 186334750., 586872. / 28879375., -36011118. / 2247378875., 0, 0 };
			static const double deltaBp[7] =
			{ 8151 / 2633750, 0, -5969249. / 559004250., 3521232. / 28879375., -846261273. / 4494757750., 4187.
					/ 36750., -1. / 25. };
			static const double aMatrix[7][7] =
			{
			{ 0, 0, 0, 0, 0, 0, 0 },
			{ 32. / 1521., 0, 0, 0, 0, 0, 0 },
			{ 4. / 169., 4 / 169, 0, 0, 0, 0, 0 },
			{ 175. / 5184., 0, 1625. / 5184., 0, 0, 0, 0 },
					{ -342497279. / 5618900760., 6827067. / 46824173., 35048741. / 102161832., -2201514. / 234120865.,
							0, 0, 0 },
					{ -7079 / 52152., 767. / 2173., 14027. / 52152., 30. / 2173., 0, 0, 0 },
					{ 4817 / 51600., 0, 388869. / 1216880., 3276. / 23575., -1142053. / 22015140., 0, 0 } };

			static const double apMatrix[7][7] =
					{
					{ 0, 0, 0, 0, 0, 0, 0 },
					{ 8. / 39., 0, 0, 0, 0, 0, 0 },
					{ 1. / 13., 3. / 13., 0, 0, 0, 0, 0 },
					{ 7385. / 6912., -9425. / 2304., 13325. / 3456., 0, 0, 0, 0 },
					{ 223324757. / 91364240., -174255393. / 18272848., 382840094. / 46824173., -39627252. / 234120865.,
							0, 0, 0 },
							{ 108475. / 36464., -9633. / 848., 7624604. / 806183., 8100. / 49979., -4568212.
									/ 19446707., 0, 0 },
							{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53.
									/ 240., 0 } };

			for (int i = r.begin(); i != r.end(); ++i)
			{
				yarg[i] = y[i] + c[4] * h * yp[i]
						+ h2
								* (aMatrix[4][0] * f1[i] + aMatrix[4][1] * f2[i] + aMatrix[4][2] * f3[i]
										+ aMatrix[4][3] * f4[i]);
				yparg[i] = yp[i]
						+ h
								* (apMatrix[4][0] * f1[i] + apMatrix[4][1] * f2[i] + apMatrix[4][2] * f3[i]
										+ apMatrix[4][3] * f4[i]);
			}
		}

	private:
		double*yarg;
		double*yparg;
		double* const f1;
		double* const f2;
		double* const f3;
		double* const f4;
		double* const y;
		double* const yp;
		double h, h2;
};
class RKNComputeArg6
{
	public:

		RKNComputeArg6(double*_yarg, double*_yparg, double*_f1, double*_f2, double*_f3, double*_f4, double*_f5,
				double*_y, double*_yp, double _h) :
				yarg(_yarg), yparg(_yparg), f1(_f1), f2(_f2), f3(_f3), f4(_f4), f5(_f5), y(_y), yp(_yp), h(_h)
		{
			h2 = _h * _h;
		}

		void operator()(const tbb::blocked_range<int> &r) const
		{
			static const double c[7] =
			{ 0.0, 8.0 / 39.0, 4.0 / 13.0, 5.0 / 6.0, 43.0 / 47.0, 1., 1. };
			static const double b[7] =
			{ 4817.0 / 51600.0, 0.0, 388869.0 / 1216880.0, 3276.0 / 23575.0, -1142053.0 / 22015140.0, 0.0, 0.0 };
			static const double bp[7] =
			{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53. / 240., 0. };
			static const double deltaB[7] =
			{ 8151. / 2633750., 0., -1377519. / 186334750., 586872. / 28879375., -36011118. / 2247378875., 0, 0 };
			static const double deltaBp[7] =
			{ 8151 / 2633750, 0, -5969249. / 559004250., 3521232. / 28879375., -846261273. / 4494757750., 4187.
					/ 36750., -1. / 25. };
			static const double aMatrix[7][7] =
			{
			{ 0, 0, 0, 0, 0, 0, 0 },
			{ 32. / 1521., 0, 0, 0, 0, 0, 0 },
			{ 4. / 169., 4 / 169, 0, 0, 0, 0, 0 },
			{ 175. / 5184., 0, 1625. / 5184., 0, 0, 0, 0 },
					{ -342497279. / 5618900760., 6827067. / 46824173., 35048741. / 102161832., -2201514. / 234120865.,
							0, 0, 0 },
					{ -7079 / 52152., 767. / 2173., 14027. / 52152., 30. / 2173., 0, 0, 0 },
					{ 4817 / 51600., 0, 388869. / 1216880., 3276. / 23575., -1142053. / 22015140., 0, 0 } };

			static const double apMatrix[7][7] =
					{
					{ 0, 0, 0, 0, 0, 0, 0 },
					{ 8. / 39., 0, 0, 0, 0, 0, 0 },
					{ 1. / 13., 3. / 13., 0, 0, 0, 0, 0 },
					{ 7385. / 6912., -9425. / 2304., 13325. / 3456., 0, 0, 0, 0 },
					{ 223324757. / 91364240., -174255393. / 18272848., 382840094. / 46824173., -39627252. / 234120865.,
							0, 0, 0 },
							{ 108475. / 36464., -9633. / 848., 7624604. / 806183., 8100. / 49979., -4568212.
									/ 19446707., 0, 0 },
							{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53.
									/ 240., 0 } };

			for (int i = r.begin(); i != r.end(); ++i)
			{
				yarg[i] = y[i] + c[5] * h * yp[i]
						+ h2
								* (aMatrix[5][0] * f1[i] + aMatrix[5][1] * f2[i] + aMatrix[5][2] * f3[i]
										+ aMatrix[5][3] * f4[i] + aMatrix[5][4] * f5[i]);
				yparg[i] = yp[i]
						+ h
								* (apMatrix[5][0] * f1[i] + apMatrix[5][1] * f2[i] + apMatrix[5][2] * f3[i]
										+ apMatrix[5][3] * f4[i] + apMatrix[5][4] * f5[i]);
			}
		}

	private:

		double*yarg;
		double*yparg;
		double* const f1;
		double* const f2;
		double* const f3;
		double* const f4;
		double* const f5;
		double* const y;
		double* const yp;
		double h, h2;
};
class RKNComputeArg7
{
	public:

		RKNComputeArg7(double*_yarg, double*_yparg, double*_f1, double*_f2, double*_f3, double*_f4, double*_f5,
				double*_f6, double*_y, double*_yp, double _h) :
				yarg(_yarg), yparg(_yparg), f1(_f1), f2(_f2), f3(_f3), f4(_f4), f5(_f5), f6(_f6), y(_y), yp(_yp), h(_h)
		{
			h2 = _h * _h;
		}

		void operator()(const tbb::blocked_range<int> &r) const
		{
			static const double c[7] =
			{ 0.0, 8.0 / 39.0, 4.0 / 13.0, 5.0 / 6.0, 43.0 / 47.0, 1., 1. };
			static const double b[7] =
			{ 4817.0 / 51600.0, 0.0, 388869.0 / 1216880.0, 3276.0 / 23575.0, -1142053.0 / 22015140.0, 0.0, 0.0 };
			static const double bp[7] =
			{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53. / 240., 0. };
			static const double deltaB[7] =
			{ 8151. / 2633750., 0., -1377519. / 186334750., 586872. / 28879375., -36011118. / 2247378875., 0, 0 };
			static const double deltaBp[7] =
			{ 8151 / 2633750, 0, -5969249. / 559004250., 3521232. / 28879375., -846261273. / 4494757750., 4187.
					/ 36750., -1. / 25. };
			static const double aMatrix[7][7] =
			{
			{ 0, 0, 0, 0, 0, 0, 0 },
			{ 32. / 1521., 0, 0, 0, 0, 0, 0 },
			{ 4. / 169., 4 / 169, 0, 0, 0, 0, 0 },
			{ 175. / 5184., 0, 1625. / 5184., 0, 0, 0, 0 },
					{ -342497279. / 5618900760., 6827067. / 46824173., 35048741. / 102161832., -2201514. / 234120865.,
							0, 0, 0 },
					{ -7079 / 52152., 767. / 2173., 14027. / 52152., 30. / 2173., 0, 0, 0 },
					{ 4817 / 51600., 0, 388869. / 1216880., 3276. / 23575., -1142053. / 22015140., 0, 0 } };

			static const double apMatrix[7][7] =
					{
					{ 0, 0, 0, 0, 0, 0, 0 },
					{ 8. / 39., 0, 0, 0, 0, 0, 0 },
					{ 1. / 13., 3. / 13., 0, 0, 0, 0, 0 },
					{ 7385. / 6912., -9425. / 2304., 13325. / 3456., 0, 0, 0, 0 },
					{ 223324757. / 91364240., -174255393. / 18272848., 382840094. / 46824173., -39627252. / 234120865.,
							0, 0, 0 },
							{ 108475. / 36464., -9633. / 848., 7624604. / 806183., 8100. / 49979., -4568212.
									/ 19446707., 0, 0 },
							{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53.
									/ 240., 0 } };

			for (int i = r.begin(); i != r.end(); ++i)
			{
				yarg[i] = y[i] + c[6] * h * yp[i]
						+ h2
								* (aMatrix[6][0] * f1[i] + aMatrix[6][2] * f3[i] + aMatrix[6][3] * f4[i]
										+ aMatrix[6][4] * f5[i] + aMatrix[6][5] * f6[i]);
				yparg[i] = yp[i]
						+ h
								* (apMatrix[6][0] * f1[i] + apMatrix[6][2] * f3[i] + apMatrix[6][3] * f4[i]
										+ apMatrix[6][4] * f5[i] + apMatrix[6][5] * f6[i]);
			}
		}

	private:
		double*yarg;
		double*yparg;
		double* const f1;
		double* const f2;
		double* const f3;
		double* const f4;
		double* const f5;
		double* const f6;
		double* const y;
		double* const yp;
		double h, h2;
};
class RKNComputeNext
{
	public:
		double erry, erryp;
		RKNComputeNext(double*_yNew, double*_ypNew, double*_f1, double*_f2, double*_f3, double*_f4, double*_f5,
				double*_f6, double*_f7, double*_y, double*_yp, double _h) :
				yNew(_yNew), ypNew(_ypNew), f1(_f1), f2(_f2), f3(_f3), f4(_f4), f5(_f5), f6(_f6), f7(_f7), y(_y), yp(
						_yp), h(_h), erry(0.0), erryp(0.0)
		{
			h2 = _h * _h;
		}
		RKNComputeNext(RKNComputeNext&b, tbb::split) :
				yNew(b.yNew), ypNew(b.ypNew), f1(b.f1), f2(b.f2), f3(b.f3), f4(b.f4), f5(b.f5), f6(b.f6), f7(b.f7), y(
						b.y), yp(b.yp), h(b.h), h2(b.h2), erry(0.0), erryp(0.0)
		{
		}

		void operator()(const tbb::blocked_range<int> &r)
		{
			static const double c[7] =
			{ 0.0, 8.0 / 39.0, 4.0 / 13.0, 5.0 / 6.0, 43.0 / 47.0, 1., 1. };
			static const double b[7] =
			{ 4817.0 / 51600.0, 0.0, 388869.0 / 1216880.0, 3276.0 / 23575.0, -1142053.0 / 22015140.0, 0.0, 0.0 };
			static const double bp[7] =
			{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53. / 240., 0. };
			static const double deltaB[7] =
			{ 8151. / 2633750., 0., -1377519. / 186334750., 586872. / 28879375., -36011118. / 2247378875., 0, 0 };
			static const double deltaBp[7] =
			{ 8151 / 2633750, 0, -5969249. / 559004250., 3521232. / 28879375., -846261273. / 4494757750., 4187.
					/ 36750., -1. / 25. };
			static const double aMatrix[7][7] =
			{
			{ 0, 0, 0, 0, 0, 0, 0 },
			{ 32. / 1521., 0, 0, 0, 0, 0, 0 },
			{ 4. / 169., 4 / 169, 0, 0, 0, 0, 0 },
			{ 175. / 5184., 0, 1625. / 5184., 0, 0, 0, 0 },
					{ -342497279. / 5618900760., 6827067. / 46824173., 35048741. / 102161832., -2201514. / 234120865.,
							0, 0, 0 },
					{ -7079 / 52152., 767. / 2173., 14027. / 52152., 30. / 2173., 0, 0, 0 },
					{ 4817 / 51600., 0, 388869. / 1216880., 3276. / 23575., -1142053. / 22015140., 0, 0 } };

			static const double apMatrix[7][7] =
					{
					{ 0, 0, 0, 0, 0, 0, 0 },
					{ 8. / 39., 0, 0, 0, 0, 0, 0 },
					{ 1. / 13., 3. / 13., 0, 0, 0, 0, 0 },
					{ 7385. / 6912., -9425. / 2304., 13325. / 3456., 0, 0, 0, 0 },
					{ 223324757. / 91364240., -174255393. / 18272848., 382840094. / 46824173., -39627252. / 234120865.,
							0, 0, 0 },
							{ 108475. / 36464., -9633. / 848., 7624604. / 806183., 8100. / 49979., -4568212.
									/ 19446707., 0, 0 },
							{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53.
									/ 240., 0 } };

			for (int i = r.begin(); i != r.end(); ++i)
			{
				double tmp1, tmp2;

				yNew[i] = y[i] + h * yp[i] + h2 * (b[0] * f1[i] + b[2] * f3[i] + b[3] * f4[i] + b[4] * f5[i]);
				ypNew[i] = yp[i] + h * (bp[0] * f1[i] + bp[2] * f3[i] + bp[3] * f4[i] + bp[4] * f5[i] + bp[5] * f6[i]);
				tmp1 = deltaB[0] * f1[i] + deltaB[2] * f3[i] + deltaB[3] * f4[i] + deltaB[4] * f5[i] + deltaB[5] * f6[i]
						+ deltaB[6] * f7[i];
				tmp2 = deltaBp[0] * f1[i] + deltaBp[2] * f3[i] + deltaBp[3] * f4[i] + deltaBp[4] * f5[i]
						+ deltaBp[5] * f6[i] + deltaBp[6] * f7[i];
				erry += tmp1 * tmp1;
				erryp += tmp2 * tmp2;
			}
		}
		void join(const RKNComputeNext&b)
		{
			erry += b.erry;
			erryp += b.erryp;
		}
	private:
		double*yNew;
		double*ypNew;
		double* const f1;
		double* const f2;
		double* const f3;
		double* const f4;
		double* const f5;
		double* const f6;
		double* const f7;
		double* const y;
		double* const yp;
		double h, h2;
};

RKNystroemDSolve::RKNystroemDSolve(FluidCurvatureRegistration *fr, int dim,
		void (FluidCurvatureRegistration::*f)(double*, double, double*, double*, int),
		int (FluidCurvatureRegistration::*out)(double, double, double*, double*, double*, double*, int) /*= NULL*/) :
		fr(fr)
{
	n = dim;
	f1 = new double[n];
	f2 = new double[n];
	f3 = new double[n];
	f4 = new double[n];
	f5 = new double[n];
	f6 = new double[n];
	yNew = new double[n];
	ypNew = new double[n];
	yarg = new double[n];
	yparg = new double[n];
	nFunctionCall = nStepCount = nAccept = 0L;
	rhs = f;
	output = out;
}
RKNystroemDSolve::~RKNystroemDSolve(void)
{
	delete[] f1;
	delete[] f2;
	delete[] f3;
	delete[] f4;
	delete[] f5;
	delete[] f6;
	delete[] yNew;
	delete[] ypNew;
	delete[] yarg;
	delete[] yparg;
}

int RKNystroemDSolve::integrate(double t0, double t1, double h, double*y, double*yp, double err, bool parallelQ)
{
	const double saveFactor = 0.9;
	double t, ddt, hnext;
	double localError, hfac, integrationLength;
	int outmessg = RKN_OUTPUT_OK;

	if ((t1 - t0) * h < 0)
		return RKN_DO_NOT_ARRIVE;
	if (n < 2)
		parallelQ = false;
	integrationLength = (double) 1.0 / (t1 - t0);
	if (NULL != output)
		((fr)->*(this->output))(t0, 0.0, y, yp, y, yp, n);
	ddt = fabs(t1 - t0);
	t = t0;
	hnext = h;
	for (;;)
	{
		h = hnext;

		if (fabs(t + 0.05 * h) <= fabs(t) || fabs(h) < RKN_TINY)
		{
			fprintf(stderr, "%s h size to small at x=%lg with h=%lg.", "RKNystroemDSolve::integrate", t, h);
			return RKN_STEP_TOO_SMALL;
		}

		if (fabs(t + h - t0) > ddt)
			h = t1 - t;

		if (parallelQ)
			localError = singleStepTBB(yNew, ypNew, t, h, y, yp);
		else
			localError = singleStep(yNew, ypNew, t, h, y, yp);
		nStepCount++;
		outmessg = RKN_OUTPUT_OK;
		if (localError < err)

		{
			if (NULL != output)
				outmessg = ((fr)->*(this->output))(t, h, y, yp, yNew, ypNew, n);
			for (int i = 0; i < n; i++)
			{
				y[i] = yNew[i];
				yp[i] = ypNew[i];
			}
			t += h;
			nAccept++;
		}

		if (RKN_OUTPUT_HALT == outmessg)
			return RKN_USER_STOP;

		hfac = saveFactor * pow(err / localError, 0.25);
		if (hfac > 5.0)
			hfac = 5.0;
		hnext = hfac * h;

		if (fabs((t - t0) * integrationLength) + RKN_TINY >= 1.0)
			break;
	}
	return RKN_OK;
}

double RKNystroemDSolve::singleStep(double*yNew, double*ypNew, double t, double h, double*y, double*yp)
{

	static const double c[7] =
	{ 0.0, 8.0 / 39.0, 4.0 / 13.0, 5.0 / 6.0, 43.0 / 47.0, 1., 1. };
	static const double b[7] =
	{ 4817.0 / 51600.0, 0.0, 388869.0 / 1216880.0, 3276.0 / 23575.0, -1142053.0 / 22015140.0, 0.0, 0.0 };
	static const double bp[7] =
	{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53. / 240., 0. };
	static const double deltaB[7] =
	{ 8151. / 2633750., 0., -1377519. / 186334750., 586872. / 28879375., -36011118. / 2247378875., 0, 0 };
	static const double deltaBp[7] =
	{ 8151 / 2633750, 0, -5969249. / 559004250., 3521232. / 28879375., -846261273. / 4494757750., 4187. / 36750., -1.
			/ 25. };
	static const double aMatrix[7][7] =
	{
	{ 0, 0, 0, 0, 0, 0, 0 },
	{ 32. / 1521., 0, 0, 0, 0, 0, 0 },
	{ 4. / 169., 4 / 169, 0, 0, 0, 0, 0 },
	{ 175. / 5184., 0, 1625. / 5184., 0, 0, 0, 0 },
	{ -342497279. / 5618900760., 6827067. / 46824173., 35048741. / 102161832., -2201514. / 234120865., 0, 0, 0 },
	{ -7079 / 52152., 767. / 2173., 14027. / 52152., 30. / 2173., 0, 0, 0 },
	{ 4817 / 51600., 0, 388869. / 1216880., 3276. / 23575., -1142053. / 22015140., 0, 0 } };

	static const double apMatrix[7][7] =
	{
	{ 0, 0, 0, 0, 0, 0, 0 },
	{ 8. / 39., 0, 0, 0, 0, 0, 0 },
	{ 1. / 13., 3. / 13., 0, 0, 0, 0, 0 },
	{ 7385. / 6912., -9425. / 2304., 13325. / 3456., 0, 0, 0, 0 },
	{ 223324757. / 91364240., -174255393. / 18272848., 382840094. / 46824173., -39627252. / 234120865., 0, 0, 0 },
	{ 108475. / 36464., -9633. / 848., 7624604. / 806183., 8100. / 49979., -4568212. / 19446707., 0, 0 },
	{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53. / 240., 0 } };

	double h2 = h * h;
	double erry = 0.0, erryp = 0.0;
	double*f7 = f2;

	((fr)->*(this->rhs))(f1, t, y, yp, n);
	for (int i = 0; i < n; i++)
	{
		yarg[i] = y[i] + c[1] * h * yp[i] + h2 * aMatrix[1][0] * f1[i];
		yparg[i] = yp[i] + h * apMatrix[1][0] * f1[i];
	}
	((fr)->*(this->rhs))(f2, t + c[1] * h, yarg, yparg, n);

	for (int i = 0; i < n; i++)
	{
		yarg[i] = y[i] + c[2] * h * yp[i] + h2 * (aMatrix[2][0] * f1[i] + aMatrix[2][1] * f2[i]);
		yparg[i] = yp[i] + h * (apMatrix[2][0] * f1[i] + apMatrix[2][1] * f2[i]);
	}
	((fr)->*(this->rhs))(f3, t + c[2] * h, yarg, yparg, n);

	for (int i = 0; i < n; i++)
	{
		yarg[i] = y[i] + c[3] * h * yp[i]
				+ h2 * (aMatrix[3][0] * f1[i] + aMatrix[3][1] * f2[i] + aMatrix[3][2] * f3[i]);
		yparg[i] = yp[i] + h * (apMatrix[3][0] * f1[i] + apMatrix[3][1] * f2[i] + apMatrix[3][2] * f3[i]);
	}
	((fr)->*(this->rhs))(f4, t + c[3] * h, yarg, yparg, n);

	for (int i = 0; i < n; i++)
	{
		yarg[i] = y[i] + c[4] * h * yp[i]
				+ h2 * (aMatrix[4][0] * f1[i] + aMatrix[4][1] * f2[i] + aMatrix[4][2] * f3[i] + aMatrix[4][3] * f4[i]);
		yparg[i] = yp[i]
				+ h
						* (apMatrix[4][0] * f1[i] + apMatrix[4][1] * f2[i] + apMatrix[4][2] * f3[i]
								+ apMatrix[4][3] * f4[i]);
	}
	((fr)->*(this->rhs))(f5, t + c[4] * h, yarg, yparg, n);

	for (int i = 0; i < n; i++)
	{
		yarg[i] = y[i] + c[5] * h * yp[i]
				+ h2
						* (aMatrix[5][0] * f1[i] + aMatrix[5][1] * f2[i] + aMatrix[5][2] * f3[i] + aMatrix[5][3] * f4[i]
								+ aMatrix[5][4] * f5[i]);
		yparg[i] = yp[i]
				+ h
						* (apMatrix[5][0] * f1[i] + apMatrix[5][1] * f2[i] + apMatrix[5][2] * f3[i]
								+ apMatrix[5][3] * f4[i] + apMatrix[5][4] * f5[i]);
	}
	((fr)->*(this->rhs))(f6, t + c[5] * h, yarg, yparg, n);

	for (int i = 0; i < n; i++)
	{
		yarg[i] = y[i] + c[6] * h * yp[i]
				+ h2
						* (aMatrix[6][0] * f1[i] + aMatrix[6][2] * f3[i] + aMatrix[6][3] * f4[i] + aMatrix[6][4] * f5[i]
								+ aMatrix[6][5] * f6[i]);
		yparg[i] = yp[i]
				+ h
						* (apMatrix[6][0] * f1[i] + apMatrix[6][2] * f3[i] + apMatrix[6][3] * f4[i]
								+ apMatrix[6][4] * f5[i] + apMatrix[6][5] * f6[i]);
	}
	((fr)->*(this->rhs))(f7, t + c[6] * h, yarg, yparg, n);

	for (int i = 0; i < n; i++)
	{
		double tmp1, tmp2;

		yNew[i] = y[i] + h * yp[i] + h2 * (b[0] * f1[i] + b[2] * f3[i] + b[3] * f4[i] + b[4] * f5[i]);
		ypNew[i] = yp[i] + h * (bp[0] * f1[i] + bp[2] * f3[i] + bp[3] * f4[i] + bp[4] * f5[i] + bp[5] * f6[i]);
		tmp1 = deltaB[0] * f1[i] + deltaB[2] * f3[i] + deltaB[3] * f4[i] + deltaB[4] * f5[i] + deltaB[5] * f6[i]
				+ deltaB[6] * f7[i];
		tmp2 = deltaBp[0] * f1[i] + deltaBp[2] * f3[i] + deltaBp[3] * f4[i] + deltaBp[4] * f5[i] + deltaBp[5] * f6[i]
				+ deltaBp[6] * f7[i];
		erry += tmp1 * tmp1;
		erryp += tmp2 * tmp2;
	}
	erry = h * h * sqrt(erry) / n;
	erryp = fabs(h) * sqrt(erryp) / n;

	nFunctionCall += 7;
	if (0.0 == erryp)
		erryp = erry;
	return sqrt(erry * erryp);

}
double RKNystroemDSolve::singleStepTBB(double*yNew, double*ypNew, double t, double h, double*y, double*yp,
		int grainSize /*= 2048*/)
{

	static const double c[7] =
	{ 0.0, 8.0 / 39.0, 4.0 / 13.0, 5.0 / 6.0, 43.0 / 47.0, 1., 1. };
	static const double b[7] =
	{ 4817.0 / 51600.0, 0.0, 388869.0 / 1216880.0, 3276.0 / 23575.0, -1142053.0 / 22015140.0, 0.0, 0.0 };
	static const double bp[7] =
	{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53. / 240., 0. };
	static const double deltaB[7] =
	{ 8151. / 2633750., 0., -1377519. / 186334750., 586872. / 28879375., -36011118. / 2247378875., 0, 0 };
	static const double deltaBp[7] =
	{ 8151 / 2633750, 0, -5969249. / 559004250., 3521232. / 28879375., -846261273. / 4494757750., 4187. / 36750., -1.
			/ 25. };
	static const double aMatrix[7][7] =
	{
	{ 0, 0, 0, 0, 0, 0, 0 },
	{ 32. / 1521., 0, 0, 0, 0, 0, 0 },
	{ 4. / 169., 4 / 169, 0, 0, 0, 0, 0 },
	{ 175. / 5184., 0, 1625. / 5184., 0, 0, 0, 0 },
	{ -342497279. / 5618900760., 6827067. / 46824173., 35048741. / 102161832., -2201514. / 234120865., 0, 0, 0 },
	{ -7079 / 52152., 767. / 2173., 14027. / 52152., 30. / 2173., 0, 0, 0 },
	{ 4817 / 51600., 0, 388869. / 1216880., 3276. / 23575., -1142053. / 22015140., 0, 0 } };

	static const double apMatrix[7][7] =
	{
	{ 0, 0, 0, 0, 0, 0, 0 },
	{ 8. / 39., 0, 0, 0, 0, 0, 0 },
	{ 1. / 13., 3. / 13., 0, 0, 0, 0, 0 },
	{ 7385. / 6912., -9425. / 2304., 13325. / 3456., 0, 0, 0, 0 },
	{ 223324757. / 91364240., -174255393. / 18272848., 382840094. / 46824173., -39627252. / 234120865., 0, 0, 0 },
	{ 108475. / 36464., -9633. / 848., 7624604. / 806183., 8100. / 49979., -4568212. / 19446707., 0, 0 },
	{ 4817. / 51600., 0., 1685099. / 3650640., 19656. / 23575., -53676491. / 88060560., 53. / 240., 0 } };

	double h2 = h * h;
	double erry = 0.0, erryp = 0.0;
	double*f7 = f2;

	if (grainSize > n)
		grainSize = n / 8;
	((fr)->*(this->rhs))(f1, t, y, yp, n);

	tbb::parallel_for(tbb::blocked_range<int>(0, n, grainSize), RKNComputeArg2(yarg, yparg, f1, y, yp, h));
	((fr)->*(this->rhs))(f2, t + c[1] * h, yarg, yparg, n);
	tbb::parallel_for(tbb::blocked_range<int>(0, n, grainSize), RKNComputeArg3(yarg, yparg, f1, f2, y, yp, h));
	((fr)->*(this->rhs))(f3, t + c[2] * h, yarg, yparg, n);
	tbb::parallel_for(tbb::blocked_range<int>(0, n, grainSize), RKNComputeArg4(yarg, yparg, f1, f2, f3, y, yp, h));
	((fr)->*(this->rhs))(f4, t + c[3] * h, yarg, yparg, n);
	tbb::parallel_for(tbb::blocked_range<int>(0, n, grainSize), RKNComputeArg5(yarg, yparg, f1, f2, f3, f4, y, yp, h));
	((fr)->*(this->rhs))(f5, t + c[4] * h, yarg, yparg, n);
	tbb::parallel_for(tbb::blocked_range<int>(0, n, grainSize),
			RKNComputeArg6(yarg, yparg, f1, f2, f3, f4, f5, y, yp, h));

	((fr)->*(this->rhs))(f6, t + c[5] * h, yarg, yparg, n);
	tbb::parallel_for(tbb::blocked_range<int>(0, n, grainSize),
			RKNComputeArg7(yarg, yparg, f1, f2, f3, f4, f5, f6, y, yp, h));
	((fr)->*(this->rhs))(f7, t + c[6] * h, yarg, yparg, n);

	RKNComputeNext next(yNew, ypNew, f1, f2, f3, f4, f5, f6, f7, y, yp, h);

	tbb::parallel_reduce(tbb::blocked_range<int>(0, n, grainSize), next);

	erry = h * h * sqrt(next.erry) / n;
	erryp = fabs(h) * sqrt(next.erryp) / n;

	nFunctionCall += 7;
	if (0.0 == erryp)
		erryp = erry;
	return sqrt(erry * erryp);

}
