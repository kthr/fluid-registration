/*
 * GaussFilterY.hpp
 *
 *  Created on: Jul 11, 2012
 *      Author: kthierbach
 */

#ifndef GAUSSFILTERY_HPP_
#define GAUSSFILTERY_HPP_

template<class T>
class GaussFilterY
{
	public:
		GaussFilterY(Image<T> *_outimg, double _C, const double*_c);
		void operator()(const tbb::blocked_range<int> &r) const;
	private:
		const int ny, channelNo;
		Image<T> *outimg;
		const double *c;
		const double C;
};

#endif /* GAUSSFILTERY_HPP_ */
