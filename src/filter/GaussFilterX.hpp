/*
 * GaussFilterX.hpp
 *
 *  Created on: Jul 11, 2012
 *      Author: kthierbach
 */

#ifndef GAUSSFILTERX_HPP_
#define GAUSSFILTERX_HPP_

template<class T>
class GaussFilterX
{
	public:
		GaussFilterX(Image<T> *_outimg, const Image<T> *_img, double _C, const double *_c);
		void operator()(const tbb::blocked_range<int> &r) const;
	private:
		int nx, channelNo;
		Image<T> const*img;
		Image<T> *outimg;
		const double*c;
		const double C;

};

#endif /* GAUSSFILTERX_HPP_ */
