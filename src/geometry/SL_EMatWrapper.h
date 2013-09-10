/*
 * SL_EMatWrapper.h
 *
 *  Created on: 2011-7-2
 *      Author: Danping Zou
 */

#ifndef SL_EMATWRAPPER_H_
#define SL_EMATWRAPPER_H_
#include "math/SL_LinAlgWarper.h"
#include "geometry/SL_5point.h"
#include "matching/SL_Matching.h"
bool estimateEMat(
		const double* K1,
		const double* K2,
		const Mat_d& surfPts1,
		const Mat_d& surfPts2,
		const Matching& matches,
		double E[9],
		int nRansac,
		double maxEpiErr,
		int minInlierNum);
bool estimateEMat(
		const double* K1,
		const double* K2,
		const Mat_d& surfPts1,
		const Mat_d& surfPts2,
		double E[9],
		int nRansac,
		double maxEpiErr,
		int minInlierNum);
#endif /* SL_EMATWRAPPER_H_ */
