/*
 * SL_EssentialMatWrapper.cpp
 *
 *  Created on: 2011-7-2
 *      Author: Danping Zou
 */

#include "SL_EMatWrapper.h"
#include "SL_Distortion.h"

bool estimateEMat(const double* K1,
		const double* K2,
		const Mat_d& pts1,
		const Mat_d& pts2,
		const Matching& matches,
		double E[9],
		int nRansac,
		double maxEpiErr,
		int minInlierNum) {
	int nm = matches.num;
	Mat_d pts1m(nm, 2), pts2m(nm, 2);
	Mat_d pts1mNorm(nm, 2), pts2mNorm(nm, 2);

	//get the matched points
	for (int i = 0; i < nm; i++) {
		int idx1 = matches[i].idx1;
		int idx2 = matches[i].idx2;

		pts1m.data[2 * i] = pts1.data[2 * idx1];
		pts1m.data[2 * i + 1] = pts1.data[2 * idx1 + 1];

		pts2m.data[2 * i] = pts2[2 * idx2];
		pts2m.data[2 * i + 1] = pts2[2 * idx2 + 1];
	}

	double iK1[9], iK2[9];
	mat33Inv(K1, iK1);
	mat33Inv(K2, iK2);

	//get the normalized coordinates
	normPoints(iK1, nm, pts1m.data, pts1mNorm.data);
	normPoints(iK2, nm, pts2m.data, pts2mNorm.data);

	int inlierNum = findEMatRansac(iK2, iK1, nm, pts2m.data, pts1m.data, pts2mNorm.data, pts1mNorm.data, E, nRansac,
			maxEpiErr);

	if (inlierNum < minInlierNum)
		return false;
	return true;
}

bool estimateEMat(const double* K1, const double* K2, const Mat_d& pts1, const Mat_d& pts2, double E[9], int nRansac, double maxEpiErr, int minInlierNum) {
	assert( pts1.m > 0 && pts1.m == pts2.m);
	int nm = pts1.m;

	Mat_d pts1Norm(nm, 2), pts2Norm(nm, 2);

	double iK1[9], iK2[9];
	mat33Inv(K1, iK1);
	mat33Inv(K2, iK2);

	//get the normalized coordinates
	normPoints(iK1, nm, pts1.data, pts1Norm.data);
	normPoints(iK2, nm, pts2.data, pts2Norm.data);

	int inlierNum = findEMatRansac(iK2, iK1, nm, pts1.data, pts2.data, pts1Norm.data, pts2Norm.data, E, nRansac, maxEpiErr);

	if (inlierNum < minInlierNum)
		return false;
	return true;
}
