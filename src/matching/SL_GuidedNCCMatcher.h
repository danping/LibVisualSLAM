/*
 * SL_GuidedNCCMatcher.h
 *
 *  Created on: 2011-3-30
 *      Author: Danping Zou
 */

#ifndef SL_GUIDEDNCCMATCHER_H_
#define SL_GUIDEDNCCMATCHER_H_
#include "SL_Matching.h"
#include "math/SL_Matrix.h"

void getSeedDisparities(const Mat_d& seed1, const Mat_d& seed2, Mat_d& seed1Disp, Mat_d& seed2Disp);
void getNearestSeeds(const Mat_d& corners, const Mat_d& seed, Mat_i& seedId);

void getDisparityMat(const Mat_d& corners1,
		const Mat_d& corners2,
		const Mat_d& seed1,
		const Mat_d& seed2,
		double dispMax,
		Mat_d& dispMat);

void greedyNCCMatch(const Mat_d& nccMat, Matching& matches);
int greedyNCCMatch(const Mat_d& nccMat, Mat_i& jInd);
void greedyGuidedNCCMatch(const Mat_d& nccMat, const Mat_d& dispMat, Matching& matches);
void lapGuidedNCCMatch(const Mat_d& nccMat, const Mat_d& dispMat, Matching& matches, double largeWeight = 1e+5);
#endif /* SL_GUIDEDNCCMATCHER_H_ */
