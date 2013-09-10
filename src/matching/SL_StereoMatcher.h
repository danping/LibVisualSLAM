/*
 * StereoMatcher.h
 *
 *  Created on: 2011-1-7
 *      Author: Danping Zou
 */

#ifndef SL_STEREOMATCHER_H_
#define SL_STEREOMATCHER_H_
#include "math/SL_Matrix.h"
#include "SL_Matching.h"

class StereoMatcherParam {
public:
	double epiDistThres;
	double descDistThres;
	double dispDistThres;

	StereoMatcherParam() {
		epiDistThres = 2;
		descDistThres = 0.4;
		dispDistThres = 50;
	}
};

class StereoMatcher : public StereoMatcherParam {
protected:
	//fundamental matrix
	const double* F;
	Mat_d epiDistMat;
	Mat_d descDistMat;
protected:
	double _getDisplacementError(const double* pts1 , const double* pts2 , const double avgDx , const double avgDy);
	void _computeEpipolarDistMat(const Mat_d& pts1 , const Mat_d& pts2 , Mat_d& distMat);
	void _computeEpipolarDisplacementMat(
			const Mat_d& pts1 ,
			const Mat_d& pts2 ,
			double avgDx ,
			double avgDy ,
			double distThes ,
			Mat_d& distMat);
	void _computeDescDistMat(const Mat_d& epiMat , const Mat_d& ds1 , const Mat_d& ds2 , Mat_d& distMat);
	int _findBestMatchRow(Mat_d& distMat , int iRow);
	int _findBestMatchCol(Mat_d& distMat , int iCol);
public:
	StereoMatcher() :
		F(0) {
	}
	StereoMatcher(const double* FMat);
	virtual ~StereoMatcher();

public:
	//use only epipolar constraint
	void match(const Mat_d& pts1 , const Mat_d& pts2 , Matching& matches);
	//use epipolar constraint + descriptors
	void match(const Mat_d& pts1 , const Mat_d& pts2 , const Mat_d& desc1 , const Mat_d& desc2 , Matching& matches);

	//use epipolar constraint + descriptors + average displacement
	void match(
			const Mat_d& pts1 ,
			const Mat_d& pts2 ,
			const Mat_d& desc1 ,
			const Mat_d& desc2 ,
			Matching& matches ,
			double avgDx ,
			double avgDy ,
			double dispThres);

	//remove the ourliers that have far 
	void removeOutlier(const Mat_d& pts1 , const Mat_d& pts2 , Matching& matches);
	void removeOutlier(
			const Mat_d& pts1 ,
			const Mat_d& pts2 ,
			Matching& matches ,
			double avgDx ,
			double avgDy ,
			double distThres);
};

#endif /* SL_STEREOMATCHER_H_ */
