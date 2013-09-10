/*
 * SL_EMatEst.h
 *
 *  Created on: 2011-7-22
 *      Author: zou
 */

#ifndef SL_EMATEST_H_
#define SL_EMATEST_H_

#include "math/SL_LinAlg.h"
#include "math/SL_Matrix.h"
#include "imgproc/SL_Image.h"
#include "geometry/SL_Distortion.h"

#include <cassert>
#include <vector>

#include "opencv2/opencv.hpp"
#define FF_RANSAC		CV_FM_RANSAC
#define FF_LMEDS		CV_FM_LMEDS
/**
 * find a fundamental matrix that 
 * pts1^T F pts2 = 0
 */
bool findFMatrix(const Mat_d& pts1, const Mat_d& pts2, Mat_d& F, double maxEpiErr, int method);
int getInlierFMatches(const Mat_d& pts1, const Mat_d& pts2, const double F[9], double maxEpiErr, Mat_uc& inlierFlags) ;
bool refineFMatrix(const Mat_d& pts1, const Mat_d& pts2, const Mat_uc& inlierFlags, Mat_d& F);

class CalibTwoCamParam {
public:
	double hessianThreshold;
	double surfRatio;
	double surfMaxDiff;
	double maxEpiErr;
	int nRansac;
	int minInlierNum;
public:
	CalibTwoCamParam() :
			hessianThreshold(100), surfRatio(0.83), surfMaxDiff(0.6), maxEpiErr(6.0), nRansac(1000), minInlierNum(20) {
	}
};
class CalibTwoCam: public CalibTwoCamParam {
public:
	CalibTwoCam();
	virtual ~CalibTwoCam();
public:
	//intrinsic matrices of two cameras
	const double* K1, *K2;
	//distortion parameters 
	const double* kc1, *kc2;

	//inverses of intrinsic matrices and distortion parameters
	double iK1[9], iK2[9];
	double ikc1[7], ikc2[7];

	//matched feature points
	Mat_d orgPts1, orgPts2;
	Mat_d normPts1, normPts2;
	Mat_uc inlierFlag;

	//estimated fundamental matrix
	double F[9];
	//estimated essential matrix
	double E[9];

	void setIntrinParam(const double* _K1, const double* _K2) {
		K1 = _K1;
		K2 = _K2;
		mat33Inv(K1, iK1);
		mat33Inv(K2, iK2);
	}

	void setDistorParam(int IW, int IH, const double* _kc1, const double* _kc2) {
		assert(K1 && K2);
		kc1 = _kc1;
		kc2 = _kc2;
		invDistorParam(IW, IH, iK1, kc1, ikc1);
		invDistorParam(IW, IH, iK2, kc2, ikc2);
	}
	void matchSURFPoints(const ImgG& img1, const ImgG& img2);
	void setMatchedPoints(const Mat_d& pts1, const Mat_d& pts2);
	void getInlierInd(std::vector<int>& ind);
	void getInlierFlag(Mat_uc& inFlag){
		inFlag.cloneFrom(inlierFlag);
	}
	
	int estimateEMatOld(double maxEpiErr = 3.0, int method = FF_RANSAC);
	int estimateEMat(double maxEpiErr = 6.0);
	void outputInlierPoints(Mat_d& pts1, Mat_d& pts2);
	void outputInlierNormPoints(Mat_d& normPts1, Mat_d& normPts2);
	void outputRTs(Mat_d& R1, Mat_d& t1, Mat_d& R2, Mat_d& t2, bool all = false);
};
#endif /* SL_EMATEST_H_ */
