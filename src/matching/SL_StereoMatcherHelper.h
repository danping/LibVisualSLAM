/*
 * SL_StereoMatcher.h
 *
 *  Created on: 2011-1-13
 *      Author: Danping Zou
 */

#ifndef SL_STEREOMACHERHELPER_H_
#define SL_STEREOMACHERHELPER_H_
#include "math/SL_Matrix.h"
#include "math/SL_LinAlg.h"
#include "geometry/SL_FundamentalMatrix.h"

#include "matching/SL_Matching.h"
#include "imgproc/SL_Image.h"

int getMatchedPts(const Matching& matches,
				const Mat_d& pts1,
				const Mat_d& pts2,
				Mat_d& pts1Matched,
				Mat_d& pts2Matched);

/*compute distance matrix for a set of 2D points*/
void getDistMat(const Mat_d& pts, Mat_d& distMat);

/*compute epipolar error matrix for image points in two cameras*/
void getEpiMat(const double* F, const Mat_d& pts1, const Mat_d& pts2, double epiDistThres, Mat_d& epiMat);

void getEpiASDMat(const double* F,
		const Mat_d& pts1,
		const Mat_d& pts2,
		const Mat_d& blks1,
		const Mat_d& blks2,
		double epiMax,
		double asdMax,
		Mat_d& epiMat,
		Mat_d& asdMat,
		double wNone = -1);

void getImageBlock(const ImgG& img, double x, double y, Mat_d& blk, int hw = 4);
void getImageBlock(const ImgG& img, double x, double y, double* blkData, int hw);
void getImageBlocks(const ImgG& img, const Mat_d& pts, Mat_d& blks, int hW = 4, double scale = 0.3);

double computeASD(int nDim, const double* blk1, const double* blk2);
double computeDescDist(int nDim, const double* ds1, const double* ds2);

/*search the nearest point of (x0,y0)*/
int searchNearestPoint(const Mat_d& pts, double x0, double y0);
int searchNearestPoint(const Mat_d& pts, double x0, double y0, double maxDist, double* dist);
//using Mahalanobis distance
int searchNearestPoint(const Mat_d& pts, double x0, double y0, double var[4], double maxDist, double* dist);

int findBestMatch(const Mat_d& desc1,
		const Mat_d& desc2,
		int idx1,
		int idx2,
		const Mat_d& desc,
		const Mat_d& pts,
		double x0,
		double y0,
		double maxDistThes,
		double maxDescThres);

void getFlaged2DPoints(const Mat_d& pts, const Mat_uc& flag, Mat_d& ptsFlaged);
#endif /* SL_STEREOMACHERHELPER_H_ */
