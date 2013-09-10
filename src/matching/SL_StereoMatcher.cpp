/*
 * FeatureMatcher.cpp
 *
 *  Created on: 2011-1-7
 *      Author: Danping Zou
 */

#include "SL_StereoMatcher.h"
#include "SL_StereoMatcherHelper.h"
#include "math/SL_LinAlg.h"
#include "geometry/SL_Geometry.h"
#include <cfloat>
#include <cmath>
StereoMatcher::StereoMatcher(const double* FMat) :
	F(FMat) {

}

StereoMatcher::~StereoMatcher() {
}

double StereoMatcher::_getDisplacementError(
		const double* pts1 ,
		const double* pts2 ,
		const double avgDx ,
		const double avgDy) {
	double dx = pts1[0] - pts2[0];
	double dy = pts1[1] - pts2[1];
	return dist2(dx, dy, avgDx, avgDy);
}
void StereoMatcher::_computeEpipolarDistMat(const Mat_d& pts1 , const Mat_d& pts2 , Mat_d& distMat) {
	getEpiMat(F, pts1, pts2, epiDistThres, distMat);
}
void StereoMatcher::_computeEpipolarDisplacementMat(
		const Mat_d& pts1 ,
		const Mat_d& pts2 ,
		double avgDx ,
		double avgDy ,
		double distThes ,
		Mat_d& distMat) {

	int M = pts1.rows;
	int N = pts2.rows;
	distMat.resize(M, N);

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			double epiErr = epipolarError(F, pts1.data + 2 * i, pts2.data + 2 * j);
			if (epiErr > epiDistThres) {
				distMat.data[i * N + j] = -1;
				continue;
			} else {
				double d = _getDisplacementError(pts1.data + 2 * i, pts2.data + 2 * j, avgDx, avgDy);
				if (d > distThes) {
					distMat.data[i * N + j] = -1;
					continue;
				}
				distMat.data[i * N + j] = epiErr;
			}
		}
	}
}

void StereoMatcher::_computeDescDistMat(const Mat_d& epiMat , const Mat_d& ds1 , const Mat_d& ds2 , Mat_d& distMat) {

	int M = ds1.rows;
	int N = ds2.rows;

	distMat.resize(M, N);
	assert(ds1.cols == ds2.cols);

	int nDim = ds1.cols;

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			if (epiMat.data[i * N + j] >= 0) {
				//compute the descriptor distance
				double desc = computeDescDist(nDim, ds1.data + nDim * i, ds2.data + nDim * j);
				distMat.data[i * N + j] = desc < descDistThres ? desc : -1;
			} else
				distMat.data[i * N + j] = -1;
		}
	}
}

int StereoMatcher::_findBestMatchRow(Mat_d& distMat , int iRow) {
	int N = distMat.cols;
	double* prow = distMat.data + iRow * N;

	int iMin = -1;
	double minVal = DBL_MAX;
	for (int i = 0; i < N; i++) {
		if (prow[i] < 0)
			continue;
		if (prow[i] < minVal) {
			iMin = i;
			minVal = prow[i];
		}
	}
	return iMin;
}
int StereoMatcher::_findBestMatchCol(Mat_d& distMat , int iCol) {
	int N = distMat.cols;
	int M = distMat.rows;
	double* pcol = distMat.data + iCol;

	int iMin = -1;
	double minVal = 0;
	for (int i = 0; i < M; i++) {
		if (pcol[i * N] < 0)
			continue;
		if (pcol[i * N] > minVal) {
			iMin = i;
			minVal = pcol[i * N];
		}
	}
	return iMin;
}

void StereoMatcher::match(const Mat_d& pts1 , const Mat_d& pts2 , Matching& matches) {
	_computeEpipolarDistMat(pts1, pts2, epiDistMat);

	int M = pts1.rows;
	int N = pts2.rows;

	Mat_i RBestMatch(M, 1);
	Mat_i CBestMatch(1, N);

	for (int i = 0; i < M; i++) {
		RBestMatch.data[i] = _findBestMatchRow(epiDistMat, i);
	}
	for (int i = 0; i < N; i++) {
		CBestMatch.data[i] = _findBestMatchCol(epiDistMat, i);
	}

	if (M <= N) {
		matches.reserve(M);
		for (int i = 0; i < M; i++) {
			int j = RBestMatch.data[i];
			if (j >= 0 && CBestMatch.data[j] == i) {
				matches.add(i, j, epiDistMat.data[i * N + j]);
			}
		}
	} else {
		matches.reserve(N);
		for (int j = 0; j < N; j++) {
			int i = CBestMatch.data[j];
			if (i >= 0 && RBestMatch.data[i] == j) {
				matches.add(i, j, epiDistMat.data[i * N + j]);
			}
		}
	}
}
void StereoMatcher::match(
		const Mat_d& pts1 ,
		const Mat_d& pts2 ,
		const Mat_d& desc1 ,
		const Mat_d& desc2 ,
		Matching& matches) {

	if (pts1.rows != desc1.rows || pts2.rows != desc2.rows)
		repErr("StereoMatcher::match - the number of points and descriptors should be the same!");

	_computeEpipolarDistMat(pts1, pts2, epiDistMat);
	_computeDescDistMat(epiDistMat, desc1, desc2, descDistMat);

	int M = pts1.rows;
	int N = pts2.rows;

	Mat_i RBestMatch(M, 1);
	Mat_i CBestMatch(1, N);

	for (int i = 0; i < M; i++) {
		RBestMatch.data[i] = _findBestMatchRow(descDistMat, i);
	}
	for (int i = 0; i < N; i++) {
		CBestMatch.data[i] = _findBestMatchCol(descDistMat, i);
	}

	if (M <= N) {
		matches.reserve(M);
		for (int i = 0; i < M; i++) {
			int j = RBestMatch.data[i];
			if (j >= 0 && CBestMatch.data[j] == i) {
				matches.add(i, j, descDistMat.data[i * N + j]);
			}
		}
	} else {
		matches.reserve(N);
		for (int j = 0; j < N; j++) {
			int i = CBestMatch.data[j];
			if (i >= 0 && RBestMatch.data[i] == j) {
				matches.add(i, j, epiDistMat.data[i * N + j]);
			}
		}
	}
}

void StereoMatcher::match(
		const Mat_d& pts1 ,
		const Mat_d& pts2 ,
		const Mat_d& desc1 ,
		const Mat_d& desc2 ,
		Matching& matches ,
		double avgDx ,
		double avgDy ,
		double dispThres) {

	if (pts1.rows != desc1.rows || pts2.rows != desc2.rows)
		repErr("StereoMatcher::match - the number of points and descriptors should be the same!");

	_computeEpipolarDisplacementMat(pts1, pts2, avgDx, avgDy, dispThres, epiDistMat);
	_computeDescDistMat(epiDistMat, desc1, desc2, descDistMat);

	int M = pts1.rows;
	int N = pts2.rows;

	Mat_i RBestMatch(M, 1);
	Mat_i CBestMatch(1, N);

	for (int i = 0; i < M; i++) {
		RBestMatch.data[i] = _findBestMatchRow(descDistMat, i);
	}
	for (int i = 0; i < N; i++) {
		CBestMatch.data[i] = _findBestMatchCol(descDistMat, i);
	}

	if (M <= N) {
		matches.reserve(M);
		for (int i = 0; i < M; i++) {
			int j = RBestMatch.data[i];
			if (j >= 0 && CBestMatch.data[j] == i) {
				matches.add(i, j, descDistMat.data[i * N + j]);
			}
		}
	} else {
		matches.reserve(N);
		for (int j = 0; j < N; j++) {
			int i = CBestMatch.data[j];
			if (i >= 0 && RBestMatch.data[i] == j) {
				matches.add(i, j, epiDistMat.data[i * N + j]);
			}
		}
	}
}
void StereoMatcher::removeOutlier(const Mat_d& pts1 , const Mat_d& pts2 , Matching& matches) {
	Matching tmpMatches;
	double avgDx = 0;
	double avgDy = 0;
	for (int i = 0; i < matches.num; i++) {
		int idx1 = matches[i].idx1;
		int idx2 = matches[i].idx2;

		double x1 = pts1.data[2 * idx1];
		double y1 = pts1.data[2 * idx1 + 1];

		double x2 = pts2.data[2 * idx2];
		double y2 = pts2.data[2 * idx2 + 1];

		avgDx = x2 - x1;
		avgDy = y2 - y1;
	}

	avgDx /= matches.num;
	avgDy /= matches.num;

	tmpMatches.reserve(matches.num);
	for (int i = 0; i < matches.num; i++) {
		int idx1 = matches[i].idx1;
		int idx2 = matches[i].idx2;

		double x1 = pts1.data[2 * idx1];
		double y1 = pts1.data[2 * idx1 + 1];

		double x2 = pts2.data[2 * idx2];
		double y2 = pts2.data[2 * idx2 + 1];

		double dx = x2 - x1;
		double dy = y2 - y1;

		double dist = dist2(dx, dy, avgDx, avgDy);
		if (dist < dispDistThres) {
			tmpMatches.add(idx1, idx2, matches[i].dist);
		}
	}
	//test
	logInfo("%d matches removed", matches.num - tmpMatches.num);
	matches.clear();
	matches.clone(tmpMatches);
}

