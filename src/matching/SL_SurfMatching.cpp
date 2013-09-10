/*
 * SL_SurfMatching.cpp
 *
 *  Created on: 2011-2-3
 *      Author: Danping Zou
 */

#include "SL_SurfMatching.h"
#include "math/SL_LinAlg.h"
#include "geometry/SL_FundamentalMatrix.h"
#include <cmath>
#include <cfloat>
#if CV_MINOR_VERSION > 3
#include "opencv2/nonfree/features2d.hpp"
#endif

float computeSurfDescDist(const float* desc0, const float* desc1, int dimDesc) {
	float dist = 0;
	for (int i = 0; i < dimDesc; i++) {
		dist += (desc0[i] - desc1[i]) * (desc0[i] - desc1[i]);
	}
	return sqrt(dist);
}
void matchSurf(const Mat_f& kdesc0, const Mat_f& kdesc1, Matching& matches,
		float ratio) {

	assert(kdesc0.rows > 0);
	assert(kdesc0.cols == kdesc1.cols);
	int dimDesc = kdesc0.cols;

	matches.clear();
	matches.reserve(std::max(kdesc0.rows, kdesc1.rows));
	float dist, d1, d2;
	for (int i = 0; i < kdesc0.rows; i++) {
		d1 = d2 = FLT_MAX;
		int jMin = -1;
		for (int j = 0; j < kdesc1.rows; j++) {
			dist = computeSurfDescDist(kdesc0.data + i * dimDesc,
					kdesc1.data + j * dimDesc, dimDesc);
			if (dist < d1) // if this feature matches better than current best
					{
				d2 = d1;
				d1 = dist;
				jMin = j;

			} else if (dist < d2) // this feature matches better than second best
					{
				d2 = dist;
			}
		}
		// If match has a d1:d2 ratio < 0.65 ipoints are a match
		if (d1 / d2 < ratio) {
			matches.add(i, jMin, d1);
		}
	}
}
void matchSurf(int dimDesc, std::vector<float>& desc0,
		std::vector<float>& desc1, Matching& matches, float ratio,
		float maxDist) {
	int len0 = (int) desc0.size();
	int len1 = (int) desc1.size();
	assert(len0 % dimDesc == 0 && len1 %dimDesc == 0);

	int npts0 = len0 / dimDesc;
	int npts1 = len1 / dimDesc;

	matches.clear();
	matches.reserve(npts0 > npts1 ? npts0 : npts1);

	float dist, d1, d2;
	for (int i = 0; i < npts0; i++) {
		d1 = d2 = FLT_MAX;
		int jMin = -1;
		for (int j = 0; j < npts1; j++) {
			dist = computeSurfDescDist(&desc0[0] + i * dimDesc,
					&desc1[0] + j * dimDesc, dimDesc);
			if (dist < d1) {
				d2 = d1;
				d1 = dist;
				jMin = j;
			} else if (dist < d2)
				d2 = dist;
		}
		if (d1 < maxDist && d1 / d2 < ratio) {
			matches.add(i, jMin, d1);
		}
	}
}
int refineMatchedPoints(const Mat_d& pts1, const Mat_d& pts2, Matching& matches,
		Matching& newMatches, double ratio) {
	double ud[2] = { 0, 0 };
	Mat_d d(matches.num, 2);

	int num = matches.num;
	for (int i = 0; i < num; i++) {
		int idx1 = matches[i].idx1;
		int idx2 = matches[i].idx2;

		d.data[2 * i] = pts2.data[2 * idx2] - pts1.data[2 * idx1];
		d.data[2 * i + 1] = pts2.data[2 * idx2 + 1] - pts1.data[2 * idx1 + 1];
		ud[0] += d.data[2 * i];
		ud[1] += d.data[2 * i + 1];
	}

	ud[0] /= num;
	ud[1] /= num;

	double cov[4] = { 0, 0, 0, 0 };
	for (int i = 0; i < num; i++) {
		double dx = d.data[2 * i] - ud[0];
		double dy = d.data[2 * i + 1] - ud[1];

		cov[0] += dx * dx;
		cov[1] += dx * dy;
		cov[3] += dy * dy;
	}

	double s = ratio * ratio;
	cov[0] /= num / s;
	cov[1] /= num / s;
	cov[3] /= num / s;
	cov[2] = cov[1];

	double icov[4];
	mat22Inv(cov, icov);

	newMatches.clear();
	newMatches.reserve(num);

	for (int i = 0; i < num; i++) {
		double dx = d.data[2 * i] - ud[0];
		double dy = d.data[2 * i + 1] - ud[1];
		double dist = dx * dx * icov[0] + 2 * dx * dy * icov[1]
				+ dy * dy * icov[2];
		if (dist < 1.0) {
			newMatches.add(matches[i].idx1, matches[i].idx2, 0);
		}
	}
	return newMatches.num;
}

int detectSURFPoints(const ImgG& img, Mat_d& surfPts,
		std::vector<float>& surfDesc, double hessianThreshold) {

	KpVec surfPtsVec;
	cv::SURF surf(hessianThreshold, 4, 2, false);
	cv::Mat cvImg(img.rows, img.cols, CV_8UC1, img.data);
	surf(cvImg, cv::Mat(), surfPtsVec, surfDesc);
	KpVec2Mat(surfPtsVec, surfPts);
	return surf.descriptorSize();
}
