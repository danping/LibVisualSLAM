/*
 * SL_StereoMacherHelper.cpp
 *
 *  Created on: 2011-1-13
 *      Author: Danping Zou
 */

#include "SL_StereoMatcherHelper.h"
#include <cmath>
#include <cfloat>
#include <cassert>
#include <opencv2/opencv.hpp>

#include "tools/cvHelper.h"
#include "geometry/SL_Geometry.h"
int getMatchedPts(const Matching& matches, const Mat_d& pts1, const Mat_d& pts2,
		Mat_d& pts1Matched, Mat_d& pts2Matched) {

	int num = matches.num;
	if (num == 0)
		repErr("no matched points!");

	pts1Matched.resize(num, 2);
	pts2Matched.resize(num, 2);

	int k = 0;
	for (int i = 0; i < num; i++) {
		int idx1 = matches[i].idx1;
		int idx2 = matches[i].idx2;
		if (idx1 < 0 || idx2 < 0)
			continue;

		double x1 = pts1.data[2 * idx1];
		double y1 = pts1.data[2 * idx1 + 1];

		double x2 = pts2.data[2 * idx2];
		double y2 = pts2.data[2 * idx2 + 1];

		pts1Matched.data[2 * k] = x1;
		pts1Matched.data[2 * k + 1] = y1;

		pts2Matched.data[2 * k] = x2;
		pts2Matched.data[2 * k + 1] = y2;

		k++;
	}
	return k;
}
void getDistMat(const Mat_d& pts, Mat_d& distMat) {
	int nPts = pts.rows;
	distMat.resize(nPts, nPts);
	for (int i = 0; i < nPts; i++) {
		distMat.data[i * nPts + i] = 0;
		for (int j = 1; j < nPts; j++) {
			double dx = pts.data[2 * i] - pts.data[2 * j];
			double dy = pts.data[2 * i + 1] - pts.data[2 * j + 1];
			distMat.data[i * nPts + j] = sqrt(dx * dx + dy * dy);
			distMat.data[j * nPts + i] = distMat.data[i * nPts + j];
		}
	}
}
void getEpiMat(const double* F, const Mat_d& pts1, const Mat_d& pts2,
		double epiDistThres, Mat_d& epiMat) {
	int M = pts1.rows;
	int N = pts2.rows;
	epiMat.resize(M, N);

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			double epiErr = epipolarError(F, pts1.data + 2 * i,
					pts2.data + 2 * j);
			epiMat.data[i * N + j] = epiErr > epiDistThres ? -1 : epiErr;
		}
	}
}
void getEpiASDMat(const double* F, const Mat_d& pts1, const Mat_d& pts2,
		const Mat_d& blks1, const Mat_d& blks2, double epiMax, double asdMax,
		Mat_d& epiMat, Mat_d& asdMat, double wNone) {

	int M = pts1.rows;
	int N = pts2.rows;

	assert(blks1.rows == M && blks2.rows == N);
	assert(blks1.cols == blks2.cols);

	epiMat.resize(M, N);
	asdMat.resize(M, N);

	int L = blks1.cols;

	for (int i = 0; i < M; i++) {
		double* pEpiData = epiMat.data + i * N;
		double* pAsdData = asdMat.data + i * N;

		for (int j = 0; j < N; j++) {
			double epiErr = epipolarError(F, pts2.data + 2 * i,
					pts1.data + 2 * j);
			if (epiErr <= epiMax) {
				//compute the ASD
				double s = 0;
				for (int k = 0; k < L; k++) {
					double d = blks1.data[i * L + k] - blks2.data[j * L + k];
					s += d * d;
				}
				s /= L;

				if (s <= asdMax) {
					*pEpiData = epiErr;
					*pAsdData = s;
				} else {
					*pEpiData = wNone;
					*pAsdData = wNone;
				}
			} else {
				*pEpiData = wNone;
				*pAsdData = wNone;
			}
			pEpiData++;
			pAsdData++;
		}
	}
}

double computeASD(int nDim, const double* blk1, const double* blk2) {
	double s = 0;
	for (int i = 0; i < nDim; i++) {
		double d = blk1[i] - blk2[i];
		s += d * d;
	}
	return s / nDim;
}
double computeDescDist(int nDim, const double* ds1, const double* ds2) {
	double s = 0;
	for (int i = 0; i < nDim; i++) {
		double d = ds1[i] - ds2[i];
		s += fabs(d);
	}
	return s / nDim;
}
void getImageBlock(const ImgG& img, double x, double y, Mat_d& blk, int hw) {
	blk.resize((2 * hw + 1), (2 * hw + 1));
	getImageBlock(img, x, y, blk.data, hw);
}
void getImageBlock(const ImgG& img, double x, double y, double* blkData,
		int hw) {
	int W = 2 * hw + 1;
	int H = 2 * hw + 1;
	int len = W * H;

	cv::Mat cvImg(img.rows, img.cols, CV_8UC1, img.data);
	cv::Size patchSize(W, H);
	cv::Mat patch;
	cv::getRectSubPix(cvImg, patchSize, cv::Point2d(x, y), patch);
	for (int i = 0; i < len; i++) {
		blkData[i] = patch.data[i];
	}
}

void getImageBlocks(const ImgG& img, const Mat_d& pts, Mat_d& blks, int hW,
		double scale) {
	if (pts.empty())
		repErr("getImageBlocks - no point is provided");

	cv::Size patchSize(2 * hW + 1, 2 * hW + 1);
	int npts = pts.rows;
	int len = (2 * hW + 1) * (2 * hW + 1);
	blks.resize(npts, len);

	cv::Mat imgRef(img.rows, img.cols, CV_8UC1, img.data);

	if (scale == 1.0) {
		for (int i = 0; i < npts; i++) {
			double x = pts.data[2 * i];
			double y = pts.data[2 * i + 1];

			cv::Mat patch;
			cv::getRectSubPix(imgRef, patchSize, cv::Point2f(x, y), patch);
			//cv::Scalar avg = cv::mean(patch);
			CvMat patchRef = patch;

			for (int j = 0; j < len; j++)
				//desc.data[i * len + j] = (patchRef.data.ptr[j] - avg(0));
				blks.data[i * len + j] = patchRef.data.ptr[j];
		}
	} else {
		cv::Mat tmpImg;
		cv::resize(imgRef, tmpImg, cv::Size(), scale, scale);
		for (int i = 0; i < npts; i++) {
			double x = pts.data[2 * i] * scale;
			double y = pts.data[2 * i + 1] * scale;

			cv::Mat patch;
			cv::getRectSubPix(tmpImg, patchSize, cv::Point2f(x, y), patch);
			//cv::Scalar avg = cv::mean(patch);
			CvMat patchRef = patch;

			for (int j = 0; j < len; j++)
				//desc.data[i * len + j] = (patchRef.data.ptr[j] - avg(0));
				blks.data[i * len + j] = patchRef.data.ptr[j];
		}
	}

}
int searchNearestPoint(const Mat_d& pts, double x0, double y0) {
	assert(pts.rows > 0);
	int num = pts.rows;

	double dMin = DBL_MAX;
	int iMin = -1;
	for (int i = 0; i < num; i++) {
		double x = pts[2 * i];
		double y = pts[2 * i + 1];

		double d = dist2(x0, y0, x, y);
		if (d < dMin) {
			iMin = i;
			dMin = d;
		}
	}
	return iMin;
}
int searchNearestPoint(const Mat_d& pts, double x0, double y0, double maxDist,
		double* dist) {
	int npts = pts.rows;
	double dMin = maxDist;
	int iMin = -1;
	for (int i = 0; i < npts; i++) {
		double x = pts.data[2 * i];
		double y = pts.data[2 * i + 1];

		double dx = x - x0;
		double dy = y - y0;
		double d = sqrt(dx * dx + dy * dy);
		if (d < dMin) {
			dMin = d;
			iMin = i;
		}
	}
	if (dist)
		*dist = dMin;
	return iMin;
}

int searchNearestPoint(const Mat_d& pts, double x0, double y0, double var[4],
		double maxDist, double* dist) {
	double ivar[4];
	mat22Inv(var, ivar);

	int npts = pts.rows;
	double dMin = maxDist;
	int iMin = -1;
	for (int i = 0; i < npts; i++) {
		double x = pts.data[2 * i];
		double y = pts.data[2 * i + 1];

		double d = mahaDist2(x, y, x0, y0, ivar);
		if (d < dMin) {
			dMin = d;
			iMin = i;
		}
	}
	if (dist)
		*dist = dMin;
	return iMin;
}
int findBestMatch(const Mat_d& desc1, const Mat_d& desc2, int idx1, int idx2,
		const Mat_d& desc, const Mat_d& pts, double x0, double y0,
		double maxDistThes, double maxDescThres) {
	int npts = pts.rows;
	int dimDesc = desc1.cols;

	const double* pDesc1 = desc1.data + idx1 * dimDesc;
	const double* pDesc2 = desc2.data + idx2 * dimDesc;
	double minDesc = maxDescThres;
	int imin = -1;
	for (int i = 0; i < npts; i++) {
		double x = pts.data[2 * i];
		double y = pts.data[2 * i + 1];

		double dx = x - x0;
		double dy = y - y0;
		double d = sqrt(dx * dx + dy * dy);
		if (d > maxDistThes)
			continue;

		double distDesc = computeDescDist(dimDesc, pDesc1,
				desc.data + i * dimDesc);
		distDesc += computeDescDist(dimDesc, pDesc2, desc.data + i * dimDesc);
		if (distDesc * 0.5 < minDesc) {
			minDesc = distDesc;
			imin = i;
		}
	}
	return imin;
}

void getFlaged2DPoints(const Mat_d& pts, const Mat_uc& flag, Mat_d& ptsFlaged) {
	if (pts.rows != flag.rows)
		repErr("getFlagedPoints - error");
	ptsFlaged.resize(pts.rows, pts.cols);
	int k = 0, nPts = pts.rows;
	for (int i = 0; i < nPts; i++) {
		if (flag.data[i]) {
			ptsFlaged.data[2 * k] = pts.data[2 * i];
			ptsFlaged.data[2 * k + 1] = pts.data[2 * i + 1];
			k++;
		}
	}
	ptsFlaged.rows = k;
}
