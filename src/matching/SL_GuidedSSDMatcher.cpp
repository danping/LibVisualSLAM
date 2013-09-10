/*
 * SL_GuidedSSDMatcher.cpp
 *
 *  Created on: 2011-1-25
 *      Author: Danping Zou
 */

#include "SL_GuidedSSDMatcher.h"
#include "SL_StereoMatcherHelper.h"
#include "lap/LAP.h"
#include <cassert>
GuidedSSDMatcher::GuidedSSDMatcher() {
//	m_pPtFlag[0] = 0;
//	m_pPtFlag[1] = 0;
}

GuidedSSDMatcher::~GuidedSSDMatcher() {
}
void GuidedSSDMatcher::_computeDescDistMat(
		const Mat_d& epiMat ,
		const Mat_d& desc1 ,
		const Mat_d& desc2 ,
		Mat_d& descMat) {

	assert(desc1.cols == desc2.cols);
	int nDim = desc1.cols;
	int M = epiMat.rows;
	int N = epiMat.cols;
	
	descMat.resize(epiMat.rows, epiMat.cols);

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			if (epiMat.data[i * N + j] < 0) {
				descMat.data[i * N + j] = -1;
			} else {
				double desc = computeDescDist(nDim, desc1.data + nDim * i, desc2.data + nDim * j);
				descMat.data[i * N + j] = desc < grayLevelThres ? desc : -1;
			}
		}
	}
}
void GuidedSSDMatcher::_computeDispDistMat(
		const Mat_d& epiMat ,
		const Mat_d& pts1 ,
		const Mat_d& pts2 ,
		Mat_d& dispMat) {

	Mat_d seed1Disp, seed2Disp;
	_getSeedDisparities(m_seedPts[0], m_seedPts[1], seed1Disp, seed2Disp);

	int M = pts1.rows;
	int N = pts2.rows;

	dispMat.resize(M, N);

	double thres = dispDistThres;//*dispDistThres;

	for (int i = 0; i < M; i++) {
		double x1 = pts1.data[2 * i];
		double y1 = pts2.data[2 * i + 1];

		int iNearestSeed = searchNearestPoint(m_seedPts[0], x1, y1);

		double predDx1 = seed1Disp.data[2 * iNearestSeed];
		double predDy1 = seed1Disp.data[2 * iNearestSeed + 1];

		for (int j = 0; j < N; j++) {
			if (epiMat.data[i * N + j] < 0) {
				dispMat.data[i * N + j] = -1;
			} else {
				double x2 = pts2.data[2 * j];
				double y2 = pts2.data[2 * j + 1];

				double dx = x2 - x1;
				double dy = y2 - y1;

				int jNearestSeed = searchNearestPoint(m_seedPts[1], x2, y2);

				double predDx2 = seed2Disp.data[2 * jNearestSeed];
				double predDy2 = seed2Disp.data[2 * jNearestSeed + 1];

				double d1 = dist2(predDx1, predDy1, dx, dy);
				double d2 = dist2(predDx2, predDy2, -dx, -dy);

				double d = d1 > d2 ? d1 : d2;//d1*d2;

				dispMat.data[i * N + j] = d > thres ? -1 : d;
			}
		}
	}
}

void GuidedSSDMatcher::_getSeedDisparities(
		const Mat_d& seed1 ,
		const Mat_d& seed2 ,
		Mat_d& seed1Disp ,
		Mat_d& seed2Disp) {
	int numMatches = seed1.rows;

	seed1Disp.resize(numMatches, 2);
	seed2Disp.resize(numMatches, 2);

	for (int i = 0; i < numMatches; i++) {
		double x1 = seed1.data[i * 2];
		double y1 = seed1.data[i * 2 + 1];
		double x2 = seed2.data[i * 2];
		double y2 = seed2.data[i * 2 + 1];

		double dx = x2 - x1;
		double dy = y2 - y1;

		seed1Disp.data[2 * i] = dx;
		seed1Disp.data[2 * i + 1] = dy;

		seed2Disp.data[2 * i] = -dx;
		seed2Disp.data[2 * i + 1] = -dy;
	}
}

void GuidedSSDMatcher::_sumWeightMat2(
		const Mat_d& mat1 ,
		const Mat_d& mat2 ,
		double w1 ,
		double w2 ,
		Mat_d& weightMat) {
	int M = mat1.rows;
	int N = mat1.cols;
	weightMat.resize(M, N);

	int len = M * N;

	const double* pMat1 = mat1.data;
	const double* pMat2 = mat2.data;
	double* pWeightMat = weightMat.data;
	for (int i = 0; i < len; i++) {
		if (pMat1[0] < 0 || pMat2[0] < 0)
			pWeightMat[0] = w_none;
		else
			pWeightMat[0] = w1 * pMat1[0] + w2 * pMat2[0];

		pMat1++;
		pMat2++;
		pWeightMat++;
	}
}
void GuidedSSDMatcher::_sumWeightMat3(
		const Mat_d& mat1 ,
		const Mat_d& mat2 ,
		const Mat_d& mat3 ,
		double w1 ,
		double w2 ,
		double w3 ,
		Mat_d& weightMat) {
	int M = mat1.rows;
	int N = mat1.cols;
	weightMat.resize(M, N);

	int len = M * N;

	const double* pMat1 = mat1.data;
	const double* pMat2 = mat2.data;
	const double* pMat3 = mat3.data;
	double* pWeightMat = weightMat.data;
	for (int i = 0; i < len; i++) {
		if (pMat1[0] < 0 || pMat2[0] < 0 || pMat3[0] < 0)
			pWeightMat[0] = w_none;
		else
			pWeightMat[0] = w1 * pMat1[0] + w2 * pMat2[0] + w3 * pMat3[0];

		pMat1++;
		pMat2++;
		pMat3++;
		pWeightMat++;
	}
}
void GuidedSSDMatcher::_getValidRowCol(Mat_d& weightMat , Mat_c& rowFlag , Mat_c& colFlag) {
	int M = weightMat.rows;
	int N = weightMat.cols;

	rowFlag.resize(M, 1);
	colFlag.resize(N, 1);

	for (int i = 0; i < M; i++) {
		bool allMinus = true;
		for (int j = 0; j < N; j++) {
			if (weightMat.data[i * N + j] >= 0) {
				allMinus = false;
				break;
			}
		}
		rowFlag.data[i] = allMinus ? 0 : 1;
	}

	for (int j = 0; j < N; j++) {
		bool allMinus = true;
		for (int i = 0; i < N; i++) {
			if (weightMat.data[i * N + j] >= 0) {
				allMinus = false;
				break;
			}
		}
		colFlag.data[j] = allMinus ? 0 : 1;
	}
}
int GuidedSSDMatcher::lapMatch(Matching& matches , int type) {
	if (type == TYPE_EPI) {
		if (m_seedPts[0] == 0) {
			getEpiMat(F, *m_pPoints[0], *m_pPoints[1], epiDistThres, m_weightMat);
		} else {
			getEpiMat(F, *m_pPoints[0], *m_pPoints[1], epiDistThres, m_epiMat);
			_computeDispDistMat(m_epiMat, *m_pPoints[0], *m_pPoints[1], m_dispMat);
			_sumWeightMat2(m_epiMat, m_dispMat, w_epi / epiSigma, w_disp / dispSigma, m_weightMat);
		}
	} else if (type == TYPE_EPI_DESC) {
		if (m_seedPts[0] == 0) {
			getEpiMat(F, *m_pPoints[0], *m_pPoints[1], epiDistThres, m_epiMat);
			_computeDescDistMat(m_epiMat, *m_pImgBlk[0], *m_pImgBlk[1], m_descMat);
			_sumWeightMat2(m_epiMat, m_descMat, w_epi / epiSigma, w_desc / graySigma, m_weightMat);
		} else {
			getEpiMat(F, *m_pPoints[0], *m_pPoints[1], epiDistThres, m_epiMat);
			_computeDescDistMat(m_epiMat, *m_pImgBlk[0], *m_pImgBlk[1], m_descMat);
			_computeDispDistMat(m_descMat, *m_pPoints[0], *m_pPoints[1], m_dispMat);
			_sumWeightMat3(
					m_epiMat,
					m_descMat,
					m_dispMat,
					w_epi / epiSigma,
					w_desc / graySigma,
					w_disp / dispSigma,
					m_weightMat);
		}
	}
	int M = m_weightMat.rows;
	int N = m_weightMat.cols;
	//solve the linear assignment problem
	if (M == N) {

		Mat_c rowFlag, colFlag;
		m_ind[0].resize(M, 1);
		m_ind[1].resize(N, 1);
		m_ind[0].fill(-1);
		m_ind[1].fill(-1);

		_getValidRowCol(m_weightMat, rowFlag, colFlag);
		if (lap(m_weightMat.data, M, N, m_ind[0].data, m_ind[1].data, rowFlag.data, colFlag.data) < 0)
			repErr("cannot solve this assignment problem\n");

		matches.reserve(M);
		for (int i = 0; i < M; i++) {
			int j = m_ind[0].data[i];
			if (j >= 0 && m_weightMat.data[i * N + j] != w_none && m_weightMat.data[i * N + j] != w_dummy)
				matches.add(i, j, 0);
		}
		return matches.num;
	}

	if (M < N) {
		Mat_d weightMat(N, N);
		for (int i = 0; i < M; i++)
			memcpy(weightMat.data + i * N, m_weightMat.data + i * N, sizeof(double) * N);
		for (int i = M; i < N; i++)
			for (int j = 0; j < N; j++)
				weightMat.data[i * N + j] = w_dummy;

		Mat_c rowFlag, colFlag;
		rowFlag.resize(N, 1);
		rowFlag.fill(1);
		colFlag.resize(N, 1);
		colFlag.fill(1);
		m_ind[0].resize(N, 1);
		m_ind[1].resize(N, 1);
		m_ind[0].fill(-1);
		m_ind[1].fill(-1);

		if (lap(weightMat.data, N, N, m_ind[0].data, m_ind[1].data, rowFlag.data, colFlag.data) < 0)
			repErr("cannot solve this assignment problem\n");

		matches.reserve(N);
		for (int i = 0; i < N; i++) {
			int j = m_ind[0].data[i];
			if (j >= 0 && weightMat.data[i * N + j] != w_none && weightMat.data[i * N + j] != w_dummy)
				matches.add(i, j, 0);
		}

	} else {
		Mat_d weightMat(M, M);
		for (int i = 0; i < M; i++) {
			memcpy(weightMat.data + i * M, m_weightMat.data + i * N, sizeof(double) * N);
			for (int j = N; j < M; j++)
				weightMat.data[i * M + j] = w_dummy;
		}

		//test
//		LogFile file("/home/Danping Zou/weightMat.txt");
//		print(file.fp,weightMat);
//		file.close();
		Mat_c rowFlag, colFlag;
		rowFlag.resize(M, 1);
		rowFlag.fill(1);
		colFlag.resize(M, 1);
		colFlag.fill(1);
		m_ind[0].resize(M, 1);
		m_ind[1].resize(M, 1);
		m_ind[0].fill(-1);
		m_ind[1].fill(-1);

		if (lap(weightMat.data, M, M, m_ind[0].data, m_ind[1].data, rowFlag.data, colFlag.data) < 0)
			repErr("cannot solve this assignment problem\n");

		matches.reserve(M);
		for (int i = 0; i < M; i++) {
			int j = m_ind[0].data[i];
			//test
			//printf("(i,j):%d - %d, weight:%lf\n",i,j,weightMat.data[i*M + j]);
			if (j >= 0 && weightMat.data[i * M + j] != w_none && weightMat.data[i * M + j] != w_dummy)
				matches.add(i, j, 0);
		}
	}
	return matches.num;
}

int GuidedSSDMatcher::lapMatch(const double *tF , const Mat_d & pts1 , const Mat_d & pts2 , Matching & matches) {
	setFMat(tF);
	m_pPoints[0] = &pts1;
	m_pPoints[1] = &pts2;

	return lapMatch(matches, TYPE_EPI);
}

int GuidedSSDMatcher::lapMatch(
		const double *tF ,
		const Mat_d & pts1 ,
		const Mat_d & pts2 ,
		const Mat_d & desc1 ,
		const Mat_d & desc2 ,
		Matching & matches) {

	setFMat(tF);
	m_pPoints[0] = &pts1;
	m_pPoints[1] = &pts2;
	m_pImgBlk[0] = &desc1;
	m_pImgBlk[1] = &desc2;

	return lapMatch(matches, TYPE_EPI_DESC);
}

void GuidedSSDMatcher::removeOutlier(const Mat_d& pts1 , const Mat_d& pts2 , Matching& matches) {
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

	removeOutlier(pts1, pts2, matches, avgDx, avgDy, dispDistThres);
}
void GuidedSSDMatcher::removeOutlier(
		const Mat_d& pts1 ,
		const Mat_d& pts2 ,
		Matching& matches ,
		double avgDx ,
		double avgDy ,
		double distThres) {

	Matching tmpMatches;
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
		double dist = dist2(avgDx, avgDy, dx, dy);
		if (dist < distThres) {
			tmpMatches.add(idx1, idx2, matches[i].dist);
		}
	}
	//test
	logInfo("%d matches removed", matches.num - tmpMatches.num);
	matches.clear();
	matches.clone(tmpMatches);
}
