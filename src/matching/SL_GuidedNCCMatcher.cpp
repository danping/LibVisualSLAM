/*
 * SL_GuidedNCCMatcher.cpp
 *
 *  Created on: 2011-3-30
 *      Author: Danping Zou
 */

#include "SL_GuidedNCCMatcher.h"
#include "SL_StereoMatcherHelper.h"
#include "geometry/SL_Geometry.h"
void getSeedDisparities(const Mat_d& seed1, const Mat_d& seed2, Mat_d& seed1Disp, Mat_d& seed2Disp) {
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

void getNearestSeeds(const Mat_d& corners, const Mat_d& seed, Mat_i& seedId) {
	int N = corners.rows;
	seedId.resize(N, 1);
	for (int i = 0; i < N; i++) {
		double x = corners.data[2 * i];
		double y = corners.data[2 * i + 1];
		int iNearest = searchNearestPoint(seed, x, y);
		seedId.data[i] = iNearest;
	}
}

void getDisparityMat(const Mat_d& corners1, const Mat_d& corners2, const Mat_d& seed1, const Mat_d& seed2, double dispMax, Mat_d& dispMat) {

	Mat_d seed1Disp, seed2Disp;
	getSeedDisparities(seed1, seed2, seed1Disp, seed2Disp);

	Mat_i seedId1, seedId2;
	getNearestSeeds(corners1, seed1, seedId1);
	getNearestSeeds(corners2, seed2, seedId2);

	int M = corners1.rows;
	int N = corners2.rows;

	dispMat.resize(M, N);

#ifdef USE_OPENMP
	#pragma omp parallel for num_threads(SLAM_CPU_CORE_NUM)
#endif
	
	for (int i = 0; i < M; i++) {
		double x1 = corners1.data[2 * i];
		double y1 = corners1.data[2 * i + 1];

		int iNearestSeed = seedId1[i];

		double predDx1 = seed1Disp.data[2 * iNearestSeed];
		double predDy1 = seed1Disp.data[2 * iNearestSeed + 1];

		for (int j = 0; j < N; j++) {
			double x2 = corners2.data[2 * j];
			double y2 = corners2.data[2 * j + 1];

			double dx = x2 - x1;
			double dy = y2 - y1;

			int jNearestSeed = seedId2[j];

			double predDx2 = seed2Disp.data[2 * jNearestSeed];
			double predDy2 = seed2Disp.data[2 * jNearestSeed + 1];

			double d1 = dist2(predDx1, predDy1, dx, dy);
			double d2 = dist2(predDx2, predDy2, -dx, -dy);

			double d = d1 < d2 ? d1 : d2; //d1*d2;

			dispMat.data[i * N + j] = d > dispMax ? -1 : d;
		}
	}
}
void greedyNCCMatch(const Mat_d& nccMat, Matching& matches) {
	int M = nccMat.rows;
	int N = nccMat.cols;

	matches.clear();
	if (M < N) {
		matches.reserve(M);
		for (int i = 0; i < M; i++) {
			double maxNcc = -1;
			int maxJ = -1;
			double* pNcc = nccMat.data + i * N;
			for (int j = 0; j < N; j++) {
				if (pNcc[j] >= 0 && pNcc[j] > maxNcc) {
					maxNcc = pNcc[j];
					maxJ = j;
				}
			}
			if (maxJ >= 0) {
				matches.add(i, maxJ, maxNcc);
			}
		}

	} else {
		matches.reserve(N);
		for (int j = 0; j < N; j++) {
			double maxNcc = -1;
			int maxI = -1;
			double* pNcc = nccMat.data + j;
			for (int i = 0; i < M; i++) {
				if (pNcc[i * N] >= 0 && pNcc[i * N] > maxNcc) {
					maxNcc = pNcc[i * N];
					maxI = i;
				}
			}
			if (maxI >= 0) {
				matches.add(maxI, j, maxNcc);
			}
		}
	}
}

int greedyNCCMatch(const Mat_d& nccMat, Mat_i& jInd) {
	int M = nccMat.rows;
	int N = nccMat.cols;

	jInd.resize(M, 1);
	jInd.fill(-1);
	int nMatched = 0;
	if (M < N) {
		for (int i = 0; i < M; i++) {
			double maxNcc = -1;
			int maxJ = -1;
			double* pNcc = nccMat.data + i * N;
			for (int j = 0; j < N; j++) {
				if (pNcc[j] >= 0 && pNcc[j] > maxNcc) {
					maxNcc = pNcc[j];
					maxJ = j;
				}
			}
			if (maxJ >= 0) {
				jInd.data[i] = maxJ;
				nMatched++;
			}
		}
	} else {
		for (int j = 0; j < N; j++) {
			double maxNcc = -1;
			int maxI = -1;
			double* pNcc = nccMat.data + j;
			for (int i = 0; i < M; i++) {
				if (pNcc[i * N] >= 0 && pNcc[i * N] > maxNcc) {
					maxNcc = pNcc[i * N];
					maxI = i;
				}
			}
			if (maxI >= 0) {
				jInd.data[maxI] = j;
				nMatched++;
			}
		}
	}
	return nMatched;
}
void greedyGuidedNCCMatch(const Mat_d& nccMat, const Mat_d& dispMat, Matching& matches) {
	int M = nccMat.rows;
	int N = nccMat.cols;

	matches.clear();

	assert(dispMat.rows == M && dispMat.cols == N);
	if (M < N) {
		int W = M + 1;
		Mat_i JMat(N, W);
		JMat.fill(0);

		matches.reserve(M);
		for (int i = 0; i < M; i++) {
			double maxNcc = -1;
			int maxJ = -1;
			double* pNcc = nccMat.data + i * N;
			double* pDisp = dispMat.data + i * N;
			for (int j = 0; j < N; j++) {
				if (pNcc[j] >= 0 && pDisp[j] >= 0 && pNcc[j] > maxNcc) {
					maxNcc = pNcc[j];
					maxJ = j;
				}
			}
			if (maxJ >= 0) {
				//record repeated assignment
				int * pJMat = JMat.data + maxJ * W;
				int k = pJMat[0];
				pJMat[k + 1] = i;
				pJMat[0] = k + 1;
			}
		}

		//check the repeated assignment
		for (int j = 0; j < N; j++) {
			int KNum = JMat.data[j * W];
			if (KNum > 1) {
				int maxI = -1;
				double maxNcc = -1;
				for (int k = 0; k < KNum; k++) {
					int i = JMat.data[j * W + k + 1];
					double ncc = nccMat.data[i * N + j];
					if (ncc > maxNcc) {
						maxNcc = ncc;
						maxI = i;
					}
				}
				//find the best assignment
				matches.add(maxI, j, nccMat.data[maxI * N + j]);
			} else if (KNum == 1) {
				int i = JMat.data[j * W + 1];
				matches.add(i, j, nccMat.data[i * N + j]);
			}
		}
	} else {
		int W = N + 1;
		Mat_i IMat(M, W);
		IMat.fill(0);

		matches.reserve(M);
		for (int j = 0; j < N; j++) {
			double maxNcc = -1;
			int maxI = -1;
			double* pNcc = nccMat.data + j;
			double* pDisp = dispMat.data + j;
			for (int i = 0; i < M; i++) {
				if (pNcc[i * N] >= 0 && pDisp[i * N] >= 0 && pNcc[i * N] > maxNcc) {
					maxNcc = pNcc[i * N];
					maxI = i;
				}
			}
			if (maxI >= 0) {
				//record repeated assignment
				int *pIMat = IMat.data + maxI * W;
				int k = pIMat[0];
				pIMat[k + 1] = j;
				pIMat[0] = k + 1;
			}
		}

		//check the repeated assignment
		for (int i = 0; i < M; i++) {
			int KNum = IMat.data[i * W];
			if (KNum > 1) {
				int maxJ = -1;
				double maxNcc = -1;
				for (int k = 0; k < KNum; k++) {
					int j = IMat.data[i * W + k + 1];
					double ncc = nccMat.data[i * N + j];
					if (ncc > maxNcc) {
						maxNcc = ncc;
						maxJ = j;
					}
				}
				//find the best assignment
				matches.add(i, maxJ, nccMat.data[i * N + maxJ]);
			} else if (KNum == 1) {
				int j = IMat.data[i * W + 1];
				matches.add(i, j, nccMat.data[i * N + j]);
			}
		}
	}
}

void lapGuidedNCCMatch(const Mat_d& nccMat, const Mat_d& dispMat, Matching& matches, double largeWeight) {
	//TODO - lapGuidedNCCMatch
}
