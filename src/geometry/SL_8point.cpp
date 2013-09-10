/*
 * SL_8point.cpp
 *
 *  Created on: 2010-11-14
 *      Author: Danping Zou
 */

#include "SL_error.h"
#include "math/SL_LinAlg.h"

#include "SL_8point.h"
#include "SL_5point.h"
#include "SL_Triangulate.h"
#include "SL_FundamentalMatrix.h"

#include <cfloat>
#include <cmath>
//have bugs 
void computeFMat8Pt(int n , const double* a , const double* b , double F[9]) {
	if (n < 7)
		repErr("more than 7 points are required.");
	double *A = new double[n * 9];
	double* S = n > 9 ? new double[n * 9] : new double[81];
	double* VT = new double[81];

	int i;
	for (i = 0; i < n; ++i) {
		double ax = a[2 * i];
		double ay = a[2 * i + 1];
		double az = 1;

		double bx = b[2 * i];
		double by = b[2 * i + 1];
		double bz = 1;

		A[i * 9] = ax * bx;
		A[i * 9 + 1] = ay * bx;
		A[i * 9 + 2] = az * bx;
		A[i * 9 + 3] = ax * by;
		A[i * 9 + 4] = ay * by;
		A[i * 9 + 5] = az * by;
		A[i * 9 + 6] = ax * bz;
		A[i * 9 + 7] = ay * bz;
		A[i * 9 + 8] = az * bz;
	}
	dgesvdFor(n, 9, A, S, VT);
	memcpy(F, VT + 8 * 9, sizeof(double) * 9);
	delete[] A;
	delete[] S;
	delete[] VT;
}
int findFMatRansac(
		const int n ,
		const double* a ,
		const double* b ,
		double* F ,
		int iterMaxNum ,
		double epiErrThes ,
		double* epiErr) {
	int k, indices[8];
	double pts1[16];
	double pts2[16];
	double tF[9];

	double errMin = DBL_MAX;
	int inlierMax = 0;

	for (k = 0; k < iterMaxNum; ++k) {
		randChoose(n, indices, 8);
		int i;
		for (i = 0; i < 8; ++i) {
			pts1[i * 2] = a[indices[i] * 2];
			pts1[i * 2 + 1] = a[indices[i] * 2 + 1];

			pts2[i * 2] = b[indices[i] * 2];
			pts2[i * 2 + 1] = b[indices[i] * 2 + 1];
		}
		computeFMat8Pt(8, pts2, pts1, tF);
		double err;
		int inlier = evaluateEMat(n, a, b, tF, epiErrThes, err);
		if (inlier > inlierMax || (inlier == inlierMax && err < errMin)) {
			inlierMax = inlier;
			errMin = err;
			memcpy(F, tF, sizeof(double) * 9);
		}
	}

	if (epiErr)
		*epiErr = errMin;
	return inlierMax;
}

