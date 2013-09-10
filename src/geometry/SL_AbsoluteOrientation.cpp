/*
 * SL_AbsoluteOrientation.cpp
 *
 *  Created on: 2010-12-3
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#include "SL_AbsoluteOrientation.h"

#include "SL_error.h"
#include "SL_Quaternion.h"
#include "SL_RigidTransform.h"

#include "math/SL_LinAlg.h"
#include "geometry/SL_Geometry.h"
#include "extern/lapack.h"

#include <cassert>

static void aoNormCoord(const int npts, const double pts[], double pts_norm[], double* ptsc = 0) {
	int i;
	double xc = 0;
	double yc = 0;
	double zc = 0;
	const double *p = pts;
	for (i = 0; i < npts; ++i) {
		xc += p[0];
		yc += p[1];
		zc += p[2];
		p += 3;
	}

	xc /= npts;
	yc /= npts;
	zc /= npts;

	double *q = pts_norm;
	p = pts;
	for (i = 0; i < npts; ++i) {
		q[0] = p[0] - xc;
		q[1] = p[1] - yc;
		q[2] = p[2] - zc;
		p += 3;
		q += 3;
	}
	if (ptsc != 0) {
		ptsc[0] = xc;
		ptsc[1] = yc;
		ptsc[2] = zc;
	}
}
void aoGetMatrix(const int npts, const double r1[], const double r2[], double A[]) {
	double Sxx = 0, Sxy = 0, Sxz = 0, Syx = 0, Syy = 0, Syz = 0, Szx = 0, Szy = 0, Szz = 0;
	int i;
	const double *pr1 = r1;
	const double *pr2 = r2;
	for (i = 0; i < npts; ++i) {
		Sxx += pr1[0] * pr2[0];
		Sxy += pr1[0] * pr2[1];
		Sxz += pr1[0] * pr2[2];

		Syx += pr1[1] * pr2[0];
		Syy += pr1[1] * pr2[1];
		Syz += pr1[1] * pr2[2];

		Szx += pr1[2] * pr2[0];
		Szy += pr1[2] * pr2[1];
		Szz += pr1[2] * pr2[2];

		pr1 += 3;
		pr2 += 3;
	}

	A[0] = Sxx + Syy + Szz;
	A[1] = Syz - Szy;
	A[2] = Szx - Sxz;
	A[3] = Sxy - Syx;
	A[4] = Syz - Szy;
	A[5] = Sxx - Syy - Szz;
	A[6] = Sxy + Syx;
	A[7] = Szx + Sxz;
	A[8] = Szx - Sxz;
	A[9] = Sxy + Syx;
	A[10] = -Sxx + Syy - Szz;
	A[11] = Syz + Szy;
	A[12] = Sxy - Syx;
	A[13] = Szx + Sxz;
	A[14] = Syz + Szy;
	A[15] = -Sxx - Syy + Szz;
}
bool aoMaxEig(double A[], double v[]) {
	char jobvl = 'N';
	char jobvr = 'V';
	int n = 4;
	int lda = n;
	double wr[4], wi[4];
	int ldvl = 1;
	double vr[16];
	int ldvr = 4;
	double work[16];
	int lwork = 16;
	int info;

	//call lapack function
	dgeev_(&jobvl, &jobvr, &n, A, &lda, wr, wi, 0, &ldvl, vr, &ldvr, work, &lwork, &info);

	if (info < 0)
		warn("Error in call to dgeev (argument %d was invalid\n", -info);
	else if (info > 0)
		warn("Error: not all eigenvalues have converged\n");

	int i, max_i = -1;
	double max_val = 0;
	for (i = 0; i < 4; ++i) {
		if (wi[i] == 0.0) {
			if (wr[i] > max_val) {
				max_val = wr[i];
				max_i = i;
			}
		}
	}
	if (max_i < 0)
		return false;

	double* pvr = vr + 4 * max_i;
	v[0] = pvr[0];
	v[1] = pvr[1];
	v[2] = pvr[2];
	v[3] = pvr[3];
	return true;
}
bool absOrient(const int npts, const double pts1[], const double pts2[], double R[], double t[]) {
	assert(npts >= 3);

	double A[16], v[4];

	double* r1 = new double[npts * 3];
	double* r2 = new double[npts * 3];

	double pc1[3], pc2[3];
	aoNormCoord(npts, pts1, r1, pc1);
	aoNormCoord(npts, pts2, r2, pc2);

	aoGetMatrix(npts, r1, r2, A);
	bool res = aoMaxEig(A, v);
	if (!res)
		return false;

	quat2mat(v, R);
	mat33ProdVec(R, pc1, pc2, t, -1, 1);
	delete[] r1;
	delete[] r2;
	return true;
}
bool absOrient3(const double M1[],
		const double M2[],
		const double M3[],
		const double M1_[],
		const double M2_[],
		const double M3_[],
		double R[],
		double t[]) {
	double M[9], M_[9];
	M[0] = M1[0];
	M[1] = M1[1];
	M[2] = M1[2];

	M[3] = M2[0];
	M[4] = M2[1];
	M[5] = M2[2];

	M[6] = M3[0];
	M[7] = M3[1];
	M[8] = M3[2];

	M_[0] = M1_[0];
	M_[1] = M1_[1];
	M_[2] = M1_[2];

	M_[3] = M2_[0];
	M_[4] = M2_[1];
	M_[5] = M2_[2];

	M_[6] = M3_[0];
	M_[7] = M3_[1];
	M_[8] = M3_[2];

	return absOrient(3, M, M_, R, t);

}
static void aoNormCoordWeighted(const int npts,
		const double w[],
		const double pts[],
		double ptsWNorm[],
		double* wcenter = 0) {
	int i;
	double xc = 0;
	double yc = 0;
	double zc = 0;

	const double *p = pts;
	double sw = 0;
	for (i = 0; i < npts; ++i) {
		xc += w[i] * p[0];
		yc += w[i] * p[1];
		zc += w[i] * p[2];
		sw += w[i];
		p += 3;
	}

	xc /= sw;
	yc /= sw;
	zc /= sw;

	double *q = ptsWNorm;
	p = pts;
	for (i = 0; i < npts; ++i) {
		q[0] = p[0] - xc;
		q[1] = p[1] - yc;
		q[2] = p[2] - zc;
		p += 3;
		q += 3;
	}
	if (wcenter != 0) {
		wcenter[0] = xc;
		wcenter[1] = yc;
		wcenter[2] = zc;
	}
}
void aoGetMatrixWeighted(const int npts, const double w[], const double r1[], const double r2[], double M[]) {
	const double *pr1 = r1;
	const double *pr2 = r2;
	memset(M, 0, sizeof(double) * 9);
	for (int i = 0; i < npts; ++i) {
		M[0] += w[i] * pr1[0] * pr2[0];
		M[1] += w[i] * pr1[0] * pr2[1];
		M[2] += w[i] * pr1[0] * pr2[2];

		M[3] += w[i] * pr1[1] * pr2[0];
		M[4] += w[i] * pr1[1] * pr2[1];
		M[5] += w[i] * pr1[1] * pr2[2];

		M[6] += w[i] * pr1[2] * pr2[0];
		M[7] += w[i] * pr1[2] * pr2[1];
		M[8] += w[i] * pr1[2] * pr2[2];

		pr1 += 3;
		pr2 += 3;
	}
}
void aoGetOptRotation(const double M[], double R[]) {
	double U[9], S[9], VT[9];
	dgesvdFor(3, 3, M, U, S, VT);
	if (mat33Det(U) * mat33Det(VT) < 0) {
		VT[6] = -VT[6];
		VT[7] = -VT[7];
		VT[8] = -VT[8];
	}
	mat33AB(U, VT, R);
	assert(mat33Det(R) > 0);
}
void absOrientWeighted(const int npts,
		const double w[],
		const double pts1[],
		const double pts2[],
		double R[],
		double t[]) {
	assert(npts >= 3);

	double M[9];
	double* r1 = new double[npts * 3];
	double* r2 = new double[npts * 3];
	double pc1[3], pc2[3];

	aoNormCoordWeighted(npts, w, pts1, r1, pc1);
	aoNormCoordWeighted(npts, w, pts2, r2, pc2);
	aoGetMatrixWeighted(npts, w, r2, r1, M);
	aoGetOptRotation(M, R);
	mat33ProdVec(R, pc1, pc2, t, -1, 1);
	delete[] r1;
	delete[] r2;
}
void absOrientRobustTukey(const int npts, const double pts1[], const double pts2[], double R[], double t[], int maxIter) {
	assert(npts >= 3);
	double* w = new double[npts];
	double* d = new double[npts];

	for (int i = 0; i < npts; i++)
		w[i] = 1.0;

	double ad = 0, sv = 0, th = 0;
	for (int k = 0; k < maxIter; k++) {
		absOrientWeighted(npts, w, pts1, pts2, R, t);
		if (k == 0) {
			//get average residual
			ad = 0;
			for (int i = 0; i < npts; i++) {
				double M[3];
				rigidTrans(R, t, pts1 + 3 * i, M);
				d[i] = dist3(M, pts2 + 3 * i);
				ad += d[i];
			}
			ad /= npts;

			//get standard deviation;
			sv = 0;
			for (int i = 0; i < npts; i++) {
				double v = d[i] - ad;
				sv += v * v;
			}

			sv /= (npts - 1);
			sv = sqrt(sv);

			th = ad + 3 * sv;
		}
		//logInfo("-------(%g,%g) --- \n", ad, sv);
		for (int i = 0; i < npts; i++) {
			if (k > 0) {
				double M[3];
				rigidTrans(R, t, pts1 + 3 * i, M);
				d[i] = dist3(M, pts2 + 3 * i);
			}

			if (d[i] > th) {
				w[i] = 0;
			} else {
				double dr = d[i] / th;
				w[i] = (1 - dr * dr);
				w[i] *= w[i];
			}
			//test
			//logInfo("[%d]:%g %g\n", i, d[i], w[i]);
		}
	}
	delete[] w;
	delete[] d;
}

void absOrientRobustL1(const int npts, const double pts1[], const double pts2[], double R[], double t[], int maxIter) {
	assert(npts >= 3);
	double* w = new double[npts];

	for (int i = 0; i < npts; i++)
		w[i] = 1.0;

	for (int k = 0; k < maxIter; k++) {
		absOrientWeighted(npts, w, pts1, pts2, R, t);
		//get the residual and recompute the weights
		logInfo("-------\n");
		for (int i = 0; i < npts; i++) {
			double M[3];
			rigidTrans(R, t, pts1 + 3 * i, M);
			double d = dist3(M, pts2 + 3 * i);
			w[i] = 1.0 / (d + 1e-6);
			//test
			//logInfo("[%d]:%g %g\n", i, d, w[i]);
		}
	}
	delete[] w;
}
//int main() {
//	//		double pts1Org[9] = {
//	//				35.1659507062997,	83.0828627896291,	58.5264091152724,
//	//				54.9723608291140,	91.7193663829810,	28.5839018820374,
//	//				75.7200229110721,	75.3729094278495,	38.0445846975357
//	//		};
//	//		double pts2Org[9] = {
//	//				67.5373770698565,	6.25152185979056,	-57.3791975698678,
//	//				82.0162888895948,	-5.65584432344439,	-25.5670527384040,
//	//				87.0437762940195,	20.9761123604317,	-18.3117160973476,
//	//		};
//
//	//	double pts1Org[9] = { 52.8533, 16.5649, 60.1982, 26.2971, 65.4079, 68.9215, 74.8152, 45.0542, 8.3821, };
//	////
//	//	double pts2Org[9] = { 38.4520, 48.8583, 14.2344, 54.4161, 70.3194, -35.2786, 46.0438, -6.7249, -14.6056, };
//
//	//1.
//	//	double pts1Org[18] = { 189.085068511010, 109.998901532673, 284.436826892585, 154.858127091123, 67.7816369112954,
//	//			351.973306892472, 147.114143128385, 226.960925200777, 325.118434048043 };
//	//
//	//	double pts2Org[18] = { 232.265142680240, 311.168659014886, 11.2862256610414, 187.314293697254, 374.310683043443,
//	//			-27.5398359359440, 321.853312859641, 304.585607392789, -83.7258257535582 };
//	//	double R[9], t[3];
//	//
//	//	absOrient(3, pts1Org, pts2Org, R, t);
//	//
//	//	printMat(3, 3, R);
//	//	printMat(3, 1, t);
//
//	double R[9] = { 0.448027, 0.88651, 0.115638, 0.295684, -0.269001, 0.916629, 0.843707, -0.376483, -0.382647 };
//	double t[3] = { 17.1431, 24.1258, 2.00525 };
//
//	int nPts = 1000;
//	Mat_d pts1Org(nPts, 3), pts1Trans(nPts, 3);
//
//	for (int i = 0; i < nPts; i++) {
//		pts1Org.data[3 * i] = rand() % 1000;
//		pts1Org.data[3 * i + 1] = rand() % 1000;
//		pts1Org.data[3 * i + 2] = rand() % 1000;
//		rigidTrans(R, t, pts1Org.data + 3 * i, pts1Trans.data + 3 * i);
//		if (i > 500) {
//			pts1Trans.data[3 * i] += 1000;
//		}
//	}
//
//	double R1[9], t1[3];
//	//1.
//	//absOrient(nPts, pts1Org.data, pts1Trans.data, R1, t1);
//
//	//2.
//	//	Mat_d w(nPts, 1);
//	//	w.fill(1.0);
//	//	absOrientWeighted(nPts, w.data, pts1Org.data, pts1Trans.data, R1, t1);
//
//	//3.
//	absOrientRobustTukey(nPts, pts1Org.data, pts1Trans.data, R1, t1, 10);
//
//	printMat(3, 3, R1);
//	printMat(3, 1, t1);
//
//	//get the residual
//	for (int i = 0; i < nPts; i++) {
//		double M[3];
//		rigidTrans(R1, t1, pts1Org.data + 3 * i, M);
//		logInfo("[%d]:%g\n", i, dist3(M, pts1Trans.data + 3 * i));
//	}
//
//	return 0;
//}
