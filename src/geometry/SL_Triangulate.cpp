/*
 * SL_Triangulate.cpp
 *
 *  Created on: 2010-11-12
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#include "SL_Triangulate.h"
#include "SL_Distortion.h"
#include "SL_Geometry.h"
#include "math/SL_LinAlg.h"
#ifdef USE_CPPMINPACK
#include "cminpack.h"
#else
#include "extern/levmar/levmar.h"
#endif

#include <cmath>

void project(const double *R, const double *t, const double *M, double* m) {
	double zn = R[6] * M[0] + R[7] * M[1] + R[8] * M[2] + t[2];
	m[0] = (R[0] * M[0] + R[1] * M[1] + R[2] * M[2] + t[0]) / zn;
	m[1] = (R[3] * M[0] + R[4] * M[1] + R[5] * M[2] + t[1]) / zn;
}

void project(const double* K, const double* R, const double* t, const double* M,
		double* m) {
	double zn = R[6] * M[0] + R[7] * M[1] + R[8] * M[2] + t[2];
	double xn = (R[0] * M[0] + R[1] * M[1] + R[2] * M[2] + t[0]) / zn;
	double yn = (R[3] * M[0] + R[4] * M[1] + R[5] * M[2] + t[1]) / zn;

	double z = xn * K[6] + yn * K[7] + K[8];
	m[0] = (xn * K[0] + yn * K[1] + K[2]) / z;
	m[1] = (xn * K[3] + yn * K[4] + K[5]) / z;
}
void project(const double* K, const double* kc, const double* R,
		const double* t, const double *M, double* m) {
	double zn = R[6] * M[0] + R[7] * M[1] + R[8] * M[2] + t[2];
	double xn = (R[0] * M[0] + R[1] * M[1] + R[2] * M[2] + t[0]) / zn;
	double yn = (R[3] * M[0] + R[4] * M[1] + R[5] * M[2] + t[1]) / zn;

	double r = sqrt(xn * xn + yn * yn);
	double rr = r * r;
	double rrrr = rr * rr;
	double rrrrrr = rrrr * rr;

	double factor_a = (1 + kc[0] * rr + kc[1] * rrrr + kc[4] * rrrrrr);
	double factor_bx = 2 * kc[2] * xn * yn + kc[3] * (rr + 2 * xn * xn);
	double factor_by = kc[2] * (rr + 2 * yn * yn) + 2 * kc[3] * xn * yn;

	double xnd = xn * factor_a + factor_bx;
	double ynd = yn * factor_a + factor_by;

	double z = K[6] * xnd + K[7] * ynd + K[8];
	m[0] = (K[0] * xnd + K[1] * ynd + K[2]) / z;
	m[1] = (K[3] * xnd + K[4] * ynd + K[5]) / z;
}
void project(const double* K, const double* kc, const double* R,
		const double* t, int npts, const double* Ms, double* ms) {
	int i;
	for (i = 0; i < npts; ++i) {
		project(K, kc, R, t, Ms + 3 * i, ms + 2 * i);
	}
}
void project(const double* R, const double* t, int npts, const double* Ms,
		double* ms) {
	for (int i = 0; i < npts; ++i) {
		project(R, t, Ms + 3 * i, ms + 2 * i);
	}
}
void project(const double* K, const double* R, const double* t, int npts,
		const double* Ms, double* ms) {
	for (int i = 0; i < npts; ++i) {
		project(K, R, t, Ms + 3 * i, ms + 2 * i);
	}
}
double reprojErrorSingle(const double * K, const double* kc, const double* R,
		const double* t, const double* M, const double* m) {
	double ms[2];
	project(K, kc, R, t, M, ms);
	double dx = m[0] - ms[0];
	double dy = m[1] - ms[1];
	return sqrt(dx * dx + dy * dy);
}

double reprojErrorSingle(const double* K, const double* R, const double* t,
		const double* M, const double* m) {
	double ms[2];
	project(K, R, t, M, ms);
	double dx = m[0] - ms[0];
	double dy = m[1] - ms[1];
	return sqrt(dx * dx + dy * dy);
}

double reprojErrCov(const double* K, const double* R, const double* t,
		const double* M, const double cov[9], const double m[2]) {
	double var[4];
	double ivar[4];
	double m0[2];
	project(K, R, t, M, m0);
	getProjectionCovMat(K, R, t, M, cov, var, 1.0);
	mat22Inv(var, ivar);
	double dx = m0[0] - m[0];
	double dy = m0[1] - m[1];

	double s = ivar[0] * dx * dx + 2 * ivar[1] * dx * dy + ivar[3] * dy * dy;
	return sqrt(s);
}

double reprojError2(const double* K, const double* kc, const double* R,
		const double* t, int npts, const double* Ms, const double* ms0) {
	double ms[2];
	double err = 0;
	int i;
	for (i = 0; i < npts; ++i) {
		project(K, kc, R, t, Ms + 3 * i, ms);
		double dx = ms0[2 * i] - ms[0];
		double dy = ms0[2 * i + 1] - ms[1];
		err += dx * dx + dy * dy;
	}
	return err;
}

double reprojError(const double* K, const double* kc, const double* R,
		const double* t, int npts, const double* Ms, const double* ms0) {
	return sqrt(reprojError2(K, kc, R, t, npts, Ms, ms0) / npts);
}
double reprojErrorEach(const double* K, const double* R, const double * t,
		int nPts, const double* Ms, const double* ms0, double* err) {
	double rm[2];
	double s = 0;
	int i;
	for (i = 0; i < nPts; ++i) {
		project(K, R, t, Ms + 3 * i, rm);
		double dx = ms0[2 * i] - rm[0];
		double dy = ms0[2 * i + 1] - rm[1];
		err[i] = sqrt(dx * dx + dy * dy);
		s += err[i];
	}
	return s / nPts;
}
double reprojError2(const double* K, const double *R, const double *t, int npts,
		const double* Ms, const double* ms0) {
	double rm[2];
	double err = 0;
	int i;
	for (i = 0; i < npts; ++i) {
		project(K, R, t, Ms + 3 * i, rm);
		double dx = ms0[2 * i] - rm[0];
		double dy = ms0[2 * i + 1] - rm[1];
		err += dx * dx + dy * dy;
	}
	return err;
}
double reprojError(const double* K, const double *R, const double *t, int npts,
		const double* Ms, const double* ms0) {
	return sqrt(reprojError2(K, R, t, npts, Ms, ms0) / npts);
}
double reprojError2(const double* R, double *t, int npts, const double * Ms,
		const double* ms0) {
	double ms[2];
	double err = 0;
	int i;
	for (i = 0; i < npts; ++i) {
		project(R, t, Ms + 3 * i, ms);
		double dx = ms0[2 * i] - ms[0];
		double dy = ms0[2 * i + 1] - ms[1];
		err += dx * dx + dy * dy;
	}
	return err;
}
double reprojError(const double* R, double *t, int npts, const double * Ms,
		const double* ms0) {
	return sqrt(reprojError2(R, t, npts, Ms, ms0) / npts);
}

double checkReprojError(const double* K, const double* kc, const double* R,
		const double* t, int npts, const double* Ms, const double* ms0,
		char* flags, double thres) {
	double ms[2];
	thres *= thres;
	int i;
	double err_sum = 0;
	for (i = 0; i < npts; ++i) {
		project(K, kc, R, t, Ms + 3 * i, ms);
		double dx = ms0[2 * i] - ms[0];
		double dy = ms0[2 * i + 1] - ms[1];
		double err = dx * dx + dy * dy;
		if (err > thres)
			flags[i] = 0;
		else {
			flags[i] = 1;
			err_sum += err;
		}
	}
	return (err_sum / npts);
}
void computePerspectiveJacobian(const double* K, const double* R,
		const double* t, const double* X, double* fjac) {
	double Q[9];
	double T[3];
	mat33AB(K, R, Q);
	mat33ProdVec(K, t, T);

	double a1 = (T[0] + Q[0] * X[0] + Q[1] * X[1] + Q[2] * X[2]);
	double a2 = (T[1] + Q[3] * X[0] + Q[4] * X[1] + Q[5] * X[2]);
	double a3 = (T[2] + Q[6] * X[0] + Q[7] * X[1] + Q[8] * X[2]);
	double a3a3 = a3 * a3;

	fjac[0] = Q[0] / a3 - (Q[6] * a1) / a3a3;
	fjac[1] = Q[1] / a3 - (Q[7] * a1) / a3a3;
	fjac[2] = Q[2] / a3 - (Q[8] * a1) / a3a3;
	fjac[3] = Q[3] / a3 - (Q[6] * a2) / a3a3;
	fjac[4] = Q[4] / a3 - (Q[7] * a2) / a3a3;
	fjac[5] = Q[5] / a3 - (Q[8] * a2) / a3a3;
}
void computePerspectiveJacobian(const double* R, const double* t,
		const double* X, double* fjac) {
	double a1 = (t[0] + R[0] * X[0] + R[1] * X[1] + R[2] * X[2]);
	double a2 = (t[1] + R[3] * X[0] + R[4] * X[1] + R[5] * X[2]);
	double a3 = (t[2] + R[6] * X[0] + R[7] * X[1] + R[8] * X[2]);
	double a3a3 = a3 * a3;
	fjac[0] = R[0] / a3 - (R[6] * a1) / a3a3;
	fjac[1] = R[1] / a3 - (R[7] * a1) / a3a3;
	fjac[2] = R[2] / a3 - (R[8] * a1) / a3a3;
	fjac[3] = R[3] / a3 - (R[6] * a2) / a3a3;
	fjac[4] = R[4] / a3 - (R[7] * a2) / a3a3;
	fjac[5] = R[5] / a3 - (R[8] * a2) / a3a3;
}

static void residual2View(double* X, double* fvec, int m, int n, void* adata) {
	CamProjParam * param = (CamProjParam*) adata;
	if (param[0].K == 0) {
		//without intrinsic matrix
		/* compute the reprojection error */
		double p[2], q[2];
		project(param[0].R, param[0].t, X, p);
		project(param[1].R, param[1].t, X, q);

		fvec[0] = param[0].m[0] - p[0];
		fvec[1] = param[0].m[1] - p[1];
		fvec[2] = param[1].m[0] - q[0];
		fvec[3] = param[1].m[1] - q[1];
	} else {
		/* compute the reprojection error */
		double p[2], q[2];
		project(param[0].K, param[0].R, param[0].t, X, p);
		project(param[1].K, param[1].R, param[1].t, X, q);

		fvec[0] = param[0].m[0] - p[0];
		fvec[1] = param[0].m[1] - p[1];
		fvec[2] = param[1].m[0] - q[0];
		fvec[3] = param[1].m[1] - q[1];
	}
}

static void residualMultiView(double* X, double* fvec, int m, int n,
		void* adata) {
	CamProjParam * param = (CamProjParam*) adata;
	if (!param->K) {
		int nview = n / 2;
		double* pfvec = fvec;
		for (int i = 0; i < nview; i++) {
			double rm[2];
			project(param[i].R, param[i].t, X, rm);
			pfvec[0] = param[i].m[0] - rm[0];
			pfvec[1] = param[i].m[1] - rm[1];
			pfvec += 2;
		}
	} else {
		int nview = n / 2;
		double* pfvec = fvec;
		for (int i = 0; i < nview; i++) {
			double rm[2];
			project(param[i].K, param[i].R, param[i].t, X, rm);
			pfvec[0] = param[i].m[0] - rm[0];
			pfvec[1] = param[i].m[1] - rm[1];
			pfvec += 2;
		}
	}
}
//static void residualMultiViewJacobi(double *X, double *fjac, int m, int n, void *adata) {
//	CamProjParam* param = (CamProjParam*) adata;
//	int nview = n / 2;
//	if (!param->K) {
//		for (int i = 0; i < nview; i++) {
//			computePerspectiveJacobian(param[i].R, param[i].t, X, fjac);
//			fjac += 6;
//		}
//	} else {
//		double* pfjac = fjac;
//		for (int i = 0; i < nview; i++) {
//			computePerspectiveJacobian(param[i].K, param[i].R, param[i].t, X, fjac);
//			pfjac += 6;
//		}
//	}
//}

void binTriangulateLin(const double* R0, const double *t0, const double *R1,
		const double *t1, const double *p, const double *q, double *x) {
	double A[12];
	double b[4];

	A[0] = R0[0] - p[0] * R0[6];
	A[1] = R0[1] - p[0] * R0[7];
	A[2] = R0[2] - p[0] * R0[8];

	A[3] = R0[3] - p[1] * R0[6];
	A[4] = R0[4] - p[1] * R0[7];
	A[5] = R0[5] - p[1] * R0[8];

	A[6] = R1[0] - q[0] * R1[6];
	A[7] = R1[1] - q[0] * R1[7];
	A[8] = R1[2] - q[0] * R1[8];

	A[9] = R1[3] - q[1] * R1[6];
	A[10] = R1[4] - q[1] * R1[7];
	A[11] = R1[5] - q[1] * R1[8];

	b[0] = t0[2] * p[0] - t0[0];
	b[1] = t0[2] * p[1] - t0[1];
	b[2] = t1[2] * q[0] - t1[0];
	b[3] = t1[2] * q[1] - t1[1];

	/* Find the least squares result */
	dgelsyFor(4, 3, 1, A, b, x);
}

void binTriangulateLin(const double K0[9], const double R0[9],
		const double t0[3], const double K1[9], const double R1[9],
		const double t1[3], const double p[2], const double q[2], double x[3]) {

	double iK0[9], iK1[9], m0[2], m1[2];
	mat33Inv(K0, iK0);
	mat33Inv(K1, iK1);
	normPoint(iK0, p, m0);
	normPoint(iK1, q, m1);
	binTriangulateLin(R0, t0, R1, t1, m0, m1, x);
}
void binTriangulate(const double R0[9], const double t0[3], const double R1[9],
		const double t1[3], const double m0[2], const double m1[2],
		double x[3]) {

	double A[12];
	double b[4];

	A[0] = R0[0] - m0[0] * R0[6];
	A[1] = R0[1] - m0[0] * R0[7];
	A[2] = R0[2] - m0[0] * R0[8];

	A[3] = R0[3] - m0[1] * R0[6];
	A[4] = R0[4] - m0[1] * R0[7];
	A[5] = R0[5] - m0[1] * R0[8];

	A[6] = R1[0] - m1[0] * R1[6];
	A[7] = R1[1] - m1[0] * R1[7];
	A[8] = R1[2] - m1[0] * R1[8];

	A[9] = R1[3] - m1[1] * R1[6];
	A[10] = R1[4] - m1[1] * R1[7];
	A[11] = R1[5] - m1[1] * R1[8];

	b[0] = t0[2] * m0[0] - t0[0];
	b[1] = t0[2] * m0[1] - t0[1];
	b[2] = t1[2] * m1[0] - t1[0];
	b[3] = t1[2] * m1[1] - t1[1];

	/* Find the least squares result */
	dgelsyFor(4, 3, 1, A, b, x);

	/* Run a non-linear optimization to refine the result */
	CamProjParam param[2];
	param[0].set(0, R0, t0, m0);
	param[1].set(0, R1, t1, m1);

	const int maxIter = 5;
	double retInfo[LM_INFO_SZ];
	dlevmar_dif(residual2View, x, 0, 3, 4, maxIter, 0, retInfo, 0, 0, param);
}

void binTriangulate(const double K0[9], const double R0[9], const double t0[3],
		const double K1[9], const double R1[9], const double t1[3],
		const double m0[2], const double m1[2], double x[3]) {

	double Q0[9], q0[3], Q1[9], q1[3];

	mat33AB(K0, R0, Q0);
	mat33ProdVec(K0, t0, q0);

	mat33AB(K1, R1, Q1);
	mat33ProdVec(K1, t1, q1);

	double A[12];
	double b[4];

	A[0] = Q0[0] - m0[0] * Q0[6];
	A[1] = Q0[1] - m0[0] * Q0[7];
	A[2] = Q0[2] - m0[0] * Q0[8];

	A[3] = Q0[3] - m0[1] * Q0[6];
	A[4] = Q0[4] - m0[1] * Q0[7];
	A[5] = Q0[5] - m0[1] * Q0[8];

	A[6] = Q1[0] - m1[0] * Q1[6];
	A[7] = Q1[1] - m1[0] * Q1[7];
	A[8] = Q1[2] - m1[0] * Q1[8];

	A[9] = Q1[3] - m1[1] * Q1[6];
	A[10] = Q1[4] - m1[1] * Q1[7];
	A[11] = Q1[5] - m1[1] * Q1[8];

	b[0] = q0[2] * m0[0] - q0[0];
	b[1] = q0[2] * m0[1] - q0[1];
	b[2] = q1[2] * m1[0] - q1[0];
	b[3] = q1[2] * m1[1] - q1[1];

	/* Find the least squares result */
	dgelsyFor(4, 3, 1, A, b, x);

	/* Run a non-linear optimization to refine the result */
	CamProjParam param[2];
	param[0].set(K0, R0, t0, m0);
	param[1].set(K1, R1, t1, m1);

	const int maxIter = 5;
	double retInfo[LM_INFO_SZ];
	dlevmar_dif(residual2View, x, 0, 3, 4, maxIter, 0, retInfo, 0, 0, param);

}

void binTriangulatePoints(const double R0[9], const double t0[3],
		const double R1[9], const double t1[3], int npts, const double p[],
		const double q[], double x[]) {
	for (int i = 0; i < npts; i++)
		binTriangulate(R0, t0, R1, t1, p + 2 * i, q + 2 * i, x + 3 * i);
}

void binTriangulatePoints(const double K0[9], const double R0[9],
		const double t0[3], const double K1[9], const double R1[9],
		const double t1[3], int npts, const double p[], const double q[],
		double x[]) {
	for (int i = 0; i < npts; i++)
		binTriangulate(K0, R0, t0, K1, R1, t1, p + 2 * i, q + 2 * i, x + 3 * i);
}

void triangulateMultiView(int nview, const double* Rs, const double* ts,
		const double nms[], double x[]) {

	double * A = new double[nview * 2 * 3];
	double * b = new double[nview * 2];

	double* pA = A;
	double* pb = b;
	const double* pR = Rs;
	const double* pT = ts;

	for (int i = 0; i < nview; i++) {
		double u = nms[2 * i];
		double v = nms[2 * i + 1];
		pA[0] = u * pR[6] - pR[0];
		pA[1] = u * pR[7] - pR[1];
		pA[2] = u * pR[8] - pR[2];

		pA[3] = v * pR[6] - pR[3];
		pA[4] = v * pR[7] - pR[4];
		pA[5] = v * pR[8] - pR[5];

		pb[0] = pT[0] - u * pT[2];
		pb[1] = pT[1] - v * pT[2];

		pA += 6;
		pb += 2;
		pR += 9;
		pT += 3;
	}
	/* Find the least squares result */
	dgelsyFor(nview * 2, 3, 1, A, b, x);
	CamProjParam* param = new CamProjParam[nview];
	for (int i = 0; i < nview; i++) {
		param[i].set(0, Rs + i * 9, ts + i * 3, nms + 2 * i);
	}

	const int maxIter = 5;
	double retInfo[LM_INFO_SZ * 2];
	dlevmar_dif(residualMultiView, x, 0, 3, nview * 2, maxIter, 0, retInfo, 0,
			0, param);

	delete[] A;
	delete[] b;
	delete[] param;
}

void triangulateMultiView(int nview, const double* Rs[], const double* ts[],
		const double* ms[], double x[]) {
	double * A = new double[nview * 2 * 3];
	double * b = new double[nview * 2];
	double* pA = A;
	double* pb = b;

	for (int i = 0; i < nview; i++) {
		double u = ms[i][0];
		double v = ms[i][1];

		pA[0] = u * Rs[i][6] - Rs[i][0];
		pA[1] = u * Rs[i][7] - Rs[i][1];
		pA[2] = u * Rs[i][8] - Rs[i][2];

		pA[3] = v * Rs[i][6] - Rs[i][3];
		pA[4] = v * Rs[i][7] - Rs[i][4];
		pA[5] = v * Rs[i][8] - Rs[i][5];

		pb[0] = ts[i][0] - u * ts[i][2];
		pb[1] = ts[i][1] - v * ts[i][2];

		pA += 6;
		pb += 2;
	}

	/* Find the least squares result */
	dgelsyFor(nview * 2, 3, 1, A, b, x);
	CamProjParam* param = new CamProjParam[nview];
	for (int i = 0; i < nview; i++) {
		param[i].set(0, Rs[i], ts[i], ms[i]);
	}

	const int maxIter = 5;
	double retInfo[LM_INFO_SZ];
	dlevmar_dif(residualMultiView, x, 0, 3, nview * 2, maxIter, 0, retInfo, 0,
			0, param);

	delete[] A;
	delete[] b;
	delete[] param;
}
void triangulateMultiViewLin(int nview, const double* R, const double* t,
		const double pts2d[], double x[]) {
	double * A = new double[nview * 2 * 3];
	double * b = new double[nview * 2];

	double* pA = A;
	double* pb = b;
	const double* pR = R;
	const double* pT = t;
	for (int i = 0; i < nview; i++) {
		double u = pts2d[2 * i];
		double v = pts2d[2 * i + 1];
		pA[0] = u * pR[6] - pR[0];
		pA[1] = u * pR[7] - pR[1];
		pA[2] = u * pR[8] - pR[2];

		pA[3] = v * pR[6] - pR[3];
		pA[4] = v * pR[7] - pR[4];
		pA[5] = v * pR[8] - pR[5];

		pb[0] = pT[0] - u * pT[2];
		pb[1] = pT[1] - v * pT[2];

		pA += 6;
		pb += 2;
		pR += 9;
		pT += 3;
	}

	/* Find the least squares result */
	dgelsyFor(nview * 2, 3, 1, A, b, x);
	delete[] A;
	delete[] b;
}
void refineTriangulation(int nview, const double* K, const double* R,
		const double* t, const double pts2d[], double x[], double* err) {

	CamProjParam* param = new CamProjParam[nview];
	for (int i = 0; i < nview; i++) {
		param[i].set(K == 0 ? 0 : K + i * 9, R + i * 9, t + i * 3,
				pts2d + 2 * i);
	}

	const int maxIter = 5;
	double retInfo[LM_INFO_SZ];
	dlevmar_dif(residualMultiView, x, 0, 3, nview * 2, maxIter, 0, retInfo, 0,
			0, param);

	if (err) {
		*err = 0;
		const double* pK = K;
		const double* pR = R;
		const double* pT = t;
		double pt[2];
		for (int i = 0; i < nview; i++) {
			if (pK)
				project(pK, pR, pT, x, pt);
			else
				project(pR, pT, x, pt);
			double dx = pts2d[2 * i] - pt[0];
			double dy = pts2d[2 * i + 1] - pt[1];
			*err += dx * dx + dy * dy;
			pK += 9;
			pR += 9;
			pT += 3;
		}
		*err = sqrt(*err / nview);
	}
	delete[] param;
}
void getTriangulateCovMat(int nview, const double* K, const double* R,
		const double* t, const double x[3], double cov[9], double sigma) {
	double * J = new double[nview * 12];
	for (int i = 0; i < nview; i++) {
		computePerspectiveJacobian(K + 9 * i, R + 9 * i, t + 3 * i, x,
				J + 6 * i);
	}
	double tmp[9];
	matATB(nview * 2, 3, nview * 2, 3, J, J, tmp);
	matInv(3, tmp, cov);

	sigma *= sigma;

	cov[0] *= sigma;
	cov[1] *= sigma;
	cov[2] *= sigma;
	cov[3] *= sigma;
	cov[4] *= sigma;
	cov[5] *= sigma;
	cov[6] *= sigma;
	cov[7] *= sigma;
	cov[8] *= sigma;

	delete[] J;
}

void getTriangulateCovMat(int nview, const double** Ks, const double** Rs,
		const double** ts, const double x[3], double cov[9], double sigma) {
	double * J = new double[nview * 6];
	for (int i = 0; i < nview; i++) {
		computePerspectiveJacobian(Ks[i], Rs[i], ts[i], x, J + 6 * i);
	}
	double tmp[9];
	matATB(nview * 2, 3, nview * 2, 3, J, J, tmp);
	matInv(3, tmp, cov);

	sigma *= sigma;

	cov[0] *= sigma;
	cov[1] *= sigma;
	cov[2] *= sigma;
	cov[3] *= sigma;
	cov[4] *= sigma;
	cov[5] *= sigma;
	cov[6] *= sigma;
	cov[7] *= sigma;
	cov[8] *= sigma;

	delete[] J;
}
void getBinTriangulateCovMat(const double* K1, const double* R1,
		const double* t1, const double* K2, const double* R2, const double* t2,
		double x[3], double cov[9], double sigma) {
	double J[12];
	computePerspectiveJacobian(K1, R1, t1, x, J);
	computePerspectiveJacobian(K2, R2, t2, x, J + 6);
	double tmp[9];
	matATB(4, 3, 4, 3, J, J, tmp);
	matInv(3, tmp, cov);

	sigma *= sigma;

	cov[0] *= sigma;
	cov[1] *= sigma;
	cov[2] *= sigma;
	cov[3] *= sigma;
	cov[4] *= sigma;
	cov[5] *= sigma;
	cov[6] *= sigma;
	cov[7] *= sigma;
	cov[8] *= sigma;

}

void getProjectionCovMat(const double* K, const double* R, const double* t,
		const double x[3], const double cov[9], double var[4], double sigma) {
	double J[6];
	computePerspectiveJacobian(K, R, t, x, J);
	double tmp[6];
	matABT(3, 3, 2, 3, cov, J, tmp);
	matAB(2, 3, 3, 2, J, tmp, var);
	var[0] += sigma * sigma;
	var[3] += sigma * sigma;
}
bool isAtCameraBack(const double* R, const double* t, const double x[]) {
	double q[3];
	mat33ProdVec(R, x, t, q, 1.0, 1.0);
	if (q[2] < 0)
		return true;
	return false;
}

void seqTriangulate(const double K[9], const double R[9], const double t[3],
		const double m[2], double M[3], double cov[3], double sigma) {

	double rm[2], dm[2];
	project(K, R, t, M, rm);
	dm[0] = m[0] - rm[0];
	dm[1] = m[1] - rm[1];
	//compute Jacobian:
	double J[6];
	computePerspectiveJacobian(K, R, t, M, J);

	//compute 2D covariance:
	double JC[6], S[4], iS[4];

	//JC = J*cov
	//S = JC*J'
	// S = J*cov*J'
	matAB(2, 3, 3, 3, J, cov, JC);
	matABT(2, 3, 2, 3, JC, J, S);
	double sigma2 = sigma * sigma;
	S[0] += sigma2;
	S[3] += sigma2;

	//iS = (sigma+S)^-1
	matInv(2, S, iS);

	//compute Kalman gain:
	double CJT[6], G[6];
	//CJT = cov*J'
	//G = cov*J'*iS = cov*J'*(sigma+J*cov*J')^-1
	matABT(3, 3, 2, 3, cov, J, CJT);
	matAB(3, 2, 2, 2, CJT, iS, G);

	//update 3D position
	double dM[3];
	matAB(3, 2, 2, 1, G, dm, dM);
	M[0] += dM[0];
	M[1] += dM[1];
	M[2] += dM[2];

	//update 3D covariance
	//cov - G*J^T*cov
	double GJ[9], GJC[9];
	matAB(3, 2, 2, 3, G, J, GJ);
	mat33AB(GJ, cov, GJC);
	mat33Diff(cov, GJC, cov);
}

void getCameraCenter(const double R[9], const double t[3], double org[3]) {
	mat33TransProdVec(R, t, org);
	org[0] = -org[0];
	org[1] = -org[1];
	org[2] = -org[2];
}
void getCameraCenterAxes(const double R[9], const double t[3], double org[3],
		double xp[3], double yp[3], double zp[3]) {
	getCameraCenter(R, t, org);

	xp[0] = R[0];
	xp[1] = R[1];
	xp[2] = R[2];

	yp[0] = R[3];
	yp[1] = R[4];
	yp[2] = R[5];

	zp[0] = R[6];
	zp[1] = R[7];
	zp[2] = R[8];
}
double getCameraDistance(const double R1[9], const double t1[3],
		const double R2[9], const double t2[3]) {
	double org1[3], org2[3];
	getCameraCenter(R1, t1, org1);
	getCameraCenter(R2, t2, org2);
	return dist3(org1, org2);
}

void CamCoord2WoldCoord(const double* R, const double* t, const double mc[3],
		double mw[3]) {
	double dmc[3] = { mc[0] - t[0], mc[1] - t[1], mc[2] - t[2] };
	mat33TransProdVec(R, dmc, mw);
}

void getZoomedIntrinsic(const double* K0, double* K, int w0, int h0, int w,
		int h) {
	memset(K, 0, sizeof(double) * 9);
	K[8] = 1;
	K[2] = K0[2] * w / w0;
	K[5] = K0[5] * h / h0;

	K[0] = K0[0] * w / w0;
	K[4] = K0[4] * h / h0;
}