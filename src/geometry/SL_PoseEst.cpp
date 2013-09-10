/*
 * SL_PoseEst.cpp
 *
 *  Created on: 2011-8-9
 *      Author: Danping Zou
 */

#include "SL_PoseEst.h"
#include "extern/lapack.h"
#include "SL_Triangulate.h"
#include "SL_Geometry.h"

extern "C" {
#include <atlas/cblas.h>
}
#include <pthread.h>
#include <cassert>
#include <algorithm>
#include <cfloat>
#include <cmath>

#include "math/SL_LinAlgWarper.h"
#include "SL_Quaternion.h"
inline void _perspectiveProj(const double* K, const double* R, const double* t,
		const double* M, double* m) {
	double zn = R[6] * M[0] + R[7] * M[1] + R[8] * M[2] + t[2];
	double xn = (R[0] * M[0] + R[1] * M[1] + R[2] * M[2] + t[0]) / zn;
	double yn = (R[3] * M[0] + R[4] * M[1] + R[5] * M[2] + t[1]) / zn;

	double z = xn * K[6] + yn * K[7] + K[8];
	m[0] = (xn * K[0] + yn * K[1] + K[2]) / z;
	m[1] = (xn * K[3] + yn * K[4] + K[5]) / z;
}

inline void _computeSO3ExpMap(const double w[3], double R[9]) {
	double theta = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
	if (theta == 0) {
		memset(R, 0, sizeof(double) * 9);
		R[0] = 1.0;
		R[4] = 1.0;
		R[8] = 1.0;
	} else {
		double hw[3] = { w[0] / theta, w[1] / theta, w[2] / theta };
		double st = sin(theta);
		double ct = 1 - cos(theta);

		double hw0hw0 = hw[0] * hw[0];
		double hw0hw1 = hw[0] * hw[1];
		double hw0hw2 = hw[0] * hw[2];
		double hw1hw1 = hw[1] * hw[1];
		double hw1hw2 = hw[1] * hw[2];
		double hw2hw2 = hw[2] * hw[2];
		R[0] = -ct * hw1hw1 - ct * hw2hw2 + 1;
		R[1] = ct * hw0hw1 - st * hw[2];
		R[2] = st * hw[1] + ct * hw0hw2;
		R[3] = st * hw[2] + ct * hw0hw1;
		R[4] = -ct * hw0hw0 - ct * hw2hw2 + 1;
		R[5] = ct * hw1hw2 - st * hw[0];
		R[6] = ct * hw0hw2 - st * hw[1];
		R[7] = st * hw[0] + ct * hw1hw2;
		R[8] = -ct * hw0hw0 - ct * hw1hw1 + 1;
	}
}
inline void _mat33AB(const double* A, const double *B, double * C) {
	C[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
	C[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
	C[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];

	C[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
	C[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
	C[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];

	C[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
	C[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
	C[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
}
/* x = s*A*x*/
inline void _mat22Ax(double s, const double* A, double* x) {
	x[0] = s * (A[0] * x[0] + A[1] * x[1]);
	x[1] = s * (A[2] * x[0] + A[3] * x[1]);
}
inline void _mat33Ax(double s, const double* A, double* x) {
	x[0] = s * (A[0] * x[0] + A[1] * x[1] + A[2] * x[2]);
	x[1] = s * (A[3] * x[0] + A[4] * x[1] + A[5] * x[2]);
	x[2] = s * (A[6] * x[0] + A[7] * x[1] + A[8] * x[2]);
}
inline void _mat33Inv(const double *A, double *Ainv) {
	double m1 = A[8] * A[4] - A[7] * A[5];
	double m2 = A[8] * A[1] - A[7] * A[2];
	double m3 = A[5] * A[1] - A[4] * A[2];

	double d = A[0] * m1 - A[3] * m2 + A[6] * m3;

	Ainv[0] = m1 / d;
	Ainv[1] = -m2 / d;
	Ainv[2] = m3 / d;

	Ainv[3] = (-A[8] * A[3] + A[6] * A[5]) / d;
	Ainv[4] = (A[8] * A[0] - A[6] * A[2]) / d;
	Ainv[5] = (-A[5] * A[0] + A[3] * A[2]) / d;

	Ainv[6] = (A[7] * A[3] - A[6] * A[4]) / d;
	Ainv[7] = (-A[7] * A[0] + A[6] * A[1]) / d;
	Ainv[8] = (A[4] * A[0] - A[3] * A[1]) / d;
}
void _mat22Inv(const double* A, double* Ainv) {
	double s = (A[0] * A[3] - A[1] * A[2]);
	Ainv[0] = A[3] / s;
	Ainv[1] = -A[1] / s;
	Ainv[2] = -A[2] / s;
	Ainv[3] = A[0] / s;
}
static void _disturbSO3(const double R0[9], double dw1, double dw2, double dw3,
		double R1[9]) {
	double w[3] = { dw1, dw2, dw3 };
	double dR[9];
	_computeSO3ExpMap(w, dR);
	_mat33AB(R0, dR, R1);
}

inline void _disturbVecX(const double T0[3], double rel_eps, double T1[3],
		double& abs_eps) {
	abs_eps = T0[0] * rel_eps;
	abs_eps = abs_eps < 1e-5 ? 1e-5 : abs_eps;
	T1[0] = T0[0] + abs_eps;
	T1[1] = T0[1];
	T1[2] = T0[2];
}
inline void _disturbVecY(const double T0[3], double rel_eps, double T1[3],
		double& abs_eps) {
	abs_eps = T0[1] * rel_eps;
	abs_eps = abs_eps < 1e-5 ? 1e-5 : abs_eps;
	T1[0] = T0[0];
	T1[1] = T0[1] + abs_eps;
	T1[2] = T0[2];
}
inline void _disturbVecZ(const double T0[3], double rel_eps, double T1[3],
		double& abs_eps) {
	abs_eps = T0[2] * rel_eps;
	abs_eps = abs_eps < 1e-5 ? 1e-5 : abs_eps;
	T1[0] = T0[0];
	T1[1] = T0[1];
	T1[2] = T0[2] + abs_eps;
}
inline void _computeJacobi(size_t n, double* x0, double* x1, double eps,
		double* J) {
	for (size_t i = 0; i < n; i++) {
		J[i] = (x1[i] - x0[i]) / eps;
	}
}

inline double _vecNorm(size_t n, double * x) {
	return cblas_dnrm2(n, x, 1);
}
inline void _subMatCopy(Mat_d& mat, int r, int c, int nr, int nc,
		const double* data) {
	assert( r>=0 && c >= 0 && r+nr <= mat.rows && c + nc <= mat.cols);
	double* ptr = mat.data + r * mat.cols + c;
	for (int i = 0; i < nr; i++) {
		memcpy(ptr + i * mat.cols, data + i * nc, sizeof(double) * nc);
	}
}
static void _cholDecomp33(const Mat_d& cov, Mat_d& cholCov) {
	assert(cov.rows == 3 && cov.cols == 3);
	cholCov.cloneFrom(cov);
	char ch = 'L';
	int n = 3;
	int lda = 3;
	int info = 0;
	dpotrf_(&ch, &n, cholCov.data, &lda, &info);
	if (info < 0) {
		repErr("_cholDecomp33 failed! %d-th argument is invalid!", info);
	} else if (info > 0) {
		repErr("_cholDecomp33 failed! cannot complete the decomposition");
	}
	cholCov(1, 0) = 0;
	cholCov(2, 0) = 0;
	cholCov(2, 1) = 0;
}
static void _cholDecomp22(const Mat_d& cov, Mat_d& cholCov) {
	assert(cov.rows == 2 && cov.cols == 2);
	cholCov.cloneFrom(cov);
	char ch = 'L';
	int n = 2;
	int lda = 2;
	int info = 0;
	dpotrf_(&ch, &n, cholCov.data, &lda, &info);
	if (info < 0) {
		repErr("_cholDecomp22 failed! %d-th argument is invalid!", info);
	} else if (info > 0) {
		repErr("_cholDecomp22 failed! cannot complete the decomposition");
	}
	cholCov(1, 0) = 0;
}

PoseEst::PoseEst() :
		m_wSmooth(1), m_damp(1.0) {
	// TODO Auto-generated constructor stub
}

PoseEst::~PoseEst() {
	// TODO Auto-generated destructor stub
}
void PoseEst::allocTmpVars() {
	M = 2 * (nMs + nMd) + 3 * nC;
	N = 6 * nC + 3 * nD;
	m_covChol.resize(nC, 9);
	m_fX.resize(M, 1);
	m_resid.resize(M, 1);
	m_deltaX.resize(N, 1);
	m_measIndSz.reserve((nS + nD) * 2);
	m_pts.reserve((nS + nD) * 2);

	m_ptInd.reserve(M);
	m_meas.reserve(M);
	m_sqrtw.reserve(M);

	size_t k = 0;

	m_pts.clear();
	m_meas.clear();
	m_ptInd.clear();
	m_measIndSz.size();

	//static map points
	for (size_t i = 0; i < m_staPts.size(); i++) {
		m_pts.push_back(m_staPts[i]);
		m_measIndSz.push_back(make_pair(k, m_staMeas[i].size()));
		for (size_t j = 0; j < m_staMeas[i].size(); j++) {
			m_ptInd[k] = i;
			m_meas.push_back(m_staMeas[i][j]);
			k++;
		}
	}
	//dynamic map points
	for (size_t i = 0; i < m_dynPts.size(); i++) {
		m_pts.push_back(m_dynPts[i]);
		m_measIndSz.push_back(make_pair(k, m_dynMeas[i].size()));
		for (size_t j = 0; j < m_dynMeas[i].size(); j++) {
			m_ptInd[k] = i + nS;
			m_meas.push_back(m_dynMeas[i][j]);
			k++;
		}
	}
	m_sqrtw.clear();
	size_t nM = nMs + nMd;
	for (size_t i = 0; i < nM; i++) {
		m_sqrtw.push_back(sqrt(m_meas[i].w));
	}

	//compute the cholesky decomposition of the inverse of prediction covariance
	Mat_d invCov(3, 3);
	m_cholInvPredCov.assign(nC, Mat_d());
	for (size_t c = 0; c < nC; c++) {
		_mat33Inv(m_predCov[c].data, invCov.data);
		_cholDecomp33(invCov, m_cholInvPredCov[c]);
	}
	if (m_bUseDefMeasCov) {
		for (size_t i = 0; i < nM; i++) {
			memcpy(m_meas[i].cov, m_sigma, sizeof(double) * 4);
		}
	}

	//compute the cholesky decomposition of the inverse of measurement covariances
	invCov.resize(2, 2);
	m_measCholInvCov.assign(nM, Mat_d());
	for (size_t i = 0; i < nM; i++) {
		_mat22Inv(m_meas[i].cov, invCov.data);
		_cholDecomp22(invCov, m_measCholInvCov[i]);
	}

	//	//compute the weight for trajectory smoothness
	//	assert(m_smoothness >= 0 && m_smoothness <1.0);
	//	double wp = 0;
	//	for (size_t i = 0; i < nM; i++) {
	//		wp += m_meas[i].w;
	//	}
	//	m_wSmooth = m_smoothness * wp / (1 - m_smoothness);

	//	//test
	//	printf("m_wSmooth:%g\n",m_wSmooth);

	m_Jc.resize(nM + nC, 1);
	for (size_t i = 0; i < nM; i++)
		m_Jc[i].resize(2, 6);

	for (size_t i = nC; i < nM + nC; i++)
		m_Jc[i].resize(3, 6);

	m_Jd.resize(nM + nC, 1);
	for (size_t i = nMs; i < nM; i++)
		m_Jd[i].resize(2, 3);

	m_W.resize(nC, nD);
	for (size_t i = 0; i < nC; i++) {
		for (size_t j = 0; j < nD; j++) {
			m_W.data[i * nD + j].resize(6, 3);
		}
	}

	m_Y.resize(nC, nD);
	for (size_t i = 0; i < nC; i++) {
		for (size_t j = 0; j < nD; j++) {
			m_Y.data[i * nD + j].resize(6, 3);
		}
	}
	m_U.resize(nC, 1);
	for (size_t i = 0; i < nC; i++) {
		m_U.data[i].resize(6, 6);
	}

	m_V.resize(nD, 1);
	for (size_t i = 0; i < nD; i++) {
		m_V.data[i].resize(3, 3);
	}

	m_iV.resize(nD, 1);
	for (size_t i = 0; i < nD; i++) {
		m_iV.data[i].resize(3, 3);
	}

	m_S.resize(6 * nC, 6 * nC);
	m_JRes.resize(6 * nC + 3 * nD, 1);

	tmpRs.resize(nC, Mat_d());
	tmpTs.resize(nC, Mat_d());
	for (size_t c = 0; c < nC; c++) {
		tmpRs[c].resize(3, 3);
		tmpTs[c].resize(3, 1);
	}
	tmpPts.resize(m_pts.size(), Point3d());

	double s = 0;
	for (size_t i = 0; i < m_pts.size(); i++) {
		s += m_pts[i].x * m_pts[i].x;
		s += m_pts[i].y * m_pts[i].y;
		s += m_pts[i].z * m_pts[i].z;
	}
	m_normX = sqrt(s);
	//test
	logInfo("m_normX:%lf\n", m_normX);
}
void PoseEst::computeResidual() {
	computefX(Rs, Ts, m_fX.data);
	size_t nM = m_meas.size();
	Meas2D* meas0 = &m_meas[0];
	//get reprojection error
	for (size_t i = 0; i < nM; i++) {
		m_resid.data[2 * i] = m_fX.data[2 * i] - meas0[i].x;
		m_resid.data[2 * i + 1] = m_fX.data[2 * i + 1] - meas0[i].y;
	}
	//get regularization error
	for (size_t c = 0; c < nC; c++) {
		m_resid.data[2 * nM + 3 * c] = m_fX.data[2 * nM + 3 * c]
				- m_predCamCenter[c].x;
		m_resid.data[2 * nM + 3 * c + 1] = m_fX.data[2 * nM + 3 * c + 1]
				- m_predCamCenter[c].x;
		m_resid.data[2 * nM + 3 * c + 2] = m_fX.data[2 * nM + 3 * c + 2]
				- m_predCamCenter[c].x;
	}
}
void PoseEst::computeResidual(const vector<Mat_d>& vecRs, vector<Mat_d>& vecTs,
		Mat_d& res) {
	assert(res.rows == (int)M);
	computefX(vecRs, vecTs, res.data);
	size_t nM = m_meas.size();
	Meas2D* meas0 = &m_meas[0];
	//get reprojection error
	for (size_t i = 0; i < nM; i++) {
		res.data[2 * i] = res.data[2 * i] - meas0[i].x;
		res.data[2 * i + 1] = res.data[2 * i + 1] - meas0[i].y;
	}
	//get regularization error
	for (size_t c = 0; c < nC; c++) {
		res.data[2 * nM + 3 * c] = res.data[2 * nM + 3 * c]
				- m_predCamCenter[c].x;
		res.data[2 * nM + 3 * c + 1] = res.data[2 * nM + 3 * c + 1]
				- m_predCamCenter[c].x;
		res.data[2 * nM + 3 * c + 2] = res.data[2 * nM + 3 * c + 2]
				- m_predCamCenter[c].x;
	}
}
void PoseEst::computefX(const vector<Mat_d>& vecRs, vector<Mat_d>& vecTs,
		double* fX) {
	size_t nM = m_meas.size();
	size_t* ptInd0 = &m_ptInd[0];
	Point3d* pts0 = &m_pts[0];
	Meas2D* meas0 = &m_meas[0];
	Mat_d* Ks0 = &Ks[0];
	const Mat_d* Rs0 = &vecRs[0];
	const Mat_d* Ts0 = &vecTs[0];
	//compute projections
	for (size_t i = 0; i < nM; i++) {
		size_t ptId = ptInd0[i];
		int viewId = meas0[i].viewId;
		_perspectiveProj(Ks0[viewId].data, Rs0[viewId].data, Ts0[viewId].data,
				pts0[ptId].M, fX + 2 * i);
	}
	//compute 3D distances between the camera centers and their predicted ones
	for (size_t c = 0; c < nC; c++) {
		getCameraCenter(Rs0[c].data, Ts0[c].data, fX + 2 * nM + 3 * c);
	}
}
void PoseEst::_computeJci(const vector<Mat_d>& vecRs, vector<Mat_d>& vecTs,
		MyMat<Mat_d>& Jc, int c, int k, double abs_eps) {
	assert(k >=0 && k < 6);
	assert(c >=0 && c < (int) nC);
	size_t nM = m_meas.size();
	//compute projections
	Mat_d ms(nM, 2);

	size_t* ptInd0 = &m_ptInd[0];
	Point3d* pts0 = &m_pts[0];
	Meas2D* meas0 = &m_meas[0];
	Mat_d* Ks0 = &Ks[0];
	const Mat_d* Rs0 = &vecRs[0];
	const Mat_d* Ts0 = &vecTs[0];

	for (size_t i = 0; i < nM; i++) {
		int viewId = meas0[i].viewId;
		if (viewId == c) {
			int ptId = ptInd0[i];
			double* K = Ks0[viewId].data;
			double* R = Rs0[viewId].data;
			double* t = Ts0[viewId].data;
			double* M = pts0[ptId].M;
			_perspectiveProj(K, R, t, M, ms.data + 2 * i);
		}
	}
	for (size_t i = 0; i < nM; i++) {
		if (meas0[i].viewId == c) {
			Jc.data[i].data[k] = (ms.data[2 * i] - m_fX.data[2 * i]) / abs_eps;
			Jc.data[i].data[k + 6] = (ms.data[2 * i + 1] - m_fX.data[2 * i + 1])
					/ abs_eps;
		}

	}
	//compute 3D distances between the camera centers and their predicted ones
	double org[3];
	getCameraCenter(vecRs[c].data, vecTs[c].data, org);
	Jc.data[nM + c].data[k] = (org[0] - m_fX.data[2 * nM + 3 * c]) / abs_eps;
	Jc.data[nM + c].data[k + 6] = (org[1] - m_fX.data[2 * nM + 3 * c + 1])
			/ abs_eps;
	Jc.data[nM + c].data[k + 12] = (org[2] - m_fX.data[2 * nM + 3 * c + 2])
			/ abs_eps;
}

class _JcData {
public:
	PoseEst* est;
	int c;
	int i;
};
void* _parallelComputeJci(void* param) {
	_JcData* pJcData = (_JcData*) param;
	PoseEst* pEst = pJcData->est;
	int c = pJcData->c;
	int i = pJcData->i;

	const double reps = 1e-4; //for disturbing the rotations
	const double teps = 1e-4; //for disturbing the translations
	double abs_teps;

	vector<Mat_d> R1s;
	vector<Mat_d> T1s;
	for (size_t k = 0; k < pEst->nC; k++) {
		R1s.push_back(pEst->Rs[k]);
		T1s.push_back(pEst->Ts[k]);
	}

	double R1[9], T1[3];
	_disturbSO3(pEst->Rs[c], i == 0 ? reps : 0, i == 1 ? reps : 0,
			i == 2 ? reps : 0, R1);
	memcpy(R1s[c].data, R1, sizeof(double) * 9);
	pEst->_computeJci(R1s, T1s, pEst->m_Jc, c, i, reps);
	//restore the rotation
	memcpy(R1s[c], pEst->Rs[c], sizeof(double) * 9);

	if (i == 0)
		_disturbVecX(pEst->Ts[c], teps, T1, abs_teps);
	else if (i == 1)
		_disturbVecY(pEst->Ts[c], teps, T1, abs_teps);
	else
		_disturbVecZ(pEst->Ts[c], teps, T1, abs_teps);
	memcpy(T1s[c].data, T1, sizeof(double) * 3);

	pEst->_computeJci(R1s, T1s, pEst->m_Jc, c, 3 + i, abs_teps);

	return 0;
}

void PoseEst::parallelComputeJc() {
//Notice!:be sure of the function 'computeResidual' being called first
	_JcData jcData[3];
	for (size_t c = 0; c < nC; c++) {
		for (size_t i = 0; i < 3; i++) {
			jcData[i].est = this;
			jcData[i].c = c;
			jcData[i].i = i;
		}
		pthread_t threads[3];
		for (int i = 1; i < 3; i++) {
			pthread_create(&threads[i], 0, _parallelComputeJci,
					(void*) &jcData[i]);
		}
		_parallelComputeJci(&jcData[0]);

		for (int i = 1; i < 3; i++) {
			pthread_join(threads[i], 0);
		}
	}
}
void PoseEst::computeJc() {
//Notice!:be sure of the function 'computeResidual' being called first
	vector<Mat_d> fX1s;
//disturbed rotations and translations
	vector<Mat_d> R1s;
	vector<Mat_d> T1s;

	const double reps = 1e-4; //for disturbing the rotations
	const double teps = 1e-4; //for disturbing the translations

	R1s.clear();
	T1s.clear();
	fX1s.resize(nC, Mat_d());

	for (size_t c = 0; c < nC; c++) {
		R1s.push_back(Rs[c]);
		T1s.push_back(Ts[c]);
		fX1s[c].resize(6, M);
	}
	vector<double> abs_teps;
	abs_teps.resize(3 * nC);
	//compute m_fX1s
	for (size_t c = 0; c < nC; c++) {
		double R1[9], T1[3];
		_disturbSO3(Rs[c], reps, 0, 0, R1);
		memcpy(R1s[c].data, R1, sizeof(double) * 9);
		_computeJci(R1s, T1s, m_Jc, c, 0, reps);

		_disturbSO3(Rs[c], 0, reps, 0, R1);
		memcpy(R1s[c].data, R1, sizeof(double) * 9);
		_computeJci(R1s, T1s, m_Jc, c, 1, reps);

		_disturbSO3(Rs[c], 0, 0, reps, R1);
		memcpy(R1s[c].data, R1, sizeof(double) * 9);
		_computeJci(R1s, T1s, m_Jc, c, 2, reps);

		//restore the rotation
		memcpy(R1s[c], Rs[c], sizeof(double) * 9);

		_disturbVecX(Ts[c], teps, T1, abs_teps[3 * c]);
		memcpy(T1s[c].data, T1, sizeof(double) * 3);
		_computeJci(R1s, T1s, m_Jc, c, 3, abs_teps[3 * c]);

		_disturbVecY(Ts[c], teps, T1, abs_teps[3 * c + 1]);
		memcpy(T1s[c].data, T1, sizeof(double) * 3);
		_computeJci(R1s, T1s, m_Jc, c, 4, abs_teps[3 * c + 1]);

		_disturbVecZ(Ts[c], teps, T1, abs_teps[3 * c + 2]);
		memcpy(T1s[c].data, T1, sizeof(double) * 3);
		_computeJci(R1s, T1s, m_Jc, c, 5, abs_teps[3 * c + 2]);

		//restore the translation
		memcpy(T1s[c], Ts[c], sizeof(double) * 3);
	}
}
#ifdef DBUG
void PoseEst::_testComputeJd() {
//	computeJd();
	size_t nM = nMs + nMd;
	FILE* fp = fopen("/home/tsou/Jd.txt", "w");
	for (size_t i = nMs; i < nM; i++) {
		fprintf(fp, "%g %g %g %g %g %g\n", m_Jd.data[i].data[0], m_Jd.data[i].data[1], m_Jd.data[i].data[2], m_Jd.data[i].data[3], m_Jd.data[i].data[4], m_Jd.data[i].data[5]);
		print(m_Jd.data[i]);
	}
	fclose(fp);
}
#endif
void PoseEst::computeJd() {
	size_t nT = nS + nD;
	const double peps = 1e-4;
	double abs_peps = 0;
	Meas2D* meas0 = &m_meas[0];
	Mat_d* Ks0 = &Ks[0];
	Mat_d* Rs0 = &Rs[0];
	Mat_d* Ts0 = &Ts[0];
	Mat_d* Jd0 = m_Jd.data;
	pair<size_t, size_t>* measIndSz0 = &m_measIndSz[0];
	Point3d* pts0 = &m_pts[0];
	for (size_t i = nS; i < nT; i++) {
		Point3d pt1;
		size_t measId0 = measIndSz0[i].first;
		size_t measNum = measIndSz0[i].second;

		//compute J^j_x
		_disturbVecX(pts0[i].M, peps, pt1.M, abs_peps);
		for (size_t j = 0; j < measNum; j++) {
			size_t measId = measId0 + j;
			size_t viewId = meas0[measId].viewId;
			double m[2];
			_perspectiveProj(Ks0[viewId], Rs0[viewId], Ts0[viewId], pt1.M, m);
			double ux = (m[0] - meas0[measId].x) / abs_peps;
			double vx = (m[1] - meas0[measId].y) / abs_peps;
			Jd0[measId].data[0] = ux;
			Jd0[measId].data[3] = vx;
		}

		//compute J^j_y
		_disturbVecY(pts0[i].M, peps, pt1.M, abs_peps);
		for (size_t j = 0; j < measNum; j++) {
			size_t measId = measId0 + j;
			size_t viewId = meas0[measId].viewId;
			double m[2];
			_perspectiveProj(Ks0[viewId], Rs0[viewId], Ts0[viewId], pt1.M, m);
			double uy = (m[0] - meas0[measId].x) / abs_peps;
			double vy = (m[1] - meas0[measId].y) / abs_peps;
			Jd0[measId].data[1] = uy;
			Jd0[measId].data[4] = vy;
		}

		//compute J^j_z
		_disturbVecY(pts0[i].M, peps, pt1.M, abs_peps);
		for (size_t j = 0; j < measNum; j++) {
			size_t measId = measId0 + j;
			size_t viewId = meas0[measId].viewId;
			double m[2];
			_perspectiveProj(Ks0[viewId], Rs0[viewId], Ts0[viewId], pt1.M, m);
			double uz = (m[0] - meas0[measId].x) / abs_peps;
			double vz = (m[1] - meas0[measId].y) / abs_peps;
			Jd0[measId].data[2] = uz;
			Jd0[measId].data[5] = vz;
		}
	}
}
void PoseEst::computeJacobian() {
	computeJc();
	//parallelComputeJc();
	computeJd();
}
void PoseEst::addWeights() {
	size_t nM = nMs + nMd;
	double* sqrtw0 = &m_sqrtw[0];
	//update Jc
	Mat_d newJ(nM, 12);
	for (size_t c = 0; c < nC; c++) {
		for (size_t i = 0; i < nM; i++) {
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 6, 2,
					sqrtw0[i], m_measCholInvCov[i].data, 2, m_Jc[i].data, 6, 0,
					newJ.data + 12 * i, 6);
			memcpy(m_Jc[i].data, newJ.data + 12 * i, sizeof(double) * 12);
		}
	}
	for (size_t c = 0; c < nC; c++) {
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 6, 2,
				m_wSmooth, m_cholInvPredCov[c].data, 2, m_Jc[nM + c].data, 6, 0,
				newJ.data + 12 * c, 6);
		memcpy(m_Jc[nM + c].data, newJ.data + 12 * c, sizeof(double) * 12);
	}

	//update Jd
	for (size_t i = nMs; i < nM; i++) {
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 3, 2,
				sqrtw0[i], m_measCholInvCov[i].data, 2, m_Jd[i].data, 3, 0,
				newJ.data + 6 * i, 3);
		memcpy(m_Jd[i].data, newJ.data + 6 * i, sizeof(double) * 6);
	}

	//update the residual
	for (size_t i = 0; i < nM; i++) {
		_mat22Ax(sqrtw0[i], m_measCholInvCov[i].data, m_resid.data + 2 * i);
	}
	for (size_t c = 0; c < nC; c++) {
		_mat33Ax(m_wSmooth, m_cholInvPredCov[c].data,
				m_resid.data + 2 * nM + 3 * c);
	}
}

#ifdef OLD_CODE
void PoseEst::addWeightsToNormEquationOld() {
	size_t nM = nMs + nMd;

//update Jc
	for (size_t c = 0; c < nC; c++) {
		for (size_t i = 0; i < nM; i++) {
			_mat22Ax(m_sqrtw[i], m_measCholInvCov[i].data, m_Jc_old.data + c * M + 2 * i);
			_mat22Ax(m_sqrtw[i], m_measCholInvCov[i].data, m_Jc_old.data + (c + 1) * M + 2 * i);
			_mat22Ax(m_sqrtw[i], m_measCholInvCov[i].data, m_Jc_old.data + (c + 2) * M + 2 * i);
			_mat22Ax(m_sqrtw[i], m_measCholInvCov[i].data, m_Jc_old.data + (c + 3) * M + 2 * i);
			_mat22Ax(m_sqrtw[i], m_measCholInvCov[i].data, m_Jc_old.data + (c + 4) * M + 2 * i);
			_mat22Ax(m_sqrtw[i], m_measCholInvCov[i].data, m_Jc_old.data + (c + 5) * M + 2 * i);
		}
	}
	for (int j = 0; j < m_Jc_old.rows; j++) {
		for (size_t c = 0; c < nC; c++) {
			_mat33Ax(m_wSmooth, m_cholInvPredCov[c].data, m_Jc_old.data + j * M + 2 * nM + 3 * c);
		}
	}

	//update Jd
	Mat_d newJ(nMd, 6);
	for (size_t i = 0; i < nMd; i++) {
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 3, 2, m_sqrtw[i], m_measCholInvCov[i].data, 2, m_Jd_old.data + 6 * i, 3, 0, newJ.data + 6 * i, 3);
		memcpy(m_Jd_old.data + i * 6, newJ.data + 6 * i, sizeof(double) * 6);
	}

	//update the residual
	for (size_t i = 0; i < nM; i++) {
		_mat22Ax(m_sqrtw[i], m_measCholInvCov[i].data, m_resid.data + 2 * i);
	}
	for (size_t c = 0; c < nC; c++) {
		_mat33Ax(m_wSmooth, m_cholInvPredCov[c].data, m_resid.data + 2 * nM + 3 * c);
	}
}

void _jacobiProd(int mode, int m, int n, double x[], double y[], void *UsrWrk) {
	assert(mode == 1 || mode == 2);
	PoseEst* poseEst = (PoseEst*) UsrWrk;
	Mat_d* pJc = &poseEst->m_Jc_old;
	Mat_d* pJd = &poseEst->m_Jd_old;
	const vector<size_t>& ptInd = poseEst->m_ptInd;
	const vector<pair<size_t, size_t> >& measIndSz = poseEst->m_measIndSz;
	size_t nMs = poseEst->nMs;
	size_t nS = poseEst->nS;
	size_t nD = poseEst->nD;
	if (mode == 1) {
		//compute y = J*x + y
		Mat_d v(m, 1);
		int nc = pJc->rows;
		cblas_dgemv(CblasColMajor, CblasNoTrans, m, nc, 1.0, pJc->data, m, x, 1, 0, v, 1);

		//#ifdef USE_OPENMP
		//#pragma omp parallel for num_threads(NUM_CUPCORE)
		//#endif
		for (int i = 0; i < pJd->rows; i++) {
			size_t j = ptInd[nMs + i];
			cblas_dgemv(CblasRowMajor, CblasNoTrans, 2, 3, 1.0, pJd->data + 6 * i, 3, x + nc + 3 * (j - nS), 1, 1.0, y + 2 * (nMs + i), 1);
		}
		cblas_daxpy(m, 1.0, v.data, 1, y, 1);
	} else {
		int nc = pJc->rows;
		cblas_dgemv(CblasRowMajor, CblasNoTrans, nc, m, 1.0, pJc->data, m, y, 1, 1.0, x, 1);
		for (size_t i = 0; i < nD; i++) {
			size_t j0 = measIndSz[nS + i].first;
			size_t jn = measIndSz[nS + i].second;
			cblas_dgemv(CblasRowMajor, CblasTrans, 2 * jn, 3, 1.0, pJd->data + 6 * (j0 - nMs), 3, y + 2 * j0, 1.0, 1.0, x + nc + 3 * i, 1);
		}
	}
}

/* call LSQR to solve J\DeltaX = -\epsilon*/
void PoseEst::solveNormEquationCG() {
//	Mat_d v(M, 1);
//	v.fill(1);
//	_jacobiProd(1, M, N, m_resid.data, v.data, this);
//	print(v);
//	Mat_d u(N, 1);
//	u.fill(0);
//	_jacobiProd(2, M, N, u.data, v.data, this);
//	print(u);
//	return;

	Mat_d b;
	b.cloneFrom(m_resid);

	Mat_d v(N, 1);
	Mat_d w(N, 1);

	int istop_out;
	int itn_out;
	double anorm_out;
	double acond_out;
	double rnorm_out;
	double arnorm_out;
	double xnorm_out;

	lsqr(M, N, _jacobiProd, m_damp, this, b.data, v.data, w.data, m_deltaX.data, 0, 1e-30, 1e-30, 1e+6, N / 2, 0, &istop_out, &itn_out, &anorm_out, &acond_out, &rnorm_out, &arnorm_out, &xnorm_out);

	print(m_Jc_old);
	print(m_Jd_old);
	print(m_resid);
	print(m_deltaX);
}
#endif
void PoseEst::buildNormEquation() {
//compute J'*J,  J'\epsilon
//J'*J = [U W;W' V];
//compute U
	Meas2D* meas0 = &m_meas[0];
	size_t nM = nMs + nMd;
	for (size_t c = 0; c < nC; c++) {
		m_U[c].fill(0);
		for (size_t i = 0; i < nM; i++) {
			if (size_t(meas0[i].viewId) == c) {
				const double* Jic = m_Jc.data[i].data;
				cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 6, 6, 2,
						1.0, Jic, 6, Jic, 6, 1.0, m_U[c].data, 6);
			}
		}
		const double* Jic = m_Jc.data[nM + c].data;
		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 6, 6, 3, 1.0, Jic,
				6, Jic, 6, 1.0, m_U[c].data, 6);
#ifdef DEBUG
		//test
		logInfo("U,%d:\n", c);
		print(m_U[c]);
#endif
	}
//compute V
	pair<size_t, size_t>* measIndSz0 = &m_measIndSz[0];
	for (size_t i = 0; i < nD; i++) {
		m_V.data[i].fill(0);
		size_t j0 = measIndSz0[i + nS].first;
		size_t nj = measIndSz0[i + nS].second;
		for (size_t k = 0; k < nj; k++) {
			size_t j = j0 + k;
			const double* Jid = m_Jd.data[j].data;
			cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 3, 3, 2, 1.0,
					Jid, 3, Jid, 3, 1.0, m_V.data[i].data, 3);
		}
	}
//compute W
	for (size_t i = 0; i < nC; i++) {
		for (size_t j = 0; j < nD; j++) {
			m_W.data[i * nD + j].fill(0);
		}
	}
	for (size_t i = 0; i < nD; i++) {
		size_t j0 = measIndSz0[i + nS].first;
		size_t nj = measIndSz0[i + nS].second;
		for (size_t k = 0; k < nj; k++) {
			size_t j = j0 + k;
			const double* Jid = m_Jd.data[j].data;
			int c = meas0[j].viewId;
			const double* Jic = m_Jc.data[j].data;
			cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 6, 3, 2, 1.0,
					Jic, 6, Jid, 3, 0, m_W.data[c * nD + i].data, 3);
		}
	}

#ifdef DEBUG
	//test
	for (size_t i = 0; i < nC; i++) {
		for (size_t j = 0; j < nD; j++) {
			//test
			logInfo("W:(%d,%d)\n", i, j);
			print(m_W[i * nD + j]);
		}
	}
#endif
}
void PoseEst::computeSchurSolveDeltaX(Mat_d& deltaX) {
	assert(deltaX.rows == m_deltaX.rows);
	//augument U and V by multiplying their diagonal elements by 1+m_damp
	for (size_t c = 0; c < nC; c++) {
		m_U.data[c].data[0] += m_damp * m_U.data[c].data[0];
		m_U.data[c].data[7] += m_damp * m_U.data[c].data[7];
		m_U.data[c].data[14] += m_damp * m_U.data[c].data[14];
		m_U.data[c].data[21] += m_damp * m_U.data[c].data[21];
		m_U.data[c].data[28] += m_damp * m_U.data[c].data[28];
		m_U.data[c].data[35] += m_damp * m_U.data[c].data[35];
#ifdef DEBUG
		//test
		logInfo("U*(%d):\n", c);
		print(m_U[c]);
#endif
	}

	for (size_t i = 0; i < nD; i++) {
#ifdef DEBUG
		logInfo("V*(%d)old:\n",i);
		print(m_V[i]);
#endif 
		m_V.data[i].data[0] += m_damp * m_V.data[i].data[0];
		m_V.data[i].data[4] += m_damp * m_V.data[i].data[4];
		m_V.data[i].data[8] += m_damp * m_V.data[i].data[8];
#ifdef DEBUG
		//test
		logInfo("V*(%d):\n", i);
		print(m_V.data[i]);
#endif
		_mat33Inv(m_V.data[i].data, m_iV[i].data);

	}

	//compute Y = W * inv(V)
	for (size_t c = 0; c < nC; c++) {
		for (size_t i = 0; i < nD; i++) {
			const double* W = m_W[c * nD + i].data;
			const double* iV = m_iV[i].data;
			double* Y = m_Y[c * nD + i].data;
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6, 3, 3, 1.0,
					W, 3, iV, 3, 0, Y, 3);

#ifdef DEBUG
			//test
			logInfo("Y:(%d,%d)\n", c, i);
			print(m_Y[c * nD + i]);
#endif
		}
	}
	MyMat<Mat_d> YW(nC, nC);
	//compute U-YW^T
	for (size_t c1 = 0; c1 < nC; c1++) {
		for (size_t c2 = 0; c2 < nC; c2++) {
			YW.data[c1 * nC + c2].resize(6, 6);
			YW.data[c1 * nC + c2].fill(0);
			double* pYW = YW.data[c1 * nC + c2].data;
			for (size_t i = 0; i < nD; i++) {
				const double* Y = m_Y[c1 * nD + i].data;
				const double* Wt = m_W[c2 * nD + i].data;
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 6, 6, 3,
						-1.0, Y, 3, Wt, 3, 1.0, pYW, 6);
			}
			if (c1 == c2) {
				cblas_daxpy(36, 1.0, m_U[c1].data, 1, pYW, 1);
			}
		}
	}
	//compute S
	for (size_t c1 = 0; c1 < nC; c1++) {
		for (size_t c2 = 0; c2 < nC; c2++) {
			_subMatCopy(m_S, 6 * c1, 6 * c2, 6, 6, YW[c1 * nC + c2].data);
		}
	}
#ifdef DEBUG
//	//test
	logInfo("S:\n");
	print(m_S);

	//test
	logInfo("Jr:\n");
	print(m_JRes);

	logInfo("r:\n");
	print(m_resid);
#endif

	//compute bc = \epsilon_c - Y'*\epsilon_d
	Mat_d bc(6 * nC, 1);
	for (size_t c = 0; c < nC; c++) {
		double* pbc = bc.data + 6 * c;
		memcpy(pbc, m_JRes.data + 6 * c, sizeof(double) * 6);
		for (size_t i = 0; i < nD; i++) {
			cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 3, -1.0,
					m_Y[c * nD + i].data, 3, m_JRes.data + 6 * nC + 3 * i, 1,
					1.0, pbc, 1);
		}
	}
#ifdef DEBUG
	//test
	logInfo("bc:\n");
	print(bc);
#endif

	//solve \delta c
	int n = 6 * nC;
	int nhs = 1;
	int info;

	char ch_flag[1] = {'L'};
	dposv_(ch_flag, &n, &nhs, m_S.data, &n, bc.data, &n, &info);

	if (info != 0) {
		//this could happen when there are too many dynamic points and too few static points
		if (verbose) {
			logInfo("damp:%g\n", m_damp);
			warn(
					"PoseEst::_computeSchur - cannot solve the symmetric linear equation! (info:%d)",
					info);
		}
		repErr(
				"PoseEst::_computeSchur - cannot solve the symmetric linear equation! (info:%d)",
				info);
	}

	//bc store the result
	memcpy(deltaX.data, bc.data, sizeof(double) * n);
	//back-substitution
	for (size_t i = 0; i < nD; i++) {
		double* pJRes = m_JRes.data + 6 * nC + 3 * i;
		for (size_t c = 0; c < nC; c++) {
			double* pW = m_W.data[c * nD + i].data;
			cblas_dgemv(CblasRowMajor, CblasTrans, 6, 3, -1.0, pW, 3,
					deltaX.data + 6 * c, 1, 1.0, pJRes, 1);
		}

		double* pInvV = m_iV[i].data;
		double* pxD = deltaX.data + 6 * nC + 3 * i;
		cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, pInvV, 3, pJRes, 1,
				0, pxD, 1);
	}

#ifdef DEBUG
//	//test
	logInfo("dx:\n");
	print(m_deltaX);
#endif
}
void PoseEst::updateParameters(Mat_d& deltaX) {
	//backup the old parameters
	for (size_t c = 0; c < nC; c++) {
		memcpy(tmpRs[c].data, Rs[c].data, sizeof(double) * 9);
		memcpy(tmpTs[c].data, Ts[c].data, sizeof(double) * 3);
	}
	std::copy(m_pts.begin(), m_pts.end(), tmpPts.begin());
//update the camera poses
	double R[9];
	for (size_t c = 0; c < nC; c++) {
		_disturbSO3(Rs[c].data, -deltaX.data[6 * c], -deltaX.data[6 * c + 1],
				-deltaX.data[6 * c + 2], R);
		memcpy(Rs[c].data, R, sizeof(double) * 9);
		Ts[c].data[0] -= deltaX.data[6 * c + 3];
		Ts[c].data[1] -= deltaX.data[6 * c + 4];
		Ts[c].data[2] -= deltaX.data[6 * c + 5];
	}
//update the positions of the dynamic points
	Point3d* pts0 = &m_pts[0];
	for (size_t i = 0; i < nD; i++) {
		pts0[i + nS].x -= deltaX.data[6 * nC + 3 * i];
		pts0[i + nS].y -= deltaX.data[6 * nC + 3 * i + 1];
		pts0[i + nS].z -= deltaX.data[6 * nC + 3 * i + 2];
	}
}
void PoseEst::restoreOldParameters() {
	for (size_t c = 0; c < nC; c++) {
		memcpy(Rs[c].data, tmpRs[c].data, sizeof(double) * 9);
		memcpy(Ts[c].data, tmpTs[c].data, sizeof(double) * 3);
	}
	std::copy(tmpPts.begin(), tmpPts.end(), m_pts.begin());
}
double PoseEst::getResidualNorm() {
	return _vecNorm(m_resid.rows, m_resid.data);
}
void PoseEst::setPoints(const vector<Point3d>& vecPts,
		const vector<vector<Meas2D> >& vecMeas, vector<bool>& dynFlags) {
	assert(vecPts.size() == dynFlags.size());
	assert(vecPts.size() == vecMeas.size());

	m_staPts.clear();
	m_dynPts.clear();

	m_staMeas.clear();
	m_dynMeas.clear();

	m_staPts.reserve(vecPts.size() * 2);
	m_dynPts.reserve(vecPts.size() * 2);
	m_staMeas.reserve(vecPts.size() * 2);
	m_dynMeas.reserve(vecPts.size() * 2);

	nMs = 0;
	nMd = 0;
	nS = 0;
	nD = 0;

	for (size_t i = 0; i < vecPts.size(); i++) {
		if (dynFlags[i]) {
			m_dynPts.push_back(vecPts[i]);
			m_dynMeas.push_back(vecMeas[i]);
			nMd += vecMeas[i].size();
			nD++;
		} else {
			m_staPts.push_back(vecPts[i]);
			m_staMeas.push_back(vecMeas[i]);
			nMs += vecMeas[i].size();
			nS++;
		}
	}
}
void PoseEst::setCamParams(const vector<Mat_d>& vecKs) {
	Ks.assign(vecKs.begin(), vecKs.end());
}
void PoseEst::setInitCamPoses(const vector<Mat_d>& vecRs,
		const vector<Mat_d>& vecTs) {
	Rs.assign(vecRs.begin(), vecRs.end());
	Ts.assign(vecTs.begin(), vecTs.end());

	nC = Rs.size();
}
void PoseEst::setPredCamCenters(const vector<Point3d>& vecOs,
		const vector<Mat_d>& vecCovOs) {
	m_predCamCenter.assign(vecOs.begin(), vecOs.end());
	m_predCov.assign(vecCovOs.begin(), vecCovOs.end());
	m_bUsedPred = true;
}

int PoseEst::apply() {
	allocTmpVars();
	int i;
	double r0, r_old = DBL_MAX, r_new;
	double dx_norm, jr;
	int stopReason = -1;
	for (i = 0; i < m_maxIter; i++) {
		computeResidual();
		computeJacobian();
		addWeights();

		r_new = getResidualNorm();
		if (i == 0)
			r0 = r_new;

		if (r_new < m_e3) {
			stopReason = 3;
			break;
		}
		if (abs(r_new - r_old) < m_e4) {
			stopReason = 4;
			break;
		}

		if (r_new < r_old) {
			m_damp = m_damp * 0.5 < 1 ? 1 : m_damp * 0.5;
			r_old = r_new;
		} else {
			//use old parameter
			restoreOldParameters();
			m_damp *= 2;
		}

		if (m_damp < 1e-8) {
			stopReason = 5;
			break;
		}
		try {
			buildNormEquation();
			computeJRes();
			jr = _vecNorm(m_JRes.rows, m_JRes.data);
			if (jr <= m_e1) {
				stopReason = 1;
				break;
			}
			computeSchurSolveDeltaX(m_deltaX);
			dx_norm = _vecNorm(m_deltaX.rows, m_deltaX.data);
			if (dx_norm < m_e2 * (m_normX + m_e2)) {
				stopReason = 2;
				break;
			}

		} catch (SL_Exception& e) {
			m_wSmooth *= 2;
			r_new = r_old + m_e4;
			continue;
		}
		//print(m_resid);
		updateParameters(m_deltaX);
	}
	if (verbose) {
		//test
		printf(
				"[%d]res :%lf (avg:%lf) r0:%lf - damp:%lf - wSmooth:%lf - jr:%lf - dx_norm:%lf \n",
				i, r_new, r_new * 2 / (nMs + nMd), r0, m_damp, m_wSmooth, jr,
				dx_norm);
	}

//	tic();
//	for (int i = 0; i < 300; i++) {
//		computeResidual();
//		computeJacobian();
//		addWeights();
//		buildNormEquation();
//		computeJRes();
//		computeSchurSolveDeltaX(m_deltaX);
//		updateParameters(m_deltaX);
//	}
//	toc();
//	logInfo("res:%lg\n", getResidualNorm());
	return stopReason;
}

void PoseEst::applyRobust(int maxIter, double maxErr) {
	apply();
	size_t nM = nMs + nMd;
	for (int k = 1; k < maxIter; k++) {
		int iOutlier = 0;
		for (size_t i = 0; i < nM; i++) {
			double err = sqrt(
					m_resid[2 * i] * m_resid[2 * i]
							+ m_resid[2 * i + 1] * m_resid[2 * i + 1]);
			if (err > maxErr) {
				m_meas[i].w = 0;
				m_meas[i].outlier = 1;
				iOutlier++;
			} else {
				double s = err / maxErr;
				double w = (1 - s * s);
				m_meas[i].w = w * w;
				m_meas[i].outlier = 0;
			}
		}
		if (iOutlier == 0)
			break;
		if (verbose) {
			logInfo("[%d] : number of outlier %d\n", k, iOutlier);
		}
		apply();
	}
}
void getCovMatFrom3DPoints(const vector<Point3d>& pts, Mat_d& cov) {
	double mu[3] = { 0, 0, 0 };
	for (size_t i = 0; i < pts.size(); i++) {
		mu[0] += pts[i].x;
		mu[1] += pts[i].y;
		mu[2] += pts[i].z;
	}

	mu[0] /= pts.size();
	mu[1] /= pts.size();
	mu[2] /= pts.size();

	matZeros(3, 3, cov);

	for (size_t i = 0; i < pts.size(); i++) {
		double dx = pts[i].x - mu[0];
		double dy = pts[i].y - mu[1];
		double dz = pts[i].z - mu[2];

		cov(0, 0) += dx * dx;
		cov(0, 1) += dx * dy;
		cov(0, 2) += dx * dz;

		cov(1, 0) += dy * dx;
		cov(1, 1) += dy * dy;
		cov(1, 2) += dy * dz;

		cov(2, 0) += dz * dx;
		cov(2, 1) += dz * dy;
		cov(2, 2) += dz * dz;
	}
	matScale(cov, 1.0 / pts.size());
}

void PoseEst::computeJRes() {
//compute Jc'*\epsilon
	Meas2D* meas0 = &m_meas[0];
	size_t nM = nMs + nMd;
	for (size_t c = 0; c < nC; c++) {
		double* Jr = m_JRes.data + 6 * c;
		memset(Jr, 0, sizeof(double) * 6);
		for (size_t i = 0; i < nM; i++) {
			size_t viewId = size_t(meas0[i].viewId);
			if (viewId == c) {
				cblas_dgemv(CblasRowMajor, CblasTrans, 2, 6, 1.0,
						m_Jc.data[i].data, 6, m_resid.data + 2 * i, 1, 1.0, Jr,
						1);
			}
		}
		cblas_dgemv(CblasRowMajor, CblasTrans, 3, 6, 1.0,
				m_Jc.data[nM + c].data, 6, m_resid.data + 2 * nM + 3 * c, 1,
				1.0, Jr, 1);
	}

	pair<size_t, size_t>* measIndSz0 = &m_measIndSz[0];
//compute Jd'*\epsilon
	for (size_t i = 0; i < nD; i++) {
		double* Jr = m_JRes.data + 6 * nC + 3 * i;

		size_t j0 = measIndSz0[i + nS].first;
		size_t jn = measIndSz0[i + nS].second;

		memset(Jr, 0, sizeof(double) * 3);

		for (size_t k = 0; k < jn; k++) {
			size_t j = j0 + k;
			cblas_dgemv(CblasRowMajor, CblasTrans, 2, 3, 1.0, m_Jd.data[j].data,
					3, m_resid.data + 2 * j, 1, 1.0, Jr, 1);
		}
	}
}

void getPCAFrom3DPoints(const vector<Point3d>& pts, Mat_d& PCAs) {
	Mat_d cov;
	getCovMatFrom3DPoints(pts, cov);
	double U[9], S[3], VT[9];
	dgesvdFor(3, 3, cov.data, U, S, VT);
	PCAs.resize(3, 3);
	memcpy(PCAs.data, VT, sizeof(double) * 9);
	crossProd(PCAs.data, PCAs.data + 3, PCAs.data + 6);
}

void predictCameraPose(const vector<Mat_d>& Rs, const vector<Mat_d>& Ts,
		Point3d& predPos, Mat_d& predCov, size_t localLen, double s1, double s2,
		double s3) {
	assert(Rs.size() > 2 && Ts.size() > 2);
	assert(Rs.size() == Ts.size());
	if (localLen == 0) {
		//nearest neighbor model
		double M[3];
		getCameraCenter(Rs[0].data, Ts[0].data, M);
		predPos.set(M[0], M[1], M[2]);
		matEyes(3, predCov);
		predCov(0, 0) = s1;
		predCov(1, 1) = s2;
		predCov(2, 2) = s3;
	} else {
		//smooth motion model
		vector<Point3d> pos;
		for (size_t t = 0; t < Rs.size() && t < localLen; t++) {
			Point3d M;
			getCameraCenter(Rs[t], Ts[t], M.M);
			pos.push_back(M);
		}
		//// 1. esitmate the average speed
		double v = 0;
		size_t n = pos.size();
		for (size_t t = 1; t < pos.size(); t++) {
			v += dist3(pos[t].M, pos[t - 1].M);
			n++;
		}
		v /= n;

		//velocity direction
		double vdir[3] = { (pos.front().x - pos.back().x) / n, (pos.front().y
				- pos.back().y) / n, (pos.front().z - pos.back().z) / n };
		predPos.x += vdir[0] * v;
		predPos.y += vdir[1] * v;
		predPos.z += vdir[2] * v;

		//compute variance matrix
		Mat_d PCAs;
		getPCAFrom3DPoints(pos, PCAs);
		double S[9] = { s1 * v, 0, 0, 0, s2 * v, 0, 0, 0, s3 * v };
		double tmp[9];
		matATB(3, 3, 3, 3, PCAs.data, S, tmp);
		predCov.resize(3, 3);
		matAB(3, 3, 3, 3, tmp, PCAs.data, predCov.data);
	}
}

void printPointVec(const vector<Point3d>& pts) {
	for (size_t i = 0; i < pts.size(); i++) {
		printf("%g\t%g\t%g\n", pts[i].x, pts[i].y, pts[i].z);
	}
}
//int main(int argc, char** argv) {
//	vector<Point3D> vecPts;
//	vector<vector<Meas2D> > vecMeas;
//	vector<bool> vecDynFlags;
//
//	vector<Mat_d> Ks;
//	vector<vector<Mat_d> > vecRs;
//	vector<vector<Mat_d> > vecTs;
//
////const char* filePath = "/home/zou/slam_results/11-08-09(13-26)_pose/119_result.txt";
//	const char* filePath = "/home/tsou/slam_results/11-08-09(21-02)_pose/105_result.txt";
//	ifstream file(filePath);
//	if (!file)
//		repErr("cannot open '%s' to read!\n", filePath);
//
//	boost::archive::text_iarchive ti(file);
//
//	ti >> vecPts;
//	ti >> vecMeas;
//	ti >> vecDynFlags;
//	ti >> Ks;
//	ti >> vecRs;
//	ti >> vecTs;
//
//	vector<Mat_d> Rs, Ts;
//
//	size_t numCams = Ks.size();
//	for (size_t c = 0; c < numCams; c++) {
//		Rs.push_back(vecRs[c][5]);
//		Ts.push_back(vecTs[c][5]);
//	}
////get camera pose prediction
//	vector<Point3D> predPos;
//	vector<Mat_d> predCov;
//	for (size_t c = 0; c < numCams; c++) {
//		vector<Mat_d> prevRs;
//		vector<Mat_d> prevTs;
//
//		prevRs.assign(vecRs[c].begin() + 1, vecRs[c].end());
//		prevTs.assign(vecTs[c].begin() + 1, vecTs[c].end());
//
//		Point3D curPos;
//		Mat_d curCov;
//		predictCameraPose(prevRs, prevTs, curPos, curCov, 20, 6.0, 3, 3);
//		predPos.push_back(curPos);
//		predCov.push_back(curCov);
//		//test
//		//print(curCov);
//	}
//
//	PoseEst pose;
////test only static points
//	vector<Point3D> staticPoints;
//	vector<vector<Meas2D> > staticMeas;
//	size_t k = 0;
//	for (size_t i = 0; i < vecPts.size(); i++) {
//		if (!vecDynFlags[i] && k < 200) {
//			staticPoints.push_back(vecPts[i]);
//			staticMeas.push_back(vecMeas[i]);
//			k++;
//		}
//	}
//
//	vector<bool> dynFlags;
//	dynFlags.resize(staticPoints.size(), true);
//	for (size_t i = 0; i < 30; i++)
//		dynFlags[i] = false;
//
////pose.setPoints(vecPts, vecMeas, vecDynFlags);
//	tic();
//	pose.setPoints(staticPoints, staticMeas, dynFlags);
//	pose.setCamParams(Ks);
//	pose.setInitCamPoses(Rs, Ts);
//	pose.setPredCamCenters(predPos, predCov);
//	pose.setSmoothWeight(1);
//	pose.applyRobust(2, 4.0);
//	toc();
////test
//	logInfo("total points:%d\n", staticPoints.size());
//	logInfo("finished!\n");
//
//	return 0;
//}

