/*
 * SL_PoseEst.h
 *
 *  Created on: 2011-8-9
 *      Author: Danping Zou
 */

#ifndef SL_POSEEST_H_
#define SL_POSEEST_H_
#include <vector>
#include "math/SL_Matrix.h"
#include "SL_Point.h"
#include "SL_BundleAdjust.h"
using namespace std;

#define POSEST_STOP 0
#define POSEST_CONTINUE 1

//#define DEBUG
class PoseEstParam {
protected:
	bool verbose;
	bool m_bUseDefMeasCov;
	bool m_bUsedPred;
	double m_sigma[4];

	//stop criteria
	double m_e1, m_e2, m_e3, m_e4;
	int m_maxIter;
public:
	PoseEstParam() :
			verbose(true), m_bUseDefMeasCov(true), m_bUsedPred(true), m_e1(1e-6), m_e2(1e-6), m_e3(1e-6), m_e4(1e-6), m_maxIter(30) {
		m_sigma[0] = m_sigma[3] = 1;
		m_sigma[1] = m_sigma[2] = 0;
	}
	void setDefMeasCov(double s1, double s2, double s3) {
		m_sigma[0] = s1;
		m_sigma[1] = s2;
		m_sigma[3] = s3;
		m_sigma[2] = s2;
		m_bUseDefMeasCov = true;
	}
};
class PoseEst: public PoseEstParam {
protected:
	vector<Point3d> m_staPts;
	vector<Point3d> m_dynPts;
	vector<vector<Meas2D> > m_staMeas;
	vector<vector<Meas2D> > m_dynMeas;
protected:
	vector<Mat_d> tmpRs; //initial camera poses
	vector<Mat_d> tmpTs;
	vector<Point3d> tmpPts;
public:

	vector<Mat_d> Ks; //camera intrinsic parameters
	vector<Mat_d> Rs; //initial camera poses
	vector<Mat_d> Ts;
	vector<Point3d> m_pts;

	vector<Mat_d> m_imgCov; //covariance of reprojections
	double m_wSmooth; //regularization weight
	vector<Mat_d> m_predCov; //covariance of the predictions
	vector<Point3d> m_predCamCenter; //predicted camera centers
protected:
	double m_damp;
	double m_normX;

	size_t nC, nMs, nMd, nS, nD;
	size_t M, N;
	//size:nCx9 -- Cholesky decomposition of prediction covariance
	Mat_d m_covChol;
	//size:(2*nMS+2*nMD+3*nC)x1  -- store the residual

	Mat_d m_fX;
	Mat_d m_resid;
	//size:(6*nC + 3*nD)x1 -- store the incremental of the solution
	Mat_d m_deltaX;
	//size:M[2x6]+C[3x6]
	MyMat<Mat_d> m_Jc;
	//size:Ms[]+ Md[2x3] + Cd[]
	MyMat<Mat_d> m_Jd;

	//size:nCxnD [6x3]
	MyMat<Mat_d> m_W;
	//size:nC [6x6]
	MyMat<Mat_d> m_U;
	//size:nD [3x3]
	MyMat<Mat_d> m_V;
	MyMat<Mat_d> m_iV; //inverse of V

	//size:nCxnD [6x3]
	MyMat<Mat_d> m_Y; //Y = W*inv(V)
	//size:6nCx6nC
	Mat_d m_S; //Schar complement of U
	//size:6nC+3nD
	Mat_d m_JRes; //J'*\epsilon

	//pointer to the first index for the corresponding point
	vector<pair<size_t, size_t> > m_measIndSz;

	//indices to the point and cameras for each measurement
	vector<size_t> m_ptInd;
	vector<Meas2D> m_meas;
	vector<Mat_d> m_measCholInvCov;
	vector<double> m_sqrtw;
	vector<Mat_d> m_cholInvPredCov;

	void allocTmpVars();
	void computeResidual();
	void computeResidual(const vector<Mat_d>& vecRs, vector<Mat_d>& vecTs, Mat_d& res);
	void computefX(const vector<Mat_d>& vecRs, vector<Mat_d>& vecTs, double* fx);

#ifdef OLD_CODE
	void computeJcOld();
	void computeJdOld();
	void addWeightsToNormEquationOld();
	/* call LSQR to solve J\DeltaX = -\epsilon*/
	void solveNormEquationCG();
#endif 
	double getResidualNorm();

#ifdef DEBUG
	void _testComputeJc();
	void _testComputeJd();
	void _testOutputJ();
#endif

	void _computeJci(const vector<Mat_d>& vecRs, vector<Mat_d>& vecTs, MyMat<Mat_d>& Jc, int c, int k, double abs_eps);
	void parallelComputeJc();
	void computeJc();
	void computeJd();
	void computeJacobian();
	void addWeights();
	void buildNormEquation();
	void computeJRes();
	void computeSchurSolveDeltaX(Mat_d& deltaX);
	void updateParameters(Mat_d& deltaX);
	void restoreOldParameters();
public:
	PoseEst();
	virtual ~PoseEst();
	void setPoints(const vector<Point3d>& vecPts, const vector<vector<Meas2D> >& vecMeas, vector<bool>& dynFlags);
	void setCamParams(const vector<Mat_d>& vecKs);
	void setInitCamPoses(const vector<Mat_d>& vecRs, const vector<Mat_d>& vecTs);
	void setPredCamCenters(const vector<Point3d>& vecOs, const vector<Mat_d>& vecCovOs);
	void setSmoothWeight(double w) {
		m_wSmooth = w;
	}
	int apply();
	void applyRobust(int maxIter, double maxEpiErr);
	friend void _jacobiProd(int mode, int m, int n, double x[], double y[], void *UsrWrk);
	friend void* _parallelComputeJci(void* param) ; 
};
void* _parallelComputeJci(void* param) ;
void _jacobiProd(int mode, int m, int n, double x[], double y[], void *UsrWrk);
void getCovMatFrom3DPoints(const vector<Point3d>& pts, Mat_d& cov);
/* each row of PCAs is the corresponding PCA vector*/
void getPCAFrom3DPoints(const vector<Point3d>& pts, Mat_d& PCAs);
void predictCameraPose(const vector<Mat_d>& Rs, const vector<Mat_d>& Ts, Point3d& predPos, Mat_d& cov, size_t localLen,
/*three scales for covariance computing*/
double s1, double s2, double s3);
void printPointVec(const vector<Point3d>& pts);
#endif /* SL_POSEEST_H_ */
