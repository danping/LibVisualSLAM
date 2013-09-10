/*
 * SL_Triangulate.h
 *
 *  Created on: 2010-11-12
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_TRIANGULATE_H_
#define SL_TRIANGULATE_H_

class CamProjParam {
public:
	const double* K;
	const double* R;
	const double* t;
	const double* m;

public:
	void set(const double* K_, const double* R_, const double* t_,
			const double* m_) {
		K = K_;
		R = R_;
		t = t_;
		m = m_;
	}
};

/* Project a point onto an image */
inline void project(const double* P, const double* M, double* m) {
	m[1] = P[8] * M[0] + P[9] * M[1] + P[10] * M[2] + P[11];
	m[0] = (P[0] * M[0] + P[1] * M[1] + P[2] * M[2] + P[3]) / m[1];
	m[1] = (P[4] * M[0] + P[5] * M[1] + P[6] * M[2] + P[7]) / m[1];
}
void project(const double* R, const double* t, const double* M, double* m);
void project(const double* K, const double* R, const double* t, const double* M,
		double* m);
void project(const double* K, const double* kc, const double* R,
		const double* t, const double *M, double* m);
void project(const double* K, const double* kc, const double* R,
		const double* t, int npts, const double* Ms, double* ms);
void project(const double* R, const double* t, int npts, const double* Ms,
		double* ms);
void project(const double* K, const double* R, const double* t, int npts,
		const double* Ms, double* ms);
inline void KRTtoP(const double* K, const double* R, const double* t,
		double* P) {
	P[0] = K[0] * R[0] + K[1] * R[3] + K[2] * R[6];
	P[1] = K[0] * R[1] + K[1] * R[4] + K[2] * R[7];
	P[2] = K[0] * R[2] + K[1] * R[5] + K[2] * R[8];
	P[3] = K[0] * t[0] + K[1] * t[1] + K[2] * t[2];
	P[4] = K[4] * R[3] + K[5] * R[6];
	P[5] = K[4] * R[4] + K[5] * R[7];
	P[6] = K[4] * R[5] + K[5] * R[8];
	P[7] = K[4] * t[1] + K[5] * t[2];
	P[8] = R[6];
	P[9] = R[7];
	P[10] = R[8];
	P[11] = t[2];
}

/* K  : intrinsic matrix 
 * kc : distortion (5x1)
 * R,t: camera pose
 * Ms : scene points
 * ms0 : image points
 * @return: error
 */
double reprojError2(const double* K, const double* kc, const double* R,
		const double* t, int npts, const double* Ms, const double* ms0);
double reprojError(const double* K, const double* kc, const double* R,
		const double* t, int npts, const double* Ms, const double* ms0);

double reprojError2(const double* K, const double *R, const double *t, int npts,
		const double* Ms, const double* ms0);
double reprojError(const double* K, const double *R, const double *t, int npts,
		const double* Ms, const double* ms0);

double reprojError2(const double* R, double *t, int npts, const double * Ms,
		const double* ms0);
double reprojError(const double* R, double *t, int npts, const double * Ms,
		const double* ms0);
double reprojErrorEach(const double* K, const double* R, const double * t,
		int nPts, const double* Ms, const double* ms0, double* err);
double reprojErrorSingle(const double * K, const double* kc, const double* R,
		const double* t, const double* M, const double* m);

double reprojErrorSingle(const double* K, const double* R, const double* t,
		const double* M, const double* m);

/*compute reprojection error using covariance matrix*/
double reprojErrCov(const double* K, const double* R, const double* t,
		const double* M, const double cov[9], const double m[2]);

/*
 * find outlier (with large reprojection_error > 3.0 pixel) 
 */
double checkReprojError(const double* K, const double* kc, const double* R,
		const double* t, int npts, const double* Ms, const double* ms0,
		char* flags, double thres = 3.0);

/**
 * compute the Jacobian matrix for perspective projection
 */
void computePerspectiveJacobian(const double* K, const double* R,
		const double* t, const double* X, double* fjac);
void computePerspectiveJacobian(const double* R, const double* t,
		const double* X, double* fjac);

/* triangulation algorithm with least square refinement */
void binTriangulate(const double R0[9], const double t0[3], const double R1[9],
		const double t1[3], const double m0[2], const double m1[2],
		double x[3]);

void binTriangulate(const double K0[9], const double R0[9], const double t0[3],
		const double K1[9], const double R1[9], const double t1[3],
		const double m0[2], const double m1[2], double x[3]);
/* linear triangulation algorithm*/
void binTriangulateLin(const double R0[9], const double t0[3],
		const double R1[9], const double t1[3], const double p[2],
		const double q[2], double x[3]);

void binTriangulateLin(const double K0[9], const double R0[9],
		const double t0[3], const double K1[9], const double R1[9],
		const double t1[3], const double p[2], const double q[2], double x[3]);

void binTriangulatePoints(const double R0[9], const double t0[3],
		const double R1[9], const double t1[3], int npts, const double p[],
		const double q[], double x[]);

void binTriangulatePoints(const double K0[9], const double R0[9],
		const double t0[3], const double K1[9], const double R1[9],
		const double t1[3], int npts, const double p[], const double q[],
		double x[]);

void triangulateMultiView(int nview, const double* Rs, const double* ts,
		const double nms[], double x[]);
void triangulateMultiView(int nview, const double* Rs[], const double* ts[],
		const double* nms[], double x[]);
void triangulateMultiViewLin(int nview, const double* R, const double* t,
		const double nms[], double x[]);

void refineTriangulation(int nview, const double* K, const double* R,
		const double* t, const double pts2d[], double x[], double* err);

void getTriangulateCovMat(int nview, const double* K, const double* R,
		const double* t, const double x[3], double cov[9], //covariance of 3D coordinates
		double sigma); //variance of 2D image points

void getTriangulateCovMat(int nview, const double** Ks, const double** Rs,
		const double** ts, const double x[3], double cov[9], double sigma);

void getBinTriangulateCovMat(const double* K1, const double* R1,
		const double* t1, const double* K2, const double* R2, const double* t2,
		double x[3], double cov[9], //covariance matrix of 3D coordinates
		double sigma); //covariance matrix of the 2D image point

/* get the covariance matrix of the image point */
void getProjectionCovMat(const double* K, const double* R, const double* t,
		const double x[3], const double cov[9], double var[4], double sigma);

bool isAtCameraBack(const double* R, const double* t, const double x[]);

/* Sequantial triangulation. 
 * When a new image point is coming, the 3D coordinates of 3D map point can be updated, together with the covarianc
 * K,R,t : parameters of incoming camera pose
 * m : incoming image point
 * M [in/out] : [old/updated]3D position of map point
 * cov [in/out] : [old/updated] covariance of 3D position estimate of map point
 * sigma : standard pixel variance of feature detection 
 */
void seqTriangulate(const double K[9], const double R[9], const double t[3],
		const double m[2], double M[3], double cov[3], double sigma);

/**
 * get the camera center
 */
void getCameraCenter(const double R[9], const double t[3], double org[3]);
/**
 * get the camera center togther with the three-axes
 */
void getCameraCenterAxes(const double R[9], const double t[3], double org[3],
		double ax[3], double ay[3], double az[3]);
/**
 * get distance between two camera centers
 */
double getCameraDistance(const double R1[9], const double t1[3],
		const double R2[9], const double t2[3]);

/*
 * transform camera coordinates into world coordinates
 */
void CamCoord2WoldCoord(const double* R, const double* t, const double mc[3],
		double mw[3]);

/**
 * get the zoomed intrinsic paramters 
 */
void getZoomedIntrinsic(const double* K0, double* K, int w0, int h0, int w, int h);
#endif /* SL_TRIANGULATE_H_ */
