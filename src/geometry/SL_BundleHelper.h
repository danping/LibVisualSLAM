/*
 * SL_Bundle.h
 *
 *  Created on: 2010-11-19
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_BUNDLEHELPER_H_
#define SL_BUNDLEHELPER_H_
#include <cmath>
/*
 * copy from 'eucsbademo.c' in sba-1.6.1 
 */
/* pointers to additional data, used for computed image projections and their jacobians */
class sbaGlobs {
public:
	double *rot0params; /* initial rotation parameters, combined with a local rotation parameterization */
	double *intrcalib; /* the 5 intrinsic calibration parameters in the order [fu, u0, v0, ar, skew],
	 * where ar is the aspect ratio fv/fu.
	 * Used only when calibration is fixed for all cameras;
	 * otherwise, it is null and the intrinsic parameters are
	 * included in the set of motion parameters for each camera
	 */
	int nccalib; /* number of calibration parameters that must be kept constant.
	 * 0: all parameters are free 
	 * 1: skew is fixed to its initial value, all other parameters vary (idx1.e. fu, u0, v0, ar) 
	 * 2: skew and aspect ratio are fixed to their initial values, all other parameters vary (idx1.e. fu, u0, v0)
	 * 3: meaningless
	 * 4: skew, aspect ratio and principal point are fixed to their initial values, only the focal length varies (idx1.e. fu)
	 * 5: all intrinsics are kept fixed to their initial values
	 * >5: meaningless
	 * Used only when calibration varies among cameras
	 */

	int ncdist; /* number of distortion parameters in Bouguet's model that must be kept constant.
	 * 0: all parameters are free 
	 * 1: 6th order radial distortion term (kc[4]) is fixed
	 * 2: 6th order radial distortion and one of the tangential distortion terms (kc[3]) are fixed
	 * 3: 6th order radial distortion and both tangential distortion terms (kc[3], kc[2]) are fixed [idx1.e., only 2nd & 4th order radial dist.]
	 * 4: 4th & 6th order radial distortion terms and both tangential distortion ones are fixed [idx1.e., only 2nd order radial dist.]
	 * 5: all distortion parameters are kept fixed to their initial values
	 * >5: meaningless
	 * Used only when calibration varies among cameras and distortion is to be estimated
	 */
	int cnp, pnp, mnp; /* dimensions */
	double *ptparams; /* needed only when bundle adjusting for camera parameters only */
	double *camparams; /* needed only when bundle adjusting for structure parameters only */
};

void img_projsRTS_jac_x(double *p , struct sba_crsm *idxij , int *rcidxs , int *rcsubs , double *jac , void *adata);
void img_projsRTS_x(double *p , struct sba_crsm *idxij , int *rcidxs , int *rcsubs , double *hx , void *adata);
void img_projsRT_x(double *p , struct sba_crsm *idxij , int *rcidxs , int *rcsubs , double *hx , void *adata);
void img_projsRT_jac_x(double *p , struct sba_crsm *idxij , int *rcidxs , int *rcsubs , double *jac , void *adata);
void img_projsKRTS_x(double *p , struct sba_crsm *idxij , int *rcidxs , int *rcsubs , double *hx , void *adata);
void img_projsKRTS_jac_x(double *p , struct sba_crsm *idxij , int *rcidxs , int *rcsubs , double *jac , void *adata);

#define FULLQUATSZ 4
#define _MK_QUAT_FRM_VEC(q, v){                                     \
  (q)[1]=(v)[0]; (q)[2]=(v)[1]; (q)[3]=(v)[2];                      \
  (q)[0]=sqrt(1.0 - (q)[1]*(q)[1] - (q)[2]*(q)[2]- (q)[3]*(q)[3]);  \
}

inline void quat2vec(double* q , double* v) {
	double mag = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
	double sg = q[0] >= 0.0 ? 1.0 : -1.0;
	mag = sg / mag;
	v[0] = q[1] * mag;
	v[1] = q[2] * mag;
	v[2] = q[3] * mag;
}
inline void quatMultFast(double q1[FULLQUATSZ] , double q2[FULLQUATSZ] , double p[FULLQUATSZ]) {
	double t1, t2, t3, t4, t5, t6, t7, t8, t9;
	//double t10, t11, t12;

	t1 = (q1[0] + q1[1]) * (q2[0] + q2[1]);
	t2 = (q1[3] - q1[2]) * (q2[2] - q2[3]);
	t3 = (q1[1] - q1[0]) * (q2[2] + q2[3]);
	t4 = (q1[2] + q1[3]) * (q2[1] - q2[0]);
	t5 = (q1[1] + q1[3]) * (q2[1] + q2[2]);
	t6 = (q1[1] - q1[3]) * (q2[1] - q2[2]);
	t7 = (q1[0] + q1[2]) * (q2[0] - q2[3]);
	t8 = (q1[0] - q1[2]) * (q2[0] + q2[3]);

	/* following fragment it equivalent to the one above */
	t9 = 0.5 * (t5 - t6 + t7 + t8);
	p[0] = t2 + t9 - t5;
	p[1] = t1 - t9 - t6;
	p[2] = -t3 + t9 - t8;
	p[3] = -t4 + t9 - t7;
}

/**Bundle adjustment by using multiple views
 * npts : number of scene points
 * meas : measurement vector (npts * nv *2)
 * covx : covariance of the measurement vector (npts * nv * 4)
 * vmask : visibility mask (npts * nv)
 * pts : 3D scene points  (npts * 3)
 * nv : number of views
 * Rs : rotation matrices (nv * 9)
 * ts : translation vectors (nv * 3)
 * max_iter : maximum number of iterations 
 */
void baRTS(
		int npts ,
		const double* meas ,
		const double* covx ,
		const char* vmask ,
		double* pts ,
		int nv ,
		int nvcon ,
		double* Rs ,
		double* ts ,
		int max_iter);
/* iK : inverse of intrinsic matrix 
 * k_ud : undistortion parameters (7x1)
 */
void baRTS(
		const double* K ,
		const double* kc ,
		const double* iK ,
		const double* k_ud ,
		int npts ,
		int nmeas ,
		const double* meas ,
		const double* covx ,
		const char* vmask ,
		double* pts ,
		int nv ,
		int nvcon ,
		double *Rs ,
		double* ts ,
		int max_iter);

/* npts:number of corresponding points
 * m0,m1 : image points at respective views
 * R0,t0, R1,t1: camera poses
 * M : scene points
 * R0_,t0_,R1_,t1_,M_ : updated values
 * max_iter : maximum number of iterations 
 */
void baRTS2(
		int npts ,
		const double* m0 ,
		const double* covx0 ,
		const double* m1 ,
		const double* covx1 ,
		const double* R0 ,
		const double* t0 ,
		const double* R1 ,
		const double* t1 ,
		const double* M ,
		double* outR0 ,
		double* outT0 ,
		double* outR1 ,
		double* outT1 ,
		double* outM ,
		int matIter);

/* iK : inverse of intrinsic matrix 
 * k_ud : undistortion parameters (7x1)
 */
void baRTS2(
		const double* iK ,
		const double* k_ud ,
		int npts ,
		const double* m0 ,
		const double* covx0 ,
		const double* m1 ,
		const double* covx1 ,
		const double* R0 ,
		const double* t0 ,
		const double* R1 ,
		const double* t1 ,
		const double* M ,
		double* R0_ ,
		double* t0_ ,
		double* R1_ ,
		double* t1_ ,
		double* M_ ,
		int max_iter);

#endif 
