/*
 * SL_RigidTransform.h
 *
 *  Created on: 2011-6-9
 *      Author: Danping Zou
 */

#ifndef SL_RIGIDTRANSFORM_H_
#define SL_RIGIDTRANSFORM_H_

/**
 * get incremental transform [dR dt; 0 1]*[R1 t1; 0 1]= [R t; 0 1];
 */
void getRigidTransFromTo(
		const double R1[9],
		const double t1[3],
		const double R[9],
		const double t[3],
		double dR[9],
		double dt[3]);
/**
 * incremental transform [ dR dt; 0 1]* [R1 t1; 0 1] -> [R t; 0 1];
 */
void rigidTransFromTo(
		const double R1[9],
		const double t1[3],
		const double dR[9],
		const double dt[3],
		double R[9],
		double t[3]);
/*
 * inverse transformation of (R,t) -> (R^T, -R^T*t)
 */
void invRigidTransFromTo(
		const double R[9],
		const double t[3],
		double iR[9],
		double it[3]);
/**
 * ptsNew[i] = R* ptsOld[i] + t
 */
void rigidTrans(
		const double* R,
		const double* t,
		int nPts,
		const double* ptsOld,
		double* ptsNew);
/**
 * ptNew = R * ptsOld + t
 */
void rigidTrans(
		const double* R,
		const double* t,
		const double ptOld[3],
		double ptNew[3]);

/**
 * pt2 = R^T*pt1 - R^T*t)
 */
void invRigidTrans(
		const double* R,
		const double* t,
		const double pt1[3],
		double pt2[3]);
/**
 * approximate the rotation matrix by SVD
 */
void approxRotationMat(const double *A, double* R);

/**
 * get the classic relative camera pose between two camera poses [R1 t1] to [R2 t2]:
 *  [I 0] and [dR dt]
 */
void getRelativePose(
		const double R1[9],
		const double t1[3],
		const double R2[9],
		const double t2[3],
		double dR[9],
		double dt[3]);
#endif /* SL_RIGIDTRANSFORM_H_ */
