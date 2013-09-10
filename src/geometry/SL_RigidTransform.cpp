/*
 * SL_RigidTransform.cpp
 *
 *  Created on: 2011-6-9
 *      Author: Danping Zou
 */

#include "SL_RigidTransform.h"
#include "math/SL_LinAlg.h"
void getRigidTransFromTo(
		const double R1[9],
		const double t1[3],
		const double R2[9],
		const double t2[3],
		double dR[9],
		double dt[3]) {
	//dR = R2* R1'
	mat33ABT(R2, R1, dR);
	//dt = - dR*t1 + t2;
	mat33ProdVec(dR, t1, t2, dt, -1.0, 1.0);
}
void rigidTransFromTo(
		const double R1[9],
		const double t1[3],
		const double dR[9],
		const double dt[3],
		double R2[9],
		double t2[3]) {
	//R2 = dR*R1;
	mat33AB(dR, R1, R2);
	//t2 = dR*t1 + dt
	mat33ProdVec(dR, t1, dt, t2, 1.0, 1.0);
}
//inverse transformation of (R,t) -> (R^T, -R^T*t)
void invRigidTransFromTo(
		const double R[9],
		const double t[3],
		double iR[9],
		double it[3]) {
	mat33Trans(R, iR);
	mat33TransProdVec(R, t, it);
	it[0] = -it[0];
	it[1] = -it[1];
	it[2] = -it[2];
}
void rigidTrans(
		const double* R,
		const double* t,
		int nPts,
		const double* ptsOld,
		double* ptsNew) {
	for (int i = 0; i < nPts; i++) {
		mat33ProdVec(R, ptsOld + 3 * i, t, ptsNew + 3 * i, 1.0, 1.0);
	}
}
void rigidTrans(
		const double* R,
		const double* t,
		const double ptOld[3],
		double ptNew[3]) {
	mat33ProdVec(R, ptOld, t, ptNew, 1.0, 1.0);
}
void invRigidTrans(
		const double* R,
		const double* t,
		const double pt1[3],
		double pt2[3]) {
	double iR[9], it[3];
	invRigidTransFromTo(R, t, iR, it);
	rigidTrans(iR, it, pt1, pt2);
}
void approxRotationMat(const double *A, double* R) {
	double U[9], S[3], VT[9];
	dgesvdFor(3, 3, A, U, S, VT);
	mat33AB(U, VT, R);
}

void getRelativePose(
		const double R1[9],
		const double t1[3],
		const double R2[9],
		const double t2[3],
		double dR[9],
		double dt[3]) {
	//dR = R2*R1'
	mat33ABT(R2, R1, dR);
	//dt = - dR*t1 + t2;
	mat33ProdVec(dR, t1, t2, dt, -1.0, 1.0);
}