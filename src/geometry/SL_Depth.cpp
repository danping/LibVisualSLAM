/*
 * SL_Depth.cpp
 *
 *  Created on: Feb 14, 2012
 *      Author: Danping Zou
 */

#include "SL_Depth.h"
#include "math/SL_LinAlg.h"
#include "geometry/SL_Distortion.h"
#include "geometry/SL_RigidTransform.h"
#include <cmath>
void depth2point3d(double nx, double ny, double d, double M0[3]) {
	double s = sqrt(nx * nx + ny * ny + 1);
	M0[0] = nx * d / s;
	M0[1] = ny * d / s;
	M0[2] = d / s;
}
void get3DCoord(const double* K, const double* R, const double* t, int u, int v,
		double d, double M[3]) {
	double invK[9];
	double M0[3];
	double nx, ny;
	mat33Inv(K, invK);
	normPoint(invK, (double) u, (double) v, nx, ny);
	depth2point3d(nx, ny, (K[0] * K[4]) < 0 ? -d : d, M0);
	invRigidTrans(R, t, M0, M);
}

