/*
 * SL_Quaternion.cpp
 *
 *  Created on: 2010-11-18
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#include "SL_Quaternion.h"
#include "SL_error.h"
#include "SL_Geometry.h"
#include "math/SL_LinAlg.h"
#include <cmath>
#include <cassert>

#define PI 3.14159265358979323846
double quat_len(const double* q) {
	return sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
}
void quat_norm(const double* q, double* qn) {
	double len = quat_len(q);
	qn[0] = q[0] / len;
	qn[1] = q[1] / len;
	qn[2] = q[2] / len;
	qn[3] = q[3] / len;
}
/*qj = q^-1*/
void quat_conj(const double* q, double* qj) {
	qj[0] = q[0];
	qj[1] = -q[1];
	qj[2] = -q[2];
	qj[3] = -q[3];
}
/*q = q1*q2*/
void quat_mult(const double* q1, const double* q2, double* q) {
	assert(q1 != q && q2 != q);
	q[0] = (q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3]);
	q[1] = (q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2]);
	q[2] = (q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1]);
	q[3] = (q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0]);
}

void quat_rot(const double* v, const double* q, double* v1) {
	double qv[4];
	qv[0] = 0;
	qv[1] = v[0];
	qv[2] = v[1];
	qv[3] = v[2];

	double qv1[4];
	double qj[4];
	quat_conj(q, qj);

	quat_mult(qv, qj, qv1);
	quat_mult(q, qv1, qv);

	v1[0] = qv[1];
	v1[1] = qv[2];
	v1[2] = qv[3];
}

void quat_rot(const double* v, const double* rot_v, double theta, double* v1) {
	//convert to radian
	theta = theta / 180 * PI * 0.5;

	//compute the rotation quaternion
	double qrot[4];
	double srot = sqrt(rot_v[0] * rot_v[0] + rot_v[1] * rot_v[1] + rot_v[2] * rot_v[2]);
	double sc = sin(theta) / srot;

	qrot[0] = cos(theta);
	qrot[1] = sc * rot_v[0];
	qrot[2] = sc * rot_v[1];
	qrot[3] = sc * rot_v[2];

	quat_rot(v, qrot, v1);
}
void mat2quat(const double* R, double*q) {
	double tr = R[0] + R[4] + R[8];
	if( tr > 0){
		double S = sqrt(1+ tr) * 2;
		      q[1] = ( R[7] - R[5] ) / S;
		      q[2] = ( R[2] - R[6] ) / S;
		      q[3] = ( R[3] - R[1] ) / S;
		      q[0] = 0.25 * S;	
	}else{
		if (R[0] > R[4] && R[0] > R[8]) { // Column 0: 
			double S = sqrt(1.0 + R[0] - R[4] - R[8]) * 2;
			q[1] = 0.25 * S;
			q[2] = (R[3] + R[1]) / S;
			q[3] = (R[2] + R[6]) / S;
			q[0] = (R[7] - R[5]) / S;
		} else if (R[4] > R[8] && R[4] > R[0]) { // Column 1: 
			double S = sqrt(1.0 + R[4] - R[0] - R[8]) * 2;
			q[1] = (R[3] + R[1]) / S;
			q[2] = 0.25 * S;
			q[3] = (R[7] + R[5]) / S;
			q[0] = (R[2] - R[6]) / S;
		} else if (R[8] > R[0] && R[8] > R[4]) { // Column 2:
			double S = sqrt(1.0 + R[8] - R[0] - R[4]) * 2;
			q[1] = (R[2] + R[6]) / S;
			q[2] = (R[7] + R[5]) / S;
			q[3] = 0.25 * S;
			q[0] = (R[3] - R[1]) / S;
		}
	}
}

void quat2mat(const double* q, double* R) {
	double Nq = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
	double s;
	if (Nq > 0.0)
		s = 2 / Nq;
	else
		s = 0.0;

	double X = q[1] * s;
	double Y = q[2] * s;
	double Z = q[3] * s;
	double wX = q[0] * X;
	double wY = q[0] * Y;
	double wZ = q[0] * Z;
	double xX = q[1] * X;
	double xY = q[1] * Y;
	double xZ = q[1] * Z;
	double yY = q[2] * Y;
	double yZ = q[2] * Z;
	double zZ = q[3] * Z;

	R[0] = 1.0 - (yY + zZ);
	R[1] = xY - wZ;
	R[2] = xZ + wY;
	R[3] = xY + wZ;
	R[4] = 1.0 - (xX + zZ);
	R[5] = yZ - wX;
	R[6] = xZ - wY;
	R[7] = yZ + wX;
	R[8] = 1.0 - (xX + yY);
}
