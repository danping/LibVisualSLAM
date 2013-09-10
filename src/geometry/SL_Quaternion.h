/*
 * SL_Quaternion.h
 *
 *  Created on: 2010-11-18
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_QUATERNION_H_
#define SL_QUATERNION_H_

double quat_len(const double* q);
void quat_norm(const double* q, double* qn);
/*qj = q^-1*/
void quat_conj(const double* q, double* qj);
/*q = q1*q2*/
void quat_mult(const double* q1, const double* q2, double* q);

/*rotate a point v by a quaternion q to produce a new point v1*/
void quat_rot(const double* v, const double* q, double* v1);

/*rotate a point v around the axis rot_v 
 * theta : rotation angle (in degree)
 */
void quat_rot(const double* v, const double* rot_v, double theta, double* v1);
/*convert a rotation matrix into a quaternion*/
void mat2quat(const double* R, double*q);
void quat2mat(const double* q, double* R);
#endif /* SL_QUATERNION_H_ */
