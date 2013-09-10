/*
 * SL_FundamentalMatrix.h
 *
 *  Created on: 2010-11-12
 *      Author: Danping Zou
 */

#ifndef SL_FMATRIX_H_
#define SL_FMATRIX_H_

/*Given a fundamental matrix , compute the epipolar error of the corresponding points*/
double epipolarError(const double *F, const double* m2, const double* m1);
/* Compute the epipolar error given inverse of intrinsic matrices and the essential matrix
 * input:
 * 	invKr, invKl : inverse of intrinsic matrices
 *  E : essential matrix
 *  rn , ln : normalized points on each camera 
 */
double epipolarError(const double* invKr, const double* invKl, const double* E, const double* rn, const double* ln);
double epipolarErrorAvg(int npts, const double* F, const double *m2, const double *m1 );

#endif /* SL_FMATRIX_H_ */
