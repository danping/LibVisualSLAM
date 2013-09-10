/*
 * SL_AbsoluteOrientation.h
 *
 *  Created on: 2010-12-3
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_AO_H_
#define SL_AO_H_

/*
 * absolute orientation algorithm using quaternion
 */
void aoGetMatrix(const int npts, const double r1[], const double r2[], double A[]);
bool aoMaxEig(double N[], double v[]);
bool absOrient(const int npts, const double pts1[], const double pts2[], double R[], double t[]);

/*
 * absolute orientation using only three point pairs
 */
bool absOrient3(const double M1[],
		const double M2[],
		const double M3[],
		const double M1_[],
		const double M2_[],
		const double M3_[],
		double R[],
		double t[]);

/**
 * weighted absolute orientation algorithm using svd
 */
void aoGetMatrixWeighted(const int npts, const double w[], const double r1[], const double r2[], double M[]);
void aoGetOptRotation(const double M[], double R[]);
void absOrientWeighted(const int npts,
		const double w[],
		const double pts1[],
		const double pts2[],
		double R[],
		double t[]);
/**
 * apply reweighted least square to obtain the robust estimation using Tukey estimator 
 */
void absOrientRobustTukey(const int npts,
		const double pts1[],
		const double pts2[],
		double R[],
		double t[],
		int maxIter = 10);
/**
 * apply reweighted least square to obtain the robust estimation using Tukey estimator 
 */
void absOrientRobustL1(const int npts,
		const double pts1[],
		const double pts2[],
		double R[],
		double t[],
		int maxIter = 10);
#endif /* SL_AO_H_ */
