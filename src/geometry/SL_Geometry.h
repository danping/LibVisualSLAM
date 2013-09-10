/*
 * SL_Geometry.h
 *
 *  Created on: Dec 15, 2011
 *      Author: Danping Zou
 */

#ifndef SLGEOMETRY_H_
#define SLGEOMETRY_H_
#include <cmath>
#include <vector>
#include "imgproc/SL_Image.h"
/* compute 3d distance*/
double dist3(const double* v1, const double* v2);
/* compute 2d distance*/
double dist2(const double* v1, const double* v2);
double dist2(const double x1, const double y1, const double x2,
		const double y2);

/*cross production v3 = v1xv2 */
void crossProd(const double* v1, const double* v2, double* v3);
void crossProd(double x1, double y1, double z1, double x2, double y2, double z2,
		double* v3);

void crossProd(double x1, double y1, double z1, double x2, double y2, double z2,
		double& x3, double& y3, double& z3);

void crossMat(const double* v, double* cmat);

/* compute MahalanobisDist */
inline double mahaDist2(const double* v1, const double* v2,
		const double ivar[4]) {
	double dx = v1[0] - v2[0];
	double dy = v1[1] - v2[1];
	return sqrt(
			ivar[0] * dx * dx + ivar[1] * dx * dy + ivar[2] * dx * dy
					+ ivar[3] * dy * dy);
}
inline double mahaDist2(double x1, double y1, double x2, double y2,
		const double ivar[4]) {
	double dx = x1 - x2;
	double dy = y1 - y2;
	return sqrt(
			ivar[0] * dx * dx + ivar[1] * dx * dy + ivar[2] * dx * dy
					+ ivar[3] * dy * dy);
}

/*compute inner production with covariance matrix*/
double covInnerProd(int m, const double* C, const double* v1, const double* v2);
/* C:3x3 matrix v1,v2 : 3x1 vectors */
double covInnerProd3(const double* C, const double* v1, const double* v2);

/**
 * get the absolute angle between points B-A-C
 */
double getAbsRadiansBetween(const double A[3], const double B[3],
		const double C[3]);

/**
 * convert polygon to mask image
 * 1 : in the region spanned by the polygon
 * 0 : out of the polygon
 */
void poly2Mask(const std::vector<double>& pts, ImgG& mask, uchar bk=0, uchar fg=255);
#endif /* SLGEOMETRY_H_ */
