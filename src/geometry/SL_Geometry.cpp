/*
 * SL_Geometry.cpp
 *
 *  Created on: Dec 15, 2011
 *      Author: Danping Zou
 */

#include "SL_Geometry.h"
#include "SL_ConvexHull2D.h"
#include "math/SL_LinAlg.h"

#include <opencv2/opencv.hpp>
#include <numeric>
using namespace std;
double dist3(const double* v1, const double* v2) {
	double dx = v1[0] - v2[0];
	double dy = v1[1] - v2[1];
	double dz = v1[2] - v2[2];
	return sqrt(dx * dx + dy * dy + dz * dz);
}
double dist2(const double* v1, const double* v2) {
	double dx = v1[0] - v2[0];
	double dy = v1[1] - v2[1];
	return sqrt(dx * dx + dy * dy);
}
double dist2(const double x1, const double y1, const double x2,
		const double y2) {
	double dx = x1 - x2;
	double dy = y1 - y2;
	return sqrt(dx * dx + dy * dy);
}

void crossProd(const double* v1, const double* v2, double* v3) {
	v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
	v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
	v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void crossProd(double x1, double y1, double z1, double x2, double y2, double z2,
		double* v3) {
	v3[0] = y1 * z2 - z1 * y2;
	v3[1] = z1 * x2 - x1 * z2;
	v3[2] = x1 * y2 - y1 * x2;
}

void crossProd(double x1, double y1, double z1, double x2, double y2, double z2,
		double& x3, double& y3, double& z3) {
	x3 = y1 * z2 - z1 * y2;
	y3 = z1 * x2 - x1 * z2;
	z3 = x1 * y2 - y1 * x2;
}

void crossMat(const double* v, double* cmat) {
	matFill(3, 3, cmat, 0);
	cmat[1] = -v[2];
	cmat[2] = v[1];
	cmat[5] = -v[0];
	cmat[3] = v[2];
	cmat[6] = -v[1];
	cmat[7] = v[0];
}
double covInnerProd(int m, const double* C, const double* v1,
		const double* v2) {
	double * t = new double[m];
	matAxpy(m, m, 1.0, C, v2, 0, 0, t);
	double d = inner_product(t, t + m, v1, 0.0);
	delete[] t;
	return d;
}
double covInnerProd3(const double* C, const double* v1, const double* v2) {
	double t[3];
	mat33ProdVec(C, v2, t);
	return inner_product(t, t + 3, v1, 0.0);
}

double getAbsRadiansBetween(const double A[3], const double B[3],
		const double C[3]) {
	double AB[3] = { B[0] - A[0], B[1] - A[1], B[2] - A[2] };
	double AC[3] = { C[0] - A[0], C[1] - A[1], C[2] - A[2] };

	double s = innerProd3(AB, AC, 0);
	double n1 = vec3Len(AB);
	double n2 = vec3Len(AC);

	return fabs(acos(s / (n1 * n2)));
}

void poly2Mask(const std::vector<double>& pts, ImgG& mask, uchar bk, uchar fg) {
	assert(!mask.empty());
	mask.fill(bk);
	cv::Mat cvMask(mask.m, mask.n, CV_8UC1, mask.data);
	int npts = pts.size() / 2;
	cv::Point2i* cvPts = new cv::Point2i[npts];
	for (int i = 0; i < npts; ++i) {
		cvPts[i].x = pts[2 * i];
		cvPts[i].y = pts[2 * i + 1];
	}
	cv::fillConvexPoly(cvMask, cvPts, npts,
			cv::Scalar(fg));
	
	delete[] cvPts;
}