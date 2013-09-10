/*
 * SL_ProbFuncs.cpp
 *
 *  Created on: 2011-10-20
 *      Author: Danping Zou
 */

#include "SL_ProbFuncs.h"
#include "SL_LinAlg.h"
#include "geometry/SL_Geometry.h"

static const double PI = 3.1415926535897932384626;
double normpdf2(const double x[], const double miu[], const double sigma[]) {
	double dx[2] = { x[0] - miu[0], x[1] - miu[1] };
	return normpdf2(dx, sigma);
}

double normpdf2(const double x[], const double sigma[]) {
	double inv2PI = 1.0 / (2 * PI);
	double det = mat22Det(sigma);
	det = 1.0 / sqrt(det);
	double invsigma[4];
	mat22Inv(sigma, invsigma);
	double dist = covInnerProd(2, invsigma, x, x);
	return inv2PI * det * exp(-0.5 * dist);
}

double normpdf3(const double x[], const double sigma[]) {
	double inv2PI = 1.0 / (2 * PI);
	double inv2PIsqrt = sqrt(inv2PI);
	double a = inv2PIsqrt * inv2PIsqrt * inv2PIsqrt;
	double det = mat33Det(sigma);
	det = 1.0 / sqrt(det);
	double invsigma[9];
	mat33Inv(sigma, invsigma);
	double dist = covInnerProd(3, invsigma, x, x);
	return a * det * exp(-0.5 * dist);
}
double normpdf_const(int dim, const double invsigma[]) {
	double inv2Pi = 1.0 / sqrt(2 * PI);
	double s = sqrt(matDet(dim, invsigma));
	for (int i = 0; i < dim; i++)
		s *= inv2Pi;
	return s;
}
double normpdf_expval(int dim, const double x[], const double invsigma[]) {
	double dist = covInnerProd(dim, invsigma, x, x);
	return exp(-0.5 * dist);
}
double normpdf(int dim, const double x[], const double sigma[]) {
	double* invSigma = new double[dim * dim];
	matInv(dim,sigma,invSigma);

	double con_val = normpdf_const(dim,invSigma);
	double exp_val = normpdf_expval(dim,x,invSigma);
	return con_val*exp_val;
}