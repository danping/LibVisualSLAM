/*
 * SL_p3p.cpp
 *
 *  Created on: 2010-12-2
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#include "SL_p3p.h"
#include "SL_error.h"
#include "SL_Geometry.h"
#include "SL_AbsoluteOrientation.h"
#include "quartSolver.h"

#include "math/SL_LinAlg.h"
#include <cmath>

void p3p_get_coeffs(const double K[], const double M1[], const double M2[], const double M3[], const double m1[],
		const double m2[], const double m3[], double d[], double cs[]) {
	d[0] = dist3(M1, M2);
	d[1] = dist3(M1, M3);
	d[2] = dist3(M2, M3);

	double iC[9], C[9];
	mat33ABT(K, K, iC);
	mat33Inv(iC, C);

	double _m1[3] = { m1[0], m1[1], 1 };
	double _m2[3] = { m2[0], m2[1], 1 };
	double _m3[3] = { m3[0], m3[1], 1 };

	double us1 = sqrt(innerProd3(_m1, _m1, C));
	double us2 = sqrt(innerProd3(_m2, _m2, C));
	double us3 = sqrt(innerProd3(_m3, _m3, C));

	cs[0] = innerProd3(_m2, _m1, C) / (us1 * us2);
	cs[1] = innerProd3(_m3, _m1, C) / (us1 * us3);
	cs[2] = innerProd3(_m3, _m2, C) / (us2 * us3);
}
void p3p_get_quartic(const double d[], const double cs[], double cfs[]) {
	double d12 = d[0];
	double d13 = d[1];
	double d23 = d[2];

	double cs12 = cs[0];
	double cs13 = cs[1];
	double cs23 = cs[2];

	double cs12_2 = cs12 * cs12;
	double cs12_3 = cs12_2 * cs12;
	double cs12_4 = cs12_3 * cs12;

	double cs13_2 = cs13 * cs13;
	double cs13_3 = cs13_2 * cs13;
	double cs13_4 = cs13_3 * cs13;

	double cs23_2 = cs23 * cs23;
	double cs23_3 = cs23_2 * cs23;
	double cs23_4 = cs23_3 * cs23;

	double d23_2 = d23 * d23;
	double d23_4 = d23_2 * d23_2;
	double d23_6 = d23_2 * d23_4;
	double d23_8 = d23_6 * d23_2;

	double d12_2 = d12 * d12;
	double d12_4 = d12_2 * d12_2;
	double d12_6 = d12_4 * d12_2;
	double d12_8 = d12_6 * d12_2;

	double d13_2 = d13 * d13;
	double d13_4 = d13_2 * d13_2;
	double d13_6 = d13_2 * d13_4;
	double d13_8 = d13_6 * d13_2;

	double cs12_2_cs13_2 = cs12_2 * cs13_2;
	double cs12_2_cs23_2 = cs12_2 * cs23_2;
	double cs13_2_cs23_2 = cs13_2 * cs23_2;

	double cs12_3_cs13_cs23 = cs12_3 * cs13 * cs23;
	double cs12_2_cs13_2_cs23_2 = cs12_2_cs13_2 * cs23_2;
	double cs12_cs13_3_cs23 = cs12 * cs13_3 * cs23;
	double cs12_cs13_cs23_3 = cs12 * cs13 * cs23_3;
	double cs12_cs13_cs23 = cs12 * cs13 * cs23;

	double e8 = 16 * cs12_4 - 64 * cs12_3_cs13_cs23 + 64 * cs12_2_cs13_2_cs23_2 + 32 * cs12_2_cs13_2 + 32
			* cs12_2_cs23_2 - 32 * cs12_2 - 64 * cs12_cs13_3_cs23 - 64 * cs12_cs13_cs23_3 + 64 * cs12_cs13_cs23 + 16
			* cs13_4 + 32 * cs13_2_cs23_2 - 32 * cs13_2 + 16 * cs23_4 - 32 * cs23_2 + 16;
	double e6 = -64 * cs12_4 * cs13_2 * d23_2 - 32 * cs12_4 * d13_2 + 32 * cs12_4 * d23_2 + 128 * cs12_3 * cs13_3
			* cs23 * d23_2 + 32 * cs12_3_cs13_cs23 * d12_2 + 96 * cs12_3_cs13_cs23 * d13_2 - 32 * cs12_3_cs13_cs23
			* d23_2 - 64 * cs12_2 * cs13_4 * d23_2 - 64 * cs12_2_cs13_2_cs23_2 * d12_2 - 64 * cs12_2_cs13_2_cs23_2
			* d13_2 - 128 * cs12_2_cs13_2_cs23_2 * d23_2 - 32 * cs12_2_cs13_2 * d12_2 - 32 * cs12_2_cs13_2 * d13_2
			+ 128 * cs12_2_cs13_2 * d23_2 - 32 * cs12_2_cs23_2 * d12_2 - 64 * cs12_2_cs23_2 * d13_2 + 32
			* cs12_2_cs23_2 * d23_2 + 32 * cs12_2 * d12_2 + 64 * cs12_2 * d13_2 - 64 * cs12_2 * d23_2 + 96
			* cs12_cs13_3_cs23 * d12_2 + 32 * cs12_cs13_3_cs23 * d13_2 - 32 * cs12_cs13_3_cs23 * d23_2 + 96
			* cs12_cs13_cs23_3 * d12_2 + 96 * cs12_cs13_cs23_3 * d13_2 + 32 * cs12_cs13_cs23_3 * d23_2 - 96
			* cs12_cs13_cs23 * d12_2 - 96 * cs12_cs13_cs23 * d13_2 + 32 * cs12_cs13_cs23 * d23_2 - 32 * cs13_4 * d12_2
			+ 32 * cs13_4 * d23_2 - 64 * cs13_2_cs23_2 * d12_2 - 32 * cs13_2_cs23_2 * d13_2 + 32 * cs13_2_cs23_2
			* d23_2 + 64 * cs13_2 * d12_2 + 32 * cs13_2 * d13_2 - 64 * cs13_2 * d23_2 - 32 * cs23_4 * d12_2 - 32
			* cs23_4 * d13_2 + 64 * cs23_2 * d12_2 + 64 * cs23_2 * d13_2 - 32 * cs23_2 * d23_2 - 32 * d12_2 - 32
			* d13_2 + 32 * d23_2;
	double e4 = 16 * cs12_4 * d13_4 - 32 * cs12_4 * d13_2 * d23_2 + 16 * cs12_4 * d23_4 - 32 * cs12_3_cs13_cs23 * d12_2
			* d13_2 - 32 * cs12_3_cs13_cs23 * d12_2 * d23_2 - 32 * cs12_3_cs13_cs23 * d13_4 + 64 * cs12_3_cs13_cs23
			* d13_2 * d23_2 - 32 * cs12_3_cs13_cs23 * d23_4 + 64 * cs12_2_cs13_2_cs23_2 * d12_2 * d13_2 + 64
			* cs12_2_cs13_2_cs23_2 * d12_2 * d23_2 + 64 * cs12_2_cs13_2_cs23_2 * d13_2 * d23_2 + 16 * cs12_2_cs13_2
			* d12_4 - 64 * cs12_2_cs13_2 * d12_2 * d23_2 + 16 * cs12_2_cs13_2 * d13_4 - 64 * cs12_2_cs13_2 * d13_2
			* d23_2 + 48 * cs12_2_cs13_2 * d23_4 + 16 * cs12_2_cs23_2 * d12_4 + 32 * cs12_2_cs23_2 * d12_2 * d13_2 + 48
			* cs12_2_cs23_2 * d13_4 - 64 * cs12_2_cs23_2 * d13_2 * d23_2 + 16 * cs12_2_cs23_2 * d23_4 - 8 * cs12_2
			* d12_4 - 48 * cs12_2 * d12_2 * d13_2 + 48 * cs12_2 * d12_2 * d23_2 - 40 * cs12_2 * d13_4 + 80 * cs12_2
			* d13_2 * d23_2 - 40 * cs12_2 * d23_4 - 32 * cs12_cs13_3_cs23 * d12_4 - 32 * cs12_cs13_3_cs23 * d12_2
			* d13_2 + 64 * cs12_cs13_3_cs23 * d12_2 * d23_2 - 32 * cs12_cs13_3_cs23 * d13_2 * d23_2 - 32
			* cs12_cs13_3_cs23 * d23_4 - 32 * cs12_cs13_cs23_3 * d12_4 - 128 * cs12_cs13_cs23_3 * d12_2 * d13_2 - 32
			* cs12_cs13_cs23_3 * d12_2 * d23_2 - 32 * cs12_cs13_cs23_3 * d13_4 - 32 * cs12_cs13_cs23_3 * d13_2 * d23_2
			+ 16 * cs12_cs13_cs23 * d12_4 + 160 * cs12_cs13_cs23 * d12_2 * d13_2 - 32 * cs12_cs13_cs23 * d12_2 * d23_2
			+ 16 * cs12_cs13_cs23 * d13_4 - 32 * cs12_cs13_cs23 * d13_2 * d23_2 + 16 * cs12_cs13_cs23 * d23_4 + 16
			* cs13_4 * d12_4 - 32 * cs13_4 * d12_2 * d23_2 + 16 * cs13_4 * d23_4 + 48 * cs13_2_cs23_2 * d12_4 + 32
			* cs13_2_cs23_2 * d12_2 * d13_2 - 64 * cs13_2_cs23_2 * d12_2 * d23_2 + 16 * cs13_2_cs23_2 * d13_4 + 16
			* cs13_2_cs23_2 * d23_4 - 40 * cs13_2 * d12_4 - 48 * cs13_2 * d12_2 * d13_2 + 80 * cs13_2 * d12_2 * d23_2
			- 8 * cs13_2 * d13_4 + 48 * cs13_2 * d13_2 * d23_2 - 40 * cs13_2 * d23_4 + 16 * cs23_4 * d12_4 + 64
			* cs23_4 * d12_2 * d13_2 + 16 * cs23_4 * d13_4 - 40 * cs23_2 * d12_4 - 112 * cs23_2 * d12_2 * d13_2 + 48
			* cs23_2 * d12_2 * d23_2 - 40 * cs23_2 * d13_4 + 48 * cs23_2 * d13_2 * d23_2 - 8 * cs23_2 * d23_4 + 24
			* d12_4 + 48 * d12_2 * d13_2 - 48 * d12_2 * d23_2 + 24 * d13_4 - 48 * d13_2 * d23_2 + 24 * d23_4;
	double e2 = -16 * cs12_2_cs23_2 * d12_4 * d13_2 - 16 * cs12_2_cs23_2 * d13_6 + 32 * cs12_2_cs23_2 * d13_4 * d23_2
			- 16 * cs12_2_cs23_2 * d13_2 * d23_4 + 8 * cs12_2 * d12_4 * d13_2 - 8 * cs12_2 * d12_4 * d23_2 + 16
			* cs12_2 * d12_2 * d13_4 - 32 * cs12_2 * d12_2 * d13_2 * d23_2 + 16 * cs12_2 * d12_2 * d23_4 + 8 * cs12_2
			* d13_6 - 24 * cs12_2 * d13_4 * d23_2 + 24 * cs12_2 * d13_2 * d23_4 - 8 * cs12_2 * d23_6 + 32
			* cs12_cs13_cs23_3 * d12_4 * d13_2 + 32 * cs12_cs13_cs23_3 * d12_2 * d13_4 + 32 * cs12_cs13_cs23_3 * d12_2
			* d13_2 * d23_2 + 8 * cs12_cs13_cs23 * d12_6 - 40 * cs12_cs13_cs23 * d12_4 * d13_2 - 8 * cs12_cs13_cs23
			* d12_4 * d23_2 - 40 * cs12_cs13_cs23 * d12_2 * d13_4 + 48 * cs12_cs13_cs23 * d12_2 * d13_2 * d23_2 - 8
			* cs12_cs13_cs23 * d12_2 * d23_4 + 8 * cs12_cs13_cs23 * d13_6 - 8 * cs12_cs13_cs23 * d13_4 * d23_2 - 8
			* cs12_cs13_cs23 * d13_2 * d23_4 + 8 * cs12_cs13_cs23 * d23_6 - 16 * cs13_2_cs23_2 * d12_6 + 32
			* cs13_2_cs23_2 * d12_4 * d23_2 - 16 * cs13_2_cs23_2 * d12_2 * d13_4 - 16 * cs13_2_cs23_2 * d12_2 * d23_4
			+ 8 * cs13_2 * d12_6 + 16 * cs13_2 * d12_4 * d13_2 - 24 * cs13_2 * d12_4 * d23_2 + 8 * cs13_2 * d12_2
			* d13_4 - 32 * cs13_2 * d12_2 * d13_2 * d23_2 + 24 * cs13_2 * d12_2 * d23_4 - 8 * cs13_2 * d13_4 * d23_2
			+ 16 * cs13_2 * d13_2 * d23_4 - 8 * cs13_2 * d23_6 - 32 * cs23_4 * d12_4 * d13_2 - 32 * cs23_4 * d12_2
			* d13_4 + 8 * cs23_2 * d12_6 + 56 * cs23_2 * d12_4 * d13_2 - 16 * cs23_2 * d12_4 * d23_2 + 56 * cs23_2
			* d12_2 * d13_4 - 64 * cs23_2 * d12_2 * d13_2 * d23_2 + 8 * cs23_2 * d12_2 * d23_4 + 8 * cs23_2 * d13_6
			- 16 * cs23_2 * d13_4 * d23_2 + 8 * cs23_2 * d13_2 * d23_4 - 8 * d12_6 - 24 * d12_4 * d13_2 + 24 * d12_4
			* d23_2 - 24 * d12_2 * d13_4 + 48 * d12_2 * d13_2 * d23_2 - 24 * d12_2 * d23_4 - 8 * d13_6 + 24 * d13_4
			* d23_2 - 24 * d13_2 * d23_4 + 8 * d23_6;
	double e0 = 16 * cs23_4 * d12_4 * d13_4 - 8 * cs23_2 * d12_6 * d13_2 - 16 * cs23_2 * d12_4 * d13_4 + 16 * cs23_2
			* d12_4 * d13_2 * d23_2 - 8 * cs23_2 * d12_2 * d13_6 + 16 * cs23_2 * d12_2 * d13_4 * d23_2 - 8 * cs23_2
			* d12_2 * d13_2 * d23_4 + d12_8 + 4 * d12_6 * d13_2 - 4 * d12_6 * d23_2 + 6 * d12_4 * d13_4 - 12 * d12_4
			* d13_2 * d23_2 + 6 * d12_4 * d23_4 + 4 * d12_2 * d13_6 - 12 * d12_2 * d13_4 * d23_2 + 12 * d12_2 * d13_2
			* d23_4 - 4 * d12_2 * d23_6 + d13_8 - 4 * d13_6 * d23_2 + 6 * d13_4 * d23_4 - 4 * d13_2 * d23_6 + d23_8;

	cfs[0] = e8;
	cfs[1] = e6;
	cfs[2] = e4;
	cfs[3] = e2;
	cfs[4] = e0;
}
int p3p_solve_xs(const double d[], const double cs[], const double cfs[], double x[], double err[]) {
	double a, b, c;
	double x1[4], x2[2], x3[2];

	double d12_2 = d[0] * d[0];
	double d13_2 = d[1] * d[1];
	double d23_2 = d[2] * d[2];

	int nx = 0;

	int nx1 = quartic(cfs[1] / cfs[0], cfs[2] / cfs[0], cfs[3] / cfs[0], cfs[4] / cfs[0], x1);
	if (nx1 == 0)
		return 0;
	int i, j, k, num = 0;
	for (i = 0; i < nx1; ++i) {
		if (x1[i] < 0)
			continue;
		a = sqrt(x1[i]);
		//solve x2
		int nx2 = quadratic(-2 * a * cs[0], x1[i] - d[0] * d[0], x2);
		for (j = 0; j < nx2; ++j) {
			if (x2[j] < 0)
				continue;
			b = x2[j];

			//solve x3
//			c = ((x1[idx1] - x2[idx2] * x2[idx2]) - (d13_2 - d23_2)) / (2 * (a * cs[1] - b * cs[2]));
//			//double a1 = ((b * b - c * c) - (d12_2 - d13_2)) / (2 * (b * cs[0] - c * cs[1]));
//			if (c >= 0) 
//			{
//				x[num] = a;
//				x[num + 1] = b;
//				x[num + 2] = c;
//				num += 3;
//			}

			int nx3 = quadratic(-2 * a * cs[1], x1[i] - d[1] * d[1], x3);
			for (k = 0; k < nx3; ++k) {
				if (x3[k] < 0)
					continue;
				c = x3[k];
				x[num] = a;
				x[num + 1] = b;
				x[num + 2] = c;
				num += 3;
			}
		}
	}
	nx = num / 3;
	if (nx == 0)
		return 0;

	if (err) {
		//compute the error for each solution
		for (num = 0; num < nx; ++num) {
			a = x[3 * num];
			b = x[3 * num + 1];
			c = x[3 * num + 2];
			double aa = a * a;
			double bb = b * b;
			double cc = c * c;

			//			err[3 * num] = fabs(aa + bb - d12_2 - 2 * a * b * cs[0]);
			//			err[3 * num + 1] = fabs(aa + cc - d13_2 - 2 * a * c * cs[1]);
			//			err[3 * num + 2] = fabs(bb + cc - d23_2 - 2 * b * c * cs[2]);
			err[num] = fabs(aa + bb - d12_2 - 2 * a * b * cs[0]) + fabs(aa + cc - d13_2 - 2 * a * c * cs[1]) + fabs(bb
					+ cc - d23_2 - 2 * b * c * cs[2]);
		}
	}
	return nx;
}
int p3p_select_best_two(const int nx, const double x[], const double err[], double bx[], double berr[]) {
	if (nx == 0)
		repErr("p3p_select_two_best - no solution is found.");

	int min_i1 = 0, min_i2 = 0;
	double min_err1 = err[0];
	double min_err2 = err[0];

	int i;
	for (i = 0; i < nx; ++i) {
		//test
		logInfo("%lf %lf %lf [%lf]\n",x[3*i],x[3*i+1],x[3*i+2],err[i]);
		if (err[i] < min_err1) {
			min_i2 = min_i1;
			min_i1 = i;
			min_err2 = min_err1;
			min_err1 = err[i];
		}
	}
	int num = 0;
	const double* px = x;

	if (min_i1 == min_i2) {
		px = x + 3 * min_i1;
		bx[0] = px[0];
		bx[1] = px[1];
		bx[2] = px[2];
		berr[0] = min_err1;
		num = 1;
	} else {
		px = x + 3 * min_i1;
		bx[0] = px[0];
		bx[1] = px[1];
		bx[2] = px[2];
		px = x + 3 * min_i2;
		bx[3] = px[0];
		bx[4] = px[1];
		bx[5] = px[2];
		berr[0] = min_err1;
		berr[1] = min_err2;
		num = 2;
	}
	return num;
}

void get_3dpt(const double iK[], const double m[], double r, double M[]) {
	double ux = iK[0] * m[0] + iK[1] * m[1] + iK[2];
	double uy = iK[4] * m[1] + iK[5];
	double s = sqrt(ux * ux + uy * uy + 1);
	M[0] = ux * r / s;
	M[1] = uy * r / s;
	M[2] = r / s;
}

//static double g_Ms[9];
//static double g_ms[6];
//
//static void project_residual(int *m, int *n, double *x, double *fvec, int *iflag) {
//	/* project the point into the two views */
//	double p[2], q[2];
//
//	project(global_R0, global_t0, x, p);
//	project(global_R1, global_t1, x, q);
//
//
//}

int p3p(const double K[], const double M1[], const double M2[], const double M3[], const double m1[],
		const double m2[], const double m3[], double Rs[], double ts[]) {
	double d[3], cs[3], cfs[5];
	double iK[9];
	mat33Inv(K, iK);
	p3p_get_coeffs(K, M1,M2,M3,m1,m2,m3, d, cs);
	p3p_get_quartic(d, cs, cfs);

	double ox[8 * 3], oerr[8 * 3];
	int onx = p3p_solve_xs(d, cs, cfs, ox, oerr);
	if (onx == 0)
		return 0;

	double x[6], err[2];
	int nx = p3p_select_best_two(onx, ox, oerr, x, err);
	
	double M1_[3],M2_[3],M3_[3];
	double* px = x, *R = Rs, *t = ts;
	int i;
	for (i = 0; i < nx; ++i) {
		get_3dpt(iK, m1, px[0], M1_);
		get_3dpt(iK, m2, px[1], M2_);
		get_3dpt(iK, m3, px[2], M3_);

		absOrient3(M1,M2,M3,M1_,M2_,M3_, R, t);
		
		//refine the estimation	
		px += 3;
		R += 9;
		t += 3;
	}
	return nx;
}
/*
void test1() {
	double K[9] = { 854.584441000000, 0, 337.046724000000, 0, 1136.89639300000, 271.145441000000, 0, 0, 1 };
	double M[9] = { 222.133913488772, 205.332638509029, 248.094904042753, 200.231711206703, 238.745523235575,
			240.865161032672, 243.434735268175, 204.221792275546, 219.989132454945 };
	double m[6] = { 2439.70172435968, -107.150415049024, 2935.48227616112, -511.591197673871, 2009.10995031492,
			-138.050327593076 };
	//solution for test
	double s[3] = { 367.782896097117, 369.143392631412, 358.063386754685 };
	double d[3], cs[3], cfs[5];
	int i;
	double x[8 * 3], err[8 * 3];
	double iK[9];
	mat33Inv(K, iK);

	double* px;

	logInfo("1 2 3\n");
	p3p_get_coeffs(K, M, M + 3, M + 6, m, m + 2, m + 4, d, cs);
	p3p_get_quartic(d, cs, cfs);
	int nx = p3p_solve_xs(d, cs, cfs, x, err);
	if (nx == 0)
		logInfo("no solution is found!\n");
	px = x;
	for (i = 0; i < nx; ++i) {
		//test
		logInfo("%lf %lf %lf | %lf\t%lf\t%lf\t| [%lf]\n", px[0], px[1], px[2], px[0] - s[0], px[1] - s[1], px[2] - s[2],
				err[i]);
		px += 3;
	}

	logInfo("1 3 2\n");
	p3p_get_coeffs(K, M, M + 6, M + 3, m, m + 4, m + 2, d, cs);
	p3p_get_quartic(d, cs, cfs);
	nx = p3p_solve_xs(d, cs, cfs, x, err);
	if (nx == 0)
		logInfo("no solution is found!\n");
	px = x;
	for (i = 0; i < nx; ++i) {
		//test
		logInfo("%lf %lf %lf | %lf\t%lf\t%lf\t| [%lf]\n", px[0], px[2], px[1], px[0] - s[0], px[2] - s[1], px[1] - s[2],
				err[i]);
		px += 3;
	}

	logInfo("2 1 3\n");
	p3p_get_coeffs(K, M + 3, M, M + 6, m + 2, m, m + 4, d, cs);
	p3p_get_quartic(d, cs, cfs);
	nx = p3p_solve_xs(d, cs, cfs, x, err);
	if (nx == 0)
		logInfo("no solution is found!\n");
	px = x;
	for (i = 0; i < nx; ++i) {
		//test
		logInfo("%lf %lf %lf | %lf\t%lf\t%lf\t| [%lf]\n", px[1], px[0], px[2], px[1] - s[0], px[0] - s[1], px[2] - s[2],
				err[i]);
		px += 3;
	}

	logInfo("2 3 1\n");
	p3p_get_coeffs(K, M + 3, M + 6, M, m + 2, m + 4, m, d, cs);
	p3p_get_quartic(d, cs, cfs);
	nx = p3p_solve_xs(d, cs, cfs, x, err);
	if (nx == 0)
		logInfo("no solution is found!\n");
	px = x;
	for (i = 0; i < nx; ++i) {
		//test
		logInfo("%lf %lf %lf | %lf\t%lf\t%lf\t| [%lf]\n", px[2], px[0], px[1], px[2] - s[0], px[0] - s[1], px[1] - s[2],
				err[i]);
		px += 3;
	}

	logInfo("3 1 2\n");
	p3p_get_coeffs(K, M + 6, M, M + 3, m + 4, m, m + 2, d, cs);
	p3p_get_quartic(d, cs, cfs);
	nx = p3p_solve_xs(d, cs, cfs, x, err);
	if (nx == 0)
		logInfo("no solution is found!\n");
	px = x;
	for (i = 0; i < nx; ++i) {
		//test
		logInfo("%lf %lf %lf | %lf\t%lf\t%lf\t| [%lf]\n", px[1], px[2], px[0], px[1] - s[0], px[2] - s[1], px[0] - s[2],
				err[i]);
		px += 3;
	}

	logInfo("3 2 1\n");
	p3p_get_coeffs(K, M + 6, M + 3, M, m + 4, m + 2, m, d, cs);
	p3p_get_quartic(d, cs, cfs);
	nx = p3p_solve_xs(d, cs, cfs, x, err);
	if (nx == 0)
		logInfo("no solution is found!\n");
	px = x;
	for (i = 0; i < nx; ++i) {
		//test
		logInfo("%lf %lf %lf | %lf\t%lf\t%lf\t| [%lf]\n", px[2], px[1], px[0], px[2] - s[0], px[1] - s[1], px[0] - s[2],
				err[i]);
		px += 3;
	}
}
int test2() {
	double K[9] = { 854.584441000000, 0, 337.046724000000, 0, 1136.89639300000, 271.145441000000, 0, 0, 1 };
	double M[9] = { 222.133913488772, 205.332638509029, 248.094904042753, 200.231711206703, 238.745523235575,
			240.865161032672, 243.434735268175, 204.221792275546, 219.989132454945 };
	double m[6] = { 2439.70172435968, -107.150415049024, 2935.48227616112, -511.591197673871, 2009.10995031492,
			-138.050327593076 };

	double d[3], cs[3], cfs[5];
	int i;
	double ox[8 * 3], oerr[8 * 3];

	//	tic();
	//	for (idx1 = 0; idx1 < 100000; ++idx1) {
	//		p3p_get_coeffs(K, M1, M2, M3, m1, m2, m3, d, cs);
	//		p3p_get_quartic(d, cs, cfs);
	//		int nx = p3p_solve_xs(d, cs, cfs, x, err);
	//		//printMat(nx, 3, x);
	//		//printMat(nx, 1, err);
	//	}
	//	toc();

	double iK[9];

	mat33Inv(K, iK);
	//test
	logInfo("M:\n");
	printMat(3, 3, M);
	p3p_get_coeffs(K, M + 6, M, M + 3, m + 4, m, m + 2, d, cs);
	//test
	logInfo("d:\n");
	logInfo("%lf %lf %lf\n", d[0], d[1], d[2]);
	//test
	logInfo("cs:\n");
	logInfo("%lf %lf %lf\n", cs[0], cs[1], cs[2]);

	p3p_get_quartic(d, cs, cfs);
	//test
	logInfo("[e8 e6 e4 e2 e0]:%lf %lf %lf %lf %lf\n", cfs[0], cfs[1], cfs[2], cfs[3], cfs[4]);

	int onx = p3p_solve_xs(d, cs, cfs, ox, oerr);
	if (onx == 0)
		logInfo("no solution is found!\n");

	double x[6], err[3];
	int nx = p3p_select_best_two(onx, ox, oerr, x, err);

	double* px = x;
	for (i = 0; i < nx; ++i) {
		//test
		//info("solotion:%lf %lf %lf [%lf]\n", px[0], px[1], px[2], oerr[3 * idx1], oerr[3 * idx1 + 1], oerr[3 * idx1 + 2]);
		logInfo("solotion:%lf %lf %lf [%lf]\n", px[0], px[1], px[2], err[i]);
		//		get_3dpt(iK,m,px[0],M_);
		//		get_3dpt(iK,m+2,px[1],M_+3);
		//		get_3dpt(iK,m+4,px[2],M_+6);
		//
		//		//test
		//		info("M_\n");
		//		printMat(3,3,M_);
		//		
		//		abs_orient(3,M,M_,R,t);
		//		
		////	printMat(3,3,M_);
		//		info("R\n");
		//		printMat(3,3,R);
		//		info("t\n");
		//		printMat(3,1,t);
		px += 3;
	}
	return 0;
}*/
//void test_p3p(){
//	double K[9] = { 854.584441000000, 0, 337.046724000000, 0, 1136.89639300000, 271.145441000000, 0, 0, 1 };
//	double M[9] = { 222.133913488772, 205.332638509029, 248.094904042753, 200.231711206703, 238.745523235575,
//			240.865161032672, 243.434735268175, 204.221792275546, 219.989132454945 };
//	double m[6] = { 2439.70172435968, -107.150415049024, 2935.48227616112, -511.591197673871, 2009.10995031492,
//			-138.050327593076 };
//	
//	double Rs[18];
//	double ts[6];
//	int npos = p3p(K,M,M+3,M+6,m,m+2,m+4,Rs,ts);
//	
//	int i;
//	for( i = 0; i < npos; ++i){
//		logInfo("R[%d]\n",i);
//		printMat(3,3,Rs+9*i);
//		logInfo("t[%d]\n",i);
//		printMat(1,3,ts+3*i);
//	}
//}
//int main() {
//	//test_p3p();
//	//main1();
//	test_p3p();
//	return 0;
//}
