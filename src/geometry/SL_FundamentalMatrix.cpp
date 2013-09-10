/*
 * SL_FundamentalMatrix.cpp
 *
 *  Created on: 2010-11-12
 *      Author: Danping Zou
 */

#include "SL_FundamentalMatrix.h"
#include "math/SL_LinAlg.h"
#include "SL_error.h"
#include "SL_Triangulate.h"
#include <cmath>
#define Vx(v) (v)[0]
#define Vy(v) (v)[1]
#define Vz(v) (v)[2]

/* m2^T*F*m1*/
double epipolarError(const double F[9], const double m2[2], const double m1[2]) {
	double Fl[3], Fr[3], pt;

	Fl[0] = F[0] * m1[0] + F[1] * m1[1] + F[2] * 1;
	Fl[1] = F[3] * m1[0] + F[4] * m1[1] + F[5] * 1;
	Fl[2] = F[6] * m1[0] + F[7] * m1[1] + F[8] * 1;

	Fr[0] = F[0] * m2[0] + F[3] * m2[1] + F[6] * 1;
	Fr[1] = F[1] * m2[0] + F[4] * m2[1] + F[7] * 1;
	Fr[2] = F[2] * m2[0] + F[5] * m2[1] + F[8] * 1;

	pt = m2[0] * Fl[0] + m2[1] * Fl[1] + 1 * Fl[2];

	return 0.5 * (fabs(pt) / sqrt(Fl[0] * Fl[0] + Fl[1] * Fl[1]) + fabs(pt) / sqrt(Fr[0] * Fr[0] + Fr[1] * Fr[1]));
}

double epipolarError(const double* invKr, const double* invKl, const double* E, const double* m2, const double* m1) {
	double El[3], S;
	El[0] = E[0] * m1[0] + E[1] * m1[1] + E[2] * 1;
	El[1] = E[3] * m1[0] + E[4] * m1[1] + E[5] * 1;
	El[2] = E[6] * m1[0] + E[7] * m1[1] + E[8] * 1;
	S = m2[0] * El[0] + m2[1] * El[1] + 1 * El[2];
	double A[9], B[9];
	mat33ATB(invKr, E, A);
	mat33AB(E, invKl, B);
	double l1[3], l2[3];
	mat33ProdVec(A, m1, l1);
	mat33TransProdVec(B, m2, l2);

	return 0.5 * (fabs(S) / sqrt(l1[0] * l1[0] + l1[1] * l1[1]) + fabs(S) / sqrt(l2[0] * l2[0] + l2[1] * l2[1]));
}
double epipolarErrorAvg(int n, const double* F, const double *m2,  const double *m1) {
	int i;
	double s = 0;
	for (i = 0; i < n; ++i) {
		s += epipolarError(m2 + 2 * i, F, m1 + 2 * i);
	}
	return s / n;
}
