/*
 * SL_ConvexHull2D.cpp
 *
 *  Created on: Feb 19, 2012
 *      Author: Danping Zou
 */

#include "SL_ConvexHull2D.h"
#include <cmath>

static int ccw(coord **P, int i, int j, int k) {
	coord a = P[i][0] - P[j][0], b = P[i][1] - P[j][1], c = P[k][0] - P[j][0], d = P[k][1] - P[j][1];
	return a * d - b * c <= 0; /* true if points i, j, k counterclockwise */
}

#define CMPM(c,A,B) \
	v = (*(coord**)A)[c] - (*(coord**)B)[c];\
	if (v>0) return 1;\
	if (v<0) return -1;

static int cmpl(const void *a, const void *b) {
	double v;
	CMPM(0,a,b);
	CMPM(1,b,a);
	return 0;
}

static int cmph(const void *a, const void *b) {
	return cmpl(b, a);
}

static int makeChain(coord** V, int n, int(*cmp)(const void*, const void*)) {
	int i, j, s = 1;
	coord* t;

	qsort(V, n, sizeof(coord*), cmp);
	for (i = 2; i < n; i++) {
		for (j = s; j >= 1 && ccw(V, i, j, j - 1); j--) {
		}
		s = j + 1;
		t = V[s];
		V[s] = V[i];
		V[i] = t;
	}
	return s;
}

static int ch2d(coord **P, int n) {
	int u = makeChain(P, n, cmpl); /* make lower hull */
	if (!n)
		return 0;
	P[n] = P[0];
	return u + makeChain(P + u, n - u + 1, cmph); /* make upper hull */
}

#define CONVEX_HULL_MAX_N 10000
void get2DConvexHull(const std::vector<double>& pts, std::vector<double>& cxh) {
	double points[CONVEX_HULL_MAX_N][2], *P[CONVEX_HULL_MAX_N];
	int n = 0;
	const size_t npts = pts.size() / 2;
	for (size_t i = 0; i < npts; i++) {
		points[i][0] = pts[2 * i];
		points[i][1] = pts[2 * i + 1];
		P[i] = points[i];
		n++;
	}

	int m = ch2d(P, n);
	cxh.reserve(4 * m);

	for (int i = 0; i < m; i++) {
		cxh.push_back(P[i][0]);
		cxh.push_back(P[i][1]);
	}
}
double getPolyArea(const std::vector<double>& poly) {
	size_t npts = poly.size() / 2;
	if (npts < 3)
		return 0;
	double s = 0;
	for (size_t i = 1; i < npts; i++) {
		double x0 = poly[2 * i - 2];
		double y0 = poly[2 * i - 1];
		double x1 = poly[2 * i];
		double y1 = poly[2 * i + 1];
		s += x0 * y1 - x1 * y0;
	}
	if (poly[0] != poly[2 * npts - 2] || poly[1] != poly[2 * npts - 1]) {
		s += poly[2 * npts - 2] * poly[1] - poly[2 * npts - 1] * poly[0];
	}
	return 0.5 * fabs(s);
}