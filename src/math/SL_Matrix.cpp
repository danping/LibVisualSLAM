/*
 * SL_Matrix.cpp
 *
 *  Created on: 2010-12-26
 *      Author: Danping Zou
 */

#include "math/SL_Matrix.h"
#include <cstdlib>
#include <ctime>

void matRand(int m, int n, Mat_d& mat) {
	mat.resize(m, n);
	int i, len = m * n;
	for (i = 0; i < len; ++i) {
		double rnd_val = double(rand()) / RAND_MAX - 0.5;
		mat.data[i] = rnd_val;
	}
}
void matRand(int m, int n, Mat_f& mat) {
	mat.resize(m, n);
	int i, len = m * n;
	for (i = 0; i < len; ++i) {
		float rnd_val = float(rand()) / RAND_MAX - 0.5;
		mat.data[i] = rnd_val;
	}
}
void matOnes(int m, int n, Mat_d& mat) {
	mat.resize(m, n);
	mat.fill(1);
}
void matZeros(int m, int n, Mat_d& mat) {
	mat.resize(m, n);
	mat.fill(0);
}
void matZeros(int m, int n, Mat_f& mat) {
	mat.resize(m, n);
	mat.fill(0);
}

void matZeros(int m, int n, double* mat) {
	memset(mat, 0, sizeof(double) * m * n);
}

void matEyes(int n, Mat_d& mat) {
	mat.resize(n, n);
	mat.fill(0);
	int i;
	for (i = 0; i < n; ++i) {
		mat.data[i * n + i] = 1;
	}
}
void matEyes(int n, double* mat) {
	memset(mat, 0, sizeof(double) * n * n);
	int i;
	for (i = 0; i < n; ++i) {
		mat[i * n + i] = 1;
	}
}
