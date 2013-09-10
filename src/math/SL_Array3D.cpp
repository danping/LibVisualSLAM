/*
 * SL_Array3D.cpp
 *
 *  Created on: Jul 14, 2012
 *      Author: Danping Zou
 */
#include "math/SL_Array3D.h"
#include <cstdlib>
#include <ctime>

void arr3Rand(int c, int m, int n, Arr3d& arr) {
	arr.resize(3, m, n);
	int i, len = 3 * m * n;
	for (i = 0; i < len; ++i) {
		double rnd_val = double(rand()) / RAND_MAX - 0.5;
		arr.data[i] = rnd_val;
	}
}

