/*
 * SL_Gradient.h
 *
 *  Created on: Feb 13, 2012
 *      Author: Danping Zou
 */

#ifndef SL_GRADIENT_H_
#define SL_GRADIENT_H_
#include <cassert>
template<class U, class V>
void compute_dpx(const U& I, V& Ix) {
	assert(Ix.m == I.m && Ix.n == I.n);
	int m = I.m;
	int n = I.n;
	for (int i = 0; i < m; ++i) {
		int j = 0;
		for (; j < n -1;++j) {
			Ix.data[i * n + j ] = I.data[i * n + j+1] - I.data[i * n + j];
		}
		Ix.data[i * n + j ] = 0;
	}
}
		template<class U, class V>
void compute_dmx(const U& I, V& Ix) {
	assert(Ix.m == I.m && Ix.n == I.n);
	int m = I.m;
	int n = I.n;
	for (int i = 0; i < m; ++i) {
		Ix.data[i * n] = I.data[i * n];
		int j =  1;
		for (; j < n - 1; ++j) {
			Ix.data[i * n + j] = I.data[i * n + j] - I.data[i * n + j - 1];
		}
		Ix.data[i * n + j ] = -I.data[i *n + j];
	}
}
		
template<class U, class V>
void compute_dpy(const U& I, V& Iy) {
	assert(Iy.m == I.m && Iy.n == I.n);
	int m = I.m;
	int n = I.n;
	for (int j = 0; j < n; ++j) {
		int i = 0;
		for (; i < m - 1; ++i) {
			Iy.data[i * n + j] = I.data[(i+1) * n + j] - I.data[i * n + j];
		}
		Iy.data[i * n + j] = 0;
	}
}
template<class U, class V>
void compute_dmy(const U& I, V& Iy) {
	assert(Iy.m == I.m && Iy.n == I.n);
	int m = I.m;
	int n = I.n;
	for (int j = 0; j < n; ++j) {
		Iy.data[j] = I.data[j];
		int i = 1;
		for (; i < m - 1; ++i) {
			Iy.data[i * n + j] = I.data[i * n + j] - I.data[(i - 1) * n + j];
		}
		
		Iy.data[i* n + j] = - I.data[(i-1)*n + j]; 
	}
}
template<class U, class V>
void compute_div(const U& Ix, const U& Iy, V& divI) {
	assert(Ix.m == Iy.m && Ix.n == Iy.n);
	V Ixx(Ix.m, Ix.n);
	V Iyy(Iy.m, Iy.n);
	compute_dmx(Ix, Ixx);
	compute_dmy(Iy, Iyy);
	
	int len = Ix.m * Ix.n;
	for (int i = 0; i < len; ++i)
	divI.data[i] = Ixx.data[i] + Iyy.data[i];
}
#endif /* SL_GRADIENT_H_ */
