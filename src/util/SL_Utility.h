/*
 * SL_Utility.h
 *
 *  Created on: 2011-6-14
 *      Author: Danping Zou
 */

#ifndef SL_UTILITY_H_
#define SL_UTILITY_H_
#include <algorithm>
#include <utility>
#include <cassert>
#include <functional>
#include "math/SL_Matrix.h"

template<class T, class U>
class PairLess {
public:
	bool operator ()(const std::pair<T, U>& pair1
			, const std::pair<T, U>& pair2) const {
		return pair1.first < pair2.first;
	}
};

template<class T, class U>
class PairGreater {
public:
	bool operator ()(const std::pair<T, U>& pair1
			, const std::pair<T, U>& pair2) const {
		return pair1.first > pair2.first;
	}
};

template<class T, class U>
void sortWithInd(U n, T* x, U ind[], bool less = true) {
	assert(n > 0);
	//construct a pair array
	std::pair<T, U>* xi = new std::pair<T, U>[n];
	for (U i = 0; i < n; i++) {
		xi[i].first = x[i];
		xi[i].second = i;
	}

	if (less)
		std::sort(xi, xi + n, PairLess<T, U>());
	else
		std::sort(xi, xi + n, PairGreater<T, U>());

	for (U i = 0; i < n; i++) {
		ind[i] = xi[i].second;
	}

	delete[] xi;
}
template<class T, class U>
void partialSortWithInd(U n, T*x, U ind[], U top_num, bool less = true) {
	assert(n > 0 && top_num > 0 && n >= top_num);
	//construct a pair array
	std::pair<T, U>* xi = new std::pair<T, U>[n];
	for (U i = 0; i < n; i++) {
		xi[i].first = x[i];
		xi[i].second = i;
	}

	if (less)
		std::partial_sort(xi, xi + top_num, xi + n, PairLess<T, U>());
	else
		std::partial_sort(xi, xi + top_num, xi + n, PairGreater<T, U>());

	for (U i = 0; i < n; i++) {
		ind[i] = xi[i].second;
	}

	delete[] xi;
}
template<class T, class U>
void cumSum(const U n, const T* x, T* cx) {
	T s = 0;
	for (U i = 0; i < n; i++) {
		s += x[i];
		cx[i] = s;
	}
}

template<class T, class U>
void cumSum(int m, int n, const T* mat, U* sum, bool row) {
	if (row) {
		for (int i = 0; i < m; ++i) {
			sum[i * n] = static_cast<U>(mat[i * n]);
			for (int j = 1; j < n; ++j) {
				sum[i * n + j] = sum[i * n + j - 1] + mat[i * n + j];
			}
		}
	} else {
		for (int j = 0; j < n; ++j) {
			sum[j] = static_cast<U>(mat[j]);
			for (int i = 1; i < m; ++i) {
				sum[i * n + j] = sum[(i - 1) * n + j] + mat[i * n + j];
			}
		}
	}
}

template<class T, class U>
void cumSum(const MyMat<T>& mat, MyMat<U>& sum, bool row) {
	sum.resize(mat.m, mat.n);
	cumSum(mat.m, mat.n, mat.data, sum.data, row);
}
#endif /* SL_UTILITY_H_ */
