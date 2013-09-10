/*
 * SL_BoxFilter.h
 *
 *  Created on: Dec 18, 2011
 *      Author: Danping Zou
 */

#ifndef SL_BOXFILTER_H_
#define SL_BOXFILTER_H_
#include "math/SL_Matrix.h"
#include "util/SL_Utility.h"

template<class T, class U>
void boxFilt2(const MyMat<T>& in, MyMat<U>& out, int r) {
	assert(r > 0);
	out.resize(in.m, in.n);
	for (int i = 0; i < in.m; i++) {
		for (int j = 0; j < in.n; j++) {
			U v = 0;
			for (int s = -r; s <= r; ++s) {
				for (int t = -r; t <= r; ++t) {
					int ii = i + s;
					int jj = j + t;
					if (ii < 0 || ii >= in.m || jj < 0 || jj >= in.n)
						continue;
					v += static_cast<U>(in.data[ii * in.n + jj]);
				}
			}
			out.data[i * in.n + j] = v;
		}
	}
}
template<class T, class U>
void fastBoxFilt2(int m, int n, const T* in, U* out, int r) {
	MyMat<U> tempSum(m, n);
	cumSum(m, n, in, tempSum.data, false);

	MyMat<T> tmpRow(m, n);
	for (int j = 0; j < n; ++j)
		for (int i = 0; i <= r; ++i)
			tmpRow.data[i * n + j] = tempSum.data[(i + r) * n + j];

	for (int j = 0; j < n; ++j)
		for (int i = r + 1; i <= m - r - 1; ++i)
			tmpRow.data[i * n + j] = tempSum.data[(i + r) * n + j]
					- tempSum.data[(i - r - 1) * n + j];
	for (int j = 0; j < n; ++j)
		for (int i = m - r; i < m; ++i)
			tmpRow.data[i * n + j] = tempSum.data[(m - 1) * n + j]
					- tempSum.data[(i - r - 1) * n + j];

	cumSum(m, n, tmpRow.data, tempSum.data, true);

	for (int i = 0; i < m; ++i)
		for (int j = 0; j <= r; ++j)
			out[i * n + j] = tempSum.data[i * n + j + r];

	for (int i = 0; i < m; ++i)
		for (int j = r + 1; j <= n - r - 1; ++j)
			out[i * n + j] = tempSum.data[i * n + j + r]
					- tempSum.data[i * n + j - r - 1];

	for (int i = 0; i < m; ++i)
		for (int j = n - r; j < n; ++j)
			out[i * n + j] = tempSum.data[i * n + n - 1]
					- tempSum.data[i * n + j - r - 1];
}
/* fast box filter using cumsum*/
template<class T, class U>
void fastBoxFilt2(const MyMat<T>& in, MyMat<U>& out, int r) {
	assert(r > 0);
	out.resize(in.m, in.n);
	fastBoxFilt2(in.m, in.n, in.data, out.data, r);
}
#endif /* SL_BOXFILTER_H_ */
