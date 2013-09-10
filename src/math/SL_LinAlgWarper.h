/*
 * SL_LinAlgWarper.h
 *
 *  Created on: 2010-11-15
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_LINALGWARPER_H_
#define SL_LINALGWARPER_H_
#include "math/SL_LinAlg.h"
#include "math/SL_Matrix.h"
inline void matScale(Mat_d& A, double scale) {
	matScale(A.rows, A.cols, A.data, scale, A.data);
}
inline void matAddScale(Mat_d& A, double val) {
	matAddScale(A.rows, A.cols, A.data, val);
}
inline void matTrans(const Mat_d& A, Mat_d& AT) {
	AT.resize(A.n, A.m);
	matTrans(A.m, A.n, A.data, AT.data);
}
inline void matAB(const Mat_d& A, const Mat_d& B, Mat_d& C) {
	C.resize(A.rows, B.cols);
	matAB(A.rows, A.cols, B.rows, B.cols, A.data, B.data, C.data);
}
inline void matATB(const Mat_d& A, const Mat_d& B, Mat_d& C) {
	C.resize(A.cols, B.cols);
	matATB(A.rows, A.cols, B.rows, B.cols, A.data, B.data, C.data);
}
inline void matAx(const Mat_d& A, const Mat_d& x, Mat_d& y) {
	assert(x.m == A.n && x.n == 1);
	y.resize(A.m, 1);
	matAxpy(A.m, A.n, 1.0, A.data, x.data, 0, 0, y.data);
}
inline void matAx(int m, int n, const double* A, const double* x, double* y) {
	matAxpy(m, n, 1.0, A, x, 0, 0, y);
}
inline void matQR(const Mat_d& A, Mat_d& Q, Mat_d& R) {
	assert(A.m > 0 && A.n > 0);
	Q.resize(A.m, A.m);
	Q.fill(0);
	R.resize(A.m, A.n);
	R.fill(0);
	dgeqrFor(A.m, A.n, A.data, Q.data, R.data);
}
#endif /* SL_LINALGWARPER_H_ */
