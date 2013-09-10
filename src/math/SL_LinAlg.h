/*
 * SL_LinAlg.h
 *
 *  Created on: 2010-11-6
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_LINALG_H_
#define SL_LINALG_H_
#ifdef WIN32
#include "atlas/cblas.h"
#else
#include "cblas.h"
#endif

#include <cmath>
#include <cstring>

inline void doubleArrCopy(double* dest, int offset, const double* src,
		int len) {
	memcpy(dest + offset * len, src, sizeof(double) * len);
}

double diffSqrSum(int n, const double* v1, const double* v2);
/* compute inner product */
double innerProd3(const double* v1, const double* v2, const double* W = 0);
double innerProd2(const double* v1, const double* v2, const double* W = 0);
/* compute length*/
inline double vec3Len(const double* v) {
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}
/* copy sub matrix */
void copySubMat(int m, int n, const double* A, int i1, int i2, int j1, int j2,
		double* sub_A);
/* transpose the 3x3 matrix*/
inline void mat33Trans(const double * A, double* AT) {
	AT[0] = A[0];
	AT[3] = A[1];
	AT[6] = A[2];
	AT[1] = A[3];
	AT[4] = A[4];
	AT[7] = A[5];
	AT[2] = A[6];
	AT[5] = A[7];
	AT[8] = A[8];
}

inline void mat44Trans(const double * A, double* AT) {
	AT[0] = A[0];
	AT[4] = A[1];
	AT[8] = A[2];
	AT[12] = A[3];
	AT[1] = A[4];
	AT[5] = A[5];
	AT[9] = A[6];
	AT[13] = A[7];
	AT[2] = A[8];
	AT[6] = A[9];
	AT[10] = A[10];
	AT[14] = A[11];
	AT[3] = A[12];
	AT[7] = A[13];
	AT[11] = A[14];
	AT[15] = A[15];
}
inline void mat44Trans(const float * A, float* AT) {
	AT[0] = A[0];
	AT[4] = A[1];
	AT[8] = A[2];
	AT[12] = A[3];
	AT[1] = A[4];
	AT[5] = A[5];
	AT[9] = A[6];
	AT[13] = A[7];
	AT[2] = A[8];
	AT[6] = A[9];
	AT[10] = A[10];
	AT[14] = A[11];
	AT[3] = A[12];
	AT[7] = A[13];
	AT[11] = A[14];
	AT[15] = A[15];
}

/* compute the inverse of 2x2 matrix */
void mat22Inv(const double* A, double* Ainv);
/* compute the inverse of 3x3 matrix */
void mat33Inv(const double *A, double *Ainv);
/* r = A*b ; b,r : 3x1 matrix*/
void mat33ProdVec(const double *A, const double *b, double *r);
/* r = alpha*A*b + beta*c ; A:3x3 matrix, b,c : 3x1vector*/
void mat33ProdVec(const double *A, const double *b, const double* c, double *r,
		double alpha, double beta);
/* r = A'*b; b,r : 3x1 matrix*/
void mat33TransProdVec(const double* A, const double* b, double *r);
/* C = A*B ; A, B, C : 3x3 matrix */
void mat33AB(const double* A, const double* B, double *C);
/* C = A^T*B ; A, B, C : 3x3 matrix */
void mat33ATB(const double* A, const double* B, double *C);
/* C = A*B^T; A, B, C : 3x3 matrix*/
void mat33ABT(const double* A, const double* B, double *C);

/* trace(A)*/
double mat33Tr(const double* A);

/* det(A) */
double mat22Det(const double* A);
double mat33Det(const double* A);
double mat44Det(const double* A);
double matDet(int m, const double* A);

/*check whether an array contains nan values*/
bool containNan(int n, double *a);
/* compute matrix transpose */
void matTrans(int m, int n, const double* Ar, double* Ac);
void matCol2Row(int m, int n, const double* Ac, double* Ar);
void matRow2Col(int m, int n, const double* Ar, double* Ac);

/* fill the mat with value*/
void matFill(int m, int n, double* A, double val);
/* B = A.*s */
void matScale(int m, int n, double* A, double s, double *B);
/* B = A + s*/
void matAddScale(int m, int n, double* A, double s);
void matAddScale(int m, int n, double* A, double s, double* B);
/* R = A - B */
void matDiff(int m, int n, double* A, double* B, double* R);
inline void mat33Diff(double* A, double* B, double* R) {
	R[0] = A[0] - B[0];
	R[1] = A[1] - B[1];
	R[2] = A[2] - B[2];
	R[3] = A[3] - B[3];
	R[4] = A[4] - B[4];
	R[5] = A[5] - B[5];
	R[6] = A[6] - B[6];
	R[7] = A[7] - B[7];
	R[8] = A[8] - B[8];
}
/* S = A + B */
void matSum(int m, int n, const double * A, const double * B, double* S);

/* C = A - B */
template<class T, class U>
void matSub(int m, int n, const T* A, const U* B, T* C) {
	int i, len = m * n;
	for (i = 0; i < len; i++)
		C[i] = A[i] - B[i];
}
/* S = (A + B)/2*/
void matAvg(int m, int n, const double* A, const double* B, double* S);

/* C = A.*B*/
template<class T, class U>
void matEleProd(int m, int n, const T *A, const U *B, U *C) {
	int i, len = m * n;
	for (i = 0; i < len; ++i) {
		C[i] = A[i] * B[i];
	}
}
/* C = A./B*/
template<class T, class U>
void matEleDiv(int m, int n, const T* A, const U* B, T* C) {
	int i, len = m * n;
	for (i = 0; i < len; ++i) {
		C[i] = A[i] / B[i];
	}
}
/* C = A B*/
void matAB(int Am, int An, int Bm, int Bn, const double *A, const double *B,
		double * C);
/* R = A^T B */
void matATB(int Am, int An, int Bm, int Bn, const double *A, const double *B,
		double *R);
/* R = A B^T */
void matABT(int Am, int An, int Bm, int Bn, const double *A, const double *B,
		double *C);
/* z = A*x + y*/
void matAxpy(int m, int n, double alpha, const double* A, const double* x,
		double beta, const double* y, double* z);
/* matrix inverse*/
bool matInv(int n, const double* A, double *invA);
bool matInvInplace(int n, double* A);

/* call fortran functions */
/* solve A*x = b by linear least square*/
int dgelsyFor(int m, int n, int nrhs, const double *A, double *b, double *x);
/* singular value decomposition A = u*s*v^T */
int dgesvdFor(int m, int n, const double *A, double * U, double *S, double *VT);
/* svd of A (only compute S,V^T)*/
int dgesvdFor(int m, int n, const double *A, double *S, double *VT);

/* n: the order of matrix A
 * A: matrix for which the eigenvectors/values are to be computed
 * evec: output array containing the eigenvectors (each row represents an eigen vector)
 * eval: output array containing the eigenvalues
 */
int dgeevFor(int n, const double *A, double *evec, double *eval);

/* QR decomposition*/
int dgeqrFor(int m, int n, const double* A, double* Q, double* R);

#endif /* SL_LINALG_H_ */
