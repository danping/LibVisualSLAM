/*
 * SL_LinAlg.cpp
 *
 *  Created on: 2010-11-6
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#include "SL_LinAlg.h"
#include "SL_error.h"
#include "extern/lapack.h"
#include <cassert>
#include <cmath>
#include <algorithm>

double diffSqrSum(int n, const double* v1, const double* v2) {
	double s = 0;
	for (int i = 0; i < n; i++) {
		double d = v2[i] - v1[i];
		s += d * d;
	}
	return s;
}
double innerProd3(const double* v1, const double* v2, const double* W) {
	if (W == 0)
		return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	return v1[0] * (W[0] * v2[0] + W[1] * v2[1] + W[2] * v2[2])
			+ v1[1] * (W[3] * v2[0] + W[4] * v2[1] + W[5] * v2[2])
			+ v1[2] * (W[6] * v2[0] + W[7] * v2[1] + W[8] * v2[2]);
}
double innerProd2(const double* v1, const double* v2, const double* W) {
	if (W == 0)
		return v1[0] * v2[0] + v1[1] * v2[1];
	return (v1[0] * v2[0] * W[0] + 2 * v1[0] * v2[1] * W[1]
			+ v1[1] * v2[1] * W[2]);
}
void copySubMat(
		int m,
		int n,
		const double* A,
		int i1,
		int i2,
		int j1,
		int j2,
		double* sub_A) {
	if (i1 < 0 || i2 >= m || j1 < 0 || j2 >= n || i1 > i2 || j1 > j2)
		repErr("copy_sub_mat() error!");
	int i, cols = j2 - j1 + 1, k = 0;
	for (i = i1; i <= i2; ++i, k += cols) {
		memcpy(sub_A + k, A + i * n + j1, cols * sizeof(double));
	}
}
void mat22Inv(const double* A, double* Ainv) {
	double s = (A[0] * A[3] - A[1] * A[2]);
	Ainv[0] = A[3] / s;
	Ainv[1] = -A[1] / s;
	Ainv[2] = -A[2] / s;
	Ainv[3] = A[0] / s;
}
void mat33Inv(const double *A, double *Ainv) {
	double m1 = A[8] * A[4] - A[7] * A[5];
	double m2 = A[8] * A[1] - A[7] * A[2];
	double m3 = A[5] * A[1] - A[4] * A[2];

	double d = A[0] * m1 - A[3] * m2 + A[6] * m3;

	Ainv[0] = m1 / d;
	Ainv[1] = -m2 / d;
	Ainv[2] = m3 / d;

	Ainv[3] = (-A[8] * A[3] + A[6] * A[5]) / d;
	Ainv[4] = (A[8] * A[0] - A[6] * A[2]) / d;
	Ainv[5] = (-A[5] * A[0] + A[3] * A[2]) / d;

	Ainv[6] = (A[7] * A[3] - A[6] * A[4]) / d;
	Ainv[7] = (-A[7] * A[0] + A[6] * A[1]) / d;
	Ainv[8] = (A[4] * A[0] - A[3] * A[1]) / d;
}

void mat33ProdVec(const double *A, const double *b, double* r) {
	r[0] = A[0] * b[0] + A[1] * b[1] + A[2] * b[2];
	r[1] = A[3] * b[0] + A[4] * b[1] + A[5] * b[2];
	r[2] = A[6] * b[0] + A[7] * b[1] + A[8] * b[2];
}

void mat33ProdVec(
		const double *A,
		const double *b,
		const double* c,
		double *r,
		double alpha,
		double beta) {
	r[0] = alpha * (A[0] * b[0] + A[1] * b[1] + A[2] * b[2]) + beta * c[0];
	r[1] = alpha * (A[3] * b[0] + A[4] * b[1] + A[5] * b[2]) + beta * c[1];
	r[2] = alpha * (A[6] * b[0] + A[7] * b[1] + A[8] * b[2]) + beta * c[2];
}
void mat33TransProdVec(const double* A, const double* b, double *r) {
	r[0] = A[0] * b[0] + A[3] * b[1] + A[6] * b[2];
	r[1] = A[1] * b[0] + A[4] * b[1] + A[7] * b[2];
	r[2] = A[2] * b[0] + A[5] * b[1] + A[8] * b[2];
}
void mat33AB(const double* A, const double *B, double * C) {
	C[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
	C[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
	C[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];

	C[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
	C[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
	C[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];

	C[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
	C[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
	C[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
}
void mat33ATB(const double* A, const double* B, double *C) {
	C[0] = A[0] * B[0] + A[3] * B[3] + A[6] * B[6];
	C[1] = A[0] * B[1] + A[3] * B[4] + A[6] * B[7];
	C[2] = A[0] * B[2] + A[3] * B[5] + A[6] * B[8];

	C[3] = A[1] * B[0] + A[4] * B[3] + A[7] * B[6];
	C[4] = A[1] * B[1] + A[4] * B[4] + A[7] * B[7];
	C[5] = A[1] * B[2] + A[4] * B[5] + A[7] * B[8];

	C[6] = A[2] * B[0] + A[5] * B[3] + A[8] * B[6];
	C[7] = A[2] * B[1] + A[5] * B[4] + A[8] * B[7];
	C[8] = A[2] * B[2] + A[5] * B[5] + A[8] * B[8];
}
void mat33ABT(const double* A, const double* B, double *C) {
	C[0] = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
	C[1] = A[0] * B[3] + A[1] * B[4] + A[2] * B[5];
	C[2] = A[0] * B[6] + A[1] * B[7] + A[2] * B[8];

	C[3] = A[3] * B[0] + A[4] * B[1] + A[5] * B[2];
	C[4] = A[3] * B[3] + A[4] * B[4] + A[5] * B[5];
	C[5] = A[3] * B[6] + A[4] * B[7] + A[5] * B[8];

	C[6] = A[6] * B[0] + A[7] * B[1] + A[8] * B[2];
	C[7] = A[6] * B[3] + A[7] * B[4] + A[8] * B[5];
	C[8] = A[6] * B[6] + A[7] * B[7] + A[8] * B[8];
}
double mat33Tr(const double* A) {
	return A[0] + A[4] + A[8];
}
double mat22Det(const double* A) {
	return A[0] * A[3] - A[1] * A[2];
}
double mat33Det(const double *A) {
	return A[0] * (A[4] * A[8] - A[5] * A[7])
			- A[1] * (A[3] * A[8] - A[5] * A[6])
			+ A[2] * (A[3] * A[7] - A[4] * A[6]);
}
double mat44Det(const double* A) {
	return A[0] * A[10] * A[15] * A[5] - A[0] * A[11] * A[14] * A[5]
			- A[1] * A[10] * A[15] * A[4] + A[1] * A[11] * A[14] * A[4]
			+ A[0] * A[11] * A[13] * A[6] - A[11] * A[13] * A[2] * A[4]
			- A[0] * A[10] * A[13] * A[7] - A[1] * A[11] * A[12] * A[6]
			+ A[10] * A[13] * A[3] * A[4] + A[11] * A[12] * A[2] * A[5]
			+ A[1] * A[10] * A[12] * A[7] - A[10] * A[12] * A[3] * A[5]
			- A[0] * A[15] * A[6] * A[9] + A[1] * A[15] * A[6] * A[8]
			+ A[15] * A[2] * A[4] * A[9] - A[15] * A[2] * A[5] * A[8]
			+ A[0] * A[14] * A[7] * A[9] - A[1] * A[14] * A[7] * A[8]
			- A[14] * A[3] * A[4] * A[9] + A[14] * A[3] * A[5] * A[8]
			+ A[13] * A[2] * A[7] * A[8] - A[13] * A[3] * A[6] * A[8]
			- A[12] * A[2] * A[7] * A[9] + A[12] * A[3] * A[6] * A[9];
}

double matDet(int m, const double* A) {
	double* tA = new double[m * m];
	int *ipiv = new int[m * m];

	memcpy(tA, A, sizeof(double) * m * m);
	int info;
	dgetrf_(&m, &m, tA, &m, ipiv, &info);
	double s = 1.0;
	for (int i = 0; i < m; i++) {
		s *= tA[i * m + i];
	}
	delete[] ipiv;
	delete[] tA;
	return s;
}
bool containNan(int n, double *a) {
	int i;
	for (i = 0; i < n; ++i)
		if (a[i] != a[i])
			return true;
	return false;
}

void matTrans(int m, int n, const double* Ar, double* Ac) {
	int i, j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			Ac[j * m + i] = Ar[i * n + j];
}
void matCol2Row(int m, int n, const double* Ac, double* Ar) {
	int i, j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			Ar[i * n + j] = Ac[j * m + i];
}
void matRow2Col(int m, int n, const double* Ar, double *Ac) {
	int i, j;
	for (j = 0; j < n; j++)
		for (i = 0; i < m; i++)
			Ac[j * m + i] = Ar[i * n + j];
}
void matFill(int m, int n, double* A, double val) {
	int i, len = m * n;
	for (i = 0; i < len; ++i)
		A[i] = val;
}
void matScale(int m, int n, double* A, double s, double *B) {
	int i, len = m * n;
	for (i = 0; i < len; ++i)
		B[i] = A[i] * s;
}
void matAddScale(int m, int n, double* A, double s) {
	int i, len = m * n;
	for (i = 0; i < len; i++)
		A[i] += s;
}
void matAddScale(int m, int n, double* A, double s, double* B) {
	int i, len = m * n;
	for (i = 0; i < len; i++)
		B[i] = A[i] + s;
}

void matDiff(int m, int n, double* A, double* B, double* R) {
	int i, len = m * n;
	for (i = 0; i < len; ++i)
		R[i] = A[i] - B[i];
}
void matSum(int m, int n, const double * A, const double * B, double* S) {
	int i, len = m * n;
	for (i = 0; i < len; ++i)
		S[i] = A[i] + B[i];
}
void matAvg(int m, int n, const double* A, const double* B, double* S) {
	int i, len = m * n;
	for (i = 0; i < len; ++i)
		S[i] = 0.5 * (A[i] + B[i]);
}
void matAB(
		int Am,
		int An,
		int Bm,
		int Bn,
		const double *A,
		const double *B,
		double * C) {
	if (An != Bm) {
		repErr(
				" The number of columns of A and the number of rows of B must be equal\n");
	}
	double alpha = 1.0;
	double beta = 0.0;
	const char* ch1 = "N";
	const char* ch2 = "N";
	dgemm_(const_cast<char*>(ch1), const_cast<char*>(ch2), &Bn, &Am, &Bm,
			&alpha, const_cast<double*>(B), &Bn, const_cast<double*>(A), &An,
			&beta, C, &Bn);
}
void matATB(
		int Am,
		int An,
		int Bm,
		int Bn,
		const double *A,
		const double *B,
		double *C) {
	if (Am != Bm) {
		repErr(
				"Error: the number of rows of A and the number of rows of B must be equal\n");
	}
	double alpha = 1.0;
	double beta = 0.0;
	const char* ch1 = "N";
	const char* ch2 = "T";
	dgemm_(const_cast<char*>(ch1), const_cast<char*>(ch2), &Bn, &An, &Bm,
			&alpha, const_cast<double*>(B), &Bn, const_cast<double*>(A), &An,
			&beta, C, &Bn);
}
void matABT(
		int Am,
		int An,
		int Bm,
		int Bn,
		const double *A,
		const double *B,
		double *C) {
	if (An != Bn) {
		repErr(
				"Error: the number of columns of A and the number of columns of B must be equal\n");
	}
	double alpha = 1.0;
	double beta = 0.0;
	const char* ch1 = "T";
	const char* ch2 = "N";
	dgemm_(const_cast<char*>(ch1), const_cast<char*>(ch2), &Bm, &Am, &Bn,
			&alpha, const_cast<double*>(B), &Bn, const_cast<double*>(A), &An,
			&beta, C, &Bm);
}
void matAxpy(
		int m,
		int n,
		double alpha,
		const double* A,
		const double* x,
		double beta,
		const double* y,
		double* z) {
	char trans = 'N';
	int incx = 1;
	int incy = 1;
	double* At = new double[m * n];
	matTrans(m, n, A, At);
	if (y) {
		memcpy(z, y, sizeof(double) * n);
	}
	dgemv_(&trans, &m, &n, &alpha, At, &m, const_cast<double*>(x), &incx, &beta,
			z, &incy);
	delete[] At;
}
bool matInv(int n, const double* A, double *invA) {
	double *At = new double[n * n];
	int m = n;
	int lda = n;
	int info;
	int *ipiv = new int[n];
	int lwork = n * 512;
	int i, j;
	double *work = new double[lwork];
	bool res = false;

	assert(At != NULL);
	assert(ipiv != NULL);
	assert(work != NULL);

	/* Transpose A info At like FORTRAN likes */
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			At[i * n + j] = A[j * n + i];

	/* Make calls to FORTRAN routines */
	dgetrf_(&m, &n, At, &lda, ipiv, &info);
	if (info != 0) {
		printf("[matrix_invert] Error[dgetrf]: %d\n", info);
		res = false;
		goto cleanup;
	}

	dgetri_(&n, At, &lda, ipiv, work, &lwork, &info);
	if (info != 0) {
		printf("[matrix_invert] Error[dgetri]: %d\n", info);
		res = false;
		goto cleanup;
	}

	/* Transpose back into Ainv */
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			invA[i * n + j] = At[j * n + i];

	res = true;
	cleanup: delete[] At;
	delete[] ipiv;
	delete[] work;
	return res;
}
bool matInvInplace(int n, double* A) {
	int m = n;
	int lda = n;
	int info;
	int *ipiv = new int[n]; //malloc(sizeof(int) * n);
	int lwork = n * 512;
	double *work = new double[lwork]; //malloc(sizeof(double) * lwork);

	/* Make calls to FORTRAN routines */
	dgetrf_(&m, &n, A, &lda, ipiv, &info);
	if (info != 0)
		goto exit;
	dgetri_(&n, A, &lda, ipiv, work, &lwork, &info);
	if (info != 0)
		goto exit;

	exit: delete[] ipiv;
	delete[] work;
	return info == 0;
}
int dgelsyFor(int m, int n, int nrhs, const double *A, double *b, double *x) {
	if (m < n) {
		repErr("dgelsy_for() : now only works when m (%d) >= n (%d) \n", m, n);
	}

	double* Atmp = new double[m * n];
	double * btmp = new double[m * nrhs];

	int lda = m;
	int ldb = m;

	int* jpvt = new int[n];

	double rcond = -1.0;
	int rank; /* Output */
	int lwork = -1;
	double* work = new double[1];

	int info;

	matRow2Col(m, n, A, Atmp);
	matRow2Col(m, nrhs, b, btmp);

	/* Query to find a good size for the work array */
	dgelsy_(&m, &n, &nrhs, Atmp, &lda, btmp, &ldb, jpvt, &rcond, &rank, work,
			&lwork, &info);

	lwork = (int) work[0];
	/* printf("Work size: %d\n", lwork); */
	//free(work);
	delete[] work;

	//work = (double*) malloc(sizeof(double) * lwork);
	work = new double[lwork];

	/* Make the FORTRAN call */
	dgelsy_(&m, &n, &nrhs, Atmp, &lda, btmp, &ldb, jpvt, &rcond, &rank, work,
			&lwork, &info);

	if (info != 0)
		printf("Error [%d] in call to dgelsy\n", info);

	/* Go from column- to row-major */
	matCol2Row(n, nrhs, btmp, x);

	delete[] Atmp;
	delete[] btmp;
	delete[] work;
	delete[] jpvt;

	return info;
}
//template<class T>
//inline T MAX(T a, T b) {
//	return a < b ? b : a;
//}
//template<class T>
//inline T MIN(T a, T b) {
//	return a < b ? a : b;
//}

int dgesvdFor(
		int m,
		int n,
		const double *A,
		double * U,
		double *S,
		double *VT) {

	double *AT, *V;
	double *UT;

	char jobu = 'a';
	char jobvt = 'a';

	int lda = m;
	int ldu = m;
	int ldvt = n;

	int lwork = 10
			* std::max(3 * std::min(m, n) + std::max(m, n), 5 * std::min(m, n));

	double *work;

	int info;

	/*convert row-major to column-major matrix*/
	AT = new double[m * n]; //(double *) malloc(sizeof(double) * m * n);
	matRow2Col(m, n, A, AT);

	/* Create temporary matrices for output of dgesvd */
	UT = new double[m * m];
	V = new double[n * n];
	work = new double[lwork];

	dgesvd_(&jobu, &jobvt, &m, &n, AT, &lda, S, UT, &ldu, V, &ldvt, work,
			&lwork, &info);

	if (info != 0) {
		warn("dgesvd_() : An error occurred\n");
	}

	/*convert to row-major*/
	matCol2Row(m, m, UT, U);
	matCol2Row(n, n, V, VT);

	delete[] AT;
	delete[] UT;
	delete[] V;
	delete[] work;

	return info;

}
int dgesvdFor(int m, int n, const double *A, double *S, double *VT) {

	double *AT, *V;
	double UT;

	char jobu = 'n';
	char jobvt = 'a';

	int lda = m;
	int ldu = 1;
	int ldvt = n;

	int lwork = 10
			* std::max(3 * std::min(m, n) + std::max(m, n), 5 * std::min(m, n));

	double *work;

	int info;

	/*convert row-major to column-major matrix*/
	AT = new double[m * n];
	matRow2Col(m, n, A, AT);

	/* Create temporary matrices for output of dgesvd */
	V = new double[n * n];
	work = new double[lwork];

	dgesvd_(&jobu, &jobvt, &m, &n, AT, &lda, S, &UT, &ldu, V, &ldvt, work,
			&lwork, &info);

	if (info != 0) {
		warn("dgesvd_() : An error occurred\n");
	}

	/*convert to row-major*/
	matCol2Row(n, n, V, VT);

	delete[] AT;
	delete[] V;
	delete[] work;

	return info;
}

int dgeevFor(int n, const double *A, double *evec, double *eval) {
	char jobvl = 'N'; /* Don't compute left eigenvectors */
	char jobvr = 'V'; /* Do compute right eigenvectors */
	int lda = n;
	double *Atmp = new double[n * n];
	double *wr = new double[n];
	double *wi = new double[n];
	double *vl = NULL;
	int ldvl = 1;
	double *vr = new double[n * n];
	int ldvr = n;
	int lwork;
	double *work, work_query[1];
	int info;

	int i, j, count = 0;

	/* Transpose the matrix for FORTRAN */
	matRow2Col(n, n, A, Atmp);

	/* Query dgeev for the optimal value of lwork */
	lwork = -1;
	dgeev_(&jobvl, &jobvr, &n, Atmp, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
			work_query, &lwork, &info);
	lwork = (int) work_query[0];
	work = new double[lwork];

	/* Make the call to dgeev */
	dgeev_(&jobvl, &jobvr, &n, Atmp, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work,
			&lwork, &info);

	if (info < 0)
		warn("Error in call to dgeev (argument %d was invalid\n", -info);
	else if (info > 0)
		warn("Error: not all eigenvalues have converged\n");

	/* Check that all eigenvalues are real */
	for (i = 0; i < n; i++) {
		if (wi[i] != 0.0) {
			// printf("[dgeev] Eigenvalue has non-zero imaginary part\n");
		} else {
			eval[count] = wr[i];

			for (j = 0; j < n; j++)
				evec[count * n + j] = vr[i * n + j];

			count++;
		}
	}

	/* Clean up */
	delete[] work;
	delete[] Atmp;
	delete[] wr;
	delete[] wi;
	delete[] vr;

	return count;
}
/* QR decomposition*/
int dgeqrFor(int m, int n, const double* A, double* Q, double* R) {
	assert(m > 0 && n > 0);

	int d = m > n ? m : n;
	double* _At = new double[d * d];
	double * _tau = new double[d];
	memset(_tau, 0, sizeof(double) * d);
	double work;
	int lwork = -1;
	int info;

	matTrans(m, n, A, _At);

	//query the lwork value
	dgeqrf_(&m, &n, _At, &m, _tau, &work, &lwork, &info);
	lwork = work;
	double* _work = new double[lwork];

	//QR decomposition
	dgeqrf_(&m, &n, _At, &m, _tau, _work, &lwork, &info);

	if (info < 0)
		warn("Error: '%d-th' argument had an illegal value,", -info);

	//read the upper triangular elements
	memset(R, 0, sizeof(double) * m * n);
	for (int j = 0; j < n; j++) {
		for (int i = 0; i <= j && i < m; i++) {
			R[i * n + j] = _At[j * m + i];
		}
	}
	{
		//compute the Q matrix
		//query the lwork value
		lwork = -1;
		int M = m;
		int N = m;
		dorgqr_(&M, &N, &N, _At, &m, _tau, &work, &lwork, &info);
		lwork = work;
		delete[] _work;
		_work = new double[lwork];
		dorgqr_(&M, &N, &N, _At, &m, _tau, _work, &lwork, &info);
	}
	matTrans(m, m, _At, Q);
	delete[] _At;
	delete[] _tau;
	delete[] _work;
	return info;
}

//int lmdifFor(void(*fcn)(int* , int* , double* , double* , int *) , int m , int n , double *xvec , double tol) {
//	int info;
//	int lwa = m * n + 5 * n + m;
//	int *iwa;
//	double *fvec, *wa;
//
//	if (n > m) {
//		reportError("Error: lmdif called with n > m\n");
//	}
//
//	iwa = (int *) malloc(sizeof(int) * n);
//	fvec = (double *) malloc(sizeof(double) * m);
//	wa = (double *) malloc(sizeof(double) * lwa);
//
//	lmdif1_(fcn, &m, &n, xvec, fvec, &tol, &info, iwa, wa, &lwa);
//
//#if 0
//	switch (info) {
//		case 0:
//		printf("Improper input parameters\n");
//		break;
//		case 1:
//		printf("Sum of squares tolerance reached\n");
//		break;
//		case 2:
//		printf("x is within tolerance\n");
//		break;
//		case 3:
//		printf("Sum of squares and x are within tolerance\n");
//		break;
//		case 4:
//		printf("fvec orthogonal\n");
//		break;
//		case 5:
//		printf("max function calls made\n");
//		break;
//		case 6:
//		printf("tolerance is too small (squares)\n");
//		break;
//		case 7:
//		printf("tolerance is too small (x)\n");
//		break;
//	}
//#endif
//
//	free(iwa);
//	free(fvec);
//	free(wa);
//
//	return info;
//}

//#include "SL_Print.h"
//int main(int argc, char** argv) {
//	int m =10;
//	int n = 7;
//	Mat_d A, Q(m, m), R(m, n);
//	//matRand(m, n, A);
//	matOnes(m, n, A);
//	print(A);
//
//	dgeqrFor(m, n, A, Q, R);
//	print(Q);
//	print(R);
//	return 0;
//}

//#include "SL_Print.h"
//int main(int argc, char** argv) {
//	int m = 7;
//	int n = 5;
//	Mat_d A, U(m, m), S(m, n), VT(n, n);
//	matRand(m, n, A);
//	print(A);
//
//	dgesvdFor(m, n, A, U.data, S.data, VT.data);
//	print( U);
//	print( S);
//	print( VT);
//	return 0;
//}

//#include "SL_Print.h"
//int main(int argc, char** argv) {
//	int m = 7;
//	int n = 5;
//	Mat_d A, x(n, 1), b;
//	matRand(m, n, A);
//	matRand(m, 1, b);
//	print(A);
//	print(b);
//
//	dgelsyFor(m,n,1,A,b,x);
//	print(x);
//	return 0;
//}
