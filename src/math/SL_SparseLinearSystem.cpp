/*
 * SL_SparseLinearSystem.cpp
 *
 *  Created on: 2011-6-15
 *      Author: Danping Zou
 */

#include "SL_SparseLinearSystem.h"
#include "util/SL_Utility.h"
#include "extern/csparse/cs.h"
#include <cassert>
bool sparseSolveLin(SparseMat& A, double b[], double x[]) {
	double* bb = new double[A.m];
	for (longInt i = 0; i < A.m; i++) {
		bb[i] = b[i];
	}
	int ok = cs_qrsol(0, (const cs*) &A, bb);
	for (longInt i = 0; i < A.n; i++) {
		x[i] = bb[i];
	}
	if (!ok)
		fprintf(stderr, "solveLinEqn failed!\n");

	delete[] bb;
	return ok;
}
bool sparseSolveLin(Triplets& T, double b[], double x[]) {
	assert(T.m >= T.n);

	cs* t = cs_spalloc(T.m, T.n, T.nmax, T.nnz, 1);
	for (longInt i = 0; i < T.nnz; i++) {
		cs_entry(t, T.ri[i], T.ci[i], T.nzval[i]);
	}
	cs* A = cs_compress(t);

	double* bb = new double[A->m];
	for (int i = 0; i < A->m; i++) {
		bb[i] = b[i];
	}
	int ok = cs_qrsol(0, A, bb);
	for (int i = 0; i < A->n; i++) {
		x[i] = bb[i];
	}
	if (!ok)
		fprintf(stderr, "solveLinEqn failed!\n");

	cs_spfree(t);
	cs_spfree(A);

	delete[] bb;
	return ok;
}
bool sparseSolveLin(Triplets& A, Triplets& B, double b[], double x[],
		double y[]) {
	assert(A.m == B.m);

	int m = A.m;
	int n1 = A.n;
	int n2 = B.n;
	int n = n1 + n2;

	assert(m >= n);
	cs* t = cs_spalloc(m, n, A.nmax + B.nmax, 1, 1);

	for (longInt i = 0; i < A.nnz; i++) {
		cs_entry(t, A.ri[i], A.ci[i], A.nzval[i]);
	}

	for (longInt i = 0; i < B.nnz; i++) {
		cs_entry(t, B.ri[i], B.ci[i] + n1, B.nzval[i]);
	}
	cs* C = cs_compress(t);

	double* bb = new double[m];
	for (int i = 0; i < m; i++) {
		bb[i] = b[i];
	}
	int ok = cs_qrsol(0, C, bb);
	for (int i = 0; i < n1; i++) {
		x[i] = bb[i];
	}
	for (int i = n1; i < n; i++) {
		y[i - n1] = bb[i];
	}
	if (!ok)
		fprintf(stderr, "solveLinEqn failed!\n");

	cs_spfree(t);
	cs_spfree(C);
	delete[] bb;
	return ok;
}
bool sparseSolveLin(SparseMat& A, SparseMat& B, double b[], double x[],
		double y[]) {
	assert(A.m == B.m);
	Triplets tA, tB;
	sparse2Triplets(A, tA);
	sparse2Triplets(B, tB);
	return sparseSolveLin(tA, tB, b, x, y);
}
void sparseMatMul(const SparseMat& A, const SparseMat& B, SparseMat& C) {
	assert(A.n == B.m);
	cs* sA = (cs*) &A;
	cs* sB = (cs*) &B;
	cs* sC = cs_multiply(sA, sB);
	C.reserve(sC->m, sC->n, sC->nzmax);
	memcpy(C.p, sC->p, sizeof(int) * (sC->n + 1));
	memcpy(C.i, sC->i, sizeof(int) * sC->nzmax);
	memcpy(C.x, sC->x, sizeof(double) * sC->nzmax);
	cs_spfree(sC);
}


//#include "SL_Print.h"
//int main(int argc, char** argv) {
//	Mat_d A, B, b;
//	matRand(8, 3, A);
//	matRand(3, 8, B);
//
//	SparseMat sA, sB, sC;
//	dense2Sparse(A, sA);
//	dense2Sparse(B, sB);
//	sparseMatMul(sA,sB, sC);
//	
//	Mat_d C;
//	sparse2Dense(sC,C);
//	
//	print(A);
//	print(B);
//	print(C);
//	
//	return 0;
//}

//#include "SL_Print.h"
//int main(int argc, char** argv) {
//	Mat_d A, B, b;
//	matRand(8, 3, A);
//	matRand(8, 3, B);
//	matRand(8, 1, b);
//
//	SparseMat sA, sB;
//	dense2Sparse(A, sA);
//	dense2Sparse(B, sB);
//
//	Mat_d x(A.cols, 1), y(B.cols, 1);
//	sparseSolveLin(sA, sB, b, x.data, y.data);
//	
//	print(A);
//	print(B);
//	print(b);
//	
//	print(x);
//	print(y);
//	
//	return 0;
//}

//#include "SL_Print.h"
//int main(int argc, const char** argv) {
//	Mat_d A;
//	matRand(6, 6, A);
//	for (int i = 0; i < A.m * A.n; i++) {
//		A.data[i] = A.data[i] < 0.3 ? 0 : A.data[i];
//	}
//	print(A);
//
//	SparseMat sA, sB;
//	Triplets tA;
//	dense2Sparse(A, sA);
//	dense2Sparse(A, sB);
//	sparse2Triplets(sA, tA);
//	triplets2Dense(tA, A);
//	print(A);
//	cs* sC = cs_multiply((cs*) &sA, (cs*) &sA);
//	Mat_d C;
//	sparse2Dense(*(SparseMat*) sC, C);
//	
//	cs_spfree(sC);
//	print(C);
//	return 0;
//}
