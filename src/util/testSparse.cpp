#include <cstdio>
#include "math/SL_SparseLinearSystem.h"
#include "tools/SL_WriteRead.h"
#include "tools/SL_Print.h"
//void testSolveNormEqn() {
//	int m = 10, n = 10;
//	Mat_d mat;
//	//matRand(m, n, mat);
//	//write(mat,"mat.txt");
//	readMat(mat, "mat.txt");
//	printMat(m, n, mat);
//
//	SparseMat cssMat;
//	dense2Sparse(mat, cssMat);
//	print(cssMat);
//
//	//test
//	printf("m,n:%d,%d\n", m, n);
//
//	double r, x[10];
//	sparseSolveNorm(cssMat, x, &r);
//
//	//test
//	logInfo("residual:%g\n", r);
//	printMat(1, 10, x);
//
//}
void testSolveLinEqn() {
	int m = 10, n = 10;
	Mat_d mat;
	//matRand(m, n, mat);
	//write(mat,"mat.txt");
	readMat(mat, "mat.txt");
	printMat(m, n, mat);

	SparseMat cssMat;
	dense2Sparse(mat, cssMat);
	print(cssMat);

	//test
	printf("m,n:%d,%d\n", m, n);

	double b[10];
	for (int i = 0; i < 10; i++)
		b[i] = 1.0;

	double x[10];
	sparseSolveLin(cssMat, b, x);
	printMat(1, 10, x);
}
void testSolveLinEqn2() {
	Mat_d Ad, b;
	//matRand(m, n, mat);
	//write(mat,"mat.txt");
	readMat(Ad, "/home/tsou/Ad.txt");
	readMat(b, "/home/tsou/b.txt");
	

	print(Ad);
	print(b);
	SparseMat A;
	dense2Sparse(Ad, A);

	print(A);
	const int m = A.m;
	const int n = A.n;
	//test
	printf("m,n:%d,%d\n", m, n);

	double* x = new double[n];
	memset(x, 0, sizeof(double) * n);
	sparseSolveLin(A, b, x);
	printMat(1, n, x);
	delete[] x;
}
void testTri2CSS() {
	Mat_d mat;
	matRand(10, 10, mat);
	print(mat);

	longInt len = mat.rows * mat.cols;
	longInt nnz = 0;
	for (longInt i = 0; i < len; i++) {
		mat.data[i] = mat.data[i] < 0.1 ? 0 : mat.data[i];
		if (mat.data[i] > 0)
			nnz++;
	}
	print(mat);

	longInt* ri = new longInt[nnz];
	longInt* ci = new longInt[nnz];
	double* nzval = new double[nnz];

	int k = 0;
	for (int i = 0; i < mat.rows; i++) {
		for (int j = 0; j < mat.cols; j++) {
			if (mat.data[i * mat.rows + j] > 0) {
				ci[k] = j;
				ri[k] = i;
				nzval[k] = mat.data[i * mat.rows + j];
				k++;
			}
		}
	}

	checkTriplets(10, 10, nnz, ri, ci, nzval);
	SparseMat cssMat;
	triplets2Sparse(10, 10, nnz, ri, ci, nzval, cssMat);

	//test
	print(cssMat);
	sparse2Dense(cssMat, mat);
	print(mat);
	//logInfo("sorted:'%s'\n", isTripletsSorted(10, 10, nnz, ri, ci, nzval) ? "true" : "false");

	delete[] ri;
	delete[] ci;
	delete[] nzval;
}
//int main(int argc, char** argv) {
//
//	//testTri2CSS();
//	//printMat(10,1,eigvecs);
//
//	//testSolveLinEqn();
//	//testSolveLinEqn2();
//	return 0;
//}
