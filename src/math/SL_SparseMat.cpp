/*
 * SL_SparseMat.cpp
 *
 *  Created on: 2011-6-14
 *      Author: Danping Zou
 */

#include "SL_SparseMat.h"
#include "SL_error.h"
#include <cstdio>
/**
 * convert triplets into compressed column storage
 * 
 * inputs:
 * m,n : dimension fo matrix
 * nnz : number of nonzeros
 * r,c : indicators to row and column indices
 * nzval: nonzero values
 * 
 * outputs:
 * rI : row indices
 * cp : column pointers
 */

bool checkTriplets(const longInt m, const longInt n, const longInt nnz, const longInt r[], const longInt c[], const double nzval[]) {
	for (longInt i = 0; i < nnz; i++) {
		if (r[i] < 0 || r[i] >= m || c[i] < 0 || c[i] >= n) {
			repErr("invavlid element (%l,%l) for a matrix %lx%l\n", r[i], c[i], m, n);
			return false;
		}
	}
	return true;
}
bool isTripletsSorted(const longInt m, const longInt n, const longInt nnz, const longInt r[], const longInt c[], const double nzval[]) {
	if (nnz == 0)
		return true;
	longInt preInd = c[0] * m + r[0];
	for (longInt i = 1; i < nnz; i++) {
		longInt ind = c[i] * m + r[i];
		if (ind < preInd) {
			return false;
		}
		preInd = ind;
	}
	return true;
}
/**
 * convert triplets to compressed column storage (CSS)
 */
void triplets2Sparse(const longInt m,
		const longInt n,
		const longInt nnz,
		const longInt r[],
		const longInt c[],
		const double nzval[],
		double *x,
		longInt rI[],
		longInt cp[]) {

	if (!checkTriplets(m, n, nnz, r, c, nzval)) {
		fprintf(stderr, "tri2css - invalid sparse matrix!\n");
		return;
	}
	longInt* nc = new longInt[n]; //number of elements in a column
	if (!isTripletsSorted(m, n, nnz, r, c, nzval)) {
		//sort triplets in column major order,
		// * Softwares like matlab, arpack++ require the row indices in each column to appear in ascending order, 
		// * but CSparse does not have such restriction
		double * w = new double[nnz];
		for (longInt k = 0; k < nnz; k++) {
			w[k] = c[k] * double(m) + r[k];
		}
		longInt* ind = new longInt[nnz];
		sortWithInd(nnz, w, ind);

		memset(nc, 0, sizeof(longInt) * n);
		for (longInt k = 0; k < nnz; k++) {
			assert(c[k] < n);
			nc[c[k]]++;
		}
		cumSum(n, nc, cp + 1);
		cp[0] = 0;

		for (longInt k = 0; k < nnz; k++) {
			x[k] = nzval[ind[k]];
			rI[k] = r[ind[k]];
		}
		delete[] ind;
		delete[] w;
	} else {
		memset(nc, 0, sizeof(longInt) * n);
		for (longInt k = 0; k < nnz; k++) {
			assert(c[k] < n);
			nc[c[k]]++;
		}

		cumSum(n, nc, cp + 1);
		cp[0] = 0;

		for (longInt k = 0; k < nnz; k++) {
			x[k] = nzval[k];
			rI[k] = r[k];
		}
	}
	delete[] nc;
}

void triplets2Dense(Triplets& trips, Mat_d& mat) {
	mat.resize(trips.m, trips.n);
	mat.fill(0);

	for (longInt i = 0; i < trips.nnz; i++) {
		longInt r = trips.ri[i];
		longInt c = trips.ci[i];
		mat.data[r * mat.n + c] = trips.nzval[i];
	}
}

void sparse2Triplets(const SparseMat& spMat, Triplets& trips) {
	longInt nnz = 0;
	for (longInt c = 0; c < spMat.n; c++) {
		for (longInt n = spMat.p[c]; n != spMat.p[c + 1]; n++)
			nnz++;
	}
	trips.reserve(spMat.m, spMat.n, nnz);
	for (longInt c = 0; c < spMat.n; c++) {
		for (longInt n = spMat.p[c]; n != spMat.p[c + 1]; n++) {
			longInt r = spMat.i[n];
			trips.add(r, c, spMat.x[n]);
		}
	}
}

void dense2Sparse(const Mat_d& mat, SparseMat& cssMat) {
	longInt k = 0;
	longInt len = mat.rows * mat.cols;
	for (longInt i = 0; i < len; i++)
		k += mat.data[i] != 0 ? 1 : 0;

	longInt* r = new longInt[k];
	longInt* c = new longInt[k];
	double* val = new double[k];

	k = 0;
	for (int j = 0; j < mat.cols; j++) {
		for (int i = 0; i < mat.rows; i++) {
			if (mat.data[i * mat.cols + j] != 0) {
				r[k] = (size_t)i;
				c[k] = (size_t)j;
				val[k] = mat.data[i * mat.cols + j];
				k++;
			}
		}
	}
	cssMat.reserve(mat.rows, mat.cols, k);
	triplets2Sparse(mat.rows, mat.cols, k, r, c, val, cssMat);

	delete[] r;
	delete[] c;
	delete[] val;
}
void sparse2Dense(const SparseMat& cssMat, Mat_d& mat) {
	mat.resize(cssMat.m, cssMat.n);
	mat.fill(0);

	for (longInt c = 0; c < cssMat.n; c++) {
		longInt n_s = cssMat.p[c];
		longInt n_e = cssMat.p[c + 1];

		for (longInt n = n_s; n < n_e; n++) {
			longInt r = cssMat.i[n];
			double val = cssMat.x[n];
			mat.data[r * mat.cols + c] = val;
		}
	}
}
void printCSS(const longInt m, const longInt n, const longInt nnz, double x[], longInt rI[], longInt cp[]) {
	printf("%lldx%lld\n", m, n);
	printf("number of nonzeros:%lld\n", nnz);
	printf("row indices:");
	for (longInt k = 0; k < nnz; k++) {
		printf("%lld ", rI[k]);
	}
	printf("\n");
	printf("non zero elements:");
	for (longInt k = 0; k < nnz; k++) {
		printf("%lf ", x[k]);
	}
	printf("\n");
	printf("column pointers:");
	for (longInt k = 0; k <= n; k++) {
		printf("%lld ", cp[k]);
	}
	printf("\n");
}

void sparseSubMat(const SparseMat& A, longInt r1, longInt r2, longInt c1, longInt c2, SparseMat& B) {
	assert(r1 >= 0 && r2 >= r1 && r2 < A.m);
	assert(c1 >= 0 && c2 >= c1 && c2 < A.n);

	Triplets tA, tB;
	sparse2Triplets(A, tA);
	tripletsSubMat(tA, r1, r2, c1, c2, tB);
	triplets2Sparse(tB, B);
}

void tripletsSubMat(const Triplets& tA, longInt r1, longInt r2, longInt c1, longInt c2, Triplets& tB) {
	assert(r1 >= 0 && r2 >= r1 && r2 < tA.m);
	assert(c1 >= 0 && c2 >= c1 && c2 < tA.n);

	longInt bm = r2 - r1 + 1;
	longInt bn = c2 - c1 + 1;

	//count the number of non-zeros in tA
	longInt nnz = 0;
	for (longInt i = 0; i < tA.nnz; i++) {
		if (tA.ri[i] >= r1 && tA.ri[i] <= r2 && tA.ci[i] >= c1 && tA.ci[i] <= c2)
			nnz++;
	}

	tB.reserve(bm, bn, nnz == 0 ? 1 : nnz);

	//set the triplets for the sub matrix
	for (longInt i = 0; i < tA.nnz; i++) {
		longInt r = tA.ri[i];
		longInt c = tA.ci[i];
		double val = tA.nzval[i];
		if (r >= r1 && r <= r2 && c >= c1 && c <= c2) {
			tB.add(r - r1, c - c1, val);
		}
	}
}
void sparseSplitCol(const SparseMat& T, longInt c0, SparseMat& A, SparseMat& B) {
	assert(c0 >= 0 && c0 <= T.n);
	Triplets tT, tA, tB;
	sparse2Triplets(T, tT);
	tripletsSplitCol(tT, c0, tA, tB);
	triplets2Sparse(tA, A);
	triplets2Sparse(tB, B);
}
void tripletsSplitCol(const Triplets& tT, longInt c0, Triplets& tA, Triplets& tB) {
	assert(c0 >= 0 && c0 <= tT.n);

	longInt am = tT.m;
	longInt an = c0;

	longInt bm = tT.m;
	longInt bn = tT.n - c0;

	tA.reserve(am, an, tT.nnz);
	tB.reserve(bm, bn, tT.nnz);
	for (longInt i = 0; i < tT.nnz; i++) {
		longInt r = tT.ri[i];
		longInt c = tT.ci[i];
		double val = tT.nzval[i];
		if (c < c0)
			tA.add(r, c, val);
		else
			tB.add(r, c - c0, val);
	}
}
void sparseSplitCol(const SparseMat& T, longInt c0, SparseMat& A, bool left) {
	assert(c0 >= 0 && c0 <= T.n);
	Triplets tT, tA;
	sparse2Triplets(T, tT);
	tripletsSplitCol(tT, c0, tA, left);
	triplets2Sparse(tA, A);
}
void tripletsSplitCol(const Triplets& tT, longInt c0, Triplets& tA, bool left) {
	assert(c0 >= 0 && c0 <= tT.n);
	if (left) {
		longInt am = tT.m;
		longInt an = c0;

		tA.reserve(am, an, tT.nnz);
		for (longInt i = 0; i < tT.nnz; i++) {
			longInt r = tT.ri[i];
			longInt c = tT.ci[i];
			double val = tT.nzval[i];
			if (c < c0)
				tA.add(r, c, val);
		}
	} else {
		longInt am = tT.m;
		longInt an = tT.n - c0;

		tA.reserve(am, an, tT.nnz);
		for (longInt i = 0; i < tT.nnz; i++) {
			longInt r = tT.ri[i];
			longInt c = tT.ci[i];
			double val = tT.nzval[i];
			if (c >= c0)
				tA.add(r, c - c0, val);
		}
	}
}
void sparseSplitRow(const SparseMat& T, longInt r0, SparseMat& A, bool upper) {
	assert(r0 >= 0 && r0 <= T.m);
	Triplets tT, tA;
	sparse2Triplets(T, tT);
	tripletsSplitRow(tT, r0, tA, upper);
	triplets2Sparse(tA, A);
}
void tripletsSplitRow(const Triplets& tT, longInt r0, Triplets& tA, bool upper) {
	assert(r0 >= 0 && r0 <= tT.m);
	if (upper) {
		longInt am = r0;
		longInt an = tT.n;

		tA.reserve(am, an, tT.nnz);
		for (longInt i = 0; i < tT.nnz; i++) {
			longInt r = tT.ri[i];
			longInt c = tT.ci[i];
			double val = tT.nzval[i];

			if (r < r0) {
				tA.add(r, c, val);
			}
		}
	} else {
		longInt am = tT.m - r0;
		longInt an = tT.n;

		tA.reserve(am, an, tT.nnz);
		for (longInt i = 0; i < tT.nnz; i++) {
			longInt r = tT.ri[i];
			longInt c = tT.ci[i];
			double val = tT.nzval[i];

			if (r >= r0) {
				tA.add(r - r0, c, val);
			}
		}
	}
}

void write(const Triplets& trips, const char* fmtstr, ...) {
	using namespace std;
	char buf[1024];
	GET_FMT_STR(fmtstr, buf)

	ofstream file(buf);
	if (!file)
		repErr("cannot open '%s' to write.", buf);
	
	trips.save(file);
	
	logInfo("write '%s' OK!\n", buf);

}
void read(Triplets& trips, const char* fmtstr, ...) {
	using namespace std;
	char buf[1024];
	GET_FMT_STR(fmtstr, buf)

	ifstream file(buf);
	if (!file)
		repErr("cannot open '%s' to read.", buf);
	
	trips.load(file);
	
	logInfo("read '%s' OK!\n", buf);
}

//#include "SL_Print.h"
//int main(int argc, char** argv) {
//	Mat_d A;
//	matRand(3, 10, A);
//	for (int i = 0; i < A.m * A.n; i++)
//		A.data[i] = A.data[i] < 0.3 ? 0 : A.data[i];
//
//	SparseMat sA;
//	dense2Sparse(A, sA);
//
//	SparseMat sB;
//	sparseSubMat(sA, 0, 2, 0, 3, sB);
//	Mat_d B;
//	sparse2Dense(sB, B);
//
//	print(A);
//	print(B);
//	
//	return 0;
//}
//#include "SL_Print.h"
//int main(int argc, char** argv) {
//	Mat_d A;
//	matRand(10, 5, A);
//	for (int i = 0; i < A.m * A.n; i++)
//		A.data[i] = A.data[i] < 0.3 ? 0 : A.data[i];
//
//	print(A);
//
//	SparseMat sA, sA1, sA2;
//	dense2Sparse(A, sA);
//	//sparseSplitCol(sA, 0, sA1, sA2);
//	sparseSplitRow(sA, 0, sA2, true);
//
//	Mat_d A1, A2;
//	sparse2Dense(sA1, A1);
//	sparse2Dense(sA2, A2);
//	print(A1);
//	print(A2);
//	return 0;
//}
