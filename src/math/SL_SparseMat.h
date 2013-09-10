/*
 * SL_SparseMat.h
 *
 *  Created on: 2011-6-14
 *      Author: Danping Zou
 */

#ifndef SL_SPARSEMAT_H_
#define SL_SPARSEMAT_H_
#include "SL_error.h"
#include <cstring>
#include <cstddef>
#include <iostream>
#include <fstream>
#include "util/SL_Utility.h"
#include "math/SL_Matrix.h"
typedef long long longInt;

class SparseMat {
public:
	longInt nzmax; /*maximum number of entries*/
	longInt m; /* mxn matrix*/
	longInt n;longInt* p; /*column indices, size n */
	longInt* i; /*row indices, size nzmax*/
	double* x; /*non-zero values*/
	longInt nz;/* # of entries in triplet matrix, always '-1' for compatibility with CSparse matrix*/

public:
	SparseMat() :
			nzmax(0), m(-1), n(-1), p(0), i(0), x(0), nz(-1) {
	}
	~SparseMat() {
		clear();
	}
	void clear() {
		if (x)
			delete[] x;
		if (i)
			delete[] i;
		if (p)
			delete[] p;
		nzmax = 0;
		m = n = 0;
		i = p = 0;
		x = 0;
		nz = -1;
	}
	void reserve(longInt rows, longInt cols, longInt num) {
		clear();
		nzmax = num;
		m = rows;
		n = cols;
		p = new longInt[n + 1];i
		= new longInt[num];
		x = new double[num];
	}
};

class Triplets {
public:
	longInt nmax;longInt m;longInt n;longInt* ri;longInt* ci;
	double* nzval;longInt nnz;
public:
	Triplets() :
			nmax(0), m(-1), n(-1), ri(0), ci(0), nzval(0), nnz(0) {
	}
	~Triplets() {
		clear();
	}
	void clear() {
		if (nzval)
			delete[] nzval;
		if (ri)
			delete[] ri;
		if (ci)
			delete[] ci;

		m = n = -1;
		nnz = 0;
		nzval = 0;
		ri = ci = 0;
		nmax = 0;
	}
	void reserve(longInt rows, longInt cols, longInt num) {
		clear();
		m = rows;
		n = cols;
		nmax = num;
		nzval = new double[num];
		ri = new longInt[num];
		ci = new longInt[num];
	}
	void add(longInt r, longInt c, double val) {
		if (nnz >= nmax) {
			repErr("Triplets::add - the maximun number of elements reached!");
			return;
		}
		ri[nnz] = r;
		ci[nnz] = c;
		nzval[nnz] = val;
		nnz++;
	}
public:
	//debug
	void print() {
		for (longInt i = 0; i < nnz; i++) {
			std::cout << "[" << i << "]: " << ri[i] << " " << ci[i] << " " << nzval[i] << std::endl;
		}
	}
	void save(std::ofstream& fp) const {
		using namespace std;
		fp << nmax << endl;
		fp << nnz << endl;
		fp << m << " " << n << endl;
		for (longInt i = 0; i < nnz; i++) {
			fp << ri[i] << " " << ci[i] << " " << nzval[i] << endl;
		}
	}
	void load(std::ifstream& fp) {
		using namespace std;
		clear();

		int rows, cols, num, maxnum;
		fp >> maxnum >> num >> rows >> cols;

		assert(maxnum >= num && rows > 0 && cols > 0 && num > 0);
		reserve(rows, cols, maxnum);

		for (longInt i = 0; i < num; ++i) {
			longInt r, c;
			double val;
			fp >> r >> c >> val;
			add(r, c, val);
		}
	}
};
/**
 * check whether these triplets is a valid form for sparse matrix
 */
bool checkTriplets(const longInt m, const longInt n, const longInt nnz, const longInt r[], const longInt c[], const double nzval[]);
/**
 * check whether these triplets are sorted in column major order 
 */
bool isTripletsSorted(const longInt m, const longInt n, const longInt nnz, const longInt r[], const longInt c[], const double nzval[]);
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

void triplets2Sparse(const longInt m,
		const longInt n,
		const longInt nnz,
		const longInt r[],
		const longInt c[],
		const double nzval[],
		double *x,
		longInt rI[],
		longInt cp[]);

inline void triplets2Sparse(const longInt m,
		const longInt n,
		const longInt nnz,
		const longInt r[],
		const longInt c[],
		const double nzval[],
		SparseMat& cssMat) {
	cssMat.reserve(m, n, nnz);
	triplets2Sparse(m, n, nnz, r, c, nzval, cssMat.x, cssMat.i, cssMat.p);
}
inline void triplets2Sparse(const Triplets& trips, SparseMat& cssMat) {
	triplets2Sparse(trips.m, trips.n, trips.nnz, trips.ri, trips.ci, trips.nzval, cssMat);
}
void triplets2Dense(Triplets& trips, Mat_d& mat);

void sparse2Triplets(const SparseMat& spMat, Triplets& trips);
/**
 * convert dense matrix into sparse matrix
 */
void dense2Sparse(const Mat_d& mat, SparseMat& cssMat);
/**
 * convert sparse matrix into dense one
 */
void sparse2Dense(const SparseMat& cssMat, Mat_d& mat);

void printCSS(longInt m, longInt n, longInt nnz, double x[], longInt rI[], longInt cp[]);

inline void print(const SparseMat& cssMat) {
	Triplets trips;
	sparse2Triplets(cssMat, trips);
	trips.print();
}

/* select sub matrix*/
void sparseSubMat(const SparseMat& A, longInt r1, longInt r2, longInt c1, longInt c2, SparseMat& B);
void tripletsSubMat(const Triplets& tA, longInt r1, longInt r2, longInt c1, longInt c2, Triplets& tB);
/* split into two column matrices*/
void sparseSplitCol(const SparseMat& T, longInt c0, SparseMat& A, SparseMat& B);
void tripletsSplitCol(const Triplets& tT, longInt c0, Triplets& tA, Triplets& tB);

/* split out one column matrix (left or right)*/
void sparseSplitCol(const SparseMat& T, longInt c0, SparseMat& A, bool left);
void tripletsSplitCol(const Triplets& tT, longInt c0, Triplets& tA, bool left);

/* split out one row matrix (upper or lower)*/
void sparseSplitRow(const SparseMat& T, longInt r0, SparseMat& A, bool upper);
void tripletsSplitRow(const Triplets& tT, longInt r0, Triplets& tA, bool upper);

/* I/O for triplets*/
void write(const Triplets& trips, const char* fmtstr, ...);
void read(Triplets& trips, const char* fmtstr, ...);
#endif /* SL_SPARSEMAT_H_ */
