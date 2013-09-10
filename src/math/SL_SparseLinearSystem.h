/*
 * SL_SparseLinearSystem.h
 *
 *  Created on: 2011-6-15
 *      Author: Danping Zou
 */

#ifndef SL_SPARSELINEARSYSTEM_H_
#define SL_SPARSELINEARSYSTEM_H_
#include "SL_SparseMat.h"

/**
 * solve Ax = b; A: mxn matrix  
 */
bool sparseSolveLin(SparseMat& A, double b[], double x[]);

/**
 * solve Ax = b using direct method 
 */
bool sparseSolveLin(Triplets& A, double b[], double x[]);

/**
 * solve [A B][x y]' = b
 */
bool sparseSolveLin(Triplets& A, Triplets& B, double b[], double x[], double y[]);
bool sparseSolveLin(SparseMat& A, SparseMat& B, double b[], double x[], double y[]);

/**
 * B = A*B
 */
void sparseMatMul(const SparseMat& A, const SparseMat& B, SparseMat& C);


/**
 *solve Ax = b using conjugate gradient 
 */
//bool sparseSolveLinCG(Triplets& A, double b[], double x[]);

#endif /* SL_SPARSELINEARSYSTEM_H_ */
