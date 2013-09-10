/*
 * SL_Math.h
 *
 *  Created on: Dec 22, 2011
 *      Author: Danping Zou
 */

#ifndef SL_MATH_H_
#define SL_MATH_H_
#include "SL_Matrix.h"
#include "SL_LinAlgWarper.h"
#include "SL_SparseMat.h"
#include "SL_SparseLinearSystem.h"
#include "SL_ProbFuncs.h"
//vectorization operations
template<class T, class U>
void div(T iter1, T iter2, U val) {
	for (T iter = iter1; iter != iter2; ++iter)
		*iter /= val;
}
template<class T, class U>
void mult(T iter1, T iter2, U val) {
	for (T iter = iter1; iter != iter2; ++iter)
		*iter *= val;
}
template<class T,class U>
void add(T iter1, T iter2, U val) {
	for (T iter = iter1; iter != iter2; ++iter)
		*iter += val;
}
template<class T,class U>
void sub(T iter1, T iter2, U val) {
	for (T iter = iter1; iter != iter2; ++iter)
		*iter -= val;
}

#endif /* SL_MATH_H_ */
