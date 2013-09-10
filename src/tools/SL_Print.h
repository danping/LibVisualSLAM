/*
 * SL_Helper.h
 *
 *  Created on: 2011-1-6
 *      Author: Danping Zou
 */

#ifndef SL_HELPER_H_
#define SL_HELPER_H_
#include "math/SL_Matrix.h"
#include "matching/SL_Matching.h"
#include "geometry/SL_Point.h"
#include "SL_TypeConversion.h"

void printMat(int m, int n, const double* mat);
void printMat(FILE* fp, int m, int n, const double* mat);
void print(const Mat_d& mat);
void print(FILE* fp, const Mat_d& mat);
void print(FILE* fp, const Mat_i& mat);
void print(const Mat_f& mat);
void print(const Mat_i& mat);
void printMat(FILE* fp, int m, int n, const int* mat);
void print(const Mat_ui& mat);
void print(const Mat_c& mat);
void print(const VecPoint2d& pts);
void print(const Matching& matches);

void printMat(int m, int n, const float* mat);
void printMat(int m, int n, const int * mat);
void printMat(int m, int n, const unsigned int * mat);
void printMat(int m, int n, const char* mat);
void printTransMat(int m, int n, const float* mat);
void printPts(int n, const Point2d* pts);

template<class T>
void releasePtrList(T& lst) {
	typename T::iterator iter = lst.begin();
	while (iter != lst.end()) {
		delete *iter;
		++iter;
	}
	lst.clear();
}
void VecPoint2D2Matrix(const VecPoint2d& pts, Mat_d& mat);
#endif /* SL_HELPER_H_ */
