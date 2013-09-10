/*
 * SL_TypeConversion.cpp
 *
 *  Created on: 2010-11-9
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#include "SL_TypeConversion.h"

void matDouble2Int(const Mat_d& md, Mat_i& mi) {
	int rows = md.rows;
	int cols = md.cols;
	if (mi.rows != rows || mi.cols != cols)
		mi.resize(rows, cols);

	int i, len = rows * cols;
	for (i = 0; i < len; ++i) {
		mi.data[i] = static_cast<int>(md.data[i] + 0.5);
	}
}
void matDouble2Float(const Mat_d& md, Mat_f& mf) {
	int rows = md.rows;
	int cols = md.cols;
	if (mf.rows != rows || mf.cols != cols)
		mf.resize(rows, cols);

	int i, len = rows * cols;
	for (i = 0; i < len; ++i) {
		mf.data[i] = static_cast<float>(md.data[i]);
	}
}
void matFloat2Double(const Mat_f& m0, Mat_d& m1) {
	int rows = m0.rows;
	int cols = m0.cols;
	if (m1.rows != rows || m1.cols != cols)
		m1.resize(rows, cols);

	int i, len = rows * cols;
	for (i = 0; i < len; ++i) {
		m1.data[i] = static_cast<double>(m0.data[i]);
	}
}
void matFloat2Int(const Mat_f& mf, Mat_i& mi) {
	int rows = mf.rows;
	int cols = mf.cols;
	if (mi.rows != rows || mi.cols != cols)
		mi.resize(rows, cols);

	int i, len = rows * cols;
	for (i = 0; i < len; ++i) {
		mi.data[i] = static_cast<int>(mf.data[i] + 0.5);
	}
}
void VecPoint2D2Matrix(const VecPoint2d& pts, Mat_d& mat) {
	mat.resize((int) pts.size(), 2);
	for (size_t i = 0; i < pts.size(); ++i) {
		mat.data[2 * i] = pts[i].x;
		mat.data[2 * i + 1] = pts[i].y;
	}
}