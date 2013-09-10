/*
 * SL_Matrix.h
 *
 *  Created on: 2010-11-6
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_MATRIX_H_
#define SL_MATRIX_H_
#include <algorithm>
#include <string.h>
#include "SL_error.h"
template<typename T>
class MyMat {
public:
	T* data;
	union {
		int rows;
		int m;
	};
	union {
		int cols;
		int n;
	};
public:
	typedef T value_type;
	MyMat(void);
	MyMat(int m, int n);
	MyMat(int m, int n, const T* outData);
	MyMat(const MyMat<T>& other) :
			data(0), rows(0), cols(0) {
		cloneFrom(other);
	}
	MyMat<T>& operator =(const MyMat<T>& other) {
		if (&other != this)
			cloneFrom(other);
		return *this;
	}
	~MyMat(void);

	void clear();
	void resize(int m, int n);
	void fill(T val);
	void copyFrom(const MyMat<T>& other);
	void copyFrom(const T* other);
	void cloneFrom(const MyMat<T>& other);
	void cloneFrom(const T* data, int rows, int cols);
	void copyTo(const MyMat<T>& other, int row_start, int row_end,
			int col_start, int col_end);

	void getRow(int i, MyMat<T>& row);
	void getCol(int j, MyMat<T>& col);
public:
	bool empty() const {
		return rows == 0 || cols == 0;
	}
	T meanSquare();
	T mean();
	T maxValue();
	T maxElement(int &imax, int &jmax);
	T minValue();
	T minElement(int &imin, int &jmin);
	T* get_row(int r) const {
		return data + r * cols;
	}
public:
	typedef T* PTYPE;
	typedef const T* CONST_PTYPE;
	operator PTYPE() {
		return data;
	}
	operator CONST_PTYPE() const {
		return data;
	}
	T* get_ptr() const {
		return data;
	}
	inline T& operator()(int i, int j) {
		return data[i * cols + j];
	}
	inline T operator()(int i, int j) const {
		return data[i * cols + j];
	}
};

/* two useful types*/
typedef MyMat<double> Mat_d;
typedef MyMat<float> Mat_f;
typedef MyMat<int> Mat_i;
typedef MyMat<unsigned int> Mat_ui;
typedef MyMat<char> Mat_c;
typedef MyMat<unsigned char> Mat_uc;

/* implementation */
template<typename T>
MyMat<T>::MyMat() :
		data(0), rows(0), cols(0) {
}
template<typename T>
MyMat<T>::MyMat(int m, int n) :
		data(0), rows(0), cols(0) {
	if (m <= 0 || n <= 0) {
		data = 0;
		rows = 0;
		cols = 0;
	} else {
		rows = m;
		cols = n;
		data = new T[m * n];
	}
}
template<typename T>
MyMat<T>::MyMat(int m, int n, const T* outData) :
		data(0), rows(0), cols(0) {
	if (m <= 0 || n <= 0) {
		data = 0;
		rows = 0;
		cols = 0;
	} else {
		rows = m;
		cols = n;
		data = new T[m * n];
		memcpy(data, outData, sizeof(T) * rows * cols);
	}
}

template<typename T>
MyMat<T>::~MyMat() {
	clear();
}
template<typename T>
void MyMat<T>::clear() {
	if (data)
		delete[] data;
	data = 0;
	rows = 0;
	cols = 0;
}
template<typename T>
void MyMat<T>::resize(int m, int n) {
	if (m == rows && n == cols) {
		return;
	} else {
		clear();
		if (m > 0 && n > 0) {
			rows = m;
			cols = n;
			data = new T[m * n];
		}
	}
}
template<typename T>
void MyMat<T>::fill(T val) {
	std::fill_n(data, rows * cols, val);
}

template<typename T>
void MyMat<T>::copyFrom(const MyMat<T>& other) {
	assert(!other.empty() && m == other.m && n == other.n);
	if (&other == this)
		return;
	if (rows != other.rows || cols != other.cols)
		repErr("Matrix::copy_from() : The size is not consistent.");
	int N = rows * cols;
	memcpy(data, other.data, sizeof(T) * N);
}
template<typename T>
void MyMat<T>::copyFrom(const T* other) {
	int N = rows * cols;
	memcpy(data, other, sizeof(T) * N);
}
template<typename T>
void MyMat<T>::cloneFrom(const MyMat<T>& other) {
	if (&other == this)
		return;
	resize(other.rows, other.cols);
	if (rows > 0 && cols > 0) {
		int N = rows * cols;
		memcpy(data, other.data, sizeof(T) * N);
	}
}
template<typename T>
void MyMat<T>::cloneFrom(const T* other, int rows, int cols) {
	resize(rows, cols);
	if (rows > 0 && cols > 0) {
		int N = rows * cols;
		memcpy(data, other, sizeof(T) * N);
	}
}
template<typename T>
void MyMat<T>::copyTo(const MyMat<T>& other, int row_start, int row_end,
		int col_start, int col_end) {
	assert(
			other.rows*other.cols == (row_end - row_start + 1)*(col_end - col_start + 1));
	int k = 0;
	for (int i = row_start; i <= row_end; i++) {
		for (int j = col_start; j <= col_end; j++) {
			other.data[k++] = data[i * cols + j];
		}
	}
}
template<typename T>
void MyMat<T>::getRow(int i, MyMat<T>& row) {
	assert(i >= 0 && i <= rows);
	row.resize(1, cols);
	for (int j = 0; j < cols; j++)
		row.data[j] = data[i * cols + j];
}
template<typename T>
void MyMat<T>::getCol(int j, MyMat<T>& col) {
	assert(j >= 0 && j <= cols);
	col.resize(rows, 1);
	for (int i = 0; i < rows; i++)
		col.data[i] = data[i * cols + j];
}
template<typename T>
T MyMat<T>::meanSquare() {
	int i = 0, len = rows * cols;
	T s = 0;
	for (i = 0; i < len; ++i)
		s += data[i] * data[i];
	s /= len;
	return s;
}
template<typename T>
T MyMat<T>::mean() {
	int i = 0, len = rows * cols;
	T s = 0;
	for (i = 0; i < len; ++i)
		s += data[i];
	s /= len;
	return s;
}
template<typename T>
T MyMat<T>::maxValue() {
	int imax, jmax;
	return maxElement(imax, jmax);
}
template<typename T>
T MyMat<T>::maxElement(int &imax, int &jmax) {
	if (empty())
		repErr("Matrix::get_max_element() :  Array cannot be empty.");

	int len = rows * cols;
	T* ptr = std::max_element(data, data + len);
	int offset = ptr - data;
	int r = offset / cols;
	int c = offset - r * cols;

	imax = r;
	jmax = c;

	return *ptr;
}
template<typename T>
T MyMat<T>::minValue() {
	int imin, jmin;
	return minElement(imin, jmin);
}
template<typename T>
T MyMat<T>::minElement(int &imin, int &jmin) {
	if (empty())
		repErr("Matrix::get_max_element() : Array cannot be empty.");

	int len = rows * cols;
	T* ptr = std::min_element(data, data + len);
	int offset = ptr - data;
	int r = offset / cols;
	int c = offset - r * cols;

	imin = r;
	jmin = c;

	return *ptr;
}

/* generate useful matrices*/
void matRand(int m, int n, Mat_d& mat);
void matRand(int m, int n, Mat_f& mat);
void matOnes(int m, int n, Mat_d& mat);
void matZeros(int m, int n, Mat_d& mat);
void matZeros(int m, int n, Mat_f& mat);
void matZeros(int m, int n, double* mat);
void matEyes(int n, Mat_d& mat);
void matEyes(int n, double* mat);
#endif /* SL_MATRIX_H_ */
