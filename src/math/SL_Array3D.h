/*
 * SL_Arrary3D.h
 *
 *  Created on: Jul 14, 2012
 *      Author: Danping Zou
 */

#ifndef SL_ARRARY3D_H_
#define SL_ARRARY3D_H_

#include <algorithm>

template<class T>
class MyArr3 {
public:
	T* data;
	union {
		int nch;
		int c;
	};
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
	MyArr3(void) :
			data(0), c(0), m(0), n(0) {

	}
	MyArr3(int _c, int _m, int _n) {
		if (_c <= 0 || _m <= 0 || _n <= 0) {
			data = 0;
			c = 0;
			m = 0;
			n = 0;
		} else {
			c = _c;
			m = _m;
			n = _n;
			data = new T[c * m * n];
		}
	}

	MyArr3(int _c, int _m, int _n, const T* externData) {
		if (_c <= 0 || _m <= 0 || _n <= 0) {
			c = 0;
			data = 0;
			m = 0;
			n = 0;
		} else {
			c = _c;
			m = _m;
			n = _n;
			data = new T[c * m * n];
			memcpy(data, externData, sizeof(T) * c * m * n);
		}
	}

	~MyArr3() {
		clear();
	}

	void clear() {
		if (data)
			delete[] data;
		data = 0;
		c = m = n = 0;
	}

	void resize(int _c, int _m, int _n) {
		if (_c == c && _m == m && _n == n)
			return;
		clear();
		if (_c > 0 && _m > 0 && _n > 0) {
			c = _c;
			m = _m;
			n = _n;
			data = new T[c * m * n];
		}
	}

	void fill(T val) {
		std::fill_n(data, c * m * n, val);
	}

	void cloneFrom(const MyArr3<T>& other) {
		if (&other == this)
			return;

		resize(other.c, other.m, other.n);
		if (c > 0 && m > 0 && n > 0) {
			int N = c * m * n;
			std::copy(other.data, other.data + N, data);
		}
	}
public:
	bool empty() const {
		return c == 0 || m == 0 || n == 0;
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
	T* getPtr() const {
		return data;
	}
	inline T& operator()(int c, int i, int j) {
		return data[c * m * n + i * n + j];
	}
	inline T& operator()(int c, int i, int j) const {
		return data[c * m * n + i * n + j];
	}
	inline T* getPtrChannel(int c) const {
		return data + c * m * n;
	}
	inline T* getPtrRow(int c, int i) const {
		return data + c * m * n + i * n;
	}
};

/*useful types*/
typedef MyArr3<double> Arr3d;
typedef MyArr3<float> Arr3f;
typedef MyArr3<int> Arr3i;
typedef MyArr3<unsigned int> Arr3ui;
typedef MyArr3<char> Arr3c;
typedef MyArr3<unsigned char> Arr3uc;

void arr3Rand(int c, int m, int n, Arr3d& arr);
#endif /* SL_ARRARY3D_H_ */
