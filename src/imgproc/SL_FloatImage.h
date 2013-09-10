/*
 * SL_FloatImage.h
 *
 *  Created on: Dec 22, 2011
 *      Author: Danping Zou
 */

#ifndef SL_FLOATIMAGE_H_
#define SL_FLOATIMAGE_H_
#include "math/SL_Matrix.h"
#include <algorithm>
#include <cassert>
#include <limits>
#include <cmath>

template<class T>
class FloatImgG {
public:
	union {
		int n;
		int w;
		int cols;
	};
	union {
		int m;
		int h;
		int rows;
	};
	T* data;
public:
	FloatImgG() :
			w(0), h(0), data(0) {
	}
	FloatImgG(int iw, int ih) {
		assert(iw > 0 && ih > 0);
		data = new T[iw * ih];
		w = iw;
		h = ih;
	}
	~FloatImgG() {
		clear();
	}

public:
	void clear() {
		if (data)
			delete[] data;
		data = 0;
		w = 0;
		h = 0;
	}
	void fill(T val) {
		std::fill_n(data, w * h, val);
	}
	void resize(int iw, int ih) {
		clear();
		data = new T[iw * ih];
		w = iw;
		h = ih;
	}
	bool empty() const {
		return w == 0 || h == 0;
	}

public:
	typedef T* PFLOAT;
	typedef const T* CONST_PFLOAT;
	operator PFLOAT() {
		return data;
	}
	operator CONST_PFLOAT() const {
		return data;
	}
	T* getPtr() const {
		return data;
	}
	T& operator()(int x, int y) const {
		return data[y * w + x];
	}
	int getMemLen() const {
		return w * h * sizeof(T);
	}

	void getAvgStd(int startx, int endx, int starty, int endy, T& std,
			T& avg) const {
		assert( startx >= 0 && starty >= 0);
		avg = 0;
		int num = 0;
		for (int i = starty; i <= endy && i < h; ++i) {
			for (int j = startx; j < endx && j < w; ++j) {
				avg += data[i * w + j];
				num++;
			}
		}
		avg /= num;

		std = 0;
		for (int i = starty; i <= endy && i < h; ++i) {
			for (int j = startx; j < endx && j < w; ++j) {
				T d = avg - data[i * w + j];
				std += d * d;
			}
		}
		std /= num;
		std = sqrt(std);
	}
	void normalize(T valmin, T valmax) {
		int len = w * h;
		for (int i = 0; i < len; ++i) {
			T v = data[i];
			v = v < valmin ? valmin : v;
			v = v > valmax ? valmax : v;
			data[i] = (v - valmin) / valmax;
		}
	}

	void normalize() {
		T valmin = getMinElement();
		T valmax = getMaxElement();
		normalize(valmin, valmax);
	}

	T getMaxElement() const {
		return *std::max_element(data, data + w * h);
	}
	T getMinElement() const {
		return *std::min_element(data, data + w * h);
	}

	void getMinMaxElements(T& valmin, T& valmax) const {
		valmin = getMinElement();
		valmax = getMaxElement();
	}

	void log() {
		int len = w * h;
		for (int i = 0; i < len; ++i) {
			data[i] = std::log(data[i]);
		}
	}
	void exp() {
		int len = w * h;
		for (int i = 0; i < len; ++i) {
			data[i] = std::exp(data[i]);
		}
	}
};
template<class T>
class FloatImgRGB {
public:
	union {
		int w;
		int n;
		int cols;
	};
	union {
		int h;
		int m;
		int rows;
	};
	T* data;
public:
	FloatImgRGB() :
			w(0), h(0), data(0) {
	}
	FloatImgRGB(int iw, int ih) {
		assert(iw > 0 && ih > 0);
		data = new T[iw * ih * 3];
		w = iw;
		h = ih;
	}
	~FloatImgRGB() {
		clear();
	}
public:
	void clear() {
		if (data)
			delete[] data;
		data = 0;
		w = 0;
		h = 0;
	}
	void fill(double r, double g, double b) {
		int len = w * h;
		for (int i = 0; i < len; ++i) {
			data[3 * i] = r;
			data[3 * i + 1] = g;
			data[3 * i + 2] = b;
		}
	}
	void resize(int iw, int ih) {
		clear();
		if (iw > 0 && ih > 0) {
			data = new T[iw * ih * 3];
			w = iw;
			h = ih;
		}
	}
	void split(FloatImgG<T>& R, FloatImgG<T>& G, FloatImgG<T>& B) const {
		R.resize(w, h);
		G.resize(w, h);
		B.resize(w, h);
		int len = w * h;
		for (int i = 0; i < len; ++i) {
			R.data[i] = data[3 * i];
			G.data[i] = data[3 * i + 1];
			B.data[i] = data[3 * i + 2];
		}
	}
	bool empty() const {
		return w == 0 || h == 0;
	}
public:
	typedef T* PFLOAT;
	typedef const T* CONST_PFLOAT;
	operator PFLOAT() {
		return data;
	}
	operator CONST_PFLOAT() const {
		return data;
	}
	T* get_ptr() const {
		return data;
	}
	T* operator ()(int x, int y) const {
		return data + y * w * 3 + x * 3;
	}
	void set(int x, int y, T r, T g, T b) {
		data[y * w * 3 + x * 3] = r;
		data[y * w * 3 + x * 3 + 1] = g;
		data[y * w * 3 + x * 3 + 2] = b;
	}
	int getMemLen() const {
		return h * w * 3 * sizeof(T);
	}
	void getAvgStd(int startx, int endx, int starty, int endy, double avg[3],
			double std[3]) const {
		std::fill_n(avg, 3, 0.0);
		int N = 0;
		for (int r = starty; r < endy && r < h; ++r) {
			for (int c = startx; c < endx && c < w; ++c) {
				avg[0] += data[3 * r * w + 3 * c];
				avg[1] += data[3 * r * w + 3 * c + 1];
				avg[2] += data[3 * r * w + 3 * c + 2];
				N++;
			}
		}
		avg[0] /= N;
		avg[1] /= N;
		avg[2] /= N;

		std::fill_n(std, 3, 0.0);
		for (int r = starty; r < endy && r < h; ++r) {
			for (int c = startx; c < endx && c < w; ++c) {
				double dR = data[3 * r * w + 3 * c] - avg[0];
				double dG = data[3 * r * w + 3 * c + 1] - avg[1];
				double dB = data[3 * r * w + 3 * c + 2] - avg[2];
				std[0] += dR * dR;
				std[1] += dG * dG;
				std[2] += dB * dB;
			}
		}
		std[0] /= N;
		std[1] /= N;
		std[2] /= N;

		std[0] = sqrt(std[0]);
		std[1] = sqrt(std[1]);
		std[2] = sqrt(std[2]);
	}
	void normalize(T valmin, T valmax) {
		int len = w * h;
		for (int i = 0; i < len; ++i) {
			T r = data[3 * i];
			T g = data[3 * i + 1];
			T b = data[3 * i + 2];
			r = r < valmin ? valmin : r;
			r = r > valmax ? valmax : r;
			g = g < valmin ? valmin : g;
			g = g > valmax ? valmax : g;
			b = b < valmin ? valmin : b;
			b = b > valmax ? valmax : b;

			data[3 * i] = (r - valmin) / valmax;
			data[3 * i + 1] = (g - valmin) / valmax;
			data[3 * i + 2] = (b - valmin) / valmax;
		}
	}
	void normalize() {
		double minrgb[3], maxrgb[3];
		getMinElement(minrgb);
		getMaxElement(maxrgb);

		T valmin = *std::min_element(minrgb, minrgb + 3);
		T valmax = *std::max_element(maxrgb, maxrgb + 3);

		normalize(valmin, valmax);
	}
	void getMinElement(T minrgb[3]) {
		int len = w * h;
		std::fill_n(minrgb, 3, std::numeric_limits<T>::max());
		for (int i = 0; i < len; ++i) {
			T r = data[3 * i];
			T g = data[3 * i + 1];
			T b = data[3 * i + 2];

			if (r < minrgb[0])
				minrgb[0] = r;
			if (g < minrgb[1])
				minrgb[1] = g;
			if (b < minrgb[2])
				minrgb[2] = b;
		}
	}
	void getMaxElement(T maxrgb[3]) {
		int len = w * h;
		std::fill_n(maxrgb, 3, std::numeric_limits<T>::min());
		for (int i = 0; i < len; ++i) {
			T r = data[3 * i];
			T g = data[3 * i + 1];
			T b = data[3 * i + 2];

			if (r > maxrgb[0])
				maxrgb[0] = r;
			if (g > maxrgb[1])
				maxrgb[1] = g;
			if (b > maxrgb[2])
				maxrgb[2] = b;
		}
	}

	void getMinMaxElements(T minrgb[3], T maxrgb[3]) {
		getMinElement(minrgb);
		getMaxElement(maxrgb);
	}

	void getMeanChannel(FloatImgG<T>& meanI) const {
		meanI.resize(w, h);
		int len = w * h;
		for (int i = 0; i < len; ++i)
			meanI.data[i] = (data[3 * i] + data[3 * i + 1] + data[3 * i + 2])
					/ 3.0;

	}
	void log() {
		int len = w * h;
		for (int i = 0; i < len; ++i) {
			data[3 * i] = std::log(data[3 * i]);
			data[3 * i + 1] = std::log(data[3 * i + 1]);
			data[3 * i + 2] = std::log(data[3 * i + 2]);
		}
	}
	void exp() {
		int len = w * h;
		for (int i = 0; i < len; ++i) {
			data[3 * i] = std::exp(data[3 * i]);
			data[3 * i + 1] = std::exp(data[3 * i + 1]);
			data[3 * i + 2] = std::exp(data[3 * i + 2]);
		}
	}
};

#endif /* SL_FLOATIMAGE_H_ */
