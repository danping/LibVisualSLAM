/*
 * SL_TypeConversion.h
 *
 *  Created on: 2010-11-9
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_TYPECONVERSION_H_
#define SL_TYPECONVERSION_H_
#include "math/SL_Matrix.h"
#include "imgproc/SL_Image.h"
#include "imgproc/SL_ImageOp.h"
#include "geometry/SL_Point.h"
#include <cassert>
#include <algorithm>
#include <iostream>
void matDouble2Int(const Mat_d& md, Mat_i& mi);
void matDouble2Float(const Mat_d& md, Mat_f& mf);
void matFloat2Double(const Mat_f& mf, Mat_d& md);
void matFloat2Int(const Mat_f& mf, Mat_i& mi);
void VecPoint2D2Matrix(const VecPoint2d& pts, Mat_d& mat);

template<class T>
void mat2ImgAuto(int W, int H, T* data, ImgG& gray) {
	assert(W > 0 && H > 0);
	gray.resize(W, H);

	int i, len = W * H;
	double max_v = static_cast<double>(*std::max_element(data, data + len));
	double min_v = static_cast<double>(*std::min_element(data, data + len));

	double r = max_v - min_v;
	if (r > 0) {
		for (i = 0; i < len; ++i)
			gray.data[i] = static_cast<uchar>((data[i] - min_v) / r * 255);
	} else {
		for (i = 0; i < len; ++i)
			gray.data[i] = 100;
	}
}
template<class T>
void mat2Img(int W, int H, const T* data, ImgG& gray, T minval, T maxval) {
	assert(maxval > minval);
	assert(W > 0 && H > 0);
	gray.resize(W, H);

	double r = maxval - minval;
	int i, len = W * H;
	if (r > 0) {
		for (i = 0; i < len; ++i) {
			T val = data[i];
			if (val > maxval)
				val = maxval;
			if (val < minval)
				val = minval;

			double new_val = (val - minval) / r * 255.0;

			if (new_val > 255.0)
				new_val = 255.0;
			if (new_val < 0.0)
				new_val = 0.0;
			gray.data[i] = static_cast<uchar>(new_val);
		}
	} else {
		for (int i = 0; i < len; ++i) {
			gray.data[i] = 100;
		}
	}
}
template<class T>
void mat2Img(const T& mat, ImgG& gray) {
	int W = mat.cols;
	int H = mat.rows;
	mat2ImgAuto(W, H, mat.data, gray);
}

template<class T>
void mat2Img(const MyMat<T>& mat, ImgG& gray, T minval, T maxval) {
	int W = mat.cols;
	int H = mat.rows;
	mat2Img(W, H, mat.data, gray, minval, maxval);
}
template<class T>
void mat2RGBImg(const MyMat<T>& mat, ImgRGB& rgb, T minval, T maxval) {
	ImgG G;
	mat2Img(mat, G, minval, maxval);
	gray2rgb(G, rgb);
}
template<class T>
void img2Mat(const ImgG& gray, MyMat<T>& mat, T maxval /*= 255.0*/) {
	int len = gray.w * gray.h;
	mat.resize(gray.h, gray.w);
	for (int i = 0; i < len; ++i) {
		mat.data[i] = (T) (gray.data[i]) / maxval;
	}
}
#endif /* SL_TYPECONVERSION_H_ */
