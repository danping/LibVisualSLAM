/*
 * GUI_ImageViewer.h
 *
 *  Created on: 2010-11-8
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef GUI_IMAGEVIEWER_H_
#define GUI_IMAGEVIEWER_H_

#include "imgproc/SL_Image.h"
#include "imgproc/SL_ImageOp.h"
#include "math/SL_Matrix.h"
#include "tools/SL_TypeConversion.h"
#include "opencv2/opencv.hpp"

#include <algorithm>

void imshow(const char* name, const ImgRGB_d& img, double scale = 1.0);
void imshow(const char* name, const ImgRGB& img, double scale = 1.0);
void imshow(const char* name, const ImgG& img, double scale = 1.0);
void imshow(const char* name, int w, int h, const uchar* data,
		bool gray = true);
template<class T>
void imshow(const char* name, int m, int n, const T* data, T minval,
		double maxval) {
	ImgG gray(n, m);
	mat2Img(n, m, data, gray, minval, maxval);
	imshow(name, gray);
}
template<class T>
void imshow(const char* name, int m, int n, const T* data) {
	int len = m * n;
	T minval = *std::min_element(data, data + len);
	T maxval = *std::max_element(data, data + len);
	cout << name << " [" << minval <<" - " << maxval << "]" << endl;
	imshow(name, m, n, data, minval, maxval);
}
template<class T>
void imshow(const char* name, const MyMat<T>& mat, T minval, T maxval) {
	ImgG gray(mat.n, mat.m);
	mat2Img(mat.n, mat.m, mat.data, gray, minval, maxval);
	imshow(name, gray);
}
template<class T>
void imshow(const char* name, const MyMat<T>& mat) {
	//test
	cout << name << endl;
	ImgG gray(mat.n, mat.m);
	mat2ImgAuto(mat.n, mat.m, mat.data, gray);
	imshow(name, gray);
}
void imshow(const char* name, const Arr3d& arr3);
void imshowext(const char* name, const ImgRGB& img);
void imshowext(const char* name, const ImgG& img);
void imshowext(const char* name, const uchar* data, int w, int h, bool gray);
#endif /* GUI_IMAGEVIEWER_H_ */
