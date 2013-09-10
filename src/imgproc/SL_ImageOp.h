/*
 * SL_ImageOp.h
 *
 *  Created on: 2010-11-8
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_IMAGEOP_H_
#define SL_IMAGEOP_H_
#include "opencv2/opencv.hpp"
#include "SL_Image.h"
#include "SL_error.h"
#include "math/SL_Array3D.h"
#include "string.h"

/* copy img1 to img2 */
template<typename T>
void copyImg(const T& img1, T& img2) {
	if (img1.w != img2.w || img1.h != img2.h)
		repErr("copy_img() : size is not the same.");
	int len = img1.getMemLen();
	memcpy(img2.data, img1.data, len);
}
/* clone img1 to img2*/
template<typename T>
void cloneImg(const T& img1, T& img2) {
	if (img1.w != img2.w || img1.h != img2.h) {
		img2.resize(img1.w, img1.h);
	}
	int len = img1.getMemLen();
	memcpy(img2.data, img1.data, len);
}

void rgb2gray(int w, int h, uchar* rgb, uchar* gray);
void rgb2gray(const ImgRGB& rgb_img, ImgG& gray_img);

void rgb2grayAvg(int w, int h, uchar* rgb, uchar* gray);
void rgb2grayAvg(const ImgRGB& rgb_img, ImgG& gray_img);

void gray2rgb(int w, int h, uchar* gray, uchar* rgb);
void gray2rgb(const ImgG& gray_img, ImgRGB& rgb_img);

void bgr2rgb(ImgRGB& rgb_out);
void scaleDownAvg(
		int rows,
		int cols,
		const uchar* input,
		int nrows,
		int ncols,
		uchar* output);
void scaleDownAvg(const ImgG& img1, ImgG& img2, double scale);

/*undistor the image by linear interpolation*/
void undistorImage(
		const ImgG& img1,
		ImgG& img2,
		const double* K,
		const double* d);
void undistorImage(
		const ImgRGB& img1,
		ImgRGB& img2,
		const double* K,
		const double* d);

/* get image block from a gray image */
void getImgBlock(int x0, int y0, int hw, const ImgG& img, double* blockData);

template<class RGB, class T>
void splitImageChannel(const RGB& rgb, MyMat<T>& R, MyMat<T>& G, MyMat<T>& B) {
	R.resize(rgb.m, rgb.n);
	G.resize(rgb.m, rgb.n);
	B.resize(rgb.m, rgb.n);
	int len = rgb.m * rgb.n;
	for (int i = 0; i < len; ++i) {
		R.data[i] = rgb.data[3 * i];
		G.data[i] = rgb.data[3 * i + 1];
		B.data[i] = rgb.data[3 * i + 2];
	}
}

template<class RGB, class T>
void img2Arr3(const RGB& rgb, MyArr3<T>& arr) {
	arr.resize(3, rgb.m, rgb.n);
	int len = rgb.m * rgb.n;
	for (int c = 0; c < 3; ++c) {
		T* pc = arr.getPtrChannel(c);
		for (int i = 0; i < len; ++i) {
			pc[i] = rgb.data[3 * i + c];
		}
	}
}

template<class U, class T>
void splitChannel(const MyArr3<U>& arr, MyMat<T>& r, MyMat<T>& g, MyMat<T>& b) {
	assert(arr.c == 3);
	r.resize(arr.m, arr.n);
	g.resize(arr.m, arr.n);
	b.resize(arr.m, arr.n);

	int len = arr.m * arr.n;
	T* pc = arr.getPtrChannel(0);
	for (int i = 0; i < len; ++i)
		r.data[i] = pc[i];

	pc = arr.getPtrChannel(1);
	for (int i = 0; i < len; ++i)
		g.data[i] = pc[i];

	pc = arr.getPtrChannel(2);
	for (int i = 0; i < len; ++i)
		b.data[i] = pc[i];

}

template<class RGB, class T>
void normalizeColorImage(
		const RGB& rgb,
		FloatImgRGB<T>& rgbn,
		T cmin = 0.0,
		T cmax = 255.0) {
	assert(cmax > cmin);
	rgbn.resize(rgb.w, rgb.h);
	int len = rgb.w * rgb.h;
	T range = (cmax - cmin);
	for (int i = 0; i < len; ++i) {
		rgbn.data[3 * i] = T(rgb.data[3 * i] - cmin) / range;
		rgbn.data[3 * i + 1] = T(rgb.data[3 * i + 1] - cmin) / range;
		rgbn.data[3 * i + 2] = T(rgb.data[3 * i + 2] - cmin) / range;
	}
}
template<class GRAY, class T>
void normalizeGrayImage(
		const GRAY& gray,
		FloatImgG<T>& grayn,
		T cmin = 0.0,
		T cmax = 255.0) {
	assert(cmax > cmin);
	grayn.resize(gray.w, gray.h);
	int len = gray.w * gray.h;
	T range = (cmax - cmin);
	for (int i = 0; i < len; ++i) {
		grayn.data[i] = T(gray.data[i] - cmin) / range;
	}
}

template<class RGB, class T>
void convert2FloatColorImage(const RGB& rgb, FloatImgRGB<T>& frgb) {
	assert(!rgb.empty());
	frgb.resize(rgb.w, rgb.h);
	int len = rgb.w * rgb.h;
	for (int i = 0; i < len; ++i) {
		frgb.data[3 * i] = rgb.data[3 * i];
		frgb.data[3 * i + 1] = rgb.data[3 * i + 1];
		frgb.data[3 * i + 2] = rgb.data[3 * i + 2];
	}
}

template<class T>
void convert2ColorImage(
		const FloatImgRGB<T>& frgb,
		ImgRGB& rgb,
		T minval = 0,
		T maxval = 1) {
	rgb.resize(frgb.w, frgb.h);
	int len = frgb.w * frgb.h;
	for (int i = 0; i < len; ++i) {
		double r = frgb.data[3 * i];
		double g = frgb.data[3 * i + 1];
		double b = frgb.data[3 * i + 2];

		r = r < minval ? minval : r;
		r = r > maxval ? maxval : r;

		g = g < minval ? minval : g;
		g = g > maxval ? maxval : g;

		b = b < minval ? minval : b;
		b = b > maxval ? maxval : b;

		rgb.data[3 * i] = uchar((r - minval) / maxval * 255.0);
		rgb.data[3 * i + 1] = uchar((g - minval) / maxval * 255.0);
		rgb.data[3 * i + 2] = uchar((b - minval) / maxval * 255.0);
	}
}

/*resize*/
void imresize(const ImgRGB& rgb_img, ImgRGB& new_img, int nw, int nh, int mode =
		cv::INTER_CUBIC);

void imresize(
		const ImgRGB& rgb_img,
		ImgRGB& new_img,
		double ratio = 1.0,
		int mode = cv::INTER_CUBIC);

void imresize(const ImgG& gray_img, ImgG& new_img, int nw, int nh, int mode =
		cv::INTER_CUBIC);

void imresize(
		const ImgG& gray_img,
		ImgG& new_img,
		double ratio = 1.0,
		int mode = cv::INTER_CUBIC);

void matresize(const Mat_d& in, Mat_d& out, int nm, int nn, int mode =
		cv::INTER_CUBIC);

void matresize(const Mat_d& in, Mat_d& out, int mode = cv::INTER_CUBIC);

void arr3resize(const Arr3d& in, Arr3d& out, int nm, int nn, int mode =
		cv::INTER_CUBIC);

template<class T, class U>
void imGetChannel(const T& rgbImg, int channel, MyMat<U>& mat) {
	assert(!rgbImg.empty());
	assert(channel >= 0 && channel < 3);

	mat.resize(rgbImg.m, rgbImg.n);
	int len = rgbImg.m * rgbImg.n;
	for (int i = 0; i < len; ++i)
		mat.data[i] = rgbImg.data[3 * i + channel];
}

void warp(ImgG& I_warped, const ImgG& I, const Mat_f& wx, const Mat_f& wy);

void warp(ImgRGB& I_warped, const ImgRGB& I, const Mat_f& wx, const Mat_f& wy);

void warp(ImgRGB& I_warped, const ImgRGB& I, const Mat_d& wx, const Mat_d& wy);

void warp(
		ImgRGB_d& I_warped,
		const ImgRGB_d& I,
		const Mat_f& wx,
		const Mat_f& wy);

void warp(Mat_d& I_warped, const Mat_d& I, const Mat_f& wx, const Mat_f& wy);
void warp(Mat_d& I_warped, const ImgG& I, const Mat_f& wx, const Mat_f& wy);

void flow2WarpMap(
		const Mat_f& ux,
		const Mat_f& uy,
		Mat_f& wx,
		Mat_f& wy,
		bool add = true);

#include <cmath>
void gammaEncode(ImgRGB& rgb);
void gammaEncode(ImgG& gray);
void gammaDecode(ImgRGB& rgb);
void gammaDecode(ImgG& gray);

void whiteBalance(ImgRGB& rgb, double wr, double wg, double wb);
void whiteBalanceInv(ImgRGB& rgb, double wr, double wg, double wb);

double getLightestColor(const ImgRGB& rgb, double color[3]);
#endif /* SL_IMAGEOP_H_ */

