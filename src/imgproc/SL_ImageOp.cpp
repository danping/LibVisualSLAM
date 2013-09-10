/*
 * SL_ImageOp.cpp
 *
 *  Created on: 2010-11-8
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#include "SL_ImageOp.h"
#include "SL_Image.h"
#include "SL_FloatImage.h"
#include "SL_error.h"
#include "math/SL_Matrix.h"
#include "math/SL_Array3D.h"
#include "math/SL_LinAlg.h"

#include "util/SL_Utility.h"
#include "tools/SL_TypeConversion.h"
#include <cmath>
#include <cassert>
#include <vector>

void rgb2gray(int w, int h, uchar* rgb, uchar* gray) {
	int i, len = w * h;
	for (i = 0; i < len; i += 4) {
		//gray[i] = rgb[j] * .2989 + rgb[j + 1] * 0.5870 + rgb[j + 2] * 0.1140;
		gray[0] = rgb[0] * .2989 + rgb[1] * 0.5870 + rgb[2] * 0.1140;
		gray[1] = rgb[3] * .2989 + rgb[4] * 0.5870 + rgb[5] * 0.1140;
		gray[2] = rgb[6] * .2989 + rgb[7] * 0.5870 + rgb[8] * 0.1140;
		gray[3] = rgb[9] * .2989 + rgb[10] * 0.5870 + rgb[11] * 0.1140;
		gray += 4;
		rgb += 12;
	}
	for (; i < len; i++) {
		gray[0] = rgb[0] * .2989 + rgb[1] * 0.5870 + rgb[2] * 0.1140;
		gray++;
	}
}
void rgb2grayAvg(int w, int h, uchar* rgb, uchar* gray) {
	int i, len = w * h;
	for (i = 0; i < len; ++i) {
		gray[i] = (uchar) (((double) (rgb[3 * i]) + rgb[3 * i + 1]
				+ rgb[3 * i + 2]) / 3.0);
	}
}
void gray2rgb(int w, int h, uchar* gray, uchar* rgb) {
	int i, j = 0, len = w * h;
	for (i = 0; i < len; ++i) {
		//gray[idx1] = rgb[idx2] * .2989 + rgb[idx2 + 1] * 0.5870 + rgb[idx2 + 2] * 0.1140;
		rgb[j] = gray[i];
		rgb[j + 1] = gray[i];
		rgb[j + 2] = gray[i];
		j += 3;
	}
}
void gray2rgb(const ImgG& gray_img, ImgRGB& rgb_img) {
	rgb_img.resize(gray_img.w, gray_img.h);
	gray2rgb(gray_img.w, gray_img.h, gray_img.data, rgb_img.data);
}
void rgb2gray(const ImgRGB& rgb_img, ImgG& gray_img) {
	gray_img.resize(rgb_img.w, rgb_img.h);
	rgb2gray(rgb_img.w, rgb_img.h, rgb_img.data, gray_img);
}

void rgb2grayAvg(const ImgRGB& rgb_img, ImgG& gray_img) {
	gray_img.resize(rgb_img.w, rgb_img.h);
	rgb2grayAvg(rgb_img.w, rgb_img.h, rgb_img.data, gray_img.data);
}
void bgr2rgb(ImgRGB& rgb) {
	if (rgb.empty())
		return;
	int len = rgb.m * rgb.n;
	for (int i = 0; i < len; ++i) {
		uchar v = rgb.data[3 * i];
		rgb.data[3 * i] = rgb.data[3 * i + 2];
		rgb.data[3 * i + 2] = v;
	}
}

void scaleDownAvg(
		int rows,
		int cols,
		const uchar* input,
		int nrows,
		int ncols,
		uchar* output) {
	//zooming ratio
	double ratio_r = double(nrows - 1) / (rows - 1);
	double ratio_c = double(ncols - 1) / (cols - 1);

	if (ratio_r > 1 || ratio_c > 1)
		repErr("scale_down_avg() : only scaling down is supported.");

	if (ratio_r == 1) {
		memcpy(output, input, rows * cols);
		return;
	}

	Mat_d counter(nrows, ncols);
	Mat_d val(nrows, ncols);
	counter.fill(0);
	val.fill(0);

	int i, len = rows * cols;
	for (i = 0; i < len; ++i) {
		int r = i / cols;
		int c = i - r * cols;

		//compute the coordinates in the downsized image
		int nr = int(ratio_r * r + 0.5);
		int nc = int(ratio_c * c + 0.5);

		val[nr * ncols + nc] += double(input[i]);
		counter[nr * ncols + nc]++;
	}
	//normalize
	len = nrows * ncols;
	for (i = 0; i < len; ++i)
		output[i] = static_cast<uchar>(val[i] / counter[i]);
}

void scaleDownAvg(const ImgG& img1, ImgG& img2, double scale) {
	if (scale > 1)
		repErr("scale_down_avg() : only scaling down is supported.");
	if (scale == 1) {
		cloneImg(img1, img2);
		return;
	}

	int w = img1.w;
	int h = img1.h;

	int nw = static_cast<int>(w * scale);
	int nh = static_cast<int>(h * scale);

	img2.resize(nw, nh);
	scaleDownAvg(h, w, img1, nh, nw, img2);
}

void undistorImage(
		const ImgG& img1,
		ImgG& img2,
		const double* K,
		const double* d) {
	int w = img1.w;
	int h = img1.h;
	img2.resize(w, h);
	Mat_d iK(3, 3);
	mat33Inv(K, iK.data);

	int x, y, k = 0;
	for (y = 0; y < h; ++y) {
		for (x = 0; x < w; ++x) {
			//get the normalized points
			double zn0 = x * iK[6] + y * iK[7] + iK[8];
			double xn0 = (x * iK[0] + y * iK[1] + iK[2]) / zn0;
			double yn0 = (x * iK[3] + y * iK[4] + iK[5]) / zn0;

			//add distortion
			double r = sqrt(xn0 * xn0 + yn0 * yn0);
			double rr = r * r;
			double rrr = rr * r;
			double rrrrr = rrr * rr;
			double rrrrrrr = rrrrr * rr;
			double r1 = r + d[0] * rrr + d[1] * rrrrr + d[4] * rrrrrrr;

			double factor = r1 / r;
			double xn = xn0 * factor;
			double yn = yn0 * factor;

			//go back to the image coordinates
			double zold = xn * K[6] + yn * K[7] + K[8];
			double xold = (xn * K[0] + yn * K[1] + K[2]) / zold;
			double yold = (xn * K[3] + yn * K[4] + K[5]) / zold;

			int x1 = static_cast<int>(xold);
			int y1 = static_cast<int>(yold);

			int x2 = x1 + 1;
			int y2 = y1 + 2;

			double dx = xold - x1;
			double dy = yold - y1;
			double dx1 = 1 - dx;
			double dy1 = 1 - dy;
			if (x1 < 0 || x2 >= w || y1 < 0 || y2 >= h)
				img2.data[k] = 0;

			img2.data[k] = img1.data[y1 * w + x1] * dx1 * dy1
					+ img1.data[y1 * w + x2] * dx * dy1
					+ img1.data[y2 * w + x2] * dx * dy
					+ img1.data[y2 * w + x1] * dx1 * dy;
			++k;
		}
	}
}
void undistorImage(
		const ImgRGB& img1,
		ImgRGB& img2,
		const double* K,
		const double* d) {
	int w = img1.w;
	int h = img1.h;
	img2.resize(w, h);
	Mat_d iK(3, 3);
	mat33Inv(K, iK.data);

	int x, y, k = 0;
	for (y = 0; y < h; ++y) {
		for (x = 0; x < w; ++x) {
			//get the normalized points
			double zn0 = x * iK[6] + y * iK[7] + iK[8];
			double xn0 = (x * iK[0] + y * iK[1] + iK[2]) / zn0;
			double yn0 = (x * iK[3] + y * iK[4] + iK[5]) / zn0;

			//add distortion
			double r = sqrt(xn0 * xn0 + yn0 * yn0);
			double rr = r * r;
			double rrr = rr * r;
			double rrrrr = rrr * rr;
			double rrrrrrr = rrrrr * rr;
			double r1 = r + d[0] * rrr + d[1] * rrrrr + d[4] * rrrrrrr;

			double factor = r1 / r;
			double xn = xn0 * factor;
			double yn = yn0 * factor;

			//go back to the image coordinates
			double zold = xn * K[6] + yn * K[7] + K[8];
			double xold = (xn * K[0] + yn * K[1] + K[2]) / zold;
			double yold = (xn * K[3] + yn * K[4] + K[5]) / zold;

			int x1 = static_cast<int>(xold);
			int y1 = static_cast<int>(yold);

			int x2 = x1 + 1;
			int y2 = y1 + 1;

			double dx = xold - x1;
			double dy = yold - y1;
			double dx1 = 1 - dx;
			double dy1 = 1 - dy;
			if (x1 < 0 || x2 >= w || y1 < 0 || y2 >= h) {
				img2.data[k] = 0;
				img2.data[k + 1] = 0;
				img2.data[k + 2] = 0;
				k += 3;
			}

			int l1 = y1 * 3 * w + x1 * 3;
			int l2 = y1 * 3 * w + x2 * 3;
			int l3 = y2 * 3 * w + x2 * 3;
			int l4 = y2 * 3 * w + x1 * 3;
			img2.data[k] = img1.data[l1] * dx1 * dy1 + img1.data[l2] * dx * dy1
					+ img1.data[l3] * dx * dy + img1.data[l4] * dx1 * dy;
			++k;
			++l1;
			++l2;
			++l3;
			++l4;
			img2.data[k] = img1.data[l1] * dx1 * dy1 + img1.data[l2] * dx * dy1
					+ img1.data[l3] * dx * dy + img1.data[l4] * dx1 * dy;
			++k;
			++l1;
			++l2;
			++l3;
			++l4;
			img2.data[k] = img1.data[l1] * dx1 * dy1 + img1.data[l2] * dx * dy1
					+ img1.data[l3] * dx * dy + img1.data[l4] * dx1 * dy;
			++k;
		}
	}
}
void getImgBlock(int x0, int y0, int hw, const ImgG& img, double* blockData) {
	const int W = img.cols;
	int k = 0;
	for (int y = y0 - hw; y <= y0 + hw; y++) {
		for (int x = x0 - hw; x <= x0 + hw; x++) {
			blockData[k++] = img.data[y * W + x];
		}
	}
}

#include <opencv2/opencv.hpp>
void imresize(
		const ImgRGB& rgb_img,
		ImgRGB& new_img,
		int nw,
		int nh,
		int mode) {
	assert(nw > 0 && nh > 0);
	new_img.resize(nw, nh);

	cv::Mat cvrgb(rgb_img.h, rgb_img.w, CV_8UC3, rgb_img.data);
	cv::Mat cvnew(new_img.h, new_img.w, CV_8UC3, new_img.data);

	cv::resize(cvrgb, cvnew, cvnew.size(), 0, 0, mode);
}

void imresize(const ImgRGB& rgb_img, ImgRGB& new_img, double ratio, int mode) {
	assert(ratio > 0);
	int nw = (int) (rgb_img.w * ratio);
	int nh = (int) (rgb_img.h * ratio);
	imresize(rgb_img, new_img, nw, nh, mode);
}
void imresize(const ImgG& gray_img, ImgG& new_img, int nw, int nh, int mode) {
	assert(nw > 0 && nh > 0);
	new_img.resize(nw, nh);

	cv::Mat cvrgb(gray_img.h, gray_img.w, CV_8UC1, gray_img.data);
	cv::Mat cvnew(new_img.h, new_img.w, CV_8UC1, new_img.data);

	cv::resize(cvrgb, cvnew, cvnew.size(), 0, 0, mode);
}
void imresize(const ImgG& gray_img, ImgG& new_img, double ratio, int mode) {
	assert(ratio > 0);
	int nw = (int) (gray_img.w * ratio);
	int nh = (int) (gray_img.h * ratio);
	imresize(gray_img, new_img, nw, nh, mode);
}
void matresize(const Mat_d& in, Mat_d& out, int nm, int nn, int mode) {
	assert(nm > 0 && nn > 0);
	assert(!in.empty());
	out.resize(nm, nn);

	cv::Mat cvin(in.m, in.n, CV_64FC1, in.data);
	cv::Mat cvout(out.m, out.n, CV_64FC1, out.data);

	cv::resize(cvin, cvout, cvout.size(), 0, 0, mode);
}
void matresize(const Mat_d& in, Mat_d& out, int mode) {
	assert(!in.empty() && !out.empty());
	cv::Mat cvin(in.m, in.n, CV_64FC1, in.data);
	cv::Mat cvout(out.m, out.n, CV_64FC1, out.data);

	cv::resize(cvin, cvout, cvout.size(), 0, 0, mode);
}

void arr3resize(const Arr3d& in, Arr3d& out, int nm, int nn, int mode) {
	assert(nm > 0 && nn > 0);
	assert(!in.empty());
	out.resize(in.c, nm, nn);

	for (int c = 0; c < in.c; ++c) {
		cv::Mat cvin(in.m, in.n, CV_64FC1, in.getPtrChannel(c));
		cv::Mat cvout(out.m, out.n, CV_64FC1, out.getPtrChannel(c));
		cv::resize(cvin, cvout, cvout.size(), 0, 0, mode);
	}
}
void warp(ImgG& I_warped, const ImgG& I, const Mat_f& wx, const Mat_f& wy) {
	assert(!I.empty());
	I_warped.resize(wx.n, wx.m);
	cv::Mat cvSrc(I.m, I.n, CV_8UC1, I.data);
	cv::Mat cvDst(I_warped.m, I_warped.n, CV_8UC1, I_warped.data);
	cv::Mat cvU(wx.m, wx.n, CV_32FC1, wx.data);
	cv::Mat cvV(wy.m, wy.n, CV_32FC1, wy.data);

	cv::remap(cvSrc, cvDst, cvU, cvV, cv::INTER_CUBIC, cv::BORDER_REFLECT,
			cv::Scalar(0));
}
void warp(ImgRGB& I_warped, const ImgRGB& I, const Mat_f& wx, const Mat_f& wy) {
	assert(!I.empty() && wx.m == I.m && wy.n == I.n);
	I_warped.resize(wx.n, wx.m);

	cv::Mat cvSrc(I.m, I.n, CV_8UC3, I.data);
	cv::Mat cvDst(I_warped.m, I_warped.n, CV_8UC3, I_warped.data);
	cv::Mat cvU(wx.m, wx.n, CV_32FC1, wx.data);
	cv::Mat cvV(wy.m, wy.n, CV_32FC1, wy.data);

	cv::remap(cvSrc, cvDst, cvU, cvV, cv::INTER_CUBIC, cv::BORDER_CONSTANT,
			cv::Scalar(0, 0, 0));
}

void warp(ImgRGB& I_warped, const ImgRGB& I, const Mat_d& wx, const Mat_d& wy) {
	Mat_f wxf, wyf;
	matDouble2Float(wx, wxf);
	matDouble2Float(wy, wyf);
	warp(I_warped, I, wxf, wyf);
}
void warp(
		ImgRGB_d& I_warped,
		const ImgRGB_d& I,
		const Mat_f& wx,
		const Mat_f& wy) {
	assert(!I.empty() && wx.m == I.m && wy.n == I.n);
	I_warped.resize(wx.n, wx.m);

	cv::Mat cvSrc(I.m, I.n, CV_64FC3, I.data);
	cv::Mat cvDst(I_warped.m, I_warped.n, CV_64FC3, I_warped.data);
	cv::Mat cvU(wx.m, wx.n, CV_32FC1, wx.data);
	cv::Mat cvV(wy.m, wy.n, CV_32FC1, wy.data);

	cv::remap(cvSrc, cvDst, cvU, cvV, cv::INTER_CUBIC, cv::BORDER_REFLECT,
			cv::Scalar(0, 0, 0));
}
void warp(Mat_d& I_warped, const Mat_d& I, const Mat_f& wx, const Mat_f& wy) {
	assert(!I.empty() && wx.m == I.m && wy.n == I.n);
	int m = wx.m;
	int n = wx.n;
	I_warped.resize(m, n);

	cv::Mat cvSrc(I.m, I.n, CV_64FC1, I.data);
	cv::Mat cvDst(m, n, CV_64FC1, I_warped.data);

	cv::Mat cvU(wx.m, wx.n, CV_32FC1, wx.data);
	cv::Mat cvV(wy.m, wy.n, CV_32FC1, wy.data);

	cv::remap(cvSrc, cvDst, cvU, cvV, cv::INTER_CUBIC, cv::BORDER_REFLECT,
			cv::Scalar(0, 0, 0));
}
void warp(Mat_d& I_warped, const ImgG& I, const Mat_f& wx, const Mat_f& wy) {
	ImgG Iw;
	warp(Iw, I, wx, wy);
	img2Mat(Iw, I_warped, 255.0);
}
void flow2WarpMap(
		const Mat_f& ux,
		const Mat_f& uy,
		Mat_f& wx,
		Mat_f& wy,
		bool add) {
	wx.resize(ux.m, ux.n);
	wy.resize(ux.m, ux.n);

	for (int y = 0; y < ux.m; ++y) {
		for (int x = 0; x < ux.n; ++x) {
			wx(y, x) = x + (add ? ux(y, x) : -ux(y, x));
			wy(y, x) = y + (add ? uy(y, x) : -uy(y, x));
		}
	}
}

void gammaEncode(ImgRGB& rgb) {
	assert(!rgb.empty());
	int len = rgb.w * rgb.h;
	for (int i = 0; i < len; ++i) {
		double r = rgb.data[3 * i];
		double g = rgb.data[3 * i + 1];
		double b = rgb.data[3 * i + 2];

		double nr = pow(r / 255.0, 0.45) * 255.0;
		double ng = pow(g / 255.0, 0.45) * 255.0;
		double nb = pow(b / 255.0, 0.45) * 255.0;

		rgb.data[3 * i] = (uchar) (nr);
		rgb.data[3 * i + 1] = (uchar) (ng);
		rgb.data[3 * i + 2] = (uchar) (nb);
	}
}
void gammaEncode(ImgG& gray) {
	assert(!gray.empty());
	int len = gray.w * gray.h;
	for (int i = 0; i < len; ++i) {
		double g = gray.data[i];
		double ng = pow(g / 255.0, 0.45) * 255.0;
		gray.data[i] = (uchar) (ng);
	}
}
void gammaDecode(ImgRGB& rgb) {
	assert(!rgb.empty());
	int len = rgb.w * rgb.h;
	for (int i = 0; i < len; ++i) {
		double r = rgb.data[3 * i];
		double g = rgb.data[3 * i + 1];
		double b = rgb.data[3 * i + 2];

		double nr = pow(r / 255.0, 2.2) * 255.0;
		double ng = pow(g / 255.0, 2.2) * 255.0;
		double nb = pow(b / 255.0, 2.2) * 255.0;

		rgb.data[3 * i] = (uchar) (nr);
		rgb.data[3 * i + 1] = (uchar) (ng);
		rgb.data[3 * i + 2] = (uchar) (nb);
	}
}

void gammaDecode(ImgG& gray) {
	assert(!gray.empty());
	int len = gray.w * gray.h;
	for (int i = 0; i < len; ++i) {
		double g = gray.data[i];
		double ng = pow(g / 255.0, 2.2) * 255.0;
		gray.data[i] = (uchar) (ng);
	}
}

void whiteBalance(ImgRGB& rgb, double wr, double wg, double wb) {
	assert(!rgb.empty());
	int len = 3 * rgb.w * rgb.h;
	for (int i = 0; i < len; i += 3) {
		double r = rgb.data[i];
		double g = rgb.data[i + 1];
		double b = rgb.data[i + 2];

		double nr = r / wr * 255.0;
		double ng = g / wg * 255.0;
		double nb = b / wb * 255.0;

		rgb.data[i] = nr >= 255.0 ? 255 : uchar(nr);
		rgb.data[i + 1] = ng >= 255.0 ? 255 : uchar(ng);
		rgb.data[i + 2] = nb >= 255.0 ? 255 : uchar(nb);
	}
}
void whiteBalanceInv(ImgRGB& rgb, double wr, double wg, double wb) {
	assert(!rgb.empty());
	int len = 3 * rgb.w * rgb.h;
	for (int i = 0; i < len; i += 3) {
		double r = rgb.data[i];
		double g = rgb.data[i + 1];
		double b = rgb.data[i + 2];

		double nr = r / 255.0 * wr;
		double ng = g / 255.0 * wg;
		double nb = b / 255.0 * wr;

		rgb.data[i] = nr >= 255.0 ? 255 : uchar(nr);
		rgb.data[i + 1] = ng >= 255.0 ? 255 : uchar(ng);
		rgb.data[i + 2] = nb >= 255.0 ? 255 : uchar(nb);
	}
}

double getLightestColor(const ImgRGB& rgb, double color[3]) {
	assert(!rgb.empty());
	using namespace std;
	ImgG gray;
	rgb2gray(rgb, gray);

	int n = rgb.w * rgb.h;
	int top_num = n * 0.001;
	top_num = top_num < 1 ? 1 : top_num;

	vector<int> ind;
	ind.resize(n);

	partialSortWithInd(n, gray.data, &ind[0], top_num, false);
	fill_n(color, 3, 0);

	for (int k = 0; k < top_num; ++k) {
		int i = ind[k];
		color[0] += rgb.data[3 * i];
		color[1] += rgb.data[3 * i + 1];
		color[2] += rgb.data[3 * i + 2];
	}

	color[0] /= top_num;
	color[1] /= top_num;
	color[2] /= top_num;

	return color[0] * color[0] + color[1] * color[1] + color[2] * color[2];

}

