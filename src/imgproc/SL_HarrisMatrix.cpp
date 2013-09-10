/*
 * SL_HarrisMatrix.cpp
 *
 *  Created on: Feb 13, 2012
 *      Author: Danping Zou
 */

#include "SL_HarrisMatrix.h"
#include "imgproc/SL_ImageOp.h"
#include "imgproc/SL_BoxFilter.h"
#include "imgproc/SL_Gradient.h"
#include "opencv2/opencv.hpp"

void computeHarrisMatrices(const ImgG& img, Mat_d& IxIx, Mat_d& IxIy,
		Mat_d& IyIy, double sigma) {
	assert(img.w > 0 && img.h > 0);
	Mat_d Ix, Iy;
	Ix.resize(img.m, img.n);
	Iy.resize(img.m, img.n);

	int m = img.h;
	int n = img.w;

	cv::Mat cvImg(m, n, CV_8UC1, img.data);
	cv::Mat cvIx(m, n, CV_64FC1, Ix.data);
	cv::Mat cvIy(m, n, CV_64FC1, Iy.data);

	cv::Sobel(cvImg, cvIx, CV_64FC1, 1, 0, 3, 1. / 255.0, 0);
	cv::Sobel(cvImg, cvIy, CV_64FC1, 0, 1, 3, 1. / 255.0, 0);

	Mat_d IxIx_org, IxIy_org, IyIy_org;
	IxIx_org.resize(m, n);
	IxIy_org.resize(m, n);
	IyIy_org.resize(m, n);

	IxIx.resize(m, n);
	IxIy.resize(m, n);
	IyIy.resize(m, n);
	int len = img.w * img.h;
	for (int i = 0; i < len; ++i) {
		IxIx_org.data[i] = Ix.data[i] * Ix.data[i];
		IxIy_org.data[i] = Ix.data[i] * Iy.data[i];
		IyIy_org.data[i] = Iy.data[i] * Iy.data[i];
	}

	cv::Mat cvIxIx_org(m, n, CV_64FC1, IxIx_org.data);
	cv::Mat cvIxIy_org(m, n, CV_64FC1, IxIy_org.data);
	cv::Mat cvIyIy_org(m, n, CV_64FC1, IyIy_org.data);

	cv::Mat cvIxIx(m, n, CV_64FC1, IxIx.data);
	cv::Mat cvIxIy(m, n, CV_64FC1, IxIy.data);
	cv::Mat cvIyIy(m, n, CV_64FC1, IyIy.data);

	cv::GaussianBlur(cvIxIx_org, cvIxIx, cv::Size2i(0, 0), sigma, sigma);
	cv::GaussianBlur(cvIxIy_org, cvIxIy, cv::Size2i(0, 0), sigma, sigma);
	cv::GaussianBlur(cvIyIy_org, cvIyIy, cv::Size2i(0, 0), sigma, sigma);
}

void computeEigenValues(const Mat_d& IxIx, const Mat_d& IxIy, const Mat_d& IyIy,
		Mat_d& v1, Mat_d& v2, Mat_d& lambda1, Mat_d& lambda2) {
	assert(IxIx.m > 0 && IxIx.n > 0);
	assert(
			IxIx.m == IxIy.m && IxIy.m == IyIy.m && IxIx.n == IxIy.n && IxIy.n == IyIy.n);

	int m = IxIx.m;
	int n = IxIx.n;

	int len = m * n;

	v1.resize(len, 2);
	v2.resize(len, 2);
	lambda1.resize(m, n);
	lambda2.resize(m, n);

	for (int i = 0; i < len; ++i) {
		double a = IxIx.data[i];
		double b = IxIy.data[i];
		double d = IyIy.data[i];

		double delta = sqrt((a - d) * (a - d) + 4 * b * b);

		double l1 = 0.5 * (a + d + delta);
		double l2 = 0.5 * (a + d - delta);

		lambda1.data[i] = l1;
		lambda2.data[i] = l2;

		double L1 = sqrt(b * b + (l1 - a) * (l1 - a));
		v1.data[2 * i] = b / L1;
		v1.data[2 * i + 1] = (l1 - a) / L1;

		v2.data[2 * i] = -v1.data[2 * i + 1];
		v2.data[2 * i + 1] = v1.data[2 * i];
	}
}

void computeHarrisResponse(const Mat_d& IxIx, const Mat_d& IxIy,
		const Mat_d& IyIy, Mat_d& res, double kata) {
	assert(IxIx.m > 0 && IxIx.n > 0);
	assert(
			IxIx.m == IxIy.m && IxIy.m == IyIy.m && IxIx.n == IxIy.n && IxIy.n == IyIy.n);

	int m = IxIx.m;
	int n = IxIx.n;

	res.resize(m, n);

	int len = m * n;

	for (int i = 0; i < len; ++i) {
		double det = IxIx.data[i] * IyIy.data[i] - IxIy.data[i] * IxIy.data[i];
		double tr = IxIx.data[i] + IyIy.data[i];
		res.data[i] = det - kata * tr * tr;
	}
}

void computeCornerEig(const ImgG& img, Mat_d& lambda1, Mat_d& lambda2,
		Mat_d& v1, Mat_d& v2, int block_size) {
	int m = img.h;
	int n = img.w;
	cv::Mat cvImg(m, n, CV_8UC1, img.data);

	Mat_f res(m * n, 6);
	cv::Mat cvRes(m, n, CV_32FC(6), res.data);
	cv::cornerEigenValsAndVecs(cvImg, cvRes, block_size, 3);

	lambda1.resize(m, n);
	lambda2.resize(m, n);
	v1.resize(m * n, 2);
	v2.resize(m * n, 2);

	int len = m * n;
	for (int i = 0; i < len; ++i) {
		lambda1.data[i] = res.data[6 * i];
		lambda2.data[i] = res.data[6 * i + 1];

		v1.data[2 * i] = res.data[6 * i + 2];
		v1.data[2 * i + 1] = res.data[6 * i + 3];
		v2.data[2 * i] = res.data[6 * i + 4];
		v2.data[2 * i + 1] = res.data[6 * i + 5];
	}
}