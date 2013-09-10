/*
 * GUI_ImageViewer.cpp
 *
 *  Created on: 2010-11-8
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#include "GUI_ImageViewer.h"
#include "imgproc/SL_ImageOp.h"
#include "imgproc/SL_ImageIO.h"
#include "tools/SL_TypeConversion.h"

#include <algorithm>
#include <cstdio>
using namespace std;
//#ifdef USE_OPENCV
void imshow(const char* name, const ImgRGB_d& img, double scale) {
	assert(!img.empty());
	ImgRGB rgb;
	convert2ColorImage(img, rgb, 0.0, 255.0);
	imshow(name, rgb, scale);
}
void imshow(const char* name, const ImgRGB& img, double scale) {
	assert(!img.empty());
	cv::Mat mat(img.rows, img.cols, CV_8UC3, img.data);
	cv::Mat mat_rgb(img.rows, img.cols, CV_8UC3);
	cv::cvtColor(mat, mat_rgb, CV_BGR2RGB);
	if (scale == 1.0) {
		cv::imshow(name, mat_rgb);
	} else {
		cv::Mat tmp_rgb(img.rows, img.cols, CV_8UC3);
		cv::resize(mat_rgb, tmp_rgb,
				cv::Size(img.rows * scale, img.cols * scale));
		cv::imshow(name, tmp_rgb);
	}
}
void imshow(const char* name, const ImgG& img, double scale) {
	assert(!img.empty());
	cv::Mat mat(img.rows, img.cols, CV_8UC1, img.data);
	if (scale == 1.0)
		cv::imshow(name, mat);
	else {
		cv::Mat tmp(img.rows, img.cols, CV_8UC1);
		cv::resize(mat, tmp, cv::Size(img.rows * scale, img.cols * scale));
		cv::imshow(name, tmp);
	}
	//cv::waitKey(-1);
}
void imshow(const char* name, int w, int h, const uchar* data, bool gray) {
	if (gray) {
		cv::Mat mat(h, w, CV_8UC1, const_cast<uchar*>(data));
		cv::imshow(name, mat);
		//cv::waitKey(-1);
	} else {
		cv::Mat mat(h, w, CV_8UC3, const_cast<uchar*>(data));
		cv::Mat mat_rgb(h, w, CV_8UC3);
		cv::cvtColor(mat, mat_rgb, CV_BGR2RGB);

		cv::imshow(name, mat_rgb);
		//cv::waitKey(-1);
	}
}

void imshow(const char* name, const Arr3d& arr3) {
	Mat_d r, g, b;
	splitChannel(arr3, r, g, b);
	char wndName[256];
	sprintf(wndName, "%s_ch_01", name);
	imshow(wndName, r);
	sprintf(wndName, "%s_ch_02", name);
	imshow(wndName, g);
	sprintf(wndName, "%s_ch_03", name);
	imshow(wndName, b);
}

void imshowext(const char* name, const ImgRGB& img) {
	char* file_path = tmpnam(0);
	savePPM(img, file_path);
	char str[1024];
	sprintf(str, "eog -n %s\n", file_path);
	system(str);
}
void imshowext(const char* name, const ImgG& img) {
	char* file_path = tmpnam(0);
	savePGM(img, file_path);
	char str[1024];
	sprintf(str, "eog -n %s\n", file_path);
	system(str);
}
void imshowext(const char* name, const uchar* data, int w, int h, bool gray) {
	char* file_path = tmpnam(0);
	if (gray)
		savePGM(data, w, h, file_path);
	else
		savePPM(data, w, h, file_path);
	char str[1024];
	sprintf(str, "eog -n %s\n", file_path);
	system(str);
}
