/*
 * SL_DrawCorners.h
 *
 *  Created on: 2010-11-10
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_DRAWCORNERS_H_
#define SL_DRAWCORNERS_H_
#include "matching/SL_Matching.h"
#include "imgproc/SL_Image.h"
#include "geometry/SL_Triangulate.h"
#include "geometry/SL_Geometry.h"
#include <vector>

void drawSquare(ImgRGB& img, double x, double y, uchar r, uchar g, uchar b);
void drawSquare(ImgRGB& img, double x, double y, double hw, uchar r, uchar g, uchar b);
void drawCross(ImgRGB& img, double x, double y, uchar r, uchar g, uchar b);
void drawPoint(ImgRGB& img, double x, double y, uchar r, uchar g, uchar b);
void drawLine(ImgRGB& img, double x1, double y1, double x2, double y2, uchar r,
		uchar g, uchar b);
void drawCorners(ImgRGB& img, const Mat_i& corners, uchar r = 0, uchar g = 255,
		uchar b = 0, uchar sym = 's'); //'s' - square, 'c' - cross, '.' - point
void drawCorners(ImgRGB& img, const Mat_d& corners, uchar r = 0, uchar g = 255,
		uchar b = 0, uchar sym = 's');
void drawCorners(ImgRGB& img, const Mat_f& corners, uchar r = 0, uchar g = 255,
		uchar b = 0, uchar sym = 's');

#include "opencv2/opencv.hpp"
#include "opencv2/features2d/features2d.hpp"

void drawKeyPoints(const ImgG& img, const Mat_d& corners, ImgRGB& rgb, uchar r,
		uchar g, uchar b);
void drawKeyPoints(const cv::Mat& imgIn,
		const std::vector<cv::KeyPoint>& keyPts, cv::Mat& imgOut,
		cv::Scalar color, int sz = 3);

void drawClusteredMatches(const ImgG& img1, const Mat_d& keyPts1,
		const ImgG& img2, const Mat_d& keyPts2, const Mat_i& flags,
		CvScalar color_tab[], ImgRGB& outImg, double scale = 1.0);

void drawMatching(const cv::Mat& img1, const std::vector<cv::KeyPoint>& keyPts1,
		const cv::Mat& img2, const std::vector<cv::KeyPoint>& keyPts2,
		const std::vector<cv::DMatch>& matches, cv::Mat& outImg, double scale =
				1.0, unsigned char* flag = 0);

void drawMatching(const ImgG& img1, const Mat_d& keyPts1, const ImgG& img2,
		const Mat_d& keyPts2, const Matching& matches, ImgRGB& outImg,
		double scale = 1.0, unsigned char* flag = 0);

void drawMatching(const ImgG& img1, const Mat_d& keyPts1, const ImgG& img2,
		const Mat_d& keyPts2, ImgRGB& outImg, double scale = 1.0,
		unsigned char* flag = 0);

void drawMatching(const ImgRGB& img1, const Mat_d& keyPts1, const ImgRGB& img2,
		const Mat_d& keyPts2, const Matching& matches, unsigned char* flag = 0);

void drawReprojectedPoints(const ImgG& img, const double* K, const double* R,
		const double* T, int nPts, const double* pts3d, const double* pts2d,
		ImgRGB& outImg);

void drawReprojectedPoints(const ImgG& img, const double* K, const double* kc,
		const double* R, const double* T, int nPts, const double* pts3d,
		const double* pts2d, ImgRGB& outImg);

void drawPoint(const ImgG& img, const double x, const double y, ImgRGB& outImg,
		int radius, uchar r, uchar g, uchar b, double thickness = 1);

void drawPoint(ImgRGB& outImg, const double x, const double y, int radius,
		uchar r, uchar g, uchar b, double thickness = 1);
void drawLine(ImgRGB& outImg, double l[], uchar r, uchar g, uchar b);

void drawEllipseToData(ImgRGB& img, double x, double y, double var[4], uchar r,
		uchar g, uchar b);

void draw2DMesh(ImgRGB& outImg, const Mat_d& pts, const Mat_i& tri);

void drawReprojectionError(ImgRGB& rgb, const Mat_d& pts3d, const Mat_d& pts2d,
		const double K[], const double R[], const double T[]);
#endif /* SL_DRAWCORNERS_H_ */
