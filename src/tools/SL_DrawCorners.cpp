/*
 * SL_DrawCorners.cpp
 *
 *  Created on: 2010-11-10
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#include "SL_DrawCorners.h"
#include "SL_TypeConversion.h"
#include "SL_Print.h"
#include "cvHelper.h"

void drawSquare(ImgRGB& img, double x, double y, uchar r, uchar g, uchar b) {
	drawLine(img, x - 3, y - 3, x + 3, y - 3, r, g, b);
	drawLine(img, x - 3, y + 3, x + 3, y + 3, r, g, b);
	drawLine(img, x - 3, y - 3, x - 3, y + 3, r, g, b);
	drawLine(img, x + 3, y - 3, x + 3, y + 3, r, g, b);
}
void drawSquare(ImgRGB& img, double x, double y, double hw, uchar r, uchar g,
		uchar b) {
	drawLine(img, x - hw, y - hw, x + hw, y - hw, r, g, b);
	drawLine(img, x - hw, y + hw, x + hw, y + hw, r, g, b);
	drawLine(img, x - hw, y - hw, x - hw, y + hw, r, g, b);
	drawLine(img, x + hw, y - hw, x + hw, y + hw, r, g, b);
}
void drawCross(ImgRGB& img, double x, double y, uchar r, uchar g, uchar b) {
	drawLine(img, x - 3, y, x + 3, y, r, g, b);
	drawLine(img, x , y - 3, x , y + 3, r, g, b);
}
void drawPoint(ImgRGB& img, double x, double y, uchar r, uchar g, uchar b) {
	int i = static_cast<int>(y + 0.5);
	int j = static_cast<int>(x + 0.5);

	int w = img.w * 3;
	img[i * w + 3 * j] = r;
	img[i * w + 3 * j + 1] = g;
	img[i * w + 3 * j + 2] = b;
}
void drawLine(ImgRGB& img, double x0, double y0, double x1, double y1, uchar r,
		uchar g, uchar b) {

	cv::Mat cvImg(img.h, img.w, CV_8UC3, img.data);
	cv::line(cvImg, cv::Point2d(x0, y0), cv::Point2d(x1, y1),cv::Scalar(r,g,b));
}
void drawCorners(ImgRGB& img_show, const Mat_i& corners, uchar r, uchar g,
		uchar b, uchar sym) {
	int nc = corners.rows;
	//draw corners
	int i;
	for (i = 0; i < nc; ++i) {
		int x0 = corners.data[2 * i];
		int y0 = corners.data[2 * i + 1];

		switch (sym) {
		case 'c':
			drawCross(img_show, x0, y0, r, g, b);
			break;
		case 's':
			drawSquare(img_show, x0, y0, r, g, b);
			break;
		case '.':
			drawPoint(img_show, x0, y0, r, g, b);
			break;
		}
	}
}

void drawCorners(ImgRGB& img_show, const Mat_d& corners, uchar r, uchar g,
		uchar b, uchar sym) {
	Mat_i corners_i;
	matDouble2Int(corners, corners_i);
	drawCorners(img_show, corners_i, r, g, b, sym);
}
void drawCorners(ImgRGB& img_show, const Mat_f& corners, uchar r, uchar g,
		uchar b, uchar sym) {
	Mat_i corners_i;
	matFloat2Int(corners, corners_i);
	drawCorners(img_show, corners_i, r, g, b, sym);
}

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"

void drawMatching(const cv::Mat& img1, const std::vector<cv::KeyPoint>& keyPts1,
		const cv::Mat& img2, const std::vector<cv::KeyPoint>& keyPts2,
		const std::vector<cv::DMatch>& matches, cv::Mat& outImg, double scale,
		unsigned char* flag) {
	int W0 = img1.cols > img2.cols ? img1.cols : img2.cols;
	int H0 = img1.rows + img2.rows;

	cv::Size sz(W0, H0);
	cv::Mat tmpOutImg;
	tmpOutImg.create(sz, CV_MAKETYPE(img1.depth(), 3));

	cv::Mat outImg1 = tmpOutImg(cv::Rect(0, 0, img1.cols, img1.rows));
	cv::Mat outImg2 = tmpOutImg(cv::Rect(0, img1.rows, img2.cols, img2.rows));

	if (img1.type() == CV_8U)
		cv::cvtColor(img1, outImg1, CV_GRAY2BGR);
	else
		img1.copyTo(outImg1);

	if (img2.type() == CV_8U)
		cv::cvtColor(img2, outImg2, CV_GRAY2BGR);
	else
		img2.copyTo(outImg2);

	cv::resize(tmpOutImg, outImg, cv::Size(), scale, scale);

	//cv::RNG& rng = cv::theRNG();
	cv::RNG rng;

	//draw matched key points and lines between them
	for (size_t i = 0; i < matches.size(); i++) {
		if (flag && flag[i] == 0)
			continue;
		int id1 = matches[i].queryIdx;
		int id2 = matches[i].trainIdx;
		cv::Scalar color = cv::Scalar(rng(256), rng(256), rng(256));

		cv::Point2f pt1 = keyPts1[id1].pt;
		cv::Point2f pt2 = keyPts2[id2].pt;

		pt1.x *= scale;
		pt1.y *= scale;

		pt2.x *= scale;
		pt2.y *= scale;

		pt2.y += img1.rows * scale;

		cv::circle(outImg, pt1, 6, color);
		cv::circle(outImg, pt2, 6, color);
		cv::line(outImg, pt1, pt2, color);
	}
}

void drawMatching(const ImgG& img1, const Mat_d& keyPts1, const ImgG& img2,
		const Mat_d& keyPts2, const Matching& matches, ImgRGB& outImg,
		double scale, unsigned char* flag) {

	cv::Mat tmpImg1(img1.rows, img1.cols, CV_8UC1, img1.data);
	cv::Mat tmpImg2(img2.rows, img2.cols, CV_8UC1, img2.data);

	KpVec tmpKeyPts1, tmpKeyPts2;
	mat2KpVec(keyPts1, tmpKeyPts1);
	mat2KpVec(keyPts2, tmpKeyPts2);

	DMatchVec tmpMatches;
	myMatch2CvMatch(matches, tmpMatches);

	cv::Mat cvImg;
	drawMatching(tmpImg1, tmpKeyPts1, tmpImg2, tmpKeyPts2, tmpMatches, cvImg,
			scale, flag);

	cvImg2ImgRGB(cvImg, outImg);
}
void drawClusteredMatches(const ImgG& img1, const Mat_d& keyPts1,
		const ImgG& img2, const Mat_d& keyPts2, const Mat_i& flags,
		CvScalar color_tab[], ImgRGB& outImg, double scale) {
	assert(keyPts1.rows == keyPts2.rows);
	int npts = keyPts1.rows;
	assert(npts > 0);

	int W0 = img1.cols > img2.cols ? img1.cols : img2.cols;
	int H0 = img1.rows + img2.rows;

	cv::Size sz(W0, H0);
	cv::Mat tmpOutImg;
	tmpOutImg.create(sz, CV_8UC3);

	cv::Mat outImg1 = tmpOutImg(cv::Rect(0, 0, img1.cols, img1.rows));
	cv::Mat outImg2 = tmpOutImg(cv::Rect(0, img1.rows, img2.cols, img2.rows));

	cv::Mat cvImg1(img1.rows, img1.cols, CV_8UC1, img1.data);
	cv::cvtColor(cvImg1, outImg1, CV_GRAY2BGR);

	cv::Mat cvImg2(img2.rows, img2.cols, CV_8UC1, img2.data);
	cv::cvtColor(cvImg2, outImg2, CV_GRAY2BGR);

	cv::Mat cvOutImg;
	cv::resize(tmpOutImg, cvOutImg, cv::Size(), scale, scale);

	//cv::RNG& rng = cv::theRNG();
	cv::RNG rng;

	//draw matched key points and lines between them
	for (int i = 0; i < npts; i++) {
		int iC = flags[i];
		cv::Scalar color = cv::Scalar(color_tab[iC]);

		cv::Point2f pt1 = cv::Point2f(keyPts1.data[2 * i],
				keyPts1.data[2 * i + 1]);
		cv::Point2f pt2 = cv::Point2f(keyPts2.data[2 * i],
				keyPts2.data[2 * i + 1]);

		pt1.x *= (float)scale;
		pt1.y *= (float)scale;

		pt2.x *= (float)scale;
		pt2.y *= (float)scale;

		pt2.y += img1.rows * (float)scale;

		cv::circle(cvOutImg, pt1, 6, color);
		cv::circle(cvOutImg, pt2, 6, color);
		cv::line(cvOutImg, pt1, pt2, color);
	}

	IplImage img = cvOutImg;
	outImg.resize(img.width, img.height);
	memcpy(outImg.data, img.imageData,
			3 * sizeof(char) * img.width * img.height);
}

void drawMatching(const ImgG& img1, const Mat_d& keyPts1, const ImgG& img2,
		const Mat_d& keyPts2, ImgRGB& outImg, double scale,
		unsigned char* flag) {

	int numPts = keyPts1.rows;
	Matching matches;
	matches.reserve(numPts);

	cv::RNG rng;

	for (int i = 0; i < numPts; i++) {
		matches.add(i, i, 0);
	}
	drawMatching(img1, keyPts1, img2, keyPts2, matches, outImg, scale, flag);
}

void drawMatching(const ImgRGB& img1, const Mat_d& keyPts1, const ImgRGB& img2,
		const Mat_d& keyPts2, const Matching& matches, unsigned char* flag) {

	cv::Mat cvImg1(img1.rows, img1.cols, CV_8UC3, img1.data);
	cv::Mat cvImg2(img2.rows, img2.cols, CV_8UC3, img2.data);

	cv::RNG rng;

	for (int i = 0; i < matches.num; i++) {
		if (flag && flag[i] == 0)
			continue;
		int idx1 = matches[i].idx1;
		int idx2 = matches[i].idx2;

		double x1 = keyPts1.data[2 * idx1];
		double y1 = keyPts1.data[2 * idx1 + 1];

		double x2 = keyPts2.data[2 * idx2];
		double y2 = keyPts2.data[2 * idx2 + 1];

		cv::Scalar color = cv::Scalar(rng(256), rng(256), rng(256));

		cv::circle(cvImg1, cv::Point2d(x1, y1), 6, color, 1, CV_AA);
		cv::circle(cvImg2, cv::Point2d(x2, y2), 6, color, 1, CV_AA);
	}
}
void drawKeyPoints(const ImgG& img, const Mat_d& corners, ImgRGB& rgb, uchar r,
		uchar g, uchar b) {
	int nPts = corners.rows;
	if (nPts == 0)
		return;

	cv::Mat cvImg(img.rows, img.cols, CV_8UC1, img.data);
	cv::Mat tmpImg(img.rows, img.cols, CV_8UC3);

	cv::cvtColor(cvImg, tmpImg, CV_GRAY2RGB);

	std::vector<cv::KeyPoint> cvPts;
	cvPts.reserve(nPts);
	for (int i = 0; i < nPts; i++) {
		double x = corners.data[2 * i];
		double y = corners.data[2 * i + 1];
		cv::circle(tmpImg, cv::Point(x, y), 3, cv::Scalar(r, g, b), 1, CV_AA);
	}

	CvMat cvtmpImg = tmpImg;
	rgb.resize(img.cols, img.rows);
	memcpy(rgb.data, cvtmpImg.data.ptr, 3 * img.rows * img.cols);
}

void drawKeyPoints(const cv::Mat& imgIn,
		const std::vector<cv::KeyPoint>& keyPts, cv::Mat& imgOut,
		cv::Scalar color, int sz) {
	if (imgOut.rows != imgIn.rows || imgOut.cols != imgIn.cols)
		imgOut.create(imgIn.rows, imgIn.cols, CV_8UC3);
	if (imgIn.type() == CV_8U) {
		cv::cvtColor(imgIn, imgOut, CV_GRAY2BGR);
	} else
		imgIn.copyTo(imgOut);

	size_t nPts = keyPts.size();
	for (size_t i = 0; i < nPts; i++) {
		cv::circle(imgOut, keyPts[i].pt, sz, color, 1, CV_AA);
	}
}

#include "geometry/SL_Triangulate.h"
void drawReprojectedPoints(const ImgG& img, const double* K, const double* R,
		const double* T, int nPts, const double* pts3d, const double* pts2d,
		ImgRGB& outImg) {

	outImg.resize(img.w, img.h);
	cv::Mat cvOutImg(img.rows, img.cols, CV_8UC3, outImg.data);
	cv::Mat cvImg(img.rows, img.cols, CV_8U, img.data);

	cv::cvtColor(cvImg, cvOutImg, CV_GRAY2RGB);

	Mat_d tmpPts2d(nPts, 2);
	project(K, R, T, nPts, pts3d, tmpPts2d);

	for (int i = 0; i < nPts; i++) {
		double x1 = tmpPts2d[2 * i];
		double y1 = tmpPts2d[2 * i + 1];

		double x0 = pts2d[2 * i];
		double y0 = pts2d[2 * i + 1];

		//cv::circle(cvOutImg, cv::Point2f(x0, y0), 3, cv::Scalar(0, 255, 0), 1, CV_AA);
		//cv::circle(cvOutImg, cv::Point2f(x1, y1), 1, cv::Scalar(0, 0, 255), 1, CV_AA);
		cv::line(cvOutImg, cv::Point2f(x1, y1), cv::Point2f(x0, y0),
				cv::Scalar(255, 0, 0), 1, CV_AA);
	}
}

void drawReprojectedPoints(const ImgG& img, const double* K, const double* kc,
		const double* R, const double* T, int nPts, const double* pts3d,
		const double* pts2d, ImgRGB& outImg) {

	outImg.resize(img.w, img.h);
	cv::Mat cvOutImg(img.rows, img.cols, CV_8UC3, outImg.data);
	cv::Mat cvImg(img.rows, img.cols, CV_8U, img.data);

	cv::cvtColor(cvImg, cvOutImg, CV_GRAY2RGB);

	Mat_d tmpPts2d(nPts, 2);
	project(K, kc, R, T, nPts, pts3d, tmpPts2d);

	for (int i = 0; i < nPts; i++) {
		double x1 = tmpPts2d[2 * i];
		double y1 = tmpPts2d[2 * i + 1];

		double x0 = pts2d[2 * i];
		double y0 = pts2d[2 * i + 1];

		cv::circle(cvOutImg, cv::Point2f(x1, y1), 3, cv::Scalar(0, 255, 0), 1,
				CV_AA);
		cv::circle(cvOutImg, cv::Point2f(x0, y0), 1, cv::Scalar(0, 0, 255), 1,
				CV_AA);
		cv::line(cvOutImg, cv::Point2f(x1, y1), cv::Point2f(x0, y0),
				cv::Scalar(255, 0, 0), 1, CV_AA);
	}
}

static void getEndPoints(double l[], int h, int w, double& x1, double& y1,
		double& x2, double& y2) {
	double a = l[0];
	double b = l[1];
	double c = l[2];

	if (fabs(a) <= fabs(b)) {
		x1 = 0;
		x2 = w;

		y1 = -c / b;
		y2 = -a / b * w - c / b;
	} else {
		y1 = 0;
		y2 = h;

		x1 = -c / a;
		x2 = -b / a * h - c / a;
	}
}
void drawPoint(const ImgG& img, const double x, const double y, ImgRGB& outImg,
		int radius, uchar r, uchar g, uchar b, double thickness) {
	outImg.resize(img.w, img.h);
	cv::Mat cvImg(img.h, img.w, CV_8U, img.data);
	cv::Mat cvOutImg(img.h, img.w, CV_8UC3, outImg.data);
	cv::cvtColor(cvImg, cvOutImg, CV_GRAY2RGB);
	cv::circle(cvOutImg, cv::Point2f(x, y), radius, cv::Scalar(r, g, b),
			thickness, CV_AA);
}
void drawPoint(ImgRGB& outImg, const double x, const double y, int radius,
		uchar r, uchar g, uchar b, double thickness) {
	cv::Mat cvOutImg(outImg.h, outImg.w, CV_8UC3, outImg.data);
	cv::circle(cvOutImg, cv::Point2f(x, y), radius, cv::Scalar(r, g, b),
			thickness, CV_AA);
}
void drawLine(ImgRGB& outImg, double l[], uchar r, uchar g, uchar b) {
	double x1, y1, x2, y2;
	getEndPoints(l, outImg.h - 1, outImg.w - 1, x1, y1, x2, y2);
	cv::Mat tmpImg(outImg.h, outImg.w, CV_8UC3, outImg.data);
	cv::line(tmpImg, cv::Point2f(x1, y1), cv::Point2f(x2, y2),
			cv::Scalar(r, g, b));
}

#include "math/SL_LinAlg.h"
void drawEllipseToData(ImgRGB& img, double x, double y, double var[4], uchar r,
		uchar g, uchar b) {

	double evec[4], eval[2];
	dgeevFor(2, var, evec, eval);

	//	int l1 = static_cast<int> (eval[0] + 0.5);
	//	int l2 = static_cast<int> (eval[1] + 0.5);
	double l1 = sqrt(eval[0]);
	double l2 = sqrt(eval[1]);

	double dx1 = l1 * evec[0];
	double dy1 = l1 * evec[1];
	double dx2 = l2 * evec[2];
	double dy2 = l2 * evec[3];

	cv::Mat cvOutImg(img.rows, img.cols, CV_8UC3, img.data);
	cv::line(cvOutImg, cv::Point2f(x - dx1, y - dy1),
			cv::Point2f(x + dx1, y + dy1), cv::Scalar(r, g, b), 1, CV_AA);
	cv::line(cvOutImg, cv::Point2f(x - dx2, y - dy2),
			cv::Point2f(x + dx2, y + dy2), cv::Scalar(r, g, b), 1, CV_AA);
	//get rotation angle
	double theta = atan2(evec[0], evec[1]);
	cv::ellipse(cvOutImg, cv::Point2f(x, y), cv::Size(l1, l2), theta, 0, 360,
			cv::Scalar(r, g, b), 1, CV_AA);
	//	double ivar[4];
	//	mat22Inv(var,ivar);

}
void draw2DMesh(ImgRGB& outImg, const Mat_d& pts, const Mat_i& tri) {
	assert(!outImg.empty());
	cv::Mat cvOutImg(outImg.rows, outImg.cols, CV_8UC3, outImg.data);
	int thickness = 1;
	cv::Scalar color = cv::Scalar(0, 255, 0);
	for (int i = 0; i < tri.m; i++) {
		int a = tri.data[3 * i];
		int b = tri.data[3 * i + 1];
		int c = tri.data[3 * i + 2];


		cv::line(cvOutImg, cv::Point2d(pts.data[2 * a], pts.data[2 * a + 1]),
				cv::Point2d(pts.data[2 * b], pts.data[2 * b + 1]), color,
				thickness, CV_AA, 0);
		cv::line(cvOutImg, cv::Point2d(pts.data[2 * b], pts.data[2 * b + 1]),
				cv::Point2d(pts.data[2 * c], pts.data[2 * c + 1]), color,
				thickness, CV_AA, 0);
		cv::line(cvOutImg, cv::Point2d(pts.data[2 * c], pts.data[2 * c + 1]),
				cv::Point2d(pts.data[2 * a], pts.data[2 * a + 1]), color,
				thickness, CV_AA, 0);
	}
}
void drawReprojectionError(ImgRGB& rgb, const Mat_d& pts3d, const Mat_d& pts2d,
		const double K[], const double R[], const double T[]) {
	int npts = pts3d.m;
	for (int i = 0; i < npts; i++) {
		double x = pts2d.data[2 * i];
		double y = pts2d.data[2 * i + 1];
		double rm[2];
		project(K, R, T, pts3d.data + 3 * i, rm);
		if (dist2(x, y, rm[0], rm[1]) > 5.0) {
			drawPoint(rgb, x, y, 5, 255, 0, 0, 1);
		} else
			drawPoint(rgb, x, y, 5, 0, 255, 0, 1);
		drawLine(rgb, x, y, rm[0], rm[1], 255, 0, 0);
	}
}