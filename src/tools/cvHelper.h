/*
 * cvHelper.h
 *
 *  Created on: 2011-1-6
 *      Author: Danping Zou
 */

#ifndef CVHELPER_H_
#define CVHELPER_H_

#include "math/SL_Matrix.h"
#include "imgproc/SL_Image.h"
#include "matching/SL_Matching.h"

#include "opencv2/features2d/features2d.hpp"
#include "opencv2/imgproc/imgproc.hpp"

typedef std::vector<cv::KeyPoint> KpVec;
typedef std::vector<cv::DMatch> DMatchVec;

void printCVMat(const cv::Mat& cvMat);

void mat2KpVec(const Mat_d& matPts, KpVec& vecPts);
void KpVec2Mat(const KpVec& vecPts, Mat_d& matPts);
void cvMat2MyMat(const cv::Mat& cvMat, Mat_d& mat);
void myMat2CvMat(const Mat_d& mat, cv::Mat& cvMat);
void cvMatch2MyMatch(const DMatchVec& cvMatch, Matching& myMatch);
void myMatch2CvMatch(const Matching& myMatch, DMatchVec& cvMatch);

void cvImg2ImgRGB(const cv::Mat& cvImg, ImgRGB& img);
#endif /* CVHELPER_H_ */
