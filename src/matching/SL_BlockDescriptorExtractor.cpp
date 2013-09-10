/*
 * BlockDescriptorExtractor.cpp
 *
 *  Created on: 2011-1-8
 *      Author: Danping Zou
 */

#include "SL_BlockDescriptorExtractor.h"
#include "opencv2/imgproc/imgproc.hpp"
BlockDescriptorExtractor::BlockDescriptorExtractor() {

}

BlockDescriptorExtractor::~BlockDescriptorExtractor() {

}

void BlockDescriptorExtractor::compute(const ImgG& img , const Mat_d& pts , Mat_d& desc) {
	if (pts.empty())
		repErr("BlockDescriptorExtractor::compute - no point is provided");
	
	cv::Size patchSize(2 * hW + 1, 2 * hW + 1);
	int nPts = pts.rows;
	int len = (2 * hW + 1) * (2 * hW + 1);
	desc.resize(nPts, len);

	cv::Mat imgRef(img.rows, img.cols, CV_8UC1, img.data);

	if (scale == 1.0) {
		for (int i = 0; i < nPts; i++) {
			double x = pts.data[2 * i];
			double y = pts.data[2 * i + 1];

			cv::Mat patch;
			cv::getRectSubPix(imgRef, patchSize, cv::Point2f(x, y), patch);
			//cv::Scalar avg = cv::mean(patch);
			CvMat patchRef = patch;

			for (int j = 0; j < len; j++)
				//desc.data[i * len + j] = (patchRef.data.ptr[j] - avg(0));
				desc.data[i * len + j] = patchRef.data.ptr[j];
		}
	} else {
		cv::Mat tmpImg;
		cv::resize(imgRef, tmpImg, cv::Size(), scale, scale);
		for (int i = 0; i < nPts; i++) {
			double x = pts.data[2 * i] * scale;
			double y = pts.data[2 * i + 1] * scale;

			cv::Mat patch;
			cv::getRectSubPix(tmpImg, patchSize, cv::Point2f(x, y), patch);
			//cv::Scalar avg = cv::mean(patch);
			CvMat patchRef = patch;

			for (int j = 0; j < len; j++)
				//desc.data[i * len + j] = (patchRef.data.ptr[j] - avg(0));
				desc.data[i * len + j] = patchRef.data.ptr[j];
		}
	}
}

