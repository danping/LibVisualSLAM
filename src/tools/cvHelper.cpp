/*
 * cvHelper.cpp
 *
 *  Created on: 2011-1-6
 *      Author: Danping Zou
 */
#include "cvHelper.h"
#include "SL_error.h"
void printCVMat(const cv::Mat& cvMat) {
	logInfo("===============[opencv matrix] ---------------\n");
	logInfo("%dx%d [%d]\n", cvMat.rows, cvMat.cols, cvMat.channels());
	if (cvMat.channels() == 1) {
		for (int i = 0; i < cvMat.rows; i++) {
			for (int j = 0; j < cvMat.cols; j++) {
				switch (cvMat.type()) {
				case CV_8U:
					logInfo("%d ", cvMat.at<unsigned char> (i, j));
					break;
				case CV_8S:
					logInfo("%d ", cvMat.at<char> (i, j));
					break;
				case CV_16U:
					logInfo("%d ", cvMat.at<unsigned short> (i, j));
					break;
				case CV_16S:
					logInfo("%d ", cvMat.at<short> (i, j));
					break;
				case CV_32S:
					logInfo("%d ", cvMat.at<int> (i, j));
					break;

				case CV_32F:
					logInfo("%lf ", cvMat.at<float> (i, j));
					break;
				case CV_64F:
					logInfo("%lf ", cvMat.at<double> (i, j));
					break;
				}
			}
			logInfo("\n");
		}
	}
	logInfo("-----------------------=======================\n");

}
void mat2KpVec(const Mat_d& matPts , KpVec& vecPts) {
	int numPts = matPts.rows;
	vecPts.reserve(numPts);
	for (int i = 0; i < numPts; i++) {
		double x = matPts.data[2 * i];
		double y = matPts.data[2 * i + 1];

		cv::KeyPoint pt(x, y, 1);
		vecPts.push_back(pt);
	}
}
void KpVec2Mat(const KpVec& vecPts , Mat_d& matPts) {
	int numPts = vecPts.size();
	matPts.resize(numPts, 2);
	for (int i = 0; i < numPts; i++) {
		matPts.data[2 * i] = vecPts[i].pt.x;
		matPts.data[2 * i + 1] = vecPts[i].pt.y;
	}
}
void cvMat2MyMat(const cv::Mat& cvMat , Mat_d& mat) {
	mat.resize(cvMat.rows, cvMat.cols);
	if (cvMat.type() == CV_64F) {
		CvMat tmpMat = cvMat;
		memcpy(mat.data, tmpMat.data.db, sizeof(double) * cvMat.rows * cvMat.cols);
	} else {
		CvMat tmpMat = cvMat;
		int len = cvMat.rows * cvMat.cols;
		for (int i = 0; i < len; i++) {
			mat.data[i] = tmpMat.data.fl[i];
		}
	}
}

void myMat2CvMat(const Mat_d& mat , cv::Mat& cvMat) {
	cvMat.create(mat.rows, mat.cols, CV_64F);
	CvMat tmpMat = cvMat;
	memcpy(tmpMat.data.db, mat.data, sizeof(double) * mat.rows * mat.cols);
}

void cvMatch2MyMatch(const DMatchVec& cvMatch , Matching& myMatch) {
	int numMatch = cvMatch.size();
	myMatch.reserve(numMatch);
	for (int i = 0; i < numMatch; i++) {
		myMatch.add(cvMatch[i].queryIdx, cvMatch[i].trainIdx, cvMatch[i].distance);
	}
}
void myMatch2CvMatch(const Matching& myMatch , DMatchVec& cvMatch) {
	int numMatch = myMatch.num;
	cvMatch.clear();
	cvMatch.reserve(numMatch * 2);
	for (int i = 0; i < numMatch; i++) {
		cvMatch.push_back(cv::DMatch(myMatch[i].idx1, myMatch[i].idx2, myMatch[i].dist));
	}
}
void cvImg2ImgRGB(const cv::Mat& cvImg , ImgRGB& img) {
	int M = cvImg.rows;
	int N = cvImg.cols;

	img.resize(N, M);

	CvMat tmpMat = cvImg;

	memcpy(img.data, tmpMat.data.ptr, M * N * 3);
	//	for( int i = 0; i < M; i++){
	//		for( int j = 0; j < N; j++){
	//			img.data[3*(i*N+j)] = cvImg.at<cv::Vec3b>(i,j)[0];
	//			img.data[3*(i*N+j)+1] = cvImg.at<cv::Vec3b>(i,j)[1];
	//			img.data[3*(i*N+j)+2] = cvImg.at<cv::Vec3b>(i,j)[2];
	//			                                               
	//		}
	//	}
}
//int main() {
//	Mat_d mat;
//	matRand(5, 5, mat);
//	print(mat);
//	cv::Mat cvMat;
//	myMat2CvMat(mat, cvMat);
//	printCVMat(cvMat);
//	return 0;
//}
