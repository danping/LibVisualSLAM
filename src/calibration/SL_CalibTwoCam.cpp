/*
 * SL_EMatEst.cpp
 *
 *  Created on: 2011-7-22
 *      Author: zou
 */

#include "SL_CalibTwoCam.h"
#include "matching/SL_SurfMatching.h"
#include "matching/SL_StereoMatcherHelper.h"
#include "geometry/SL_5point.h"
#include "geometry/SL_Triangulate.h"

bool findFMatrix(const Mat_d& pts1, const Mat_d& pts2, Mat_d& F, double maxEpiErr, int method) {
	assert( pts1.rows >0 && pts1.rows == pts2.rows);
	assert(method == FF_RANSAC || method == FF_LMEDS);
	int npts = pts1.rows;
	CvMat* cvPts1 = cvCreateMat(npts, 2, CV_64FC1);
	CvMat* cvPts2 = cvCreateMat(npts, 2, CV_64FC1);

	for (int i = 0; i < npts; i++) {
		cvPts1->data.db[2 * i] = pts2.data[2 * i];
		cvPts1->data.db[2 * i + 1] = pts2.data[2 * i + 1];
		cvPts2->data.db[2 * i] = pts1.data[2 * i];
		cvPts2->data.db[2 * i + 1] = pts1.data[2 * i + 1];
	}

	CvMat* cvF = cvCreateMat(3, 3, CV_64FC1);

	if (cvFindFundamentalMat(cvPts1, cvPts2, cvF, method, maxEpiErr) == 0) {
		cvReleaseMat(&cvPts1);
		cvReleaseMat(&cvPts2);
		cvReleaseMat(&cvF);
		return false;
	}

	F.resize(3, 3);
	memcpy(F.data, cvF->data.db, sizeof(double) * 9);

	cvReleaseMat(&cvPts1);
	cvReleaseMat(&cvPts2);
	cvReleaseMat(&cvF);
	return true;
}

int getInlierFMatches(const Mat_d& pts1, const Mat_d& pts2, const double F[9], double maxEpiErr, Mat_uc& inlierFlags) {
	assert( pts1.rows >0 && pts1.rows == pts2.rows);
	int npts = pts1.rows;
	inlierFlags.resize(npts, 1);

	int k = 0;
	for (int i = 0; i < npts; i++) {
		double err = epipolarError(F, pts1.data + 2 * i, pts2.data + 2 * i);
		if (err > maxEpiErr) {
			inlierFlags.data[i] = 0;
		} else {
			inlierFlags.data[i] = 1;
			k++;
		}
	}
	//test
	logInfo("getInlierFMatches - total number :%d , inlier number :%d\n", npts, k);
	return k;
}

bool refineFMatrix(const Mat_d& pts1, const Mat_d& pts2, const Mat_uc& inlierFlags, Mat_d& F) {
	assert( pts1.rows >0 && pts1.rows == pts2.rows);
	assert( pts1.rows == inlierFlags.rows);

	int npts = 0;
	for (int i = 0; i < inlierFlags.rows; i++)
		if (inlierFlags.data[i] > 0)
			npts++;

	CvMat* cvPts1 = cvCreateMat(npts, 2, CV_64FC1);
	CvMat* cvPts2 = cvCreateMat(npts, 2, CV_64FC1);

	int k = 0;
	for (int i = 0; i < pts1.rows; i++) {
		if (inlierFlags.data[i] > 0) {
			cvPts1->data.db[2 * k] = pts2.data[2 * k];
			cvPts1->data.db[2 * k + 1] = pts2.data[2 * k + 1];
			cvPts2->data.db[2 * k] = pts1.data[2 * k];
			cvPts2->data.db[2 * k + 1] = pts1.data[2 * k + 1];
			k++;
		}
	}

	CvMat* cvF = cvCreateMat(3, 3, CV_64FC1);

	if (cvFindFundamentalMat(cvPts1, cvPts2, cvF, CV_FM_8POINT) == 0) {
		cvReleaseMat(&cvPts1);
		cvReleaseMat(&cvPts2);
		cvReleaseMat(&cvF);
		return false;
	}

	F.resize(3, 3);
	memcpy(F.data, cvF->data.db, sizeof(double) * 9);

	cvReleaseMat(&cvPts1);
	cvReleaseMat(&cvPts2);
	cvReleaseMat(&cvF);
	return true;
}
CalibTwoCam::CalibTwoCam() :
		K1(0), K2(0), kc1(0), kc2(0) {
	// TODO Auto-generated constructor stub

}

CalibTwoCam::~CalibTwoCam() {
	// TODO Auto-generated destructor stub
}
void CalibTwoCam::matchSURFPoints(const ImgG& img1, const ImgG& img2) {
	assert(!img1.empty() && !img2.empty());
	assert(K1 && K2);

	Mat_d surfPts1, surfPts2;
	std::vector<float> surfDesc1, surfDesc2;

	int dim = detectSURFPoints(img1, surfPts1, surfDesc1, hessianThreshold);
	detectSURFPoints(img2, surfPts2, surfDesc2, hessianThreshold);

	Matching surfMatchesOrg;
	matchSurf(dim, surfDesc1, surfDesc2, surfMatchesOrg, (float) surfRatio, (float) surfMaxDiff);

	//get matched points;
	Mat_d pts1, pts2;
	getMatchedPts(surfMatchesOrg, surfPts1, surfPts2, pts1, pts2);

	setMatchedPoints(pts1, pts2);
}
void CalibTwoCam::setMatchedPoints(const Mat_d& pts1, const Mat_d& pts2) {
	assert(K1 && K2);
	orgPts1.cloneFrom(pts1);
	orgPts2.cloneFrom(pts2);

	normPts1.resize(orgPts1.m, orgPts1.n);
	normPts2.resize(orgPts2.m, orgPts2.n);

	if (kc1) {
		undistorNormPoints(iK1, ikc1, orgPts1.m, orgPts1.data, normPts1.data);
		imagePoints(K1, orgPts1.m, normPts1.data, orgPts1.data);
	} else
		normPoints(iK1, orgPts1.m, orgPts1.data, normPts1.data);

	if (kc2) {
		undistorNormPoints(iK2, ikc2, orgPts2.m, orgPts2.data, normPts2.data);
		imagePoints(K2, orgPts2.m, normPts2.data, orgPts2.data);
	} else
		normPoints(iK2, orgPts2.m, orgPts2.data, normPts2.data);
}
void CalibTwoCam::getInlierInd(std::vector<int>& ind) {
	for (int i = 0; i < inlierFlag.m; i++) {
		if (inlierFlag.data[i] > 0)
			ind.push_back(i);
	}
}
int CalibTwoCam::estimateEMatOld(double maxEpiErr, int method) {
	assert(iK1 && iK2);
	int n = orgPts1.rows;
	if (n < minInlierNum)
		repErr("estimateEmat - too little matches");

	//find the fundamental matrix
	Mat_d matF;
	findFMatrix(orgPts1, orgPts2, matF, maxEpiErr, method);
	int nInliers = getInlierFMatches(orgPts1, orgPts2, matF.data, maxEpiErr, inlierFlag);
	getEMat(K1, K2, matF.data, E);
	return nInliers;
}
int CalibTwoCam::estimateEMat(double maxEpiErr){
	assert(iK1 && iK2);
	int n = orgPts1.rows;
	assert(n > 0);
	if (n < minInlierNum)
		repErr("estimateEmat - too little matches");

	Mat_d normPts1(n,2),normPts2(n,2);
	normPoints(iK1,n,orgPts1,normPts1);
	normPoints(iK2,n,orgPts2,normPts2);

	findEMatRansac(iK2, iK1, n, orgPts2.data, orgPts1.data, normPts2.data,normPts1.data, E, nRansac, maxEpiErr);
	
	Mat_d matF(3, 3);
	getFMatK(K2, K1, E,matF.data);
	memcpy(F,matF.data,sizeof(double)*9);

	int nInliers = getInlierFMatches(orgPts2, orgPts1, matF, maxEpiErr, inlierFlag);
	return nInliers;
}
void CalibTwoCam::outputInlierPoints(Mat_d& pts1, Mat_d& pts2) {
	assert(!orgPts1.empty() && !orgPts2.empty());
	assert(orgPts1.m == orgPts2.m);
	assert(!inlierFlag.empty() && inlierFlag.m == orgPts1.m);
	getFlaged2DPoints(orgPts1, inlierFlag, pts1);
	getFlaged2DPoints(orgPts2, inlierFlag, pts2);
}
void CalibTwoCam::outputInlierNormPoints(Mat_d& pts1norm, Mat_d& pts2norm) {
	assert(!orgPts1.empty() && !orgPts2.empty());
	assert(orgPts1.m == orgPts2.m);
	assert(!inlierFlag.empty() && inlierFlag.m == orgPts1.m);
	getFlaged2DPoints(normPts1, inlierFlag, pts1norm);
	getFlaged2DPoints(normPts2, inlierFlag, pts2norm);
}
void CalibTwoCam::outputRTs(Mat_d& R1, Mat_d& t1, Mat_d& R2, Mat_d& t2, bool all) {
	R2.resize(3, 3);
	t2.resize(3, 1);

	matEyes(3, R1);
	matZeros(3, 1, t1);

	Mat_d normInPts1, normInPts2;
	getFlaged2DPoints(normPts1, inlierFlag, normInPts1);
	getFlaged2DPoints(normPts2, inlierFlag, normInPts2);
	int nPts = normInPts1.m;

	double Rs[36], ts[12];
	decompEMat(E, Rs, ts);
	if( !all)
		selectFacingRT(Rs, ts, nPts, normInPts1.data, normInPts2.data, R2.data, t2.data);	
	else{
		R2.resize(4,9);
		t2.resize(4,3);

		R2.copyFrom(Rs);
		t2.copyFrom(ts);
	}
}

////test finding fundamental matrix
//int main(int argc, char** argv) {
//	ImgG img1, img2;
//	imread(img1, "/home/zou/images/cam_0.pgm");
//	imread(img2, "/home/zou/images/cam_1.pgm");
//
//	Mat_d surfPts1, surfPts2;
//	std::vector<float> surfDesc1, surfDesc2;
//
//	tic();
//	double hessianThreshold = 100;
//	double surfRatio = 0.8;
//	double surfMaxDiff = 0.6;
//	int dim = detectSURFPoints(img1, surfPts1, surfDesc1, hessianThreshold);
//	detectSURFPoints(img2, surfPts2, surfDesc2, hessianThreshold);
//	toc();
//
//	Matching surfMatchesOrg;
//	matchSurf(dim, surfDesc1, surfDesc2, surfMatchesOrg, surfRatio, surfMaxDiff);
//
//	//get matched points;
//	Mat_d pts1, pts2;
//	getMatchedPts(surfMatchesOrg, surfPts1, surfPts2, pts1, pts2);
//
//	Mat_d F;
//	findFMatrix(pts1, pts2, F, FF_LMEDS);
//
//	Mat_uc inlierFlags;
//	getInlierFMatches(pts1, pts2, F, 3.0, inlierFlags);
//
//	ImgRGB outImg;
//	drawMatching(img1, pts1, img2, pts2, outImg, 1.0, inlierFlags.data);
//
//	imshow("inlier matches", outImg);
//	cv::waitKey(-1);
//	print(F);
//	return 0;
//}

//int main(int argc, char** argv) {
//	ImgG img1, img2;
//	imread(img1, "/home/zou/record_07123/0460.png");
//	imread(img2, "/home/zou/record_07123/0690.png");
//
//	CalibTwoCam ematEst;
//	double K[9] = { 575.81575, 0, 320.240, 0, 575.81575, 320.240, 0, 0, 1 };
//	ematEst.setIntrinParam(K, K);
//	ematEst.matchSURFPoints(img1, img2);
//	ematEst.estimateEMat();
//
//	ImgRGB outImg;
//	drawMatching(img1, ematEst.orgPts1, img2, ematEst.orgPts2, outImg, 1.0, ematEst.inlierFlag.data);
//	imshow("test", outImg);
//	cv::waitKey(-1);
//
//	Mat_d R1, t1, R2, t2, pts1, pts2, x;
//	ematEst.outputInlierPoints(pts1, pts2);
//	ematEst.outputRTs(R1, t1, R2, t2);
//
//	int npts = pts1.m;
//	x.resize(npts, 3);
//	binTriangulatePoints(K, R1, t1, K, R2, t2, pts1.m, pts1.data, pts2.data, x.data);
//
//	writeMat(R1, "/home/zou/R1.txt");
//	writeMat(t1, "/home/zou/t1.txt");
//	writeMat(R2, "/home/zou/R2.txt");
//	writeMat(t2, "/home/zou/t2.txt");
//	writeMat(pts1, "/home/zou/pts1.txt");
//	writeMat(pts2, "/home/zou/pts2.txt");
//	writeMat(x, "/home/zou/x.txt");
//
//	printMat(3, 3, ematEst.E);
//
//	return 0;
//}s
//int main( int argc, char** argv){
//	using namespace std;
//	
//	cv::SURF surf(4000,4,2);
//	cv::Mat rgb = cv::imread("/home/zou/woods.jpg");
//	cv::Mat gray;
//	cv::cvtColor(rgb,gray,CV_RGB2GRAY);
//	
//	vector<cv::KeyPoint> points;
//	vector<float> desc;
//	tic();
//	#pragma omp parallel for num_threads(4)
//	for( int i = 0; i < 4; i++)
//		surf(gray,cv::Mat(),points,desc);
//	toc();
//	return 0;
//}

//int main(int argc, char** argv) {
//	ImgG img1, img2;
//	imread(img1, "/home/zou/images/cam_1.pgm");
//	imread(img2, "/home/zou/images/cam_2.pgm");
//
//	imshow("test1", img1);
//	imshow("test2", img2);
//	cv::waitKey(-1);
//	CalibTwoCam ematEst;
//	double K[9] = { 613.63471, 0, 399.5, 0, 613.63471, 239.5, 0, 0, 1 };
//	ematEst.setIntrinParam(K, K);
//	ematEst.matchSURFPoints(img1, img2);
//	ematEst.estimateEMat();
//
//	ImgRGB outImg, outImg0;
//	drawMatching(img1, ematEst.orgPts1, img2, ematEst.orgPts2, outImg, 1.0, 0);
//	imshow("test0", outImg);
//	drawMatching(img1, ematEst.orgPts1, img2, ematEst.orgPts2, outImg, 1.0, ematEst.inlierFlag.data);
//	imshow("test", outImg);
//	cv::waitKey(-1);
//
//	Mat_d R1, t1, R2, t2, pts1, pts2, x;
//	ematEst.outputInlierPoints(pts1, pts2);
//	ematEst.outputRTs(R1, t1, R2, t2);
//
//	int npts = pts1.m;
//	x.resize(npts, 3);
//	binTriangulatePoints(K, R1, t1, K, R2, t2, pts1.m, pts1.data, pts2.data, x.data);
//
//	writeMat(3, 3, ematEst.E, "/home/zou/E.txt");
//	writeMat(R1, "/home/zou/R1.txt");
//	writeMat(t1, "/home/zou/t1.txt");
//	writeMat(R2, "/home/zou/R2.txt");
//	writeMat(t2, "/home/zou/t2.txt");
//	writeMat(pts1, "/home/zou/pts1.txt");
//	writeMat(pts2, "/home/zou/pts2.txt");
//	writeMat(x, "/home/zou/x.txt");
//
//	printMat(3, 3, ematEst.E);
//
//	return 0;
//}

