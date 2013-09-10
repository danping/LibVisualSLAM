/*
 * SL_5point.h
 *
 *  Created on: 2010-11-17
 *      Author: Danping Zou
 */

#ifndef SL_5POINT_H_
#define SL_5POINT_H_

/*randomly choose m elements in an n-element array*/
void randChoose(int n, int* indices, int m);
/* compute fundamental matrix 
 * F = inv(K1)^T E inv(K2)
 */
void getFMat(const double* invK1, const double* invK2, const double* E, double* F);
/* compute the essential matrix
 * E = K1^T*F*K2
 */
void getEMat(const double* K1, const double* K2, const double* F, double* E);
void getFMatK(const double* K1, const double* K2, const double *E , double* F);
void computeEpipolarLine(const double x, const double y, const double* F, double l[3]);
void computeEpipolarLine(const double* F, const double x, double y, double l[3]);

int computeEMat5Pt(int n, const double* a, const double* b, double* E);

int evaluateEMat(const int n,
		const double * a,
		const double * b,
		const double* E,
		double epiErrThes,
		double& epiErr,
		unsigned char* flag = 0);

/* find essential matrix by using RANSAC algorithm (using five point algorithm)
 *  a^T * E * b = 0
 * input:
 * 	n : number of correspondence
 *  a, b : normalized points
 *  iterMax : maximum number of RNASAC steps
 *  epiErrThres : threshold of epipolar error
 * output:
 * 	E : essential matrix
 *  epiErr : sum of epipolar error
 */
int findEMatRansac(const int n,
		const double* a,
		const double* b,
		double* E,
		int iterMax,
		double epiErrThes,
		double* epiErr = 0);

/* find essential matrix by using RANSAC algorithm (using five point algorithm)
 *  an^T*E*bn = 0
 * input:
 *	 	invK : 3x3 inverse of intrinsic matrices of two views
 * 		invKd  : 1x7 inverse distortion parameters
 * 		n : number of correspondence
 * 	    a, b : corresponding feature points
 * 		an, bn : normalized points
 * output:
 * 		E: essential matrix
 */
int findEMatRansac(const double* invKa,
		const double* invKb,
		int n,
		const double* a,
		const double* b,
		const double* an,
		const double* bn,
		double* E,
		int iterMax,
		double epiErrThres,
		double* epiErr = 0);

double findBestEMat(int n, const double* a, const double *b, const double* E, int nE, double *EBest, int& iMin);

/*get camera pose from the essential matrix
 * input :
 *	 E: Essential matrix
 * 	 n : number of point pairs
 * 	 p1,p2 : corresponding feature points
 * output:
 *   R, t : camera pose
 */
int getCamPoseFromEMat(double *E, int n, double *p1, double *p2, double *R, double *t);

/* use five point algorithm to simultaneously estimate the camera pose and point coordinates
 * input: 
 * 	invK : 3x3 inverse of intrinsic matrices of two views
 * 	invKd  : 1x7 inverse distortion parameters
 * 	nc : number of correspondences
 * 	a,b :  ncx2 array to store the corresponding points
 * 	thres : threshold for ransec algorithm
 * 	npts : number of reconstructed 3d points
 * 	pts  : npts x 3 array to store the map points
 * 
 * output:
 * 	return the reprojection error
 */
int fivePointReconstruct(const double* invK,
		const double* invKd,
		int nc,
		const double* a,
		const double* b,
		double thres,
		double* R,
		double *t,
		double* pts,
		int* inlierIdx = 0,
		int ransacIterNum = 20);

void getInlier2DPoints(int nInliers, const int* inlierIdx, const double* pts, double* inPts);

double computeMeanError(int n, const double* pts, const double* pts1);
/** generate essential matrix from the camera poses*/
void formEMat(const double* R0, const double* T0, const double* R1, const double* T1, double* E);

/* decompose the four camera matrices from the essential matrix*/
/* E = [t]xR*/
void decompEMat(const double E[9], double R[/*4x9*/], double t[/*4x3*/]);



/* select the camera poses where two cameras are facing the scene points*/
int selectFacingRT(const double Rs[/*4x9*/],
		const double ts[/*4x3*/],
		int nPts,
		const double* pts1,
		const double* pts2,
		double R1[9],
		double t1[9]);
#endif /* SL_5POINT_H_ */
