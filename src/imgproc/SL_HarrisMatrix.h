/*
 * SL_HarrisMatrix.h
 *
 *  Created on: Feb 13, 2012
 *      Author: Danping Zou
 */

#ifndef SL_HARRISMATRIX_H_
#define SL_HARRISMATRIX_H_
#include "imgproc/SL_Image.h"
/**
 * img : w x h
 * harrisMat : (w x h) x 4 matrix 
 */
void computeHarrisMatrices(const ImgG& img, Mat_d& IxIx, Mat_d& IxIy,
		Mat_d& IyIy, double sigma);

void computeEigenValues(const Mat_d& IxIx, const Mat_d& IxIy, const Mat_d& IyIy,
		Mat_d& v1, Mat_d& v2, Mat_d& lambda1, Mat_d& lambda2);

void computeHarrisResponse(const Mat_d& IxIx, const Mat_d& IxIy,
		const Mat_d& IyIy, Mat_d& res, double kata);

void computeCornerEig(const ImgG& img, Mat_d& lambda1, Mat_d& lambda2,
		Mat_d& v1, Mat_d& v2, int block_size);

#endif /* SL_HARRISMATRIX_H_ */
