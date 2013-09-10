/*
 * SL_SurfMatching.h
 *
 *  Created on: 2011-2-3
 *      Author: Danping Zou
 */

#ifndef SL_SURFMATCHING_H_
#define SL_SURFMATCHING_H_
#include "SL_Matching.h"

#include "math/SL_Matrix.h"
#include "imgproc/SL_Image.h"
#include "tools/cvHelper.h"

float computeSurfDescDist(const float* desc0, const float* desc1, int dimDesc);
void matchSurf(const Mat_f& kdesc0, const Mat_f& kdesc1, Matching& matches, float ratio);

void matchSurf(int dimDesc, std::vector<float>& desc1, std::vector<float>& desc2, Matching& matches, float ratio, float maxDist);
/**
 * refine the matching results by remove the correspondences with different directions from the majority
 */
int refineMatchedPoints(const Mat_d& pts1, const Mat_d& pts2, Matching& matches, Matching& newMatches, double ratio);
/**
 * detect SURF feature points on the image
 * @return the dimension of the descriptor
 */
int detectSURFPoints(const ImgG& img, Mat_d& surfPts, std::vector<float>& surfDesc, double hessianThreshold = 2000);
#endif /* SL_SURFMATCHING_H_ */
