/*
 * SL_Distortion.h
 *
 *  Created on: 2010-11-6
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_DISTORTION_H_
#define SL_DISTORTION_H_

/*get the max range of the normalized image plane
 * input:
 * 	w, h : image size
 * 	iK : inverse of the intrinsic matrix
 * output:
 * 	return the maximum radius 
 */
double maxRangeNormPlane(int w, int h, const double* iK);
/* compute the polynomial for inverting the distortion
 *  r : maximum range of image points
 *  kc : distortion parameters (5x1) 
 *  k_ud : undistortion polynomial (stored as a 7x1 vector)
 */
void invDistorParam(double r, const double *kc, double *k_ud);
/* compute the polynomial for inverting the distortion
 *  w, h : image size
 *  kc : distortion parameters (5x1) 
 *  k_ud : undistortion polynomial (stored as a 7x1 vector)
 */
void invDistorParam(int w, int h, const double* iK, const double* kc, double* k_ud);
/* get the normalized image points with distortion removal
 * iK : inverse of intrinsic matrix
 * k_ud : undistortion parameters
 * n : number of points
 * a : nx2 array to store the image points
 * an : nx2 array to store the normalized image points
 */
void undistorNormPoints(const double* iK, const double* k_ud, int n, const double* a, double* an);
void undistorPoint(const double* K, const double* k_ud, const double* pt, double * undistPt);

/*get the normalize image points pn = K^-1 p*/
void normPoints(const double* iK, int n, const double* a, double *an);
void normPoint(const double* iK, const double* a, double *an);
void normPoint(const double* iK, const double x, const double y, double& xNorm, double& yNorm);

void getInvK(const double* K, double* invK);
/*get the image points from the normalized ones*/
void imagePoints(const double* K, int n, const double* an, double* a);
void imagePoint(const double* K, const double* an, double* a);
void imagePoints(const double* K, const double* kc, int n, const double* an, double* a);
#endif /* SL_DISTORTION_H_ */
