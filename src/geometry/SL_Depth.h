/*
 * SL_Depth.h
 *
 *  Created on: Feb 14, 2012
 *      Author: Danping Zou
 */

#ifndef SL_DEPTH_H_
#define SL_DEPTH_H_

void depth2point3d(double nx, double ny, double d, double M0[3]);
void get3DCoord(const double* K, const double* R, const double* t, int u, int v,
		double d, double M[3]);

#endif /* SL_DEPTH_H_ */
