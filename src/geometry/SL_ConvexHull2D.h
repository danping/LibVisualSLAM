/*
 * SL_ConvexHull2D.h
 *
 *  Created on: 2011-5-10
 *      Author: Danping Zou
 */

#ifndef SL_CONVEXHULL2D_H_
#define SL_CONVEXHULL2D_H_

#include <cstdlib>
#include <cstdio>
#include <vector>
typedef double coord;
void get2DConvexHull(const std::vector<double>& pts, std::vector<double>& cxh);
double getPolyArea(const std::vector<double>& poly);
#endif /* SL_CONVEXHULL2D_H_ */
