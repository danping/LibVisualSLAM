/*
 * SL_BundleAdjust.h
 *
 *  Created on: 2011-7-12
 *      Author: Danping Zou
 */

#ifndef SL_BUNDLEADJUST_H_
#define SL_BUNDLEADJUST_H_
#include "math/SL_Matrix.h"
#include "SL_Point.h"
#include "SL_BundleHelper.h"

#include <utility>
#include <vector>

class Meas2D {
public:
	int viewId;
	union {
		struct {
			double x, y;
		};
		double m[2];
	};
	double w;
	double cov[4];
	int outlier;
public:
	Meas2D() {
		cov[0] = cov[1] = cov[2] = cov[3] = 0;
	}
	Meas2D(int camId, double fx, double fy) :
			viewId(camId), x(fx), y(fy), w(1.0), outlier(0) {
		cov[0] = cov[1] = cov[2] = cov[3] = 0;
	}
};

inline bool operator <(const Meas2D& item1, const Meas2D& item2) {
	return item1.viewId < item2.viewId;
}

void bundleAdjust(int nconcam, const std::vector<Mat_d>& Ks, /*intrinsic matrices*/
std::vector<Mat_d>& Rs, /*rotation matrices*/
std::vector<Mat_d>& ts, /*translations*/
int nconpts, std::vector<Point3d>& pts3d, /*3D points*/
std::vector<std::vector<Meas2D> >& pts2d, /*2D projections*/
int nIter = 50, double* info = 0);
/*
 * use M-estimator (Tukey)
 */
void bundleAdjustRobust(int nconcam, const std::vector<Mat_d>& Ks,
		std::vector<Mat_d>& Rs, std::vector<Mat_d>& ts, int nconpts,
		std::vector<Point3d>& pts3d, std::vector<std::vector<Meas2D> >& pts2d,
		int maxErr = 3, int nIter = 10, int nInnerIter = 50);

#endif /* SL_BUNDLEADJUST_H_ */
