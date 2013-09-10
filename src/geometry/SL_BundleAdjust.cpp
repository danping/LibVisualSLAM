/*
 * SL_BundleAdjust.cpp
 *
 *  Created on: 2011-7-12
 *      Author: Danping Zou
 */

#include "SL_BundleAdjust.h"
#include "SL_Geometry.h"
#include "SL_Triangulate.h"
#include "SL_Quaternion.h"

#include "math/SL_LinAlg.h"
#include "extern/sba-1.6/sba.h"

void bundleAdjust(int nconcam, const std::vector<Mat_d>& Ks,
		std::vector<Mat_d>& Rs, std::vector<Mat_d>& ts, int nconpts,
		std::vector<Point3d>& pts3d, std::vector<std::vector<Meas2D> >& pts2d,
		int nIter, double* info) {

	const int numCams = (int) Rs.size();
	const int numPts = (int) pts3d.size();

	const int cnp = 11;
	const int pnp = 3;
	const int mnp = 2;

	sbaGlobs glob;
	glob.cnp = cnp;
	glob.pnp = pnp;
	glob.mnp = mnp;

	Mat_d rotParam(numCams, FULLQUATSZ);
	glob.rot0params = rotParam.data;

	//set initial camera parameters
	for (int i = 0; i < numCams; ++i) {
		mat2quat(Rs[i].data, glob.rot0params + FULLQUATSZ * i);
	}

	glob.intrcalib = 0;
	glob.nccalib = 5;

	glob.camparams = 0;
	glob.ptparams = 0;

	/* call sparse LM routine */
	double opts[SBA_OPTSSZ];
	opts[0] = SBA_INIT_MU * 1E-4;
	opts[1] = SBA_STOP_THRESH;
	opts[2] = SBA_STOP_THRESH;
	opts[3] = 0; //0.05 * numMeas; // uncomment to force termination if the average reprojection error drops below 0.05
	opts[4] = 1E-16; // uncomment to force termination if the relative reduction in the RMS reprojection error drops below 1E-05

	Mat_d matParamVec(numCams * cnp + numPts * pnp, 1);
	double* paramVec = matParamVec.data;

	//set the remaining camera parameters (intrinsic parameters and translations)
	double* pParamVec = paramVec;
	for (int i = 0; i < numCams; ++i) {
		double* pKs = Ks[i].data;
		double* pTs = ts[i].data;
		pParamVec[0] = pKs[0];
		pParamVec[1] = pKs[2];
		pParamVec[2] = pKs[5];
		pParamVec[3] = pKs[4] / pKs[0];
		pParamVec[4] = pKs[1];

		pParamVec[5] = 0; //rotation
		pParamVec[6] = 0;
		pParamVec[7] = 0;
		pParamVec[8] = pTs[0]; //translation
		pParamVec[9] = pTs[1];
		pParamVec[10] = pTs[2];

		pParamVec += cnp;
	}

	//set the point parameters
	double* pParamPoints = paramVec + numCams * cnp;

	for (int i = 0; i < numPts; i++) {
		pParamPoints[0] = pts3d[i].x;
		pParamPoints[1] = pts3d[i].y;
		pParamPoints[2] = pts3d[i].z;
		pParamPoints += 3;
	}

	//set the 2d measurements and their weights
	Mat_d featPts(numPts * numCams, 2);
	Mat_c vmask(numPts, numCams);
	vmask.fill(0);

	//set weights
	Mat_d weights(numPts, numCams);

	int k = 0;
	for (int i = 0; i < numPts; i++) {
		int num2DPts = pts2d[i].size();
		for (int j = 0; j < num2DPts; j++) {
			featPts.data[2 * k] = pts2d[i][j].x;
			featPts.data[2 * k + 1] = pts2d[i][j].y;
			weights.data[k] = pts2d[i][j].w;
			k++;
			int viewId = pts2d[i][j].viewId;
			vmask.data[i * numCams + viewId] = 1;
		}
	}

	int maxIter = nIter;
	double sbaInfo[SBA_INFOSZ];
	if (sba_motstr_levmar_wx(numPts, nconpts, numCams, nconcam, vmask.data,
			paramVec, cnp, pnp, featPts.data, weights.data, mnp,
			img_projsKRTS_x, img_projsKRTS_jac_x, (void *) (&glob), maxIter, 0,
			opts, sbaInfo) == SBA_ERROR) {
		repErr("bundle adjustment failed!\n");
	}

	//test
	logInfo(
			"initial error:%lf, final error:%lf #iterations:%lf stop reason:%lf\n",
			sqrt(sbaInfo[0] / k), sqrt(sbaInfo[1] / k), sbaInfo[5], sbaInfo[6]);

	//update the camera poses
	pParamVec = paramVec;
	double dq[4], q[4];
	for (int i = 0; i < numCams; i++) {
		_MK_QUAT_FRM_VEC(dq, pParamVec + 5);
		quatMultFast(dq, glob.rot0params + i * FULLQUATSZ, q);
		quat2mat(q, Rs[i].data);
		memcpy(ts[i].data, pParamVec + 8, sizeof(double) * 3);
		pParamVec += cnp;
	}

	//update the map points
	pParamPoints = paramVec + numCams * cnp;
	for (int i = 0; i < numPts; i++) {
		memcpy(pts3d[i].M, pParamPoints, sizeof(double) * pnp);
		pParamPoints += pnp;
	}

	if (info) {
		memcpy(info, sbaInfo, sizeof(double) * SBA_INFOSZ);
	}
}
void bundleAdjustRobust(int nconcam, const std::vector<Mat_d>& Ks,
		std::vector<Mat_d>& Rs, std::vector<Mat_d>& ts, int nconpts,
		std::vector<Point3d>& pts3d, std::vector<std::vector<Meas2D> >& pts2d,
		int maxErr, int nIter, int nInnerIter) {

	double sbaInfo[SBA_INFOSZ];
	bundleAdjust(nconcam, Ks, Rs, ts, nconpts, pts3d, pts2d, nInnerIter,
			sbaInfo);
	for (int k = 1; k <= nIter; k++) {
		//compute the new weights according to the reprojection error
		for (size_t i = 0; i < pts2d.size(); i++) {
			for (size_t j = 0; j < pts2d[i].size(); j++) {
				int camId = pts2d[i][j].viewId;
				double* m = pts2d[i][j].m;

				double* K = Ks[camId].data;
				double* R = Rs[camId].data;
				double* t = ts[camId].data;

				double rm[2];
				project(K, R, t, pts3d[i].M, rm);

				double err = dist2(rm, m);

				if (err > maxErr) {
					pts2d[i][j].w = 0;
					pts2d[i][j].outlier = 1;
				} else {
					double s = err / maxErr;
					double w = (1 - s * s);
					pts2d[i][j].w = w * w;
					pts2d[i][j].outlier = 0;
				}
			}
		}
		bundleAdjust(nconcam, Ks, Rs, ts, nconpts, pts3d, pts2d, nInnerIter,
				sbaInfo);
	}
}
