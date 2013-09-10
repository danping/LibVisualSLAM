/*
 * SL_Debug.cpp
 *
 *  Created on: 2010-11-6
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#include "SL_Debug.h"
#include "SL_error.h"

void draw_reprojection_points(
		const double* K ,
		const double* kc ,
		int v_id ,
		int nv ,
		const double* Rs ,
		const double* ts ,
		int npts ,
		const double* pts ,
		const double* meas ,
		const char* vmask) {
	const double* R = Rs + v_id * 9;
	const double* t = ts + v_id * 3;

	ImgRGB img(800, 450);
	img.fill(0, 0, 0);
	Mat_d rp_n(npts, 2);
	Mat_d mms_n(npts, 2);

	Mat_d rp(npts, 2);
	Mat_d mms(npts, 2);

	Mat_i indices(npts, nv);
	indices.fill(-1);

	int k = 0;
	for (int i = 0; i < npts * nv; i++) {
		if (vmask[i] == 1)
			indices.data[i] = k++;
	}

	k = 0;
	double m[2];

	Mat_i A(npts, 1);
	for (int i = 0; i < npts; i++) {
		if (vmask[i * nv + v_id] == 0)
			continue;

		project(K, kc, R, t, pts + 3 * i, m);
		rp_n.data[2 * k] = m[0];
		rp_n.data[2 * k + 1] = m[1];

		int ind = indices[i * nv + v_id];
		const double* pmeas = meas + 2 * ind;

		mms_n.data[2 * k] = pmeas[0];
		mms_n.data[2 * k + 1] = pmeas[1];

		A.data[k] = i;
		k++;
	}

	rp_n.rows = k;
	mms_n.rows = k;
	//	ms.rows = k;
	mms.rows = k;

	//	imagePoints(K, kc, npts, ms_n, ms);
	imagePoints(K, kc, npts, mms_n, mms);
	drawCorners(img, rp_n, 0, 0, 255, 's');
	drawCorners(img, mms, 255, 0, 0, 'c');

	for (int i = 0; i < k; i++) {
		//draw_cross(img, rp_n.data[idx1 * 2], rp_n.data[idx1 * 2 + 1], 255, 255, 0);
		drawLine(img, rp_n.data[i * 2], rp_n.data[i * 2 + 1], mms.data[i * 2], mms.data[i * 2 + 1], 255, 255, 0);
		double dist = dist2(rp_n.data + 2 * i, mms.data + 2 * i);
		//test
		int id = A.data[i];
		logInfo("[id:%d](%lf,%lf,%lf) -> (%lf,%lf) - %lf - measure:(%lf,%lf)\n", id, pts[3 * id], pts[3 * id + 1], pts[3
				* id + 2], rp_n[2 * i], rp_n[2 * i + 1], dist, mms[2 * i], mms[2 * i + 1]);
	}

	char buf[256];
	sprintf(buf, "%d", v_id);
	imshow(buf, img);
}
