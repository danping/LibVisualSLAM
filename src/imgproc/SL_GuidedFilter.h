/*
 * SL_GuidedFilter.h
 *
 *  Created on: Dec 18, 2011
 *      Author: Danping Zou
 */

#ifndef SL_GUIDEDFILTER_H_
#define SL_GUIDEDFILTER_H_
#include "SL_Image.h"
#include "SL_ImageOp.h"
#include "SL_BoxFilter.h"
#include "math/SL_Matrix.h"
#include "math/SL_LinAlgWarper.h"
#include "tools/SL_TypeConversion.h"
#include <cassert>

template<class T, class U>
void guidedFilt(const MyMat<T>& I, const MyMat<U>& p, MyMat<U>& q, int r,
		double eps) {
	assert(I.m == p.m && I.n == p.n && I.m > 0 && I.n > 0);
	assert(r > 0);
	q.resize(p.m, p.n);

	Mat_i ones(p.m, p.n), num(p.m, p.n);
	ones.fill(1);
	boxFilt2(ones, num, r);

	/* compute I*p*/
	MyMat<U> Ip(p.m, p.n);
	matEleProd(p.m, p.n, I.data, p.data, Ip.data);

	MyMat<U> mI, mp, mIp;
	/* compute mean I,p, Ip*/
	fastBoxFilt2(I, mI, r);
	fastBoxFilt2(p, mp, r);
	fastBoxFilt2(Ip, mIp, r);

	matEleDiv(p.m, p.n, mI.data, num.data, mI.data);
	matEleDiv(p.m, p.n, mp.data, num.data, mp.data);
	matEleDiv(p.m, p.n, mIp.data, num.data, mIp.data);

	/* compute var I*/
	MyMat<U> dI(p.m, p.n);
	matSub(p.m, p.n, I.data, mI.data, dI.data);
	matEleProd(p.m, p.n, dI.data, dI.data, dI.data);

	MyMat<U> sI;
	fastBoxFilt2(dI, sI, r);
	matEleDiv(p.m, p.n, sI.data, num.data, sI.data);

	MyMat<U> a(p.m, p.n), b(p.m, p.n);
	q.resize(p.m, p.n);

	int len = p.m * p.n;
	for (int i = 0; i < len; ++i) {
		a.data[i] = (mIp.data[i] - mI.data[i] * mp.data[i])
				/ (sI.data[i] + eps);
		b.data[i] = mp.data[i] - a.data[i] * mI.data[i];
	}

	MyMat<U> ma, mb;
	fastBoxFilt2(a, ma, r);
	fastBoxFilt2(b, mb, r);

	matEleDiv(p.m, p.n, ma.data, num.data, ma.data);
	matEleDiv(p.m, p.n, mb.data, num.data, mb.data);

	for (int i = 0; i < len; ++i)
		q.data[i] = a.data[i] * I.data[i] + b.data[i];

}
template<class RGB, class U>
void guidedFiltColor(const RGB& rgb, const U* p, U* q, int r, double eps) {
	int h = rgb.h;
	int w = rgb.w;
	Mat_i ones(h, w), num(h, w);
	ones.fill(1);
	fastBoxFilt2(ones, num, r);

	MyMat<U> I[3];
	splitImageChannel(rgb, I[0], I[1], I[2]);

	/* compute I*p */
	MyMat<U> Ip[3];
	for (int i = 0; i < 3; ++i) {
		Ip[i].resize(h, w);
		matEleProd(h, w, I[i].data, p, Ip[i].data);
	}

	MyMat<U> mI[3], mIp[3], mp(h, w);

	/* compute means of I, Ip, and p */
	for (int i = 0; i < 3; ++i) {
		fastBoxFilt2(I[i], mI[i], r);
		fastBoxFilt2(Ip[i], mIp[i], r);
		matEleDiv(h, w, mI[i].data, num.data, mI[i].data);
		matEleDiv(h, w, mIp[i].data, num.data, mIp[i].data);
	}
	fastBoxFilt2(h, w, p, mp.data, r);
	matEleDiv(h, w, mp.data, num.data, mp.data);

	int len = h * w;

	MyMat<U> covIp(len, 3);
	for (int i = 0; i < len; ++i) {
		covIp.data[i * 3] = mIp[0].data[i] - mI[0].data[i] * mp.data[i];
		covIp.data[i * 3 + 1] = mIp[1].data[i] - mI[1].data[i] * mp.data[i];
		covIp.data[i * 3 + 2] = mIp[2].data[i] - mI[2].data[i] * mp.data[i];
	}

	/* compute Sigma*/
	MyMat<U> dI[3];
	for (int i = 0; i < 3; ++i) {
		dI[i].resize(h, w);
		matSub(h, w, I[i].data, mI[i].data, dI[i].data);
	}
	int ind[6][2] =
			{ { 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 1 }, { 1, 2 }, { 2, 2 } };

	MyMat<U> II[6];
	for (int i = 0; i < 6; ++i) {
		II[i].resize(h, w);
		matEleProd(h, w, I[ind[i][0]].data, I[ind[i][1]].data, II[i].data);
	}

	MyMat<U> mII[6];
	for (int i = 0; i < 6; ++i) {
		mII[i].resize(h, w);
		fastBoxFilt2(II[i], mII[i], r);
		matEleDiv(h, w, mII[i].data, num.data, mII[i].data);
	}

	for (int i = 0; i < 6; ++i) {
		for (int k = 0; k < len; ++k)
			mII[i].data[k] -= mI[ind[i][0]].data[k] * mI[ind[i][1]].data[k];
	}

	MyMat<U> A[3], mA[3], B(h, w), mB(h, w);
	for (int i = 0; i < 3; ++i) {
		A[i].resize(h, w);
	}

	for (int i = 0; i < len; ++i) {
		double sigma[9] = { mII[0].data[i] + eps, mII[1].data[i],
				mII[2].data[i], mII[1].data[i], mII[3].data[i] + eps,
				mII[4].data[i], mII[2].data[i], mII[4].data[i], mII[5].data[i]
						+ eps };

		double isigma[9], a[3];
		mat33Inv(sigma, isigma);
		matAB(3, 3, 3, 1, isigma, covIp.data + 3 * i, a);

		B.data[i] = mp.data[i] - a[0] * mI[0].data[i] - a[1] * mI[1].data[i]
				- a[2] * mI[2].data[i];

		A[0].data[i] = a[0];
		A[1].data[i] = a[1];
		A[2].data[i] = a[2];
	}

	for (int i = 0; i < 3; ++i) 
		fastBoxFilt2(A[i], mA[i], r);
	fastBoxFilt2(B, mB, r);

	for (int i = 0; i < len; ++i) {
		q[i] = mA[0].data[i] * I[0].data[i] + mA[1].data[i] * I[1].data[i]
				+ mA[2].data[i] * I[2].data[i] + mB.data[i];
	}
	matEleDiv(h, w, q, num.data, q);
}

template<class RGB, class U>
void guidedFiltColor(const RGB& rgb, const MyMat<U>& p, MyMat<U>& q, int r,
		double eps) {
	assert(rgb.m == p.m && rgb.n == p.n && rgb.m > 0 && rgb.n > 0);
	assert(r > 0);
	int w = rgb.w;
	int h = rgb.h;
	q.resize(h, w);
	guidedFiltColor(rgb, p.data, q.data, r, eps);
}
#endif /* SL_GUIDEDFILTER_H_ */
