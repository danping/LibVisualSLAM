/*
 * SL_Image.cpp
 *
 *  Created on: 2010-11-8
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#include "SL_Image.h"
#include <algorithm>
ImgRGB::ImgRGB() :
		w(0), h(0), data(0) {

}
ImgRGB::ImgRGB(int w_, int h_) :
		w(0), h(0), data(0) {
	if (h_ > 0 && w_ > 0) {
		w = w_;
		h = h_;
		data = new uchar[w * h * 3];
	}
}
ImgRGB::~ImgRGB() {
	clear();
}

void ImgRGB::clear() {
	if (data) {
		delete[] data;
		data = 0;
	}
	w = 0;
	h = 0;
}
void ImgRGB::fill(uchar r, uchar g, uchar b) {
	int i, len = w * h;
	uchar* pdata = data;
	for (i = 0; i < len; i++) {
		pdata[0] = r;
		pdata[1] = g;
		pdata[2] = b;
		pdata += 3;
	}
}
void ImgRGB::resize(int wd, int ht) {
	clear();
	if (ht > 0 && wd > 0) {
		w = wd;
		h = ht;
		data = new uchar[h * w * 3];
	}
}
void ImgRGB::split(ImgG& R, ImgG& G, ImgG& B) const {
	R.resize(w, h);
	G.resize(w, h);
	B.resize(w, h);

	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			R.data[i * w + j] = data[3 * i * w + 3 * j];
			G.data[i * w + j] = data[3 * i * w + 3 * j + 1];
			B.data[i * w + j] = data[3 * i * w + 3 * j + 2];
		}
	}
}
#include "tools/SL_Print.h"
void ImgRGB::getStdAvg(int startx, int endx, int start_row, int end_row,
		double avg[3], double std[3]) const {
	std::fill_n(avg, 3, 0.0);
	for (int r = start_row; r < end_row && r < h; ++r) {
		for (int c = startx; c < endx && c < w; ++c) {
			avg[0] += data[3 * r * w + 3 * c];
			avg[1] += data[3 * r * w + 3 * c + 1];
			avg[2] += data[3 * r * w + 3 * c + 2];
		}
	}
	int N = (end_row - start_row) * w;
	avg[0] /= N;
	avg[1] /= N;
	avg[2] /= N;
	
	std::fill_n(std, 3, 0.0);
	for (int r = start_row; r < end_row && r < h; ++r) {
		for (int c = startx;  c < endx && c < w; ++c) {
			double dR = data[3 * r * w + 3 * c] - avg[0];
			double dG = data[3 * r * w + 3 * c + 1] - avg[1];
			double dB = data[3 * r * w + 3 * c + 2] - avg[2];
			std[0] += dR * dR;
			std[1] += dG * dG;
			std[2] += dB * dB;
		}
	}
	std[0] /= N;
	std[1] /= N;
	std[2] /= N;
}
ImgG::ImgG() :
		w(0), h(0), data(0) {

}
ImgG::ImgG(int w_, int h_) :
		w(w_), h(h_), data(0) {
	if (h_ > 0 && w_ > 0) {
		w = w_;
		h = h_;
		data = new uchar[h * w];
	} else {
		w = 0;
		h = 0;
		data = 0;
	}
}
ImgG::~ImgG() {
	clear();
}
void ImgG::clear() {
	if (data) {
		delete[] data;
		data = 0;
	}
	w = 0;
	h = 0;
}
void ImgG::resize(int wd, int ht) {
	if (wd == w && ht == h)
		return;
	clear();
	if (ht > 0 && wd > 0) {
		w = wd;
		h = ht;
		data = 0;
		data = new uchar[h * w];
	}
}