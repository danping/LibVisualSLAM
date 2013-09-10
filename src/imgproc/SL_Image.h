/*
 * SL_Image.h
 *
 *  Created on: 2010-11-8
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_IMAGE_H_
#define SL_IMAGE_H_

typedef unsigned char uchar;
#include <cstring>
class ImgG {
public:
	union {
		int n;
		int w;
		int cols;
	};
	union {
		int m;
		int h;
		int rows;
	};
	uchar* data;
public:
	ImgG();
	ImgG(int w, int h);
	~ImgG();

public:
	void clear();
	void fill(uchar val) {
		memset(data, val, w * h * sizeof(uchar));
	}
	void resize(int w, int h);
	bool empty() const {
		return w == 0 || h == 0;
	}

public:
	typedef uchar* PUCHAR;
	typedef const uchar* CONST_PUCHAR;
	operator PUCHAR() {
		return data;
	}
	operator CONST_PUCHAR() const {
		return data;
	}
	uchar* getPtr() const {
		return data;
	}
	uchar& operator()(int y, int x) const {
		return data[y * w + x];
	}
	int getMemLen() const {
		return w * h;
	}
};

class ImgRGB {
public:
	union {
		int w;
		int n;
		int cols;
	};
	union {
		int h;
		int m;
		int rows;
	};
	uchar* data;
public:
	ImgRGB();
	ImgRGB(int w, int h);
	~ImgRGB();

public:
	void clear();
	void fill(uchar r, uchar g, uchar b);
	void resize(int w, int h);
	void split(ImgG& R, ImgG& G, ImgG& B) const;
	bool empty() const {
		return w == 0 || h == 0;
	}
public:
	typedef uchar* PUCHAR;
	typedef const uchar* CONST_PUCHAR;
	operator PUCHAR() {
		return data;
	}
	operator CONST_PUCHAR() const {
		return data;
	}
	uchar* get_ptr() const {
		return data;
	}
	uchar* operator ()(int y, int x) const {
		return data + y * w * 3 + x * 3;
	}
	int getMemLen() const {
		return h * w * 3;
	}
	void getStdAvg(
			int startx,
			int endx,
			int starty,
			int endy,
			double avg[3],
			double std[3]) const;

	uchar* getRowPtr(int y) const {
		return data + y * w * 3;
	}
};

#include "SL_FloatImage.h"
typedef FloatImgG<double> ImgG_d;
typedef FloatImgG<float> ImgG_f;
typedef FloatImgRGB<double> ImgRGB_d;
typedef FloatImgRGB<float> ImgRGB_f;
#endif /* SL_IMAGE_H_ */
