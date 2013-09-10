/*
 * SL_MatchKeyPoint.h
 *
 *  Created on: 2010-11-9
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_MATCHING
#define SL_MATCHING

#include <vector>
#include <cassert>
#include "math/SL_Matrix.h"

class MatchItem {
public:
	int idx1, idx2;
	double dist;
};

class Matching {
protected:
	int _max_num;
public:
	MatchItem* data;
	int num;
public:
	Matching() :
			_max_num(0), data(0), num(0) {
	}
	~Matching() {
		clear();
	}
	bool empty() {
		return num == 0;
	}
	void clear() {
		if (data)
			delete[] data;
		_max_num = 0;
		num = 0;
		data = 0;
	}
	void reset() {
		num = 0;
		data = 0;
	}
	void clone(const Matching& other) {
		if (&other != this) {
			reserve(other.num);
			num = other.num;
			memcpy(data, other.data, sizeof(MatchItem) * num);
		}
	}
	void reserve(int n) {
		clear();
		if (n > 0) {
			data = new MatchItem[n];
			_max_num = n;
		}
	}
	void add(int i, int j, double dist) {
		assert(num < _max_num);
		data[num].idx1 = i;
		data[num].idx2 = j;
		data[num].dist = dist;
		num++;
	}
	MatchItem& operator[](int i) {
		return data[i];
	}
	MatchItem operator[](int i) const {
		return data[i];
	}
};
#endif /* SL_MATCHING */
