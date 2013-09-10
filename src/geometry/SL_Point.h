/*
 * point.h
 *
 *  Created on: 2010-11-12
 *      Author: Danping Zou
 */

#ifndef POINT_H_
#define POINT_H_

#include <vector>
typedef long long longInt;
using namespace std;

template<class T>
class Point2 {
public:
	union {
		struct {
			T x, y;
		};
		T m[2];
	};
public:
	Point2(const Point2& other) :
			x(other.x), y(other.y) {
	}
	Point2() {
	}
	Point2(T x0, T y0) {
		set(x0, y0);
	}
	Point2& operator =(const Point2& other) {
		if (&other != this) {
			x = other.x;
			y = other.y;
		}
		return *this;
	}
	void set(double x0, double y0) {
		m[0] = x0;
		m[1] = y0;
	}
};

typedef Point2<double> Point2d;
typedef Point2<float> Point2f;
typedef Point2<int> Point2i;

template<class T>
inline bool Point2CompareLess(const Point2<T>& pt1, const Point2<T>& pt2) {
	if (pt1.y < pt2.y)
		return true;
	if (pt1.y == pt2.y) {
		if (pt1.x <= pt2.x)
			return true;
	}
	return false;
}

template<class T>
class Point3 {
public:
	union {
		struct {
			T x, y, z;
		};
		T M[3];
	};
public:
	Point3() {
	}
	Point3(const Point3& other) :
			x(other.x), y(other.y), z(other.z){
	}
	Point3(T x0, T y0, T z0) {
		set(x0, y0, z0);
	}

	Point3& operator =(const Point3& other) {
		if (&other == this) {
			return *this;
		}
		x = other.x;
		y = other.y;
		z = other.z;
		return *this;
	}

	void set(T x0, T y0, T z0) {
		M[0] = x0;
		M[1] = y0;
		M[2] = z0;
	}
};

typedef Point3<double> Point3d;
typedef Point3<float> Point3f;
typedef Point3<int> Point3i;

class Point3dId: public Point3d {
public:
	longInt id;
public:
	Point3dId() :
			id(-1) {
	}
	Point3dId(double x0, double y0, double z0, int id0) :
			Point3(x0, y0, z0), id(id0) {
	}
};
typedef vector<Point2d> VecPoint2d;
typedef vector<Point3d> VecPoint3d;

#endif /* POINT_H_ */
