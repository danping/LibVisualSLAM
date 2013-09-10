/*
 * SL_WriteRead.h
 *
 *  Created on: 2010-11-13
 *      Author: Danping Zou
 */

#ifndef SL_WRITEREAD_H_
#define SL_WRITEREAD_H_
#include "SL_error.h"
#include "geometry/SL_Point.h"
#include "math/SL_Matrix.h"

#include "matching/SL_Matching.h"

#include <string>
#include <fstream>
#include <vector>
#include <sstream>

/* points IO */
void write(const VecPoint2d& pts, const char* fmtstr, ...);
void read(VecPoint2d& pts, const char* fmtstr, ...);

template<class Mat>
void writeMat(const Mat& mat, const char* fmtstr, ...) {
	char filePath[1024];
	GET_FMT_STR(fmtstr, filePath);
	using namespace std;
	ofstream file(filePath);
	if (!file)
		repErr("write - cannot open '%s'!\n", filePath);
	file.precision(24);
	for (int i = 0; i < mat.rows; i++) {
		for (int j = 0; j < mat.cols; j++) {
			file << mat.data[i * mat.cols + j] << " ";
		}
		file << endl;
	}
	file.close();
	logInfo("save '%s' [ok]\n", filePath);
}
template<class Mat>
void readMat(Mat& mat, const char* fmtstr, ...) {
	char filePath[1024];
	GET_FMT_STR(fmtstr, filePath);
	using namespace std;

	ifstream file(filePath);
	if (!file)
		repErr("read - cannot open '%s'!\n", filePath);

	typedef typename Mat::value_type T;
	vector<vector<T> > data;
	string line;

	int m = 0, n = 0;
	while (!getline(file, line).eof()) {
		istringstream reader(line);
		vector<T> lineData;
		int k = 0;
		while (!reader.eof()) {
			T val;
			reader >> val;
			if (!reader.fail()) {
				lineData.push_back(val);
				k++;
			} else
				break;
		}
		data.push_back(lineData);
		if (k > n) {
			n = k;
		}
		m++;
	}
	file.close();

	mat.resize(m, n);
	mat.fill(0);

	for (size_t i = 0; i < data.size(); i++) {
		for (size_t j = 0; j < data[i].size(); j++) {
			mat.data[i * n + j] = data[i][j];
		}
	}
	
	logInfo("read '%s' OK!\n",filePath);
}

/* array IO*/
void writeMat(int m, int n, const double* data, const char* fmstr, ...);
/* matching IO*/
void write(const Matching& matches, const char* fmstr, ...);
void read(Matching& matches, const char* fmtstr, ...);
/* read intrinsic parameters */
bool readIntrinParam(const char* fname, double* K); //K : 3x3 matrix
/* read intrinsic parameters and distortion */
bool readIntrinDistParam(const char* fname, double* K, double* D); // D : 1x5 vector
/* read K,R,t*/
bool readKRT(double* K, double* R, double* t, const char* fmstr, ...);

/*write K,R,t*/
void writeKRT(double* K, double* R, double* t, const char* fmstr, ...);

#include "geometry/SL_BundleAdjust.h"
#include "geometry/SL_Point.h"
void write(const vector<vector<Meas2D> >& vec_meas, const char* fmtstr, ...);
void read(vector<vector<Meas2D> >& vec_meas, const char* fmtstr, ...);

void write(const vector<Point3d>& vec_pts, const char* fmtstr, ...);
void read(vector<Point3d>& vec_pts, const char* fmtstr, ...);
void write(const vector<Mat_d>& vec_mat, const char* fmtstr, ...);
void read(vector<Mat_d>& vec_mat, const char* fmtstr, ...);
#endif /* SL_WRITEREAD_H_ */
