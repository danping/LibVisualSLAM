/*
 * SL_Helper.cpp
 *
 *  Created on: 2011-1-6
 *      Author: Danping Zou
 */
#include "SL_Print.h"
#include "SL_error.h"

void printMat(int m , int n , const double* mat) {
	int i, j;
	logInfo("-------------+++++++++++++\n");
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			logInfo("%.12g ", mat[i * n + j]);
		}
		logInfo("\n");
	}
	logInfo("+++++++++++++-------------\n");
}

void printMat(FILE* fp , int m , int n , const double* mat) {
	int i, j;
	fprintf(fp, "-------------+++++++++++++\n");
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			fprintf(fp, "%g ", mat[i * n + j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "+++++++++++++-------------\n");

}
void printMat(FILE* fp , int m , int n , const int* mat) {
	int i, j;
	fprintf(fp, "-------------+++++++++++++\n");
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			fprintf(fp, "%d ", mat[i * n + j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "+++++++++++++-------------\n");

}
void print(const Mat_d& mat) {
	printMat(mat.rows, mat.cols, mat.data);
}
void print(FILE* fp, const Mat_d& mat){
	printMat(fp, mat.rows, mat.cols, mat.data);
}
void print(FILE* fp, const Mat_i& mat){
	printMat(fp, mat.rows, mat.cols, mat.data);
}
void print(const Mat_f& mat) {
	printMat(mat.rows, mat.cols, mat.data);
}
void print(const Mat_i& mat) {
	printMat(mat.rows, mat.cols, mat.data);
}
void print(const Mat_ui& mat) {
	printMat(mat.rows, mat.cols, mat.data);
}
void print(const Mat_c& mat) {
	printMat(mat.rows, mat.cols, mat.data);
}
void print(const VecPoint2d& pts) {
	size_t n = pts.size();
	logInfo("num:%d\n", n);
	for (size_t i = 0; i < n; ++i)
		logInfo("%g\t%g\n", pts[i].x, pts[i].y);

}
void print(const Matching& matches) {
	int n = matches.num;
	logInfo("matches number:%d\n", n);
	for (int i = 0; i < n; i++) {
		logInfo("  %d -- %d (%g)\n", matches.data[i].idx1, matches.data[i].idx2, matches.data[i].dist);
	}

}

void printPts(int n , const Point2d* pts) {
	logInfo("num:%d\n", n);
	int i;
	for (i = 0; i < n; ++i)
		logInfo("%g\t%g\n", pts[i].m[0], pts[i].m[1]);
}
void printMat(int m , int n , const float* mat) {
	int i, j;
	logInfo("-------------+++++++++++++\n");
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			logInfo("%f ", mat[i * n + j]);
		}
		logInfo("\n");
	}
	logInfo("+++++++++++++-------------\n");
}
void printTransMat(int m , int n , const float* mat) {
	int i, j;
	logInfo("-------------+++++++++++++\n");
	for (j = 0; j < n; ++j) {
		for (i = 0; i < m; ++i) {
			logInfo("%f ", mat[i * n + j]);
		}
		logInfo("\n");
	}
	logInfo("+++++++++++++-------------\n");
}
void printMat(int m , int n , const int * mat) {
	int i, j;
	logInfo("-------------+++++++++++++\n");
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			logInfo("%d ", mat[i * n + j]);
		}
		logInfo("\n");
	}
	logInfo("+++++++++++++-------------\n");
}
void printMat(int m , int n , const unsigned int * mat) {
	int i, j;
	logInfo("-------------+++++++++++++\n");
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			logInfo("%d ", mat[i * n + j]);
		}
		logInfo("\n");
	}
	logInfo("+++++++++++++-------------\n");
}
void printMat(int m , int n , const char* mat) {
	int i, j;
	logInfo("-------------+++++++++++++\n");
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			logInfo("%d ", mat[i * n + j]);
		}
		logInfo("\n");
	}
	logInfo("+++++++++++++-------------\n");
}
