/*
* SL_WriteRead.cpp
*
*  Created on: 2010-11-13
*      Author: Danping Zou
*/

#include "SL_WriteRead.h"
#include "SL_error.h"
#include "SL_Print.h"
#include <ctype.h>
#include <cstdio>

#define MAXSTRLEN  2048 /* 2K */
/* get rid of the rest of a line upto \n or EOF */
#define SKIP_LINE(f){                                                       \
	char buf[MAXSTRLEN];                                                        \
	while(!feof(f))                                                           \
	if(!fgets(buf, MAXSTRLEN-1, f) || buf[strlen(buf)-1]=='\n') break;      \
}
void write(const VecPoint2d& pts, const char* fmtstr, ...) {
	char buf[1024];
	GET_FMT_STR(fmtstr, buf)

		FILE* fp = fopen(buf, "w");
	if (!fp)
		repErr("cannot open %s to write.", buf);

	/*write the number of points*/
	fprintf(fp, "%lu", pts.size());
	for (size_t i = 0; i < pts.size(); ++i) {
		fprintf(fp, " %lf %lf", pts[i].x, pts[i].y);
	}
	fprintf(fp, "\n");
	fclose(fp);
}
void read(VecPoint2d& pts, const char* fmtstr, ...) {
	char buf[1024];
	GET_FMT_STR(fmtstr, buf)

		FILE* fp = fopen(buf, "r");
	if (!fp)
		repErr("cannot open %s to read.", buf);

	pts.clear();

	/*read the number of points*/
	int i, npt;
	fscanf(fp, "%d", &npt);
	pts.resize(npt);

	for (i = 0; i < npt; ++i) {
		double x, y;
		fscanf(fp, "%lf%lf", &x, &y);
		pts.push_back(Point2d(x, y));
	}
	fclose(fp);
}
void writeMat(int m, int n, const double* data, const char* fmtstr, ...) {
	char buf[1024];
	GET_FMT_STR(fmtstr, buf)

		FILE* fp = fopen(buf, "w");
	if (!fp)
		repErr("cannot open %s to write.", buf);
	/*write the dimension*/
	int i, j;
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			fprintf(fp, "%g ", data[i * n + j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	logInfo("save '%s' [ok]\n", buf);
}

void write(const Matching& matches, const char* fmtstr, ...) {
	char buf[1024];
	GET_FMT_STR(fmtstr, buf)

		FILE* fp = fopen(buf, "w");
	if (!fp)
		repErr("cannot open %s to write.", buf);
	/*write the number of matches*/
	int nMatches = matches.num;
	fprintf(fp, "%d\n ", nMatches);
	for (int i = 0; i < nMatches; i++) {
		fprintf(fp, "%d %d %lf\n", matches[i].idx1, matches[i].idx2,
			matches[i].dist);
	}
	fclose(fp);
	logInfo("save '%s' [ok]\n", buf);
}
void read(Matching& matches, const char* fmtstr, ...) {
	char buf[1024];
	GET_FMT_STR(fmtstr, buf)

		FILE* fp = fopen(buf, "r");
	if (!fp)
		repErr("cannot open %s to read.", buf);
	int nMatches;
	fscanf(fp, "%d", &nMatches);
	matches.reserve(nMatches);

	for (int i = 0; i < nMatches; i++) {
		int idx1, idx2;
		double dist;
		fscanf(fp, "%d%d%lf", &idx1, &idx2, &dist);
		matches.add(idx1, idx2, dist);
	}
	fclose(fp);
	logInfo("load '%s' [ok]\n", buf);
}
/* reads (from "fp") "nvals" doubles into "vals".
* Returns number of doubles actually read, EOF on end of file, EOF-1 on error
*/
static int _readNDouble(FILE *fp, double *vals, int nvals) {
	register int i;
	int n, j;

	for (i = n = 0; i < nvals; ++i) {
		j = fscanf(fp, "%lf", vals + i);
		if (j == EOF)
			return EOF;

		if (j != 1 || ferror(fp))
			return EOF - 1;

		n += j;
	}

	return n;
}
/* reads (from "fp") "nvals" doubles without storing them.
* Returns EOF on end of file, EOF-1 on error
*/
static int _skipNDouble(FILE *fp, int nvals) {
	register int i;
	int j;
	for (i = 0; i < nvals; ++i) {
		j = fscanf(fp, "%*f");
		if (j == EOF)
			return EOF;
		if (ferror(fp))
			return EOF - 1;
	}
	return nvals;
}

/* reads (from "fp") "nvals" ints into "vals".
* Returns number of ints actually read, EOF on end of file, EOF-1 on error
*/
static int _readNInts(FILE *fp, int *vals, int nvals) {
	register int i;
	int n, j;

	for (i = n = 0; i < nvals; ++i) {
		j = fscanf(fp, "%d", vals + i);
		if (j == EOF)
			return EOF;
		if (j != 1 || ferror(fp))
			return EOF - 1;
		n += j;
	}
	return n;
}

bool readIntrinParam(const char* fname, double* K) {
	FILE *fp;
	int i, ch = EOF;

	if ((fp = fopen(fname, "r")) == NULL) {
		repErr("cannot open file %s, exiting\n", fname);
	}

	while (!feof(fp) && (ch = fgetc(fp)) == '#') /* skip comments */
		SKIP_LINE(fp);

	if (feof(fp)) {
		K[0] = K[1] = K[2] = K[3] = K[4] = K[5] = K[6] = K[7] = K[8] = 0.0;
		fclose(fp);
		return false;
	}

	ungetc(ch, fp);

	for (i = 0; i < 3; i++) {
		fscanf(fp, "%lf%lf%lf\n", K, K + 1, K + 2);
		K += 3;
	}

	fclose(fp);
	return true;
}

bool readIntrinDistParam(const char* fname, double* K, double* D) {
	FILE *fp;
	int i, ch = EOF;

	if ((fp = fopen(fname, "r+")) == NULL) {
		repErr("cannot open file %s, exiting\n", fname);
	}

	while (!feof(fp) && (ch = fgetc(fp)) == '#') /* skip comments */
		SKIP_LINE(fp);

	if (feof(fp)) {
		K[0] = K[1] = K[2] = K[3] = K[4] = K[5] = K[6] = K[7] = K[8] = 0.0;
		D[0] = D[1] = D[2] = D[3] = D[4] = 0;
		fclose(fp);
		return false;
	}

	ungetc(ch, fp);

	for (i = 0; i < 9; i++) {
		fscanf(fp, "%lf", K);
		K++;
	}
	while (!feof(fp) && isspace(ch = fgetc(fp)))
		;
	ungetc(ch, fp);
	while (!feof(fp) && (ch = fgetc(fp)) == '#') /* skip comments */
		SKIP_LINE(fp);

	if (feof(fp)) {
		D[0] = D[1] = D[2] = D[3] = D[4] = 0;
		return false;
	}
	D[0] = 1;
	D[1] = D[2] = D[3] = D[4] = 0;

	ungetc(ch, fp);

	for (i = 0; i < 5; i++) {
		char buf[256];
		fscanf(fp, "%s", buf);
		D[i] = atof(buf);
	}
	fclose(fp);

	return true;
}

bool readKRT(double* K, double* R, double* t, const char* fmtstr, ...) {
	char fname[1024];
	GET_FMT_STR(fmtstr, fname)

		FILE *fp;
	int i, ch = EOF;

	if ((fp = fopen(fname, "r+")) == NULL) {
		repErr("cannot open file %s, exiting\n", fname);
	}

	while (!feof(fp) && (ch = fgetc(fp)) == '#') /* skip comments */
		SKIP_LINE(fp);

	if (feof(fp)) {
		fclose(fp);
		return false;
	}

	ungetc(ch, fp);

	for (i = 0; i < 9; i++) {
		fscanf(fp, "%lf", K);
		K++;
	}
	while (!feof(fp) && isspace(ch = fgetc(fp)))
		;

	ungetc(ch, fp);
	while (!feof(fp) && (ch = fgetc(fp)) == '#') /* skip comments */
		SKIP_LINE(fp);

	if (feof(fp)) {
		return false;
	}

	ungetc(ch, fp);

	for (i = 0; i < 9; i++) {
		char buf[256];
		fscanf(fp, "%s", buf);
		R[i] = atof(buf);
	}

	while (!feof(fp) && isspace(ch = fgetc(fp)))
		;
	ungetc(ch, fp);
	while (!feof(fp) && (ch = fgetc(fp)) == '#') /* skip comments */
		SKIP_LINE(fp);

	if (feof(fp)) {
		return false;
	}

	ungetc(ch, fp);
	for (i = 0; i < 3; i++) {
		char buf[256];
		fscanf(fp, "%s", buf);
		t[i] = atof(buf);
	}
	fclose(fp);
	return true;
}

void writeKRT(double* K, double* R, double* t, const char* fmtstr, ...) {
	char fname[1024];
	GET_FMT_STR(fmtstr, fname)

		ofstream file(fname);
	if (!file)
		repErr("cannot open file '%s' to write.", fname);

	file << "#intrinsic matrix:" << endl;
	file << K[0] << " " << K[1] << " " << K[2] << endl;
	file << K[3] << " " << K[4] << " " << K[5] << endl;
	file << K[6] << " " << K[7] << " " << K[8] << endl;

	file << "#rotation matrix:" << endl;
	file << R[0] << " " << R[1] << " " << R[2] << endl;
	file << R[3] << " " << R[4] << " " << R[5] << endl;
	file << R[6] << " " << R[7] << " " << R[8] << endl;

	file << "#translation:" << endl;
	file << t[0] << endl;
	file << t[1] << endl;
	file << t[2] << endl;
}

void write(const vector<vector<Meas2D> >& vec_meas, const char* fmtstr, ...) {
	char fname[1024];
	GET_FMT_STR(fmtstr, fname)

		ofstream file(fname);
	if (!file)
		repErr("cannot open file '%s' to write.", fname);

	file << vec_meas.size() << endl;

	for (size_t i = 0; i < vec_meas.size(); ++i) {
		file << vec_meas[i].size() << endl;
		for (size_t j = 0; j < vec_meas[i].size(); ++j) {
			file << vec_meas[i][j].viewId << " " << vec_meas[i][j].x << " "
				<< vec_meas[i][j].y << " " << vec_meas[i][j].w << " ";
			file << vec_meas[i][j].cov[0] << " " << vec_meas[i][j].cov[1] << " "
				<< vec_meas[i][j].cov[2] << " " << vec_meas[i][j].cov[3]
			<< " " << vec_meas[i][j].outlier << endl;
		}
	}
	logInfo("save '%s' [ok]\n", fname);
}
void read(vector<vector<Meas2D> >& vec_meas, const char* fmtstr, ...) {
	char filePath[1024];
	GET_FMT_STR(fmtstr, filePath);
	using namespace std;

	ifstream file(filePath);
	if (!file)
		repErr("read - cannot open '%s'!\n", filePath);

	size_t nmeas;
	file >> nmeas;

	vec_meas.clear();
	for (size_t i = 0; i < nmeas; ++i) {
		size_t npts;
		file >> npts;
		vec_meas.push_back(vector<Meas2D>());
		vec_meas.back().resize(npts);

		for (size_t j = 0; j < npts; ++j) {
			int viewId;
			double x, y, w, cov[4];
			int outlier;

			file >> viewId >> x >> y >> w >> cov[0] >> cov[1] >> cov[2]
			>> cov[3] >> outlier;

			vec_meas.back()[j].viewId = viewId;
			vec_meas.back()[j].x = x;
			vec_meas.back()[j].y = y;
			vec_meas.back()[j].w = w;
			vec_meas.back()[j].cov[0] = cov[0];
			vec_meas.back()[j].cov[1] = cov[1];
			vec_meas.back()[j].cov[2] = cov[2];
			vec_meas.back()[j].cov[3] = cov[3];
			vec_meas.back()[j].outlier = outlier;
		}
	}
}
void write(const vector<Point3d>& vec_pts, const char* fmtstr, ...) {
	char fname[1024];
	GET_FMT_STR(fmtstr, fname)

		ofstream file(fname);
	if (!file)
		repErr("cannot open file '%s' to write.", fname);

	file << vec_pts.size() << endl;
	for (size_t i = 0; i < vec_pts.size(); ++i) {
		file << vec_pts[i].x << " " << vec_pts[i].y << " " << vec_pts[i].z
			<< endl;
	}
	logInfo("save '%s' [ok]\n", fname);
}
void read(vector<Point3d>& vec_pts, const char* fmtstr, ...) {
	char filePath[1024];
	GET_FMT_STR(fmtstr, filePath);
	using namespace std;

	ifstream file(filePath);
	if (!file)
		repErr("read - cannot open '%s'!\n", filePath);

	size_t npts;
	file >> npts;

	vec_pts.clear();
	for (size_t i = 0; i < npts; ++i) {
		double x, y, z;
		file >> x >> y >> z;
		vec_pts.push_back(Point3d(x, y, z));
	}

}
void write(const vector<Mat_d>& vec_mat, const char* fmtstr, ...) {
	char fname[1024];
	GET_FMT_STR(fmtstr, fname)

		ofstream file(fname);
	if (!file)
		repErr("cannot open file '%s' to write.", fname);

	file << vec_mat.size() << endl;
	for (size_t i = 0; i < vec_mat.size(); ++i) {
		file << vec_mat[i].m << " " << vec_mat[i].n << endl;
		int len = vec_mat[i].m * vec_mat[i].n;
		for (int j = 0; j < len; ++j) {
			file << vec_mat[i].data[j] << " ";
		}
		file << endl;
	}
}
void read(vector<Mat_d>& vec_mat, const char* fmtstr, ...) {
	char filePath[1024];
	GET_FMT_STR(fmtstr, filePath);
	using namespace std;

	ifstream file(filePath);
	if (!file)
		repErr("read - cannot open '%s'!\n", filePath);

	size_t nmat;
	file >> nmat;
	vec_mat.clear();
	for (size_t i = 0; i < nmat; ++i) {
		int m, n;
		file >> m >> n;
		vec_mat.push_back(Mat_d(m, n));
		int len = m * n;
		for (int j = 0; j < len; ++j) {
			file >> vec_mat.back().data[j];
		}
	}
}
