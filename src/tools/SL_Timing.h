/*
 * SL_Timing.h
 *
 *  Created on: 2011-5-18
 *      Author: Danping Zou
 */

#ifndef SL_TIMING_H_
#define SL_TIMING_H_
#include "SL_error.h"
#include <vector>
#include <cstdio>

#define TIMING_DATA_LEN 4

class Timing {
public:
	int f;
	double t;
	double data[TIMING_DATA_LEN];
public:
	Timing();
	Timing(int frame, double tm) :
			f(frame), t(tm) {
		memset(data, 0, sizeof(double) * TIMING_DATA_LEN);
	}
	Timing(const Timing& other) :
			f(other.f), t(other.t) {
		memcpy(data, other.data, sizeof(double) * TIMING_DATA_LEN);
	}
	Timing& operator =(const Timing& other) {
		if (&other != this) {
			f = other.f;
			t = other.t;
			memcpy(data, other.data, sizeof(double) * TIMING_DATA_LEN);
		}
		return *this;
	}
	~Timing();
};

class Timings {
public:
	std::vector<Timing> tms;
public:
	void clear() {
		tms.clear();
	}
	Timing& add(int f, double t) {
		tms.push_back(Timing(f, t));
		return tms.back();
	}
	void save(const char* fmtstr, ...) {
		char buf[1024];
		GET_FMT_STR(fmtstr,buf)

		FILE * fp = fopen(buf, "w");
		if (!fp)
		repErr("Timings::save - cannot open the file!\n");

		for (size_t i = 0; i < tms.size(); i++) {
			fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", tms[i].f, tms[i].t, tms[i].data[0], tms[i].data[1],
					tms[i].data[2], tms[i].data[3]);
		}
		fclose(fp);
	}
};
#endif /* SL_TIMING_H_ */
