/*
 * SL_Tictoc.h
 *
 *  Created on: 2010-11-9
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_TICTOC_H_
#define SL_TICTOC_H_

#include <ctime>
#ifndef WIN32
	#include <sys/time.h>
#else
	#include "winsock2.h"
	#include <windows.h>
	int gettimeofday(struct timeval *tp, void *tzp);
#endif

#include <cstdio>
//for measurement of timing
class TimeMeasurer {
protected:
	struct timeval tm_s, tm_e;
public:
	TimeMeasurer(void) {
	}
	~TimeMeasurer(void) {
	}
public:
	void tic() {
		gettimeofday(&tm_s, 0);
	}
	double toc() {
		gettimeofday(&tm_e, 0);
		double tm = get_pass_time();
		return tm;
	}
	double get_pass_time() {
		double tm = (tm_e.tv_sec - tm_s.tv_sec) * 1000.0 + (tm_e.tv_usec - tm_s.tv_usec) / 1000.0;
		return tm;
	}
};
void tic();
double toc();
#endif /* SL_TICTOC_H_ */
