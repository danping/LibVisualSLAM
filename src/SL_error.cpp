/*
 * SL_error.cpp
 *
 *  Created on: 2010-11-5
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#include "SL_error.h"
SL_Exception::SL_Exception(void) {
}
SL_Exception::SL_Exception(const char* str) {
	strcpy(err_str, str);
}
SL_Exception::~SL_Exception(void) {

}
void warn(const char* fmtstr, ...) {
	char buf[1024];
	GET_FMT_STR(fmtstr, buf)
	fprintf(stderr, buf);
}
void repErr(const char* fmtstr, ...) {
	char buf[1024];
	GET_FMT_STR(fmtstr, buf)
	logInfo("%s\n", buf);
	throw SL_Exception(buf);
}

void logInfo(const char* fmtstr, ...) {
	char buf[1024];
	GET_FMT_STR(fmtstr, buf)
	//can be replaced by other methods to output the information
#ifdef WIN32
	//TRACE1("%s",buf);
	printf("%s", buf);
#else
	printf("%s", buf);
#endif
}
void LogFile::open(const char* fmtstr, ...){
	char buf[1024];
	GET_FMT_STR(fmtstr, buf);
	fp = fopen(buf, "w+");
}
LogFile::LogFile(const char* fmtstr, ...) {
	char buf[1024];
	GET_FMT_STR(fmtstr, buf);
	fp = fopen(buf, "w+");
}
void LogFile::print(const char* fmtstr, ...) {
	char buf[1024];
	GET_FMT_STR(fmtstr, buf);
	fprintf(fp, "%s", buf);
}
