/*
 * SL_error.h
 *
 *  Created on: 2010-11-5
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_ERROR_H_
#define SL_ERROR_H_

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>

#define GET_FMT_STR(fmstr,buf) \
	va_list args;\
	va_start(args,fmtstr);\
	vsprintf(buf,fmstr,args);\
	va_end(args);

class SL_Exception {
protected:
	char err_str[1024];
public:
	SL_Exception(void);
	SL_Exception(const char* str);
	~SL_Exception(void);
	const char* what() const {
		return err_str;
	}
};
//report a warning without interrupting the execution of the program
void warn(const char* fmtstr, ...);
//report an error exception
void repErr(const char* fmtstr, ...);
//just log out information
void logInfo(const char* fmtstr, ...);


class LogFile {
public:
	FILE* fp;
public:
	LogFile(){}
	LogFile(const char* fmtstr, ...);
	~LogFile() {
		if (fp)
			fclose(fp);
	}
	void print(const char* fmtstr, ...);
	void open(const char* fmtstr, ...);
	void close() {
		fclose(fp);
		fp = 0;
	}
};
#endif /* SL_ERROR_H_ */
