/*
 * SL_USBCamReader.h
 *
 *  Created on: 2012-4-10
 *      Author: elezoud
 */

#ifndef SL_USBCAMREADER_H_
#define SL_USBCAMREADER_H_
#include "SL_VideoReader.h"
#include "SL_Tictoc.h"

class USBCamReader:public VideoReader {
protected:
	CvCapture* videoCap;
	TimeMeasurer _tm;
	uint32_t _timstamp; 
public:
	int camid;
public:
	USBCamReader() :
		videoCap(0), camid(-1){
			avi = false;
	}
	virtual ~USBCamReader();
	virtual void open();
	virtual void grabFrame();
	virtual void readCurFrameRGB(unsigned char* imgdata);
	virtual void readCurFrameGray(unsigned char* grayImgData);
	virtual void readCurFrame(unsigned char* rgbdata, unsigned char* graydata);
	virtual uint32_t getTimeStamp();
};
#endif /* SL_USBCAMREADER_H_ */
