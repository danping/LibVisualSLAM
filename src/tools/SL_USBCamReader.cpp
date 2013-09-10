/*
 * SL_USBCamReader.cpp
 *
 *  Created on: 2012-4-10
 *      Author: elezoud
 */

#include "SL_USBCamReader.h"
#include "SL_error.h"

void USBCamReader::open() {
	if (videoCap)
		cvReleaseCapture(&videoCap);
	videoCap = cvCaptureFromCAM(camid);
	cvSetCaptureProperty(videoCap, CV_CAP_PROP_FRAME_WIDTH, 640);
	cvSetCaptureProperty(videoCap, CV_CAP_PROP_FRAME_HEIGHT, 480);
	cvSetCaptureProperty(videoCap, CV_CAP_PROP_FPS, 30);
	if (!videoCap) {
		repErr("ERROR: Fail to detect the USB camera \n");
	}
	_w = (int)cvGetCaptureProperty(videoCap, CV_CAP_PROP_FRAME_WIDTH);
	_h = (int)cvGetCaptureProperty(videoCap, CV_CAP_PROP_FRAME_HEIGHT);
	if (_w <=0 || _h <= 0)
		repErr("ERROR: Fail to open the USB camera \n");

	_tm.tic();
}
USBCamReader::~USBCamReader(){
	if(videoCap)
		cvReleaseCapture(&videoCap);
}

void USBCamReader::grabFrame() {
	assert(videoCap);
	cvGrabFrame(videoCap);
}
void USBCamReader::readCurFrameRGB(unsigned char* imgdata) {
	assert(videoCap);
	IplImage* img = cvRetrieveFrame(videoCap);
	memcpy(imgdata, img->imageData, _w * _h * 3);
}
void USBCamReader::readCurFrameGray(unsigned char* grayImgData) {
	assert(videoCap);
	IplImage* img = cvRetrieveFrame(videoCap);
	cv::Mat rawFrame(img);
	cv::Mat videoFrame(_h, _w, CV_8UC1, grayImgData);
	cv::cvtColor(rawFrame, videoFrame, CV_RGB2GRAY);
}
void USBCamReader::readCurFrame(unsigned char* rgbdata,
		unsigned char* graydata) {
	assert(videoCap);
	IplImage* img = cvRetrieveFrame(videoCap);
	cv::Mat rawFrame(img);

	cv::Mat rgbImg(_h, _w, CV_8UC3, rgbdata);
	cv::cvtColor(rawFrame, rgbImg, CV_BGR2RGB);

	cv::Mat videoFrame(_h, _w, CV_8UC1, graydata);
	cv::cvtColor(rawFrame, videoFrame, CV_RGB2GRAY);
}
uint32_t USBCamReader::getTimeStamp(){
	double ts = _tm.get_pass_time();
	return (uint32_t)(ts+0.5);
}