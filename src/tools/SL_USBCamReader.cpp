/*
 * SL_USBCamReader.cpp
 *
 *  Created on: 2012-4-10
 *      Author: elezoud
 */

#include "SL_USBCamReader.h"
#include "SL_error.h"

void USBCamReader::open() {
#if(CV_MAJOR_VERSION <= 2)                
	if (videoCap)
		cvReleaseCapture(&videoCap);
    
	videoCap = cvCaptureFromCAM(camid);
    if (!videoCap) {
		repErr("ERROR: Fail to detect the USB camera %d\n", camid);
	}
    
	cvSetCaptureProperty(videoCap, CV_CAP_PROP_FRAME_WIDTH, 640);
	cvSetCaptureProperty(videoCap, CV_CAP_PROP_FRAME_HEIGHT, 480);
	cvSetCaptureProperty(videoCap, CV_CAP_PROP_FPS, 30);
	_w = (int)cvGetCaptureProperty(videoCap, CV_CAP_PROP_FRAME_WIDTH);
	_h = (int)cvGetCaptureProperty(videoCap, CV_CAP_PROP_FRAME_HEIGHT);
#else
    if( !videoCap.open(camid)){
        repErr("ERROR: Fail to detect the USB camera %d \n", camid);
    }
    videoCap.set(CV_CAP_PROP_FRAME_WIDTH,640);
    videoCap.set(CV_CAP_PROP_FRAME_WIDTH,480);
    videoCap.set(CV_CAP_PROP_FPS,30);
    
    _w = (int) videoCap.get(CV_CAP_PROP_FRAME_WIDTH);
    _w = (int) videoCap.get(CV_CAP_PROP_FRAME_HEIGHT);
#endif
            
    if (_w <=0 || _h <= 0)
		repErr("ERROR: Fail to open the USB camera \n");

	_tm.tic();
}
USBCamReader::~USBCamReader(){
#if(CV_MAJOR_VERSION <= 2)                
	if(videoCap)
		cvReleaseCapture(&videoCap);
#endif
}

void USBCamReader::grabFrame() {
#if(CV_MAJOR_VERSION <= 2)                    
	assert(videoCap);
	cvGrabFrame(videoCap);
#else
    videoCap.grab();
#endif
}
void USBCamReader::readCurFrameRGB(unsigned char* imgdata) {
#if(CV_MAJOR_VERSION <= 2)                    
	assert(videoCap);
	IplImage* img = cvRetrieveFrame(videoCap);
	memcpy(imgdata, img->imageData, _w * _h * 3);
#else
    cv::Mat cvImg(_w,_h,CV_8UC3, imgdata);
    videoCap.retrieve(cvImg);
#endif
}
void USBCamReader::readCurFrameGray(unsigned char* grayImgData) {
#if(CV_MAJOR_VERSION <= 2)                
	assert(videoCap);
	IplImage* img = cvRetrieveFrame(videoCap);
	cv::Mat rawFrame(img);
#else
    cv::Mat rawFrame;
    videoCap.retrieve(rawFrame);    
#endif
	cv::Mat videoFrame(_h, _w, CV_8UC1, grayImgData);
	cv::cvtColor(rawFrame, videoFrame, CV_RGB2GRAY);
}
void USBCamReader::readCurFrame(unsigned char* rgbdata,
		unsigned char* graydata) {
#if(CV_MAJOR_VERSION <= 2)                    
	assert(videoCap);
	IplImage* img = cvRetrieveFrame(videoCap);
	cv::Mat rawFrame(img);
#else
    cv::Mat rawFrame;
    videoCap.retrieve(rawFrame);
#endif

	cv::Mat rgbImg(_h, _w, CV_8UC3, rgbdata);
	cv::cvtColor(rawFrame, rgbImg, CV_BGR2RGB);

	cv::Mat videoFrame(_h, _w, CV_8UC1, graydata);
	cv::cvtColor(rawFrame, videoFrame, CV_RGB2GRAY);
}
uint32_t USBCamReader::getTimeStamp(){
	double ts = _tm.get_pass_time();
	return (uint32_t)(ts+0.5);
}
