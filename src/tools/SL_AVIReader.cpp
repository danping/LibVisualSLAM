/*
 * SL_AVIReader.cpp
 *
 *  Created on: 2012-4-10
 *      Author: elezoud
 */
#include "SL_AVIReader.h"
#include "SL_error.h"
void AVIReader::open() {
#if(CV_MAJOR_VERSION <= 2)
	if (videoCap)
		cvReleaseCapture(&videoCap);
	videoCap = cvCaptureFromFile(filePath.c_str());
	if (!videoCap)
		repErr("Cannot open %s\n", filePath.c_str());
	_w = (int) cvGetCaptureProperty(videoCap, CV_CAP_PROP_FRAME_WIDTH);
	_h = (int) cvGetCaptureProperty(videoCap, CV_CAP_PROP_FRAME_HEIGHT);
#else
    videoCap.release();
    if( !videoCap.open(filePath))
        repErr("Cannot open %s\n", filePath.c_str());
        
    _w = videoCap.get(CV_CAP_PROP_FRAME_WIDTH);
    _h = videoCap.get(CV_CAP_PROP_FRAME_HEIGHT);
#endif

	assert(_w > 0 && _h > 0);
	_timestamp = 0;
}

void AVIReader::grabFrame() {
#if(CV_MAJOR_VERSION <= 2)    
    assert(videoCap);
	cvGrabFrame(videoCap);
#else
    videoCap.grab();
#endif
}
//skip: setting skip frame numbers can only be used with video
void AVIReader::skip(int nSkippedFrame) {
#if(CV_MAJOR_VERSION <= 2)    
	cvSetCaptureProperty(videoCap, CV_CAP_PROP_POS_FRAMES, nSkippedFrame);
#else
    videoCap.set(CV_CAP_PROP_POS_FRAMES, nSkippedFrame);
#endif
}
void AVIReader::grabReadFrame(int Frame, unsigned char* imgdata) {
	skip(Frame);
	grabFrame();
	readCurFrameRGB(imgdata);
}
void AVIReader::readCurFrameRGB(unsigned char* imgdata) {
#if(CV_MAJOR_VERSION <= 2)    
    assert(videoCap);
	IplImage* img = cvRetrieveFrame(videoCap);
	memcpy(imgdata, img->imageData, _w * _h * 3);
#else
    cv::Mat cvImg(_h,_w, CV_8UC3,imgdata);
    videoCap.retrieve(cvImg);
#endif
}
void AVIReader::readCurFrameGray(unsigned char* grayImgData) {
#if(CV_MAJOR_VERSION <= 2)        
	assert(videoCap);
	IplImage* img = cvRetrieveFrame(videoCap);
	cv::Mat rawFrame(img);
	cv::Mat videoFrame(_h, _w, CV_8UC1, grayImgData);
#else
    cv::Mat rawFrame;
    cv::Mat videoFrame(_h, _w, CV_8UC1, grayImgData);    
    videoCap.retrieve(rawFrame);
#endif
    if(rawFrame.type() == CV_8UC1)
		rawFrame.copyTo(videoFrame);
	else
		cv::cvtColor(rawFrame, videoFrame, CV_RGB2GRAY);
}
void AVIReader::readCurFrame(unsigned char* rgbdata,
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
int AVIReader::getTotalFrame() {
#if(CV_MAJOR_VERSION <= 2)            
	return (int)cvGetCaptureProperty(videoCap, CV_CAP_PROP_FRAME_COUNT);
#else
    return (int) videoCap.get(CV_CAP_PROP_FRAME_COUNT);
#endif
}

void AVIReader::releaseCamera(){
#if(CV_MAJOR_VERSION <= 2)                
	cvReleaseCapture(&videoCap);
#endif
}

uint32_t AVIReader::getTimeStamp(){
	return _timestamp++;
}
