/*
 * SL_ImageIO.h
 *
 *  Created on: Dec 22, 2011
 *      Author: Danping Zou
 */

#ifndef SL_IMAGEIO_H_
#define SL_IMAGEIO_H_
#include "SL_Image.h"
/* save as ppm (color image)*/
void savePPM(const ImgRGB& rgb_img, const char* fmtstr, ...);
void savePPM(const uchar* data, int w, int h, const char* fmtstr, ...);

/* save as pgm (gray image)*/
void savePGM(const ImgG& gray_img, const char* fmtstr, ...);
void savePGM(const uchar* data, int w, int h, const char* fmtstr, ...);

/* read pgm files*/
void loadPGM(ImgG& gray_img, const char* fmtstr, ...);
void loadPPM(ImgRGB& rgb_img, const char* fmtstr, ...);

/*read image*/
void imread(ImgG& gray_img, const char* fmtstr, ...);
void imread(ImgRGB& rgb_img, const char* fmtstr, ...);
void imwrite(const ImgG& gray_img, const char* fmtstr, ...);
void imwrite(const ImgRGB& rgb_img, const char* fmtstr, ...);

#endif /* SL_IMAGEIO_H_ */
