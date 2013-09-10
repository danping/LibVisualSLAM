/*
 * SL_ImageIO.cpp
 *
 *  Created on: Dec 22, 2011
 *      Author: Danping Zou
 */

#include "SL_ImageIO.h"

void savePPM(const uchar* data, int w, int h, const char* fmtstr, ...) {
	assert(data);
	char file_path[1024];
	GET_FMT_STR(fmtstr, file_path)

	FILE *pFile;
	int y;

	// Open file
	pFile = fopen(file_path, "wb");
	if (pFile == NULL)
		repErr("save_ppm() : cannot open file %s.", file_path);

	// Write header
	fprintf(pFile, "P6\n%d %d\n255\n", w, h);

	// Write pixel data
	for (y = 0; y < h; y++)
		fwrite(data + y * w * 3, 1, w * 3, pFile);

	// Close file
	fclose(pFile);

}
void savePPM(const ImgRGB& rgb_img, const char* fmtstr, ...) {
	if (rgb_img.empty())
		repErr("save_ppm() : empty image.");

	char file_path[1024];
	GET_FMT_STR(fmtstr, file_path)

	FILE *pFile;
	int y;

	// Open file
	pFile = fopen(file_path, "wb");
	if (pFile == NULL)
		repErr("save_ppm() : cannot open file %s.", file_path);

	int w = rgb_img.w;
	int h = rgb_img.h;

	// Write header
	fprintf(pFile, "P6\n%d %d\n255\n", w, h);

	int l = w*h;
	unsigned char* img_data = new unsigned char[3*l];
	for (int i = 0; i < l; ++i) {
		img_data[3*i] = rgb_img.data[3*i+2];
		img_data[3*i+1] = rgb_img.data[3*i+1];
		img_data[3*i+2] = rgb_img.data[3*i];
	}
	// Write pixel data
	for (y = 0; y < h; y++)
		fwrite(img_data + y * w * 3, 1, w * 3, pFile);

	// Close file
	fclose(pFile);
	
	delete[] img_data;
}
void savePGM(const uchar* data, int w, int h, const char* fmtstr, ...) {
	assert(data);
	char file_path[1024];
	GET_FMT_STR(fmtstr, file_path)

	FILE *pFile;
	int y;

	// Open file
	pFile = fopen(file_path, "wb");
	if (pFile == NULL)
		repErr("save_ppm() : cannot open file %s.", file_path);

	// Write header
	fprintf(pFile, "P5\n%d %d\n255\n", w, h);

	// Write pixel data
	for (y = 0; y < h; y++)
		fwrite(data + y * w, 1, w, pFile);

	// Close file
	fclose(pFile);
}
void savePGM(const ImgG& gray_img, const char* fmtstr, ...) {
	assert(gray_img);
	if (gray_img.empty())
		repErr("save_ppm() : empty image.");

	char file_path[1024];
	GET_FMT_STR(fmtstr, file_path)

	FILE *pFile;
	int y;

	// Open file
	pFile = fopen(file_path, "wb");
	if (pFile == NULL)
		repErr("save_ppm() : cannot open file %s.", file_path);

	int w = gray_img.w;
	int h = gray_img.h;

	// Write header
	fprintf(pFile, "P5\n%d %d\n255\n", w, h);

	// Write pixel data
	for (y = 0; y < h; y++)
		fwrite(gray_img.data + y * w, 1, w, pFile);

	// Close file
	fclose(pFile);
	logInfo("save '%s' [ok]\n", file_path);
}

void loadPGM(ImgG& gray_img, const char* fmtstr, ...) {
	char file_path[1024];
	GET_FMT_STR(fmtstr, file_path)

	const int PPMREADBUFLEN = 1024;
	char buf[PPMREADBUFLEN], *t;
	int w, h, d;
	int r;

	FILE* pf = fopen(file_path, "rb");
	if (pf == 0)
		repErr("read_pgm() - cannot open the file:%s\n", file_path);

	t = fgets(buf, PPMREADBUFLEN, pf);
	if ((t == NULL) || (strncmp(buf, "P5\n", 3) != 0))
		repErr("read_pgm() - P5 is required!");
	do { /* Px formats can have # comments after first line */
		t = fgets(buf, PPMREADBUFLEN, pf);
		if (t == 0)
			repErr("read_pgm() - error");
	} while (strncmp(buf, "#", 1) == 0);
	r = sscanf(buf, "%u %u", &w, &h);
	if (r < 2)
		repErr("read_pgm() - error");

	r = fscanf(pf, "%u\n", &d);
	if ((r < 1) || (d != 255))
		repErr("read_pgm() - error");

	gray_img.resize(w, h);
	size_t len = w * h;

	size_t rd = fread(gray_img.data, sizeof(char), len, pf);
	if (rd < len) {
		gray_img.clear();
		repErr("read_pgm() - '%s' failed!", file_path);
	}
	logInfo("load '%s' [ok]\n", file_path);
}
void loadPPM(ImgRGB& rgb_img, const char* fmtstr, ...) {
	char file_path[1024];
	GET_FMT_STR(fmtstr, file_path);

	char buff[16];
	FILE *fp;
	int c, rgb_comp_color;
	//open PPM file for reading
	fp = fopen(file_path, "rb");
	if (!fp) {
		fprintf(stderr, "Unable to open file '%s'\n", file_path);
		exit(1);
	}

	//read image format
	if (!fgets(buff, sizeof(buff), fp)) {
		perror(file_path);
		exit(1);
	}

	//check the image format
	if (buff[0] != 'P' || buff[1] != '6') {
		fprintf(stderr, "Invalid image format (must be 'P6')\n");
		exit(1);
	}

	//check for comments
	c = getc(fp);
	while (c == '#') {
		while (getc(fp) != '\n')
			;
		c = getc(fp);
	}

	ungetc(c, fp);
	//read image size information
	int w, h;
	if (fscanf(fp, "%d %d", &w, &h) != 2) {
		fprintf(stderr, "Invalid image size (error loading '%s')\n", file_path);
		exit(1);
	}

	//read rgb component
	if (fscanf(fp, "%d", &rgb_comp_color) != 1) {
		fprintf(stderr, "Invalid rgb component (error loading '%s')\n",
				file_path);
		exit(1);
	}

	//check rgb component depth
	if (rgb_comp_color != 255) {
		fprintf(stderr, "'%s' does not have 8-bits components\n", file_path);
		exit(1);
	}

	while (fgetc(fp) != '\n')
		;
	//memory allocation for pixel data
	rgb_img.resize(w, h);

	//read pixel data from file
	if (fread(rgb_img.data, 3 * w, h, fp) != (size_t) h) {
		fprintf(stderr, "Error loading image '%s'\n", file_path);
		exit(1);
	}

	//swap the r,b channel
	int l = w * h;
	for (int i = 0; i < l; ++i) {
		uchar t = rgb_img.data[3 * i];
		rgb_img.data[3 * i] = rgb_img.data[3 * i + 2];
		rgb_img.data[3 * i + 2] = t;
	}
	fclose(fp);
}
//use opencv to read images
#include <opencv2/opencv.hpp>
void imread(ImgG& gray_img, const char* fmtstr, ...) {
	char file_path[1024];
	GET_FMT_STR(fmtstr, file_path);
	cv::Mat cvImg = cv::imread(file_path, 0);
	gray_img.resize(cvImg.cols, cvImg.rows);
	memcpy(gray_img.data, cvImg.data, sizeof(uchar) * cvImg.cols * cvImg.rows);
}
void imread(ImgRGB& rgb_img, const char* fmtstr, ...) {
	using namespace std;
	char file_path[1024];
	GET_FMT_STR(fmtstr, file_path);
	//test
	cout <<"imread: " << file_path << endl;
	cv::Mat cvImg = cv::imread(file_path, 1);
	cv::Mat cvBGR(rgb_img.h, rgb_img.w, CV_8UC3);
	cv::cvtColor(cvImg, cvBGR, cv::COLOR_RGB2BGR);
	
	rgb_img.resize(cvImg.cols, cvImg.rows);
	memcpy(rgb_img.data, cvBGR.data,
			sizeof(uchar) * cvBGR.cols * cvBGR.rows * 3);
}
void imwrite(const ImgG& gray_img, const char* fmtstr, ...) {
	char file_path[1024];
	GET_FMT_STR(fmtstr, file_path);
	cv::Mat cvImg(gray_img.h, gray_img.w, CV_8UC1, gray_img.data);
	cv::imwrite(file_path, cvImg);
}
void imwrite(const ImgRGB& rgb_img, const char* fmtstr, ...) {
	char file_path[1024];
	GET_FMT_STR(fmtstr, file_path);
	cv::Mat cvImg(rgb_img.h, rgb_img.w, CV_8UC3, rgb_img.data);
	cv::Mat cvBGR(rgb_img.h, rgb_img.w, CV_8UC3);
	cv::cvtColor(cvImg, cvBGR, cv::COLOR_RGB2BGR);
	cv::imwrite(file_path, cvBGR);
}