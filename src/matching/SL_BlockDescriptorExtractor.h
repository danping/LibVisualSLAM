/*
 * BlockDescriptorExtractor.h
 *
 *  Created on: 2011-1-8
 *      Author: Danping Zou
 */

#ifndef BLOCKDESCRIPTOREXTRACTOR_H_
#define BLOCKDESCRIPTOREXTRACTOR_H_
#include "math/SL_Matrix.h"
#include "imgproc/SL_Image.h"

class BlockDescriptorExtractorParam {
public:
	int hW;
	double scale;
public:
	BlockDescriptorExtractorParam() {
		hW = 4;
		scale = 0.3;
	}
};
class BlockDescriptorExtractor : public BlockDescriptorExtractorParam {
public:
	BlockDescriptorExtractor();
	~BlockDescriptorExtractor();
public:
	void compute(const ImgG& img , const Mat_d& pts , Mat_d& desc);
};

//the above will be replaced by the `getImageBlocks' in `SL_StereoMatcherHelper.h'
#endif /* BLOCKDESCRIPTOREXTRACTOR_H_ */
