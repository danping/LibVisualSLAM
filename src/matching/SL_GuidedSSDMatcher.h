/*
 * SL_GuidedSSDMatcher.h
 *
 *  Created on: 2011-1-25
 *      Author: Danping Zou
 */

#ifndef SL_GUIDEDSTEREOMATCHER_H_
#define SL_GUIDEDSTEREOMATCHER_H_
#include "SL_Matching.h"
#include "SL_StereoMatcherHelper.h"

#include "math/SL_Matrix.h"
#include "math/SL_LinAlg.h"
#include "geometry/SL_Geometry.h"

class GuidedSSDMatcherParam {
public:
	double epiDistThres; 	//threshold of epipolar error
	double grayLevelThres; 	//threshold of descriptor
	double dispDistThres; 	//threshold of disparity

	double epiSigma;
	double graySigma;
	double dispSigma;
	//weights of different terms
	double w_epi; //epipolar error
	double w_desc; //descriptor error
	double w_disp; //displacement error
	double w_none; //penalty for no correspondence
	double w_dummy; 
public:
	GuidedSSDMatcherParam() {
		epiDistThres = 4.5;
		grayLevelThres = 60;
		//descDistThres = 0.8/64;
		dispDistThres = 120;

		epiSigma = 1;
		graySigma = 2.0;
		//descSigma = 0.06*0.1/64;
		dispSigma = 40;

		w_epi = 1.0;
		w_desc = 1.0;
		w_disp = 1.0;
		w_none = 100000;
		w_dummy = 20000;
	}
};
class GuidedSSDMatcher : public GuidedSSDMatcherParam {
public:
	enum {
		TYPE_EPI = 0, //epipolar constraint
		TYPE_EPI_DESC = 1
	//epipolar constraint + descriptor
	};
	const Mat_d* m_pPoints[2];
	const Mat_d* m_pImgBlk[2];
//	const Mat_uc* m_pPtFlag[2];

	double F[9];

	Mat_d m_seedPts[2];
	Mat_d speed1Disp;
	Mat_d speed2Disp;

	Mat_d m_epiMat;
	Mat_d m_descMat;
	Mat_d m_dispMat;
	Mat_d m_weightMat;

	Mat_i m_ind[2];
public:
	GuidedSSDMatcher();
	virtual ~GuidedSSDMatcher();
protected:
	void _computeDescDistMat(const Mat_d& epiMat , const Mat_d& desc1 , const Mat_d& desc2 , Mat_d& descMat);
	void _computeDispDistMat(const Mat_d& epiMat , const Mat_d& pts1 , const Mat_d& pts2 , Mat_d& dispMat);
	void _getSeedDisparities(const Mat_d& seed1 , const Mat_d& seed2 , Mat_d& speed1Disp , Mat_d& speed2Disp);
	void _sumWeightMat2(const Mat_d& mat1 , const Mat_d& mat2 , double w1 , double w2 , Mat_d& weightMat);
	void _sumWeightMat3(const Mat_d& mat1 ,	const Mat_d& mat2 ,	const Mat_d& mat3 ,	double w1 ,	double w2 ,	double w3 ,	Mat_d& weightMat);
	void _getValidRowCol(Mat_d& weigthMat , Mat_c& rowFlag , Mat_c& colFlag);
public:
	void setSeeds(const Mat_d& seedPts1 , const Mat_d& seedPts2 , const Matching& matches) {
		getMatchedPts(matches, seedPts1, seedPts2, m_seedPts[0], m_seedPts[1]);
	}
	void setSeedPts(const Mat_d& seedPts1 , const Mat_d& seedPts2) {
		m_seedPts[0].cloneFrom(seedPts1);
		m_seedPts[1].cloneFrom(seedPts2);
	}

	void setFMat(const double* tF) {
		memcpy(F, tF, sizeof(double) * 9);
	}
	void setPoints(const Mat_d& pts1 , const Mat_d& gray1 , const Mat_d& pts2 , const Mat_d& gray2) {
		if (pts1.rows != gray1.rows || pts2.rows != gray2.rows)
			repErr("GuidedStereoMatcher::setPoints - error!");

		m_pPoints[0] = &pts1;
		m_pPoints[1] = &pts2;
		m_pImgBlk[0] = &gray1;
		m_pImgBlk[1] = &gray2;
	}
	int lapMatch(Matching& matches , int type = TYPE_EPI_DESC);
	//use only epipolar constraint
	int lapMatch(const double* tF , const Mat_d& pts1 , const Mat_d& pts2 , Matching& matches);
	//use epipolar constraint + descriptors
	int lapMatch(
			const double* tF ,
			const Mat_d& pts1 ,
			const Mat_d& pts2 ,
			const Mat_d& desc1 ,
			const Mat_d& desc2 ,
			Matching& matches);
	void removeOutlier(const Mat_d& pts1 , const Mat_d& pts2 , Matching& matches);
	void removeOutlier(
			const Mat_d& pts1 ,
			const Mat_d& pts2 ,
			Matching& matches ,
			double avgDx ,
			double avgDy ,
			double distThres);
};
#endif /* SL_GUIDEDSTEREOMATCHER_H_ */
