/*
 * SL_Debug.h
 *
 *  Created on: 2010-11-6
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_DEBUG_H_
#define SL_DEBUG_H_
#include "math/SL_Matrix.h"
#include "math/SL_LinAlg.h"

#include "geometry/SL_Geometry.h"
#include "geometry/SL_Point.h"
#include "geometry/SL_Distortion.h"
#include "geometry/SL_Triangulate.h"

#include "imgproc/SL_Image.h"

#include "SL_Print.h"
#include "SL_WriteRead.h"
#include "SL_DrawCorners.h"
#include "SL_Tictoc.h"

#include "GUI_ImageViewer.h"
template<class T>
int check_same_mat(T& mat1, T& mat2) {
	if (mat1.rows != mat2.rows || mat1.cols != mat2.cols)
		return -1;
	int len = mat1.rows * mat1.cols;
	for (int i = 0; i < len; ++i) {
		if (mat1.data[i] != mat2.data[i]) {
			int y = i / mat1.cols;
			int x = i - y * mat1.cols;
			logInfo("(%d,%d)%d != %d\n", x, y, mat1.data[i], mat2.data[i]);
			return -i;
		}
	}
	return 1;
}
void draw_reprojection_points(const double* K,
		const double* kc,
		int v_id,
		int nv,
		const double* Rs,
		const double* ts,
		int npts,
		const double* pts,
		const double* meas,
		const char* vmask);

#endif /* SL_DEBUG_H_ */
