/*
 * SL_8point.h
 *
 *  Created on: 2010-11-14
 *      Author: Danping Zou
 */

#ifndef SL_8POINT_H_
#define SL_8POINT_H_
void computeFMat8Pt(int n , const double* a , const double* b , double* E);
int findFMatRansac(
		const int n ,
		const double* a ,
		const double* b ,
		double* F ,
		int iterMaxNum ,
		double epiErrThes ,
		double* epiErr);
#endif /* SL_8POINT_H_ */
