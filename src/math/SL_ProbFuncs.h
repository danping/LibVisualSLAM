/*
 * SL_ProbFuncs.h
 *
 *  Created on: 2011-10-20
 *      Author: Danping Zou
 */

#ifndef PROBFUNCS_H_
#define PROBFUNCS_H_
/* 
 * compute normal distribution 
 */
double normpdf2(const double x[], const double miu[], const double sigma[]);
double normpdf2(const double x[], const double sigma[]);
double normpdf3(const double x[], const double miu[], const double sigma[]);
double normpdf3(const double x[], const double sigma[]);

/* compute (2*pi)^{-dim/2} det(\sigma)^{1/2}*/
double normpdf_const(int dim, const double invsigma[]);
/* compute exp(-0.5 x'*invsigma*x)*/
double normpdf_expval(int dim, const double x[], const double invsigma[]);

double normpdf(int dim, const double x[], const double sigma[]);
#endif /* PROBFUNCS_H_ */
