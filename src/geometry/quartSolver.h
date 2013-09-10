/*
 * quartic.h
 *
 *  Created on: 2010-12-2
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef QUARTIC_H_
#define QUARTIC_H_

int quadratic(double b, double c, double rts[]);
int cubic(double p, double q, double r, double v3[]);
void cubnewton(double p, double q, double r, int n3, double v3[]);
int quartic(double a, double b, double c, double d, double rts[]);
int descartes(double a, double b, double c, double d, double rts[]);
int ferrari(double a, double b, double c, double d, double rts[]);
int neumark(double a, double b, double c, double d, double rts[]);
int yacfraid(double a, double b, double c, double d, double rts[]);
int chris(double a, double b, double c, double d, double rts[]);
#endif /* QUARTIC_H_ */
