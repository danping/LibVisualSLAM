/*
 * SL_p3p.h
 *
 *  Created on: 2010-12-2
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#ifndef SL_P3P_H_
#define SL_P3P_H_

/*three point pose estimation*/

/* get coeffs for equation: dij^2 = x_i^2 + x_j^2 - 2 x_i*x_j*cos(theta_ij)
 * input:
 * 	K : intrinsic matrix of the camera
 * 	M1,M2,M3 : three scene points. Each one is 3x1 vector
 * 	m1,m2,m2 : projections of the three points (2x1 vectors)
 * output:
 * 	d[3] : store d12 d13 d23
 * 	cs[3] : store cos(theta_12),cos(theta_13),cos(theta_23)
 */
void p3p_get_coeffs(const double K[], const double M1[], const double M2[], const double M3[], const double m1[],
		const double m2[], const double m3[], double d[], double cs[]);
/* compute the quartic equation for x_1^2:
 * a*x_1^8 + b *x_1^6 + c*x_1^4 + d*x_1^2 + e = 0
 * input:
 * 	d[3] : d12,d13,d23
 *  cs[3] : cos(theta_12),cos(theta_13),cos(theta_23)
 * output:
 *  cfs[5] : a,b,c,d,e
 */
void p3p_get_quartic(const double d[], const double cs[], double cfs[]);

/**
 * solve the equations to get x1,x2,x3
 * d12^2 = x1^2 + x2^2 - 2*x1*x2*cos(theta_12)
 * d13^2 = x1^2 + x2^2 - 2*x1*x3*cos(theta_13)
 * d23^2 = x2^2 + x3^2 - 2*x2*x3*cos(theta_23)
 * input:
 * 	d[3] : d12,d13,d23
 *  cs[3] : cos(theta_12),cos(theta_13),cos(theta_23)
 * output:
 *  x[] : x1,x2,x3
 *  err[] : error
 *  @return : number of solutions
 */
int p3p_solve_xs(const double d[], const double cs[], const double cfs[], double x[], double err[]);

/** 
 * select two best solutions and refine it 
 */
int p3p_select_best_two(const int nx, const double x[], const double err[], double bx[], double berr[]);

/**
 * three point algorithm to estimate the camera pose based on the following papers
 *  
 * @article{quan2002linear,
  title={{Linear n-point camera pose determination}},
  author={Quan, L. and Lan, Z.},
  journal={Pattern Analysis and Machine Intelligence, IEEE Transactions on},
  volume={21},
  number={8},
  pages={774--780},
  issn={0162-8828},
  year={2002},
  publisher={IEEE}
}
 * @article{horn1988closed,
  title={{Closed-form solution of absolute orientation using orthonormal matrices}},
  author={Horn, B.K.P. and Hilden, H.M. and Negahdaripour, S.},
  journal={Journal of the Optical Society of America A},
  volume={5},
  number={7},
  pages={1127--1135},
  issn={1084-7529},
  year={1988},
  publisher={Citeseer}
}
 *  
 */
int p3p(const double K[], const double M1[], const double M2[], const double M3[], const double m1[],
		const double m2[], const double m3[], double Rs[], double ts[]);

#endif /* SL_P3P_H_ */
