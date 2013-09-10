/*
 * LAP.h
 *
 *  Created on: 2009-8-14
 *      Author: Danping Zou
 */

#ifndef LAP_H_
#define LAP_H_

void printMatrix(double W[] , int m , int n , double* u = 0 , double* v = 0);
void printMatrix(double W[] , int m , int n , bool* rf , bool* cf);
void printResult(
		double W[] ,
		int m ,
		int n ,
		int* ar ,
		int* ac ,
		double* u = 0 ,
		double* v = 0 ,
		bool* rf = 0 ,
		bool *cf = 0);

int getProblematicRow();
/**
 * solve the lap by using successive shortest path(SSP) method
 */
double lapSSP(double W[] , int m , //mxn matrix  (m < n)
		int n ,
		int x[] , //columns assigned to rows
		int y[] , //rows assigned to columns
		double u[] , //dual row variables
		double v[] , //dual column variables
		double d[] , //shortest path lengths
		int frow[] , //free rows
		int pred[] , //predecessor nodes in the shortest path
		int q[] , //queue for store columns
		//	scaned (k = 1...low -1)
		//  labeled and unscanned (k = low ... up-1)
		//  unlabeled ( k = up ... n)
		char rf[] , // flags for valid rows (true - valid)
		char cf[] // flags for valid columns (true - valid)
		);

/**
 * one step of successive shortest path algorithm 
 */
double lapSSPOneStep(
		double W[] ,
		int m ,
		int n ,
		int x[] ,
		int y[] ,
		double u[] ,
		double v[] ,
		double d[] ,
		int row ,
		int pred[] ,
		int q[] ,
		char rf[] ,
		char cf[]);

double lap(double W[] , int m , //mxn matrix  (m < n)
		int n ,
		int x[] , //columns assigned to rows (-1 with no assignment)
		int y[] , //rows assigned to columns
		char rf[] , //flags for valid rows (0 - invalid, 1 - valid)
		char cf[] //flags for valid columns ( 0 - invalid, 1 - valid);
		);

#endif /* LAP_H_ */
