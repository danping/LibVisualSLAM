/*
 * LAP.cpp
 *
 *  Created on: 2009-8-14
 *      Author: Danping
 */
#include "LAP.h"
#include <cstdlib>
#include <cassert>
#include <cfloat>
#include <cstdio>
#include <cstring>
static int problematic_row = -1;
void printMatrix(double W[] , int m , int n , double* u , double* v) {
	int i, j;
	printf("Matrix:\n");
	if (v) {
		for (i = 0; i < n; ++i)
			printf("\t %.0f ", v[i]);
	}
	printf("\n\t");
	for (i = 0; i < n; ++i)
		printf("========");
	printf("\n");
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			double val = W[i * n + j];
			if (val < 0)
				printf("\t # ");
			else
				printf("\t %.1f ", val);
		}
		if (u)
			printf("\t | %.1f", u[i]);
		printf("\n");
	}
	printf("\t");
	for (i = 0; i < n; ++i)
		printf("========");
	printf("\n");
}
void printMatrix(double W[] , int m , int n , bool* rf , bool* cf) {
	int i, j;
	printf("Matrix:\n");
	if (cf) {
		for (i = 0; i < n; ++i)
			printf("\t %c ", cf[i] > 0 ? '*' : 'o');
	}
	printf("\n\t");
	for (i = 0; i < n; ++i)
		printf("========");
	printf("\n");
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			double val = W[i * n + j];
			if (val < 0)
				printf("\t # ");
			else
				printf("\t %.1f ", val);
		}
		if (rf)
			printf("\t | %c", rf[i] ? '*' : 'o');
		printf("\n");
	}
	printf("\t");
	for (i = 0; i < n; ++i)
		printf("========");
	printf("\n");
}
void print_result(double W[] , int m , int n , int* ar , int* ac , double* u , double* v) {
	int i, j;
	printf("Result:\n");
	if (v) {
		for (i = 0; i < n; ++i)
			printf("\t %.0f ", v[i]);
	}
	printf("\n\t");
	for (i = 0; i < n; ++i)
		printf("========");
	printf("\n");
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			double val = W[i * n + j];
			if (ar && ar[i] == j)
				printf("\t[%.0f]", val);
			else {
				if (val < 0)
					printf("\t # ");
				else
					printf("\t %.0f ", val);
			}
		}
		if (u)
			printf("\t | %.0f", u[i]);
		printf("\n");
	}
	printf("\t");
	for (i = 0; i < n; ++i)
		printf("========");
	printf("\n");
}

void printResult(double W[] , int m , int n , int* ar , int* ac , double* u , double* v , bool* rf , bool* cf) {
	int i, j;
	printf("Result:\n");
	if (v) {
		for (i = 0; i < n; ++i)
			printf("\t %.0f ", v[i]);
	}
	printf("\n");
	if (cf) {
		for (i = 0; i < n; ++i)
			printf("\t %c ", cf[i] > 0 ? '*' : 'o');
	}
	printf("\n\t");
	for (i = 0; i < n; ++i)
		printf("========");
	printf("\n");
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			double val = W[i * n + j];
			if (ar && ar[i] == j)
				printf("\t[%.0f]", val);
			else {
				if (val < 0)
					printf("\t # ");
				else
					printf("\t %.0f ", val);
			}
		}
		if (rf)
			printf("\t%c", rf[i] > 0 ? '*' : 'o');
		if (u)
			printf("\t | %.0f", u[i]);
		printf("\n");
	}
	printf("\t");
	for (i = 0; i < n; ++i)
		printf("========");
	printf("\n");
}
int getProblematicRow() {
	return problematic_row;
}
double lapSSP(double W[] , //Cost matrix
		int m , // m x n  
		int n ,
		int x[] , //columns assigned to rows
		int y[] , //rows assigned to columns
		double u[] , //dual row variables
		double v[] , //dual column variables
		double d[] , //shortest path lengths
		int frow[] , //free row
		int pred[] , //predecessor nodes in the shortest path
		int q[] , //queue for store columns
		//	scaned (k = 1...low -1)
		//  labeled and unscanned (k = low ... up-1)
		//  unlabeled ( k = up ... n)
		char rf[] , // flags of valid rows (true - valid)
		char cf[] // flags of valid columns (true - valid)
) {
	int i, j, nvr, nvc, f, nf, r, c, low, up, last, uc, r1, c1, u1;
	double dist, min_dist, cost, min_cost;
	int* vr = new int[m];	//valid row
	int* vc = new int[n];	//valid column

	nvr = 0, nvc = 0;
	for (i = 0; i < m; ++i) {
		if (rf[i] > 0) {
			vr[nvr++] = i;
			x[i] = -1;
		}
	}
	for (i = 0; i < n; ++i) {
		if (cf[i] > 0) {
			vc[nvc++] = i;
			y[i] = -1;
		}
	}

	//column reduction
	if (nvr == nvc) {
		for (j = nvc - 1; j >= 0; --j) {
			c = vc[j];
			min_cost = DBL_MAX;
			int min_r = -1;
			for (i = 0; i < nvr; ++i) {
				r = vr[i];
				cost = W[r * n + c];
				if (cost >= 0 && cost < min_cost) {
					min_cost = cost;
					min_r = r;
				}
			}
			if (min_r < 0) {
				min_cost = -1;
				goto end;
			}
			v[c] = min_cost;
			if (x[min_r] < 0) {
				x[min_r] = c;
				y[c] = min_r;
			}
		}
	} else {
		for (c = 0; c < nvc; ++c) {
			v[vc[c]] = 0;
		}
		for (i = 0; i < nvr; ++i) {
			r = vr[i];
			min_cost = DBL_MAX;
			int min_c = -1;
			for (j = nvc - 1; j >= 0; --j) {
				c = vc[j];
				cost = W[r * n + c];
				if (cost >= 0 && cost < min_cost) {
					min_cost = cost;
					min_c = c;
				}
			}
			if (min_c < 0) {
				min_cost = -1;
				goto end;
			}
			if (y[min_c] < 0) {
				y[min_c] = r;
				x[r] = min_c;
			}
		}
	}
	//get free rows
	nf = 0;
	for (i = 0; i < nvr; ++i) {
		r = vr[i];
		if (x[r] < 0)
			frow[nf++] = r;
	}
	//start successively finding shortest path in the residual graph
	for (f = 0; f < nf; ++f) {
		r = frow[f];
		//initialize the shortest paths
		for (j = 0; j < nvc; ++j) {
			c = vc[j];
			cost = W[r * n + c];
			if (cost < 0)
				d[c] = DBL_MAX;
			else
				d[c] = cost - v[c]; //u[r] are set to zero initially  
			pred[c] = r;
			q[j] = c;
		}
		//start dijkstra shortest path algorithm
		low = up = 0;
		uc = -1; //unassigned column
		while (uc < 0) {
			if (low == up) {
				if (low == nvr || low == nvc) {
					problematic_row = r;
					min_cost = -1;
					goto end;
				}
				//find column with new value of minimum d
				last = low - 1;
				min_dist = d[q[up++]];
				for (j = up; j < nvc; ++j) {
					c = q[j];
					dist = d[c];
					if (dist <= min_dist) {
						if (dist < min_dist) {
							up = low;
							min_dist = dist;
						}
						q[j] = q[up];
						q[up++] = c;
					}
				}
				for (j = low; j < up; ++j) {
					c = q[j];
					if (y[c] < 0) {
						//the shortest path arrives at an unassigned column
						//go to augmentation
						uc = c;
						break;
					}
				}
			}//end if(low == up)
			if (uc < 0) {
				c1 = q[low++];
				r1 = y[c1];
				u1 = W[r1 * n + c1] - v[c1] - min_dist;
				for (j = up; j < nvc; ++j) {
					c = q[j];
					cost = W[r1 * n + c];
					if (cost < 0)
						continue;
					dist = cost - v[c] - u1;
					if (d[c] == DBL_MAX || dist < d[c]) {
						d[c] = dist; //update the distance
						pred[c] = r1; //recored the predecessor row
						if (dist == min_dist) {
							if (y[c] < 0) {
								//go to augmentation immediately
								uc = c;
								break;
							}
							q[j] = q[up];
							q[up++] = c;
						}
					}
				}
			}
		}//end while(uc < 0)
		assert(uc >= 0);
		//update column variables
		for (j = 0; j <= last; ++j) {
			c = q[j];
			v[c] = v[c] + d[c] - min_dist;
		}
		//augmentation
		do {
			r1 = pred[uc];
			y[uc] = r1;
			c = uc;
			uc = x[r1];
			x[r1] = c;
		} while (r1 != r);
	}//end for( f = 0; f < nf; ++f)

	min_cost = 0;
	for (i = 0; i < nvr; ++i) {
		r = vr[i];
		c = x[r];
		cost = W[r * n + c];
		if (cost < 0) {
			problematic_row = r;
			min_cost = -1;
			goto end;
		}
		u[r] = cost - v[c];
		min_cost += cost;
	}
	end: delete[] vr;
	delete[] vc;
	return min_cost;
}
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
		char cf[]) {

	int i, j, nvc, r, c, low, up, last, uc, r1, c1, u1;
	double dist, min_dist, cost, min_cost;
	int* vc = new int[n];
	nvc = 0;
	for (i = 0; i < n; ++i) {
		if (cf[i] > 0) {
			vc[nvc++] = i;
		}
	}
	y[x[row]] = -1;
	x[row] = -1;
	r = row;
	//initialize the shortest paths
	for (j = 0; j < nvc; ++j) {
		c = vc[j];
		cost = W[r * n + c];
		if (cost < 0)
			d[c] = DBL_MAX;
		else
			d[c] = cost - v[c] - u[r];
		pred[c] = r;
		q[j] = c;
	}
	//start dijkstra shortest path algorithm
	low = up = 0;
	uc = -1; //unassigned column
	while (uc < 0) {
		if (low == up) {
			//find column with new value of minimum d
			last = low - 1;
			min_dist = d[q[up++]];
			for (j = up; j < nvc; ++j) {
				c = q[j];
				dist = d[c];
				if (dist <= min_dist) {
					if (dist < min_dist) {
						up = low;
						min_dist = dist;
					}
					q[j] = q[up];
					q[up++] = c;
				}
			}
			for (j = low; j < up; ++j) {
				c = q[j];
				if (y[c] < 0) {
					//the shortest path arrives at an unassigned column
					//go to augmentation
					uc = c;
					break;
				}
			}
		}//end if(low == up)
		if (uc < 0) {
			c1 = q[low++];
			r1 = y[c1];
			u1 = W[r1 * n + c1] - v[c1] - min_dist;
			for (j = up; j < nvc; ++j) {
				c = q[j];
				cost = W[r1 * n + c];
				if (cost < 0)
					continue;
				dist = cost - v[c] - u1;
				if (d[c] == DBL_MAX || dist < d[c]) {
					d[c] = dist; //update the distance
					pred[c] = r1; //recored the predecessor row
					if (dist == min_dist) {
						if (y[c] < 0) {
							//go to augmentation immediately
							uc = c;
							break;
						}
						q[j] = q[up];
						q[up++] = c;
					}
				}
			}
		}
	}//end while(uc < 0)
	//update column variables
	for (j = 0; j <= last; ++j) {
		c = q[j];
		v[c] = v[c] + d[c] - min_dist;
	}
	//augmentation
	do {
		r1 = pred[uc];
		y[uc] = r1;
		c = uc;
		uc = x[r1];
		x[r1] = c;
	} while (r1 != r);

	min_cost = 0;
	for (j = 0; j < nvc; ++j) {
		c = vc[j];
		r = y[c];
		if (r < 0 || rf[r] <= 0)
			continue;
		cost = W[r * n + c];
		if (cost < 0) {
			min_cost = -1;
			goto end;
		}
		u[r] = cost - v[c];
		min_cost += cost;
	}
	end: delete[] vc;
	return min_cost;
}

double lap(double W[] , int m , //mxn matrix  (m < n)
		int n ,
		int x[] , //columns assigned to rows 
		int y[] , //rows assigned to columns
		char rf[] , //flags for valid rows (0 - invalid, 1 - valid)
		char cf[] //flags for valid columns ( 0 - invalid, 1 - valid)
) {
	//temporal variables
	double* u = new double[m];
	double* v = new double[n];
	double* d = new double[n];
	int* frow = new int[m];
	int* pred = new int[n];
	int* q = new int[n];

	double s = lapSSP(W, m, n, x, y, u, v, d, frow, pred, q, rf, cf);

	delete[] u;
	delete[] v;
	delete[] d;
	delete[] frow;
	delete[] pred;
	delete[] q;

	return s;
}

