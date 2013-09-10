/*
 * SL_Graph.h
 *
 *  Created on: 2011-7-6
 *      Author: Danping Zou
 */

#ifndef SL_GRAPH_H_
#define SL_GRAPH_H_
#include <cassert>
#include "math/SL_Matrix.h"
/**
 * find connected components defined by a symmetric matrix M = M^T
 */
template<typename T>
void findConnectedComponents(const T& M, std::vector<Mat_i>& conns) {
	assert(M.m == M.n);
	int num = M.m;
	conns.clear();
	int* VQ = new int[num];
	int* CON = new int[num];
	bool* flag = new bool[num];
	for (int i = 0; i < num; i++) {
		if (0 == flag[i]) {
			//start to grow
			int nVQ = 0, nCon = 0;
			nVQ = nCon = 0;
			CON[nCon++] = i;
			VQ[nVQ++] = i;
			flag[i] = true;
			while (nVQ > 0) {
				int j0 = VQ[nVQ - 1];
				nVQ--;
				for (int j = 0; j < num; j++) {
					if (j != j0 && !flag[j] && M.data[j0 * num + j] > 0) {
						CON[nCon++] = j;
						VQ[nVQ++] = j;
						flag[j] = true;
					}
				}
			}
			//output the connected components
			Mat_i conn(1, nCon);
			for (int k = 0; k < nCon; k++)
				conn.data[k] = CON[k];
			conns.push_back(conn);
		}

	}
	delete[] VQ;
	delete[] CON;
	delete[] flag;
}
#endif /* SL_GRAPH_H_ */
