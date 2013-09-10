/////////////////////////////////////////////////////////////////////////////////
//// 
////  CRS sparse matrices manipulation routines
////  Copyright (C) 2004-2008 Manolis Lourakis (lourakis at ics forth gr)
////  Institute of Computer Science, Foundation for Research & Technology - Hellas
////  Heraklion, Crete, Greece.
////
////  This program is free software; you can redistribute it and/or modify
////  it under the terms of the GNU General Public License as published by
////  the Free Software Foundation; either version 2 of the License, or
////  (at your option) any later version.
////
////  This program is distributed in the hope that it will be useful,
////  but WITHOUT ANY WARRANTY; without even the implied warranty of
////  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
////  GNU General Public License for more details.
////
///////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>

#include "sba.h"

static void sba_crsm_print(struct sba_crsm *sm, FILE *fp);
static void sba_crsm_build(struct sba_crsm *sm, int *m, int nr, int nc);

/* allocate a sparse CRS matrix */
void sba_crsm_alloc(struct sba_crsm *sm, int nr, int nc, int nnz) {
	int msz;

	sm->nr = nr;
	sm->nc = nc;
	sm->nnz = nnz;
	msz = 2 * nnz + nr + 1;
	sm->val = (int *) malloc(msz * sizeof(int)); /* required memory is allocated in a single step */
	if (!sm->val) {
		fprintf(stderr, "SBA: memory allocation request failed in sba_crsm_alloc() [nr=%d, nc=%d, nnz=%d]\n", nr, nc,
				nnz);
		exit(1);
	}
	sm->colidx = sm->val + nnz;
	sm->rowptr = sm->colidx + nnz;
}

/* free a sparse CRS matrix */
void sba_crsm_free(struct sba_crsm *sm) {
	sm->nr = sm->nc = sm->nnz = -1;
	free(sm->val);
	sm->val = sm->colidx = sm->rowptr = NULL;
}

static void sba_crsm_print(struct sba_crsm *sm, FILE *fp) {
	register int i;

	fprintf(fp, "matrix is %dx%d, %d non-zeros\nval: ", sm->nr, sm->nc, sm->nnz);
	for (i = 0; i < sm->nnz; ++i)
		fprintf(fp, "%d ", sm->val[i]);
	fprintf(fp, "\ncolidx: ");
	for (i = 0; i < sm->nnz; ++i)
		fprintf(fp, "%d ", sm->colidx[i]);
	fprintf(fp, "\nrowptr: ");
	for (i = 0; i <= sm->nr; ++i)
		fprintf(fp, "%d ", sm->rowptr[i]);
	fprintf(fp, "\n");
}

/* build a sparse CRS matrix from a dense one. intended to serve as an example for sm creation */
static void sba_crsm_build(struct sba_crsm *sm, int *m, int nr, int nc) {
	int nnz;
	register int i, j, k;

	/* count nonzeros */
	for (i = nnz = 0; i < nr; ++i)
		for (j = 0; j < nc; ++j)
			if (m[i * nc + j] != 0)
				++nnz;

	sba_crsm_alloc(sm, nr, nc, nnz);

	/* fill up the sm structure */
	for (i = k = 0; i < nr; ++i) {
		sm->rowptr[i] = k;
		for (j = 0; j < nc; ++j)
			if (m[i * nc + j] != 0) {
				sm->val[k] = m[i * nc + j];
				sm->colidx[k++] = j;
			}
	}
	sm->rowptr[nr] = nnz;
}

/* returns the index of the (idx1, idx2) element. No bounds checking! */
int sba_crsm_elmidx(struct sba_crsm *sm, int i, int j) {
	register int low, high, mid, diff;

	low = sm->rowptr[i];
	high = sm->rowptr[i + 1] - 1;

	/* binary search for finding the element at column idx2 */
	while (low <= high) {
		/* following early termination test seems to actually slow down the search */
		//if(idx2<sm->colidx[low] || idx2>sm->colidx[high]) return -1; /* not found */

		/* mid=low+((high-low)>>1) ensures no index overflows */
		mid = (low + high) >> 1; //(low+high)/2;
		diff = j - sm->colidx[mid];
		if (diff < 0)
			high = mid - 1;
		else if (diff > 0)
			low = mid + 1;
		else
			return mid;
	}

	return -1; /* not found */
}

/* similarly to sba_crsm_elmidx() above, returns the index of the (idx1, idx2) element using the
 * fact that the index of element (idx1, jp) was previously found to be jpidx. This can be
 * slightly faster than sba_crsm_elmidx(). No bounds checking!
 */
int sba_crsm_elmidxp(struct sba_crsm *sm, int i, int j, int jp, int jpidx) {
	register int low, high, mid, diff;

	diff = j - jp;
	if (diff > 0) {
		low = jpidx + 1;
		high = sm->rowptr[i + 1] - 1;
	} else if (diff == 0)
		return jpidx;
	else { /* diff<0 */
		low = sm->rowptr[i];
		high = jpidx - 1;
	}

	/* binary search for finding the element at column idx2 */
	while (low <= high) {
		/* following early termination test seems to actually slow down the search */
		//if(idx2<sm->colidx[low] || idx2>sm->colidx[high]) return -1; /* not found */

		/* mid=low+((high-low)>>1) ensures no index overflows */
		mid = (low + high) >> 1; //(low+high)/2;
		diff = j - sm->colidx[mid];
		if (diff < 0)
			high = mid - 1;
		else if (diff > 0)
			low = mid + 1;
		else
			return mid;
	}

	return -1; /* not found */
}

/* returns the number of nonzero elements in row idx1 and
 * fills up the vidxs and jidxs arrays with the val and column
 * indexes of the elements found, respectively.
 * vidxs and jidxs are assumed preallocated and of max. size sm->nc
 */
int sba_crsm_row_elmidxs(struct sba_crsm *sm, int i, int *vidxs, int *jidxs) {
	register int j, k;

	for (j = sm->rowptr[i], k = 0; j < sm->rowptr[i + 1]; ++j, ++k) {
		vidxs[k] = j;
		jidxs[k] = sm->colidx[j];
	}

	return k;
}

/* returns the number of nonzero elements in col idx2 and
 * fills up the vidxs and iidxs arrays with the val and row
 * indexes of the elements found, respectively.
 * vidxs and iidxs are assumed preallocated and of max. size sm->nr
 */
int sba_crsm_col_elmidxs(struct sba_crsm *sm, int j, int *vidxs, int *iidxs) {
	register int *nextrowptr = sm->rowptr + 1;
	register int i, l;
	register int low, high, mid, diff;

	for (i = l = 0; i < sm->nr; ++i) {
		low = sm->rowptr[i];
		high = nextrowptr[i] - 1;

		/* binary search attempting to find an element at column idx2 */
		while (low <= high) {
			//if(idx2<sm->colidx[low] || idx2>sm->colidx[high]) break; /* not found */

			mid = (low + high) >> 1; //(low+high)/2;
			diff = j - sm->colidx[mid];
			if (diff < 0)
				high = mid - 1;
			else if (diff > 0)
				low = mid + 1;
			else { /* found */
				vidxs[l] = mid;
				iidxs[l++] = i;
				break;
			}
		}
	}

	return l;
}

/* a more straighforward (but slower) implementation of the above function */
/***
 int sba_crsm_col_elmidxs(struct sba_crsm *sm, int idx2, int *vidxs, int *iidxs)
 {
 register int idx1, k, pyramidLevel;

 for(idx1=pyramidLevel=0; idx1<sm->nr; ++idx1)
 for(k=sm->rowptr[idx1]; k<sm->rowptr[idx1+1]; ++k)
 if(sm->colidx[k]==idx2){
 vidxs[pyramidLevel]=k;
 iidxs[pyramidLevel++]=idx1;
 }

 return pyramidLevel;
 }
 ***/

#if 0
/* returns 1 if there exists a row idx1 having columns idx2 and k,
 * idx1.e. a row idx1 s.t. elements (idx1, idx2) and (idx1, k) are nonzero;
 * 0 otherwise
 */
int sba_crsm_common_row(struct sba_crsm *sm, int idx2, int k)
{
	register int idx1, low, high, mid, diff;

	if(idx2==k) return 1;

	for(idx1=0; idx1<sm->nr; ++idx1) {
		low=sm->rowptr[idx1];
		high=sm->rowptr[idx1+1]-1;
		if(idx2<sm->colidx[low] || idx2>sm->colidx[high] || /* idx2 not found */
				k<sm->colidx[low] || k>sm->colidx[high]) /* k not found */
		continue;

		/* binary search for finding the element at column idx2 */
		while(low<=high) {
			mid=(low+high)>>1; //(low+high)/2;
			diff=idx2-sm->colidx[mid];
			if(diff<0)
			high=mid-1;
			else if(diff>0)
			low=mid+1;
			else
			goto jfound;
		}

		continue; /* idx2 not found */

		jfound:
		if(idx2>k) {
			low=sm->rowptr[idx1];
			high=mid-1;
		}
		else {
			low=mid+1;
			high=sm->rowptr[idx1+1]-1;
		}

		if(k<sm->colidx[low] || k>sm->colidx[high]) continue; /* k not found */

		/* binary search for finding the element at column k */
		while(low<=high) {
			mid=(low+high)>>1; //(low+high)/2;
			diff=k-sm->colidx[mid];
			if(diff<0)
			high=mid-1;
			else if(diff>0)
			low=mid+1;
			else /* found */
			return 1;
		}
	}

	return 0;
}
#endif

#if 0

/* sample code using the above routines */

main()
{
	int mat[7][6]= {
		{	10, 0, 0, 0, -2, 0},
		{	3, 9, 0, 0, 0, 3},
		{	0, 7, 8, 7, 0, 0},
		{	3, 0, 8, 7, 5, 0},
		{	0, 8, 0, 9, 9, 13},
		{	0, 4, 0, 0, 2, -1},
		{	3, 7, 0, 9, 2, 0}
	};

	struct sba_crsm sm;
	int idx1, idx2, k, pyramidLevel;
	int vidxs[7], /* max(6, 7) */
	jidxs[6], iidxs[7];

	sba_crsm_build(&sm, mat[0], 7, 6);
	sba_crsm_print(&sm, stdout);

	for(idx1=0; idx1<7; ++idx1) {
		for(idx2=0; idx2<6; ++idx2)
		printf("%3d ", ((k=sba_crsm_elmidx(&sm, idx1, idx2))!=-1)? sm.val[k] : 0);
		printf("\n");
	}

	for(idx1=0; idx1<7; ++idx1) {
		k=sba_crsm_row_elmidxs(&sm, idx1, vidxs, jidxs);
		printf("row %d\n", idx1);
		for(pyramidLevel=0; pyramidLevel<k; ++pyramidLevel) {
			idx2=jidxs[pyramidLevel];
			printf("%d %d  ", idx2, sm.val[vidxs[pyramidLevel]]);
		}
		printf("\n");
	}

	for(idx2=0; idx2<6; ++idx2) {
		k=sba_crsm_col_elmidxs(&sm, idx2, vidxs, iidxs);
		printf("col %d\n", idx2);
		for(pyramidLevel=0; pyramidLevel<k; ++pyramidLevel) {
			idx1=iidxs[pyramidLevel];
			printf("%d %d  ", idx1, sm.val[vidxs[pyramidLevel]]);
		}
		printf("\n");
	}

	sba_crsm_free(&sm);
}
#endif
