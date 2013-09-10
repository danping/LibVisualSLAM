#ifndef LAPACK_H_
#define LAPACK_H_

#ifdef __cplusplus 	
extern "C" {
#endif		

int dgelsy_(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int * jpvt, double *rcond, int *rank, double *work, int * lwork, int *info);

int dgetrf_(int *m, int *n, double *a, int * lda, int *ipiv, int *info);

int dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);

void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);

void dptsv_(int* n, int* nrhs, double* d, double* e, double* b, int* ldb, int* info);

void dposv_(char* uplo, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, int* info);

void dgeev_(char *jobvl, char *jobvr, int *n, double *A, int *lda, double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);

void dgeevx_(char balanc, char jobvl, char jobvr, char sense, int n, double *a, int lda, double *wr, double *wi, double *vl, int ldvl, double *vr, int ldvr, int *ilo, int *ihi, double *scale, double *abnrm, double *rcone, double *rconv, int *info);

void dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *A, int *lda, double *S, double *U, int *ldu, double *VT, int *ldvt, double *work, int *lwork, int *info);

void dgemm_(char*, char*, int*, int *, int *, double *, double *, int*, double *, int *, double *, double *, int*);

int dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

int dorgqr_(int *m, int *n, int *k, double * a, int *lda, double *tau, double *work, int *lwork, int *info);

int dgemv_(char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

#ifdef __cplusplus
}
#endif

#endif /* LAPACK_H_ */

