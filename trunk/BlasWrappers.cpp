/*
 *  BlasWrappers.cpp
 *
 *  Created by Spyros Schismenos on 07/06/2012.
 *
 */ 

// some documentation can be found here 
// https://developer.apple.com/library/mac/#documentation/Accelerate/Reference/BLAS_Ref/Reference/reference.html

#include "BlasWrappers.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_sort_vector.h"


double norm2(const gsl_vector* x)
{
	return cblas_dnrm2(x->size, x->data,x->stride);
}

void My_dgemm(const enum CBLAS_TRANSPOSE TransA,
			  const enum CBLAS_TRANSPOSE TransB, const double alpha, const gsl_matrix *A,
			  const gsl_matrix* B,
			  const double beta, gsl_matrix* C)
{
	int m,n,k;
	if (TransA == CblasNoTrans)
	{
		m = A->size1;
		k = A->size2;
	}
	else {
		m = A->size2;
		k = A->size1;
	}
	if (TransB ==CblasNoTrans)
		n = B->size2;
	else
		n = B->size1;
	cblas_dgemm(CblasRowMajor, TransA, TransB, m, n, k,alpha, A->data, A->tda, B->data, B->tda, beta, C->data, C->tda);
	
}

void My_dgemv(const enum CBLAS_TRANSPOSE TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
{
	int m,n;
	if (TransA == CblasTrans)
	{
		m = A->size1;
		n = A->size2;
	}
	else {
		m = A->size2;
		n = A->size1;
	}
	cblas_dgemv(CblasRowMajor,TransA, m, n, alpha,A->data, A->tda, x->data,x->stride, beta,y->data, y->stride);
}

void My_dsyrk(const enum CBLAS_UPLO Uplo,const enum CBLAS_TRANSPOSE Trans, double alpha, const gsl_matrix * A, double beta, gsl_matrix * C)
{
	int K;
	if (Trans == CblasNoTrans)
		K = A->size2;
	else
		K = A->size1;
	
	int N = C->size1;
	cblas_dsyrk(CblasRowMajor, Uplo, Trans, N, K, alpha, A->data, A->tda, beta, C->data, C->tda);
}

void My_dsyr(const enum CBLAS_UPLO Uplo,const enum CBLAS_TRANSPOSE Trans, double alpha, const gsl_matrix * A, double beta, gsl_matrix * C)
{
	int N = C->size1;
	cblas_dsyr(CblasRowMajor, Uplo, N, alpha, A->data, A->tda, C->data, C->tda);
}

void My_dtrsm(const enum CBLAS_SIDE Side,const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, 
			  double alpha, const gsl_matrix * A, gsl_matrix * B)
{
	int M = B->size1;
	int N = B->size2;
	cblas_dtrsm(CblasRowMajor, Side , Uplo , TransA ,Diag ,M ,N, alpha, A->data, A->tda, B->data, B->tda);
}

void My_dtrsv(const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, 
			const gsl_matrix * A, gsl_vector* b)
{
	int N = b->size;
	cblas_dtrsv(CblasRowMajor, Uplo , TransA ,Diag, N, A->data, A->tda, b->data, b->stride);
}

void My_drotg(double* a, double* b, double* c, double* s)
{
	cblas_drotg(a,b,c,s);
}

void My_drot(gsl_vector* x, gsl_vector* y, const double c, const double s)
{
	cblas_drot (x->size, x->data,x->stride,y->data,y->stride,c,s);
}

void My_dger(double alpha, const gsl_vector * x, const gsl_vector * y, gsl_matrix * A)
{
	cblas_dger(CblasRowMajor, A->size1, A->size2, alpha, x->data, x->stride, y->data, y->stride, A->data, A->tda);
}

double My_ddot(const gsl_vector * x, const gsl_vector* y)
{
	return cblas_ddot(x->size, x->data, x->stride, y->data, y->stride);
}

void orderMatrix(const gsl_matrix* x, gsl_matrix* y)
{
	int n = x->size1;
	int m = x->size2;
	gsl_vector* x_norms = gsl_vector_alloc(m);
	for	(int i =0;i<m;i++)
	{
		gsl_vector_const_view xcol = gsl_matrix_const_column(x,i);
		gsl_vector_set(x_norms, i, -norm2(&xcol.vector));
	}
	gsl_permutation* p = gsl_permutation_alloc(m);
	gsl_sort_vector_index(p, x_norms);
	for (int i=0; i<n; i++) {
		for (int j=0; j<m; j++) {
			gsl_matrix_set(y, i, j, gsl_matrix_get(x, i, gsl_permutation_get(p, j)));
		}
	}
	gsl_vector_free(x_norms);
	gsl_permutation_free(p);
}

void orderMatrix(const gsl_matrix* x, gsl_matrix* y, const gsl_matrix* M)
{
	int n = x->size1;
	int m = x->size2;
	gsl_matrix* invM = gsl_matrix_alloc(n,n);
	gsl_matrix_memcpy(invM,M);	
	int info=0;
	char lower = 'U';
	int lda = invM->tda;
	dpotrf_(&lower, &n, invM->data, &lda, &info);
	dpotri_(&lower, &n, invM->data, &lda, &info);
	for (int i=0; i<n; i++) {
		for (int j=i+1 ; j<n; j++) {
			gsl_matrix_set(invM,i,j,gsl_matrix_get(invM,j,i)) ;
		}
	}
	gsl_vector* x_ell_norms = gsl_vector_alloc(m);
	gsl_vector* temp = gsl_vector_alloc(n);
	for	(int i =0;i<m;i++)
	{
		gsl_vector_const_view xcol = gsl_matrix_const_column(x,i);
		My_dgemv(CblasNoTrans, 1.0, invM, &xcol.vector, 0.0, temp);
		gsl_vector_set(x_ell_norms, i, -My_ddot(&xcol.vector, temp));
	}
	gsl_permutation* p = gsl_permutation_alloc(m);
	gsl_sort_vector_index(p, x_ell_norms);
	for (int i=0; i<n; i++) {
		for (int j=0; j<m; j++) {
			gsl_matrix_set(y, i, j, gsl_matrix_get(x, i, gsl_permutation_get(p, j)));
		}
	}
	gsl_vector_free(x_ell_norms);
	gsl_vector_free(temp);
	gsl_matrix_free(invM);
	gsl_permutation_free(p);
	
}

double blas_sum(const gsl_vector* x)
{
	double res = 0.0;
	gsl_vector* ones = gsl_vector_alloc(x->size); gsl_vector_set_all(ones, 1.0);
	res = My_ddot(x, ones);
	gsl_vector_free(ones);
	return res;
}

double My_dasum(const gsl_vector *x)
{
	return cblas_dasum(x->size, x->data, x->stride);
}

void My_daxpy(gsl_vector* y, const gsl_vector* x, double alpha)
{
	cblas_daxpy(y->size, alpha, x->data, x->stride, y->data, y->stride);
}

void My_dscal(gsl_vector* x, const double alpha)
{
	cblas_dscal(x->size, alpha, x->data, x->stride);
}

void My_dswap(gsl_vector* x, gsl_vector* y)
{
	cblas_dswap(x->size, x->data, x->stride, y->data, y->stride);
}

void getAbsMax(const gsl_vector* x, double* maxValue, int* maxPosition)
{
	*maxPosition = cblas_idamax(x->size, x->data, x->stride);
	*maxValue = gsl_vector_get(x, *maxPosition);
}
