/*
 * BlasWrappers.h
 *
 *  Created on: 31 Oct 2012
 *      Author: spyros
 */

#ifndef BLASWRAPPERS_H_
#define BLASWRAPPERS_H_

#include <math.h>
#include <cblas.h>
#include <gsl/gsl_matrix.h>

double My_ddot(const gsl_vector * x, const gsl_vector* y);

double norm2(const gsl_vector* x);

void My_dgemm(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const double alpha, const gsl_matrix *A,
                          const gsl_matrix* B, const double beta, gsl_matrix* C);

void My_dsyrk(const enum CBLAS_UPLO Uplo,const enum CBLAS_TRANSPOSE Trans, double alpha, const gsl_matrix * A, double beta, gsl_matrix * C);

void My_dsyr(const enum CBLAS_UPLO Uplo,const enum CBLAS_TRANSPOSE Trans, double alpha, const gsl_matrix * A, double beta, gsl_matrix * C);

void My_dtrsm(const enum CBLAS_SIDE Side,const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                          double alpha, const gsl_matrix * A, gsl_matrix * B);

void My_dtrsv(const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                          const gsl_matrix * A, gsl_vector* B);

void My_drotg(double* a, double* b, double* c, double* s);

void My_drot(gsl_vector* x, gsl_vector* y, const double c, const double s);

void My_dgemv(const enum CBLAS_TRANSPOSE TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y);

void My_dger (double alpha, const gsl_vector * x, const gsl_vector * y, gsl_matrix * A);

double blas_sum(const gsl_vector* x);

double My_dasum(const gsl_vector* x);

void My_daxpy(gsl_vector* y, const gsl_vector* x, double alpha);

void My_dscal(gsl_vector* x, const double alpha);

void My_dswap(gsl_vector* x, gsl_vector* y);

void getAbsMax(const gsl_vector* x, double* maxValue, int* maxPosition);

#endif /* BLASWRAPPERS_HPP_ */
