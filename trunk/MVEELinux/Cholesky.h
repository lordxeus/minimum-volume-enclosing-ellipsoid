/*
 * Cholesky.h
 *
 *  Created on: 31 Oct 2012
 *      Author: spyros
 */

#ifndef CHOLESKY_H_
#define CHOLESKY_H_

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

//Makes R the lower triangular cholesky factor of M
void getChol(const gsl_matrix* M, gsl_matrix* R);

void cholUpdate(gsl_matrix* R, gsl_vector* x);

void cholDowndate(gsl_matrix* R, gsl_vector* x, int* info);

void updateROld(gsl_matrix* R,double* factor,gsl_vector* xj,double* tau, int* down_err);

void updateR(gsl_matrix* R,double* factor,gsl_vector* p,gsl_vector* z, double* tau);

void dpotrf(const char *in,
                   int &n,
                   double *a,
                   int &lda,
                   int *info);

void dpotri(const char *in,
                   int &n,
                   double *a,
                   int &lda,
                   int *info);


#endif /* CHOLESKY_H_ */
