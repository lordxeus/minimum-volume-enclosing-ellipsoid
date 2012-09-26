#ifndef CHOLESKY
#define CHOLESKY

/*
 *  Cholesky.h
 *
 *  Created by Spyros Schismenos on 26/09/2012.
 *  Routines related to Cholesky factorization
 *  and Cholesky update
 */

//For LAPACK
#include <Accelerate/Accelerate.h>

//some gsl stuff
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

//Makes R the lower triangular cholesky factor of M
void getChol(const gsl_matrix* M, gsl_matrix* R);

void cholUpdate(gsl_matrix* R, gsl_vector* x);

void cholDowndate(gsl_matrix* R, gsl_vector* x, int* info);

void updateROld(gsl_matrix* R,double* factor,gsl_vector* xj,double* tau, int* down_err);

void updateR(gsl_matrix* R,double* factor,gsl_vector* p,gsl_vector* z, double* tau);

#endif