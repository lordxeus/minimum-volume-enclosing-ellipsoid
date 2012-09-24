#ifndef MINVOLUME_H
#define MINVOLUME_H
/*
 *  minVolume.h
 *  LinAlg
 *
 *  Created by Spyros Schismenos on 24/06/2012.
 *
 */

#include <iostream>
#include <stdio.h>
#include <vector>
#include "Utils.h"
#include "BlasWrappers.h"

//For LAPACK
#include <Accelerate/Accelerate.h>

//some gsl stuff
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

const double SMALLE = 0.00000001;
const double VERYSMALL =  0.000000000001;

//even faster code, still contaminates x
void getVarFastFast(gsl_vector* var, gsl_matrix* x);

void cholUpdate(gsl_matrix* R, gsl_vector* x);

void cholDowndate(gsl_matrix* R, gsl_vector* x, int* info);

//This method emulates the find(u>0) of MATLAB
//posIndices has the same size as x, and has 1 if element is positive, 0 otherwise.
void findPositive(gsl_vector* posIndices, const gsl_vector* x);

// min(x(u>0)) in MATLAB
void minOnSubset(const gsl_vector* x, const gsl_vector* u, double* minValue, int* minPosition);

//this does the equivalent of output=original(indices)
//here indices is of size n and has 0 or 1
//at least one 1 is needed
//output has to be already defined to be of size sum(indices) 
void getSubVector(gsl_vector* original,const gsl_vector* indices, gsl_vector* output);

//this does the equivalent of output=original(:,indices)
//here indices is of size m and contains 0's or 1's
//at least one 1 is needed
//output has to be already defined to be of size [n,sum(indices)] 
void getSubMatrixFromColumns(const gsl_matrix* input, const gsl_vector_int* indices, gsl_matrix* output);

//norm(X) in Matlab
double norm2(const gsl_matrix* X);

double calcNormDiff(const gsl_matrix* M, const gsl_matrix* R, double factor);

//Calculates X * U * XT in an efficient way
// The symmetric matrix M's lower part is calculated.
// Cholesky should be able to cope with that...
void getXUXT(const gsl_matrix* X,const gsl_vector* u, gsl_matrix* M, int doWhole);

//Makes R the lower triangular cholesky factor of M
void getChol(const gsl_matrix* M, gsl_matrix* R);

void updateR(gsl_matrix* R,double* factor,gsl_vector* xj,double* tau, int* down_err);

void updateRMike(gsl_matrix* R,double* factor,gsl_vector* p,gsl_vector* z, double* tau, int* down_err);

void updatevar(gsl_vector* var,double* tau,double* mult,gsl_matrix* Mxj,gsl_matrix* XX);

void updateVarFast(gsl_vector* var,double* tau,double* mult,const gsl_vector* Mxj,const gsl_matrix* XX);

void updateVarFast(gsl_vector* var,double* tau,double* mult,const gsl_matrix* Mxj,const gsl_matrix* XX);

//the special case n=m
int minVol(const gsl_matrix* X,
		   gsl_vector* uu,
		   gsl_matrix* M_out,
		   int n,
		   const double tol=0.0000001,
		   const int KKY=0,
		   const int maxit=100000,
		   const bool doPrint=false,
		   const bool isInitialized=false);

int minVol(const gsl_matrix* X,
		   gsl_vector* uu,
		   gsl_matrix* M_out,
		   const double tol=0.0000001,
		   const int KKY=0,
		   const int maxit=100000,
		   const bool doPrint=false,
		   const bool isInitialized=false);

bool OptimalityTest(const gsl_matrix* X, const gsl_matrix* M,double tol, const gsl_vector* u);

void initwt(const gsl_matrix* X, gsl_vector* u);

void initwt6(const gsl_matrix* X, gsl_vector* u);

int calculateBoundAdd(const gsl_matrix* M, int k, const gsl_matrix* X, const double& optimal, double& bound);

void calculateBoundRemove(const gsl_matrix* M, int k, const gsl_vector* u, const double& optimal, double& bound);

void calculateBoundAddRemove(const gsl_matrix* M, int k , int l, const gsl_matrix* X, const gsl_vector* u, const double& optimal, double& bound);

//version that is used for when k and l are indices in X, but u is of size h<m.
void calculateBoundAddRemove(const gsl_matrix* M, int k , int l,int i, const gsl_matrix* X, const gsl_vector* u, const double& optimal, double& bound);

void getRandomMatrixMVEE(gsl_matrix* X, gsl_rng* r, double& optimalValue, int h);

#endif
