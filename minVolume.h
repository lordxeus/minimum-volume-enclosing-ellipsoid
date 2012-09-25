#ifndef MINVOLUME_H
#define MINVOLUME_H
/*
 *  minVolume.h
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
void getVar(gsl_vector* var, gsl_matrix* x);

double calcNormDiff(const gsl_matrix* M, const gsl_matrix* R, double factor);

//Calculates X * U * XT in an efficient way
// The symmetric matrix M's lower part is calculated.
// Cholesky should be able to cope with that...
void getXUXT(const gsl_matrix* X,const gsl_vector* u, gsl_matrix* M, int doWhole);

void updatevarOld(gsl_vector* var,double* tau,double* mult,gsl_matrix* Mxj,gsl_matrix* XX);

void updateVar(gsl_vector* var,double* tau,double* mult,const gsl_vector* Mxj,const gsl_matrix* XX);

void updateVar(gsl_vector* var,double* tau,double* mult,const gsl_matrix* Mxj,const gsl_matrix* XX);

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
