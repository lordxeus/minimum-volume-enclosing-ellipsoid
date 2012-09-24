#ifndef LEGACY_H
#define LEGACY_H
/*
 *  Legacy.h
 *  LinAlg
 *
 *  Created by Spyros Schismenos on 30/06/2012.
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

void initDiag(const gsl_vector* u, gsl_matrix* U);

void getSubMatrixFromColumnsAndRows(const gsl_matrix* input, const gsl_vector* CIndices, const gsl_vector* RIndices, gsl_matrix* output);
//This method calculates M = X*diag(u)*tr(X)
// LEGACY, probably should be removed
void getXUXT1(const gsl_matrix* X,const gsl_vector* u, gsl_matrix* M);

//This method calculates M = X*diag(u)*tr(X) efficiently hopefully
// only mm*(n+2*n^2)
void getXUXT2(const gsl_matrix* X,const gsl_vector* u, gsl_matrix* M);

//This method calculates M = X*diag(u)*tr(X) efficiently hopefully
// only mm*(n+2*n^2)
void getXUXT3(const gsl_matrix* X,const gsl_vector* u, gsl_matrix* M);

#endif