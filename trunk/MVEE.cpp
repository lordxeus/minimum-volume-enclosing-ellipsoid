/*
 *  MVEE.cpp
 *  LinAlg
 *
 *  Created by Spyros Schismenos on 27/06/2012.
 *
 */

#include "MVEE.h"

MVEE::MVEE(gsl_matrix* XX):
tol(0.0000001),
KKY(0),
maxit(100000),
doPrint(false)
{
	n = XX->size1;
	m = XX->size2;
	X = gsl_matrix_calloc(n,m);
	M = gsl_matrix_calloc(n,n);
	u = gsl_vector_calloc(m);
	gsl_matrix_memcpy(X,XX);
}

MVEE::MVEE(gsl_matrix* XX,const double& tol, const int& KKY, const int& maxit, const bool& doPrint)
:tol(tol),
KKY(KKY),
maxit(maxit),
doPrint(doPrint)
{
	n = XX->size1;
	m = XX->size2;
	X = gsl_matrix_alloc(n,m);
	M = gsl_matrix_alloc(n,n);
	u = gsl_vector_alloc(m);
	gsl_matrix_memcpy(X,XX);
}

MVEE::~MVEE()
{
	gsl_matrix_free(X);
	gsl_matrix_free(M);
	gsl_vector_free(u);
}

int MVEE::get_n() const {
	return n;
}

int MVEE::get_m() const {
	return m;
}

void MVEE::solve()
{
	isSolved = minVol(X,u,M,tol,KKY,maxit,doPrint);
}

void MVEE::printResults() const
{
	print(M);
}

bool MVEE::isOptimal(double tol)
{
	return OptimalityTest(X,M,tol,u);
}

void MVEE::calculateOptimalValue()
{
	gsl_matrix* LU=gsl_matrix_alloc(M->size1, M->size2);
	gsl_matrix_memcpy(LU, M);
	gsl_permutation* p = gsl_permutation_alloc(M->size1);
	int signum = 0;
	gsl_linalg_LU_decomp (LU, p, &signum);
	optimalValue = log(gsl_linalg_LU_det(LU,signum));
	gsl_matrix_free(LU);
	gsl_permutation_free(p);
}

void MVEE::get_u(gsl_vector* uu) const
{
	gsl_vector_memcpy(uu, u);
};

void MVEE::set_doPrint(bool _doPrint) {doPrint = _doPrint;}
