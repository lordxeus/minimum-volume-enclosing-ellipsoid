#ifndef UTILS
#define UTILS
/*
 *  Utils.h
 *  LinAlg
 *
 *  Created by Spyros Schismenos on 07/06/2012.
 *   
 *  This contains some basic utilities, mainly simple wrappers around commonly used functions
 */

#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <vector>

class MyCombination
{
public:
	gsl_vector_int* data;
	int n;
	int k;
	int h;
	bool isChild;
	MyCombination(int n, int k, int h, bool isChild=false);
	~MyCombination();
	int advance();
	MyCombination getHigherDim();
};

//Random number utilities
void getRandomVector(gsl_vector* M,gsl_rng * r);

void getRandomMatrix(gsl_matrix* M,gsl_rng * r);

//End of random number utilities

//Print utilities
void print(const gsl_matrix* X);

void print(const gsl_vector* x);

void print(const gsl_vector_int* x);

void printMATLAB(const gsl_matrix* X);

void printMATLAB(const gsl_vector* x);

//End of print utilities

// ------------------------------------------------------------------
// --------------vector utilities------------------------------------
void getMax(const gsl_vector* x, double* maxValue, int* maxPosition);

void getMin(const gsl_vector* x, double* minValue, int* minPosition);

//This method calculates the sum of a vector
double sum(const gsl_vector* x);

int sum(const gsl_vector_int* x);

// This method calculates sum(X.*X,1) in MATLAB
// slow code, but works with const x
void getVarOld(gsl_vector* var,const gsl_matrix* x);

// fast code, but it contaminates x
void getVarOldFast(gsl_vector* var, gsl_matrix* x);

void findThresh(const gsl_vector* x, gsl_vector* indices, double thresh);

void findThresh1or2(const gsl_vector* x1, const gsl_vector* x2, gsl_vector_int* indices, double thresh1, double thresh2);

//--------------end of vector utilities------------------------------
// ------------------------------------------------------------------

void initAct(gsl_vector_int* act);

double diffclock(clock_t clock1,clock_t clock2);

bool tolZero(const double& x1, const double& x2, const double tol);

void getNotI(const gsl_vector_int* x, gsl_vector_int* y);

double Myavg(const std::vector<double>& x);

double norm2(const gsl_matrix* X);

// min(x(u>eps)) in MATLAB
void minOnSubset(const gsl_vector* x, const gsl_vector* u, double* minValue, int* minPosition, double eps = 0.00000001);

//this does the equivalent of output=original(indices)
//here indices is of size n and has 0 or 1
//at least one 1 is needed
//output has to be already defined to be of size sum(indices) 
void getSubVector(gsl_vector* original,const gsl_vector_int* indices, gsl_vector* output);
void getSubVector(gsl_vector_int* original,const gsl_vector_int* indices, gsl_vector_int* output);

//this does the equivalent of output=original(:,indices)
//here indices is of size m and contains 0's or 1's
//at least one 1 is needed
//output has to be already defined to be of size [n,sum(indices)] 
void getSubMatrixFromColumns(const gsl_matrix* input, const gsl_vector_int* indices, gsl_matrix* output);

//This method emulates the find(u>eps) of MATLAB
//posIndices has the same size as x, and has 1 if element is positive, 0 otherwise.
void findPositive(gsl_vector* posIndices, const gsl_vector* x,double eps=0.000000000001);

#endif