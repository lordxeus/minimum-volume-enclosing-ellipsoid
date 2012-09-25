/*
 *  UtilsNoGSL.cpp
 *  LinAlg
 *
 *  Created by Spyros Schismenos on 07/08/2012.
 *  Copyright 2012 JPMorgan. All rights reserved.
 *
 */

#include "UtilsNoGSL.h"

//Print utilities
void printMatrix(const double* X, const int* XRows, const int* XCols)
{
	for (int i=0; i<(*XRows); i++) {
		for (int j=0; j<(*XCols); j++) {
		printf("%d,%d,%d:   %e\n",i,j,i*(*XRows)+j,X[i*(*XRows)+j]);
		}
	}
}

void printVector(const double* x,const int* xSize)
{
	for (int i=0; i<(*xSize); i++) {
		printf("%d: %e\n",i,x[i]);
	}
}

void printVector(const int* x,const int* xSize)
{
	for (int i=0; i<(*xSize); i++) {
		printf("%d: %d\n",i,x[i]);
	}
}

//End of print utilities

// --------------vector utilities------------------------------------
void getMax(const double* x, double* maxValue, int* maxPosition)
{
}

void getMin(const double* x, double* minValue, int* minPosition)
{
}

//This method calculates the sum of a vector
double sum(const double* x)
{
	return 0.0;
}

int sum(const int* x)
{
	return 1;
}

void initAct(int* act)
{
}
