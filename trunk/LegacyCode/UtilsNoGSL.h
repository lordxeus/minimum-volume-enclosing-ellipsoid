#ifndef UTILSNOGSL
#define UTILSNOGSL

#include <time.h>
#include <vector>

//Print utilities
void printMatrix(const double* X, const int* XRows, const int* XCols);

void printVector(const double* x, const int* xSize);

void printVector(const int* x,const int* xSize);

//End of print utilities

// ------------------------------------------------------------------
// --------------vector utilities------------------------------------
void getMax(const double* x, double* maxValue, int* maxPosition);

void getMin(const double* x, double* minValue, int* minPosition);

//This method calculates the sum of a vector
double sum(const double* x);

int sum(const int* x);

void initAct(int* act);

#endif