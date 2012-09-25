/*
 *  Utils.cpp
 *  LinAlg
 *
 *  Created by Spyros Schismenos on 07/06/2012.
 *
 */

#include "Utils.h"

// we always assume that k<=h, otherwise it doesn't make sense!
MyCombination::MyCombination(int n, int k, int h, bool isChild):n(n),k(k), h(h), isChild(isChild)
{
	data = gsl_vector_int_alloc(k);
	for (int i=0; i<k; i++) {
		gsl_vector_int_set(data, i, i);
	}
}

MyCombination::~MyCombination()
{
	gsl_vector_int_free(data);
}

int MyCombination::advance()
{
	if (!isChild) {
		for (int i=k-1; i>=0; i--) {
			if (i==k-1) {
				if (gsl_vector_int_get(data, i)<n-1+k-h) {
					gsl_vector_int_set(data, i, gsl_vector_int_get(data, i)+1);
					return 1;
				}
			}
			else {
				if (i>0) {
					if (gsl_vector_int_get(data,i)<gsl_vector_int_get(data, i+1)-1) {
						gsl_vector_int_set(data, i, gsl_vector_int_get(data, i)+1);
						for (int j=i+1; j<k; j++) {
							gsl_vector_int_set(data, j, gsl_vector_int_get(data, i)+j-i);
						}
						return 1; 
					}
				}
				if (i==0) {
					if (gsl_vector_int_get(data, i)==n-h) {
						return 0;
					}
					else {
						gsl_vector_int_set(data, 0, gsl_vector_int_get(data, 0)+1);
						for (int j=1; j<k; j++) {
							gsl_vector_int_set(data, j, gsl_vector_int_get(data, 0)+j);
						}
						return 1;
					}
				}
			}
		}
	}
	else{ //if we are a child, only move the last point!
		if (gsl_vector_int_get(data, k-1)<n-1+k-h) {
			gsl_vector_int_set(data, k-1, gsl_vector_int_get(data, k-1)+1);
			return 1;
		}
		else {
			return 0;
		}

	}
	//Should never end up here
	return 2;
}

MyCombination MyCombination::getHigherDim()
{
	//Verify that k<=n-1 (Should easily never happen)
	//Verify that data(k)<n-1 (Should make sure it doesn't happen)
	MyCombination newCombination(n, k+1,h,true);
	for (int i=0; i<k; i++) {
		gsl_vector_int_set(newCombination.data, i, gsl_vector_int_get(data, i));
	}
	gsl_vector_int_set(newCombination.data, k, gsl_vector_int_get(data, k-1)+1);
	return newCombination;
}

//Random number utilities
void getRandomMatrix(gsl_matrix* M,gsl_rng * r)
{
	int n = M->size1;
	int m = M->size2;
	for (int i = 0; i < n; i++) 
	{
		for (int j=0; j<m;j++)
		{
			//double u = gsl_rng_uniform (r);
			double u = gsl_ran_gaussian(r,1.0); 
			gsl_matrix_set(M,i,j,u);
		}
	}
}

void getRandomVector(gsl_vector* M,gsl_rng * r)
{
	int n = M->size;
	for (int i = 0; i < n; i++) 
	{
		double u = gsl_ran_gaussian(r,1.0); 
		gsl_vector_set(M,i,u);
	}
}

//End of random number utilities

// min and max
void getMax(const gsl_vector* x, double* maxValue, int* maxPosition)
{
	*maxPosition = gsl_vector_max_index(x);
	*maxValue = gsl_vector_get(x,*maxPosition);
}

void getMin(const gsl_vector* x, double* minValue, int* minPosition)
{
	*minPosition = gsl_vector_min_index(x);
	*minValue = gsl_vector_get(x,*minPosition);
}

//This method calculates the sum of a vector
double sum(const gsl_vector* x)
{
	size_t n = x->size;
	double res = 0.0;
	for (int i=0; i<n; i++) {
		res += gsl_vector_get(x,i);
	}
	return res;
}

//This method calculates the sum of a vector
int sum(const gsl_vector_int* x)
{
	size_t n = x->size;
	int res = 0;
	for (int i=0; i<n; i++) {
		res += gsl_vector_int_get(x,i);
	}
	return res;
}

void getVarOld(gsl_vector* var,const gsl_matrix* x)
{
	size_t n = x->size1;
	size_t m = x->size2;
	for (int i=0; i<m; i++) {
		for (int j=0; j<n; j++) {
			double xij = gsl_matrix_get(x,j,i);
			gsl_vector_set(var,i, gsl_vector_get(var,i)+xij*xij);
		}
	}
}

void getVarOldFast(gsl_vector* var, gsl_matrix* x)
{
	gsl_matrix_mul_elements(x,x); //x is now contaminated!
	int n = x->size1;
	gsl_vector_view x0 =  gsl_matrix_row(x,0);
	gsl_vector_memcpy(var, &x0.vector);
	for (int i=1 ;i<n;i++)
	{
		gsl_vector_view xi =  gsl_matrix_row(x,i);
		gsl_vector_add(var,&xi.vector);
	}
	
}

void findThresh(const gsl_vector* x, gsl_vector* indices, double thresh)
{
	int i; size_t n = x->size;
	for (i=0; i<n; i++) {
		if (gsl_vector_get(x,i)>thresh){
			gsl_vector_set(indices,i,1.0);
		}
	}
}

void findThresh1or2(const gsl_vector* x1, const gsl_vector* x2, gsl_vector_int* indices, double thresh1, double thresh2)
{
	size_t n = x1->size;
	int i=0;
	for (i=0; i<n; i++) {
		if ((gsl_vector_get(x1,i)>thresh1)||(gsl_vector_get(x2,i)>thresh2)){
			gsl_vector_int_set(indices,i,1);
		}
	}
}

//end of min and max


double diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
	return diffms;
} 

void print(const gsl_matrix* X)
{
    size_t n = X->size1;
	size_t m = X->size2;
	int i; int j;
	for (i=0; i<n; i++) {
		for (j=0; j<m; j++) {
			printf("(%d , %d) element :%.5e\n",i,j,gsl_matrix_get(X,i,j));
		}
	}
	printf("------------------\n");
}

void print(const gsl_vector* x)
{
	size_t n = x->size;
	int i;
	for (i=0; i<n; i++) {
		printf("(%d) element :%.5e\n",i,gsl_vector_get(x,i));
	}
	printf("------------------\n");
}

void print(const gsl_vector_int* x)
{
	size_t n = x->size;
	int i;
	for (i=0; i<n; i++) {
		printf("%d ",gsl_vector_int_get(x,i));
	}
}

void printMATLAB(const gsl_matrix* X)
{
	size_t n = X->size1;
	size_t m = X->size2;
	printf("matrix=[");
	for (int i=0; i<n; i++) {
		for (int j=0; j<m; j++) {
			printf(" %.9e ",gsl_matrix_get(X,i,j));
			if (j<m-1) {
				printf(",");
			}
		} 
		if (i<n-1) {
			printf(";");
		}
		else
			printf("]");
	}
	printf("------------------\n");
}
	

void printMATLAB(const gsl_vector* x)
{
	size_t n = x->size;
	printf("vector=[");
	for (int i=0; i<n; i++) {
		printf(" %.5e ",gsl_vector_get(x,i));
		if (i<n-1) {
			printf(";");
		}
		else
			printf("]");
		
	}
	printf("------------------\n");
}

void initAct(gsl_vector_int* act)
{
	int n = act->size;
	for (int i=0; i<n; i++) {
		gsl_vector_int_set(act,i, i);
	}
}

bool tolZero(const double& x1, const double& x2, const double tol)
{
	return	((x1-x2)<tol&&(x1-x2)>-tol);
}

void getNotI(const gsl_vector_int* x, gsl_vector_int* y)
{
	//print(x); printf("\n");
	int h =x->size;
	int m = h + y->size;
	int count = 0;
	int count2 = 0;
	for (int i=0; i<m; i++) {
		if (count<h) {
			if (gsl_vector_int_get(x, count)==i) {
				count++;
			}
			else {
				gsl_vector_int_set(y, count2, i);
				//print(y);
				//printf("\n");
				count2++;
			}
		}
		else {
			gsl_vector_int_set(y, count2, i);
			//print(y);
			//printf("\n");
			count2++;
		}

	}
}


double Myavg(const std::vector<double>& x)
{
	double res = 0.0;
	for (int i=0; i<x.size(); i++) {
		res = res + x[i];
	}
	return res/(double)x.size();
}
