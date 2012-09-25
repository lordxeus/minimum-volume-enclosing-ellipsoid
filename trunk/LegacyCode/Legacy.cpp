/*
 *  Legacy.cpp
 *  LinAlg
 *
 *  Created by Spyros Schismenos on 30/06/2012.
 *
 */

#include "Legacy.h"


void initDiag(const gsl_vector* u, gsl_matrix* U) {}
/*
{
	int i; size_t n = u->size;
	for (i=0; i<n; i++) {
		gsl_matrix_set(U,i,i,gsl_vector_get(u,i));	
	}
}
*/
void getSubMatrixFromColumnsAndRows(const gsl_matrix* input, const gsl_vector* CIndices, const gsl_vector* RIndices, gsl_matrix* output)
{}
/*
{
	size_t n = input->size1;
	size_t m = input->size2;
	int i=0; int j=0; int k_i=0; int k_j=0;
	for (i=0; i<n; i++) {
		k_j=0;
		if (gsl_vector_get(RIndices,i)>0.5) {
			for (j=0; j<m; j++) {
				if (gsl_vector_get(CIndices,j)>0.5) {
					gsl_matrix_set(output,k_i,k_j,gsl_matrix_get(input,i,j));
					k_j++;
				}
			}
			k_i++;
		}
	}
}
*/

//This method calculates M = X*diag(u)*tr(X)
//TODO: make it more efficient
void getXUXT1(const gsl_matrix* X,const gsl_vector* u, gsl_matrix* M) {}
/*
{
	gsl_matrix* U = gsl_matrix_calloc(u->size,u->size); initDiag(u,U);	
	
	gsl_vector* fu = gsl_vector_calloc(u->size); findThresh(u, fu, VERYSMALL); int lfu = sum(fu);
	
	gsl_matrix* X_fu = gsl_matrix_calloc(X->size1,lfu); getSubMatrixFromColumns(X, fu, X_fu);
	
	gsl_matrix* U_fu = gsl_matrix_calloc(lfu,lfu); getSubMatrixFromColumnsAndRows(U,fu,fu,U_fu);
	gsl_matrix* XU = gsl_matrix_calloc(X->size1,lfu); //TODO: check if there is a better multiplication for diagonal matrices
	
	My_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X_fu, U_fu , 0.0,XU);
	My_dgemm(CblasNoTrans, CblasTrans, 1.0, XU, X_fu , 0.0,M); 
	
	//now free some memory
	gsl_matrix_free(XU);
	gsl_matrix_free(X_fu);
	gsl_matrix_free(U);
	gsl_matrix_free(U_fu);
	gsl_vector_free(fu);
	
}
*/
//This method calculates M = X*diag(u)*tr(X) efficiently hopefully
// only mm*(n+2*n^2)
void getXUXT2(const gsl_matrix* X,const gsl_vector* u, gsl_matrix* M) {}
/*
{	
	size_t n = X->size1;
	size_t m = X->size2;
	size_t i,j;
	gsl_matrix* Xrow = gsl_matrix_calloc(n,1);
	for (i=0; i<m; i++) {
		if (gsl_vector_get(u,i)>VERYSMALL) {
			for (j=0; j<n; j++) {
				gsl_matrix_set(Xrow,j,0,gsl_matrix_get(X,j,i));
			}
			My_dgemm(CblasNoTrans,CblasTrans, gsl_vector_get(u,i), Xrow, Xrow, 1.0, M);
		}
	}
	gsl_matrix_free(Xrow);
	
}
*/

//This method calculates M = X*diag(u)*tr(X) efficiently hopefully
// only mm*(n+2*n^2)
void getXUXT3(const gsl_matrix* X,const gsl_vector* u, gsl_matrix* M) {}
/*
{	
	size_t n = X->size1;
	size_t m = X->size2;
	size_t i,j;
	gsl_matrix* Xrow = gsl_matrix_calloc(n,1);
	for (i=0; i<m; i++) {
		double ui = gsl_vector_get(u,i);
		if (ui>VERYSMALL) {
			for (j=0; j<n; j++) {
				gsl_matrix_set(Xrow,j,0,gsl_matrix_get(X,j,i));
			}
			My_dsyrk(CblasLower, CblasNoTrans, ui,Xrow, 1.0, M);
		}
	}
	gsl_matrix_free(Xrow);
	for (i=0; i<n; i++) {
		for (j=0; j<i; j++) {
			gsl_matrix_set(M,j,i,gsl_matrix_get(M,i,j));
		}
	}	
}
*/