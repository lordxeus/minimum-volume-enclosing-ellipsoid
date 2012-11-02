/*
 * Cholesky.cpp
 *
 *  Created on: 31 Oct 2012
 *      Author: spyros
 */

#include "BlasWrappers.h"
#include "Cholesky.h"

extern "C" void dpotrf_(const char *in,
                   int &n,
                   double *a,
                   int &lda,
                   int *info);

extern "C" void dpotri_(const char *in,
                   int &n,
                   double *a,
                   int &lda,
                   int *info);


void dpotrf(const char *in,
                   int &n,
                   double *a,
                   int &lda,
                   int *info)

{
	dpotrf_(in,n,a,lda,info);
}

void dpotri(const char *in,
                   int &n,
                   double *a,
                   int &lda,
                   int *info)

{
	dpotri_(in,n,a,lda,info);
}
void cholUpdate(gsl_matrix* R, gsl_vector* x)
{
        int n = R->size1;
        int i = 0;
        double c; double s;
        for (i=0;i<n;i++) {
                double* a = gsl_matrix_ptr(R,i,i);
                double* b = gsl_vector_ptr(x, i);
                My_drotg(a,b,&c,&s);
                if ((*a)<0.0) {
                        *a = - (*a);
                        c = - c;
                        s = - s;
                }
                if (i<n-1) {
                        gsl_vector_view Ri = gsl_matrix_column(R, i);
                        gsl_vector_view Rii = gsl_vector_subvector(&Ri.vector, i+1, n-i-1);
                        gsl_vector_view xi = gsl_vector_subvector(x, i+1, n-i-1);
                        My_drot(&Rii.vector,&xi.vector,c,s);
                }
        }
}

void cholDowndate(gsl_matrix* R, gsl_vector* x, int* info)
{
        int n = R->size1;
        int i;
        gsl_vector* c = gsl_vector_calloc(n);
        gsl_vector* s = gsl_vector_calloc(n);

        My_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, R, x);
        double aux = norm2(x);
        if (aux>1.0) {
                *info = -1;
                return;
        }
        if (aux>0.5) {
                aux = sin(acos(aux));
        }
        else {
                aux = cos(asin(aux));
        }
        for (i=n-1; i>=0; i--) {
                double c1,s1;
                My_drotg(&aux,gsl_vector_ptr(x,i),&c1,&s1);
                if (aux<0.0) {
                        aux = -aux;
                        gsl_vector_set(c,i, - c1);
                        gsl_vector_set(s,i, - s1);
                }
                else {
                        gsl_vector_set(c,i,c1);
                        gsl_vector_set(s,i,s1);
                }

        }
        for (int j=0; j<n; j++) {
                aux = 0.0;
                for (int ii=0; ii<=j; ii++) {
                        i = j-ii;
                        double temp = gsl_vector_get(c,i)*aux+gsl_vector_get(s,i)*gsl_matrix_get(R,j,i);
                        double temp2 = gsl_matrix_get(R,j,i);
                        gsl_matrix_set(R,j,i,gsl_vector_get(c,i)*temp2-gsl_vector_get(s,i)*aux);
                        aux = temp;
                }
        }
        gsl_vector_free(c);
        gsl_vector_free(s);
}

void getChol(const gsl_matrix* M, gsl_matrix* R)
{
        int n = M->size1;
        gsl_matrix_memcpy(R, M);
        int info=0;
        char lower = 'U';
        int lda = R->tda;
        dpotrf(&lower, n, R->data, lda, &info);
}

void updateROld(gsl_matrix* R,double* factor,gsl_vector* xj,double* tau, int* down_err)
{
        *factor = (*factor) / (1.0 - (*tau));
        int p=0;
        My_dscal(xj,sqrt(GSL_MAX(*tau,- *tau)*(*factor)));
        if (*tau>0.0) {
                cholUpdate(R, xj);
        }
        else {
                cholDowndate(R, xj, &p);
        }
        if (p>0) {
                *down_err = 1;
        }
        else {
                *down_err = 0;
        }
}

void updateR(gsl_matrix* R,double* factor,gsl_vector* p,gsl_vector* z, double* tau)
{
        *factor = (*factor) / (1.0 - (*tau));

        My_dscal(p,sqrt(GSL_MAX(*tau,- *tau)*(*factor)));

        My_dscal(z,sqrt(GSL_MAX(*tau,- *tau)*(*factor)));


        gsl_vector* w = gsl_vector_alloc(z->size); gsl_vector_memcpy(w, z);
        gsl_vector* s = gsl_vector_calloc(w->size+1);
        int n = w->size;

        gsl_vector_set(s, n-1, gsl_vector_get(p, n-1)*gsl_vector_get(p, n-1));
        for (int i=n-2; i>=0; i--) {
                gsl_vector_set(s, i, gsl_vector_get(s, i+1)+gsl_vector_get(p, i)*gsl_vector_get(p, i));
        }


        double a = 1.0;
        if (*tau < 0.0) {
                a = -1.0;
        }
        double sigma = a/(1.0+sqrt(1.0+a*gsl_vector_get(s, 0)));
        double q;
        double theta;
        double sigma1;
        double beta;
        double rho;

        gsl_vector* d2 = gsl_vector_alloc(n);
        gsl_vector_view d22 = gsl_matrix_diagonal(R);
        gsl_vector_memcpy(d2, &d22.vector);

        for (int j=0; j<n; j++) {

                q = gsl_pow_2(gsl_vector_get(p, j));

                theta = 1.0 + sigma * q;

                gsl_vector_set(s, j+1, gsl_vector_get(s, j)-q);

                rho =  sqrt(theta*theta+sigma*sigma*q*gsl_vector_get(s, j+1));

                beta = a * gsl_vector_get(p, j) * gsl_matrix_get(R, j, j);

                gsl_matrix_set(R, j, j, rho * gsl_matrix_get(R, j, j));

                beta = beta/gsl_matrix_get(R, j, j)/gsl_matrix_get(R, j, j);

                a = a / rho/rho;
                sigma1 = sigma* (1.0 + rho)/(rho*(theta + rho));
                sigma = sigma1;
                //for (int r = j+1; r<n; r++) {
                //      gsl_vector_set(w, r, gsl_vector_get(w, r)-gsl_vector_get(p, j)*gsl_matrix_get(R, r, j));
                //      gsl_matrix_set(R, r, j, gsl_matrix_get(R, r, j)/gsl_vector_get(d2, j)+beta*gsl_vector_get(w, r));
                //      gsl_matrix_set(R, r, j, gsl_matrix_get(R, r, j)*gsl_matrix_get(R, j, j));
                //}
                if (j<n-1) {
                        gsl_vector_view wr = gsl_vector_subvector(w, j+1, n-j-1);
                        gsl_vector_view Rr = gsl_matrix_subcolumn(R, j, j+1, n-j-1);
                        My_daxpy(&wr.vector, &Rr.vector, -gsl_vector_get(p, j));
                        My_dscal(&Rr.vector, 1.0/gsl_vector_get(d2, j));
                        My_daxpy(&Rr.vector, &wr.vector, beta);
                        My_dscal(&Rr.vector, gsl_matrix_get(R, j, j));
                }

        }

        //clean up
        gsl_vector_free(w);
        gsl_vector_free(s);
        gsl_vector_free(d2);

}


