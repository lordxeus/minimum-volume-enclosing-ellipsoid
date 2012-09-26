/*
 *  minVolume.cpp
 *
 *  Created by Spyros Schismenos on 24/06/2012.
 *
 */

#include "minVolume.h" 
#include "Cholesky.h"

void getVar(gsl_vector* var, gsl_matrix* x)
{
	gsl_matrix_mul_elements(x,x); //x is now contaminated!
	int n = x->size1;
	gsl_vector* ones = gsl_vector_alloc(n); gsl_vector_set_all(ones, 1.0);
	My_dgemv(CblasTrans, 1.0, x, ones, 0.0, var);
	gsl_vector_free(ones);
}

double calcNormDiff(const gsl_matrix* M, const gsl_matrix* R, double factor)
{
	gsl_matrix* Mtemp = gsl_matrix_alloc(M->size1,M->size2);
	gsl_matrix_memcpy(Mtemp, M);
	My_dgemm(CblasNoTrans, CblasTrans, -1.0, R, R, factor, Mtemp);
	double normDiff = norm2(Mtemp)/GSL_MAX(factor,-factor)/norm2(M);
	gsl_matrix_free(Mtemp);
	return normDiff;
}

//Calculates X * U * XT in an efficient way
// The symmetric matrix M's lower part is calculated.
// Cholesky should be able to cope with that...
void getXUXT(const gsl_matrix* X,const gsl_vector* u, gsl_matrix* M, int doWhole)
{	
	size_t n = X->size1;
	size_t m = X->size2;
	size_t i,j;
	for (i=0; i<m; i++) {
		double ui = gsl_vector_get(u,i);
		if (ui>VERYSMALL) {
			gsl_matrix_const_view Xrow = gsl_matrix_const_submatrix(X,0,i,n,1); 
			for (j=0; j<n; j++) {
			}
			My_dsyr(CblasLower, CblasNoTrans, ui,&Xrow.matrix, 1.0, M);
		}
	}
	if (doWhole==1) {
		for (i=0; i<n; i++) {
			for (j=0; j<i; j++) {
				gsl_matrix_set(M,j,i,gsl_matrix_get(M,i,j));
			}
		}
	}
}

void updateVar(gsl_vector* var,double* tau,double* mult,const gsl_vector* Mxj,const gsl_matrix* XX)
{	
	int m = XX->size2;
	gsl_vector* temp = gsl_vector_alloc(m);
	double tauModified = 1.0/(1.0-*tau);
	My_dgemv(CblasTrans,1.0,XX,Mxj,0.0,temp);
	gsl_vector_mul(temp, temp);
	//gsl_vector_scale(temp,*mult);
	//gsl_vector_sub(var,temp);
	My_daxpy(var, temp, -(*mult));
	//gsl_vector_scale(var, tauModified);
	My_dscal(var, tauModified);
	gsl_vector_free(temp);
}

int minVol(const gsl_matrix* X,
		   gsl_vector* u_out,
		   gsl_matrix* M_out,
		   int n,
		   const double tol,
		   const int KKY,
		   const int maxit,
		   const bool doPrint,
		   const bool isInitialized)

{
	gsl_vector* u1 = gsl_vector_alloc(n);
	gsl_vector_set_all(u1, 1.0/n);
	gsl_matrix* M1 = gsl_matrix_calloc(n,n); getXUXT(X,u1,M1,1); 
	gsl_matrix_memcpy(M_out, M1);
	gsl_matrix_free(M1);
	gsl_vector_free(u1);
	return GSL_SUCCESS;
}

int minVol(const gsl_matrix* X,
		   gsl_vector* u_out,
		   gsl_matrix* M_out,
		   const double tol,
		   const int KKY,
		   const int maxit,
		   const bool doPrint,
		   const bool isInitialized)
{
	int numRecomputes = 0;
	//Finds the minimum-volume ellipsoid containing the columns of the matrix X using the
	// Fedorov-Wynn-Franke-Wolfe method, with Wolfe-Atwood away steps if KKY=0
	
	// The algorithm also uses the method of Harman and Pronzato to eliminate points that are found
	//to be inessential.
	//The algorithm returns an ellipsoid providing a (1+tol)*n-rounding of the convex hull of the columns of X in the 
	//n-space. 
	//Set tol = esp/n to get an (1+eps) approximation of the minimum volume ellipsoid.
	
	//On output
	//u determines the optimal weights on the m columns of X
	//M = XUXT determines the shape of the ellipsoid xT*inv(M)*x <=n

	int n = X->size1; int m = X->size2;

	// Handle the special case m = n
	if (m==n) {
		return minVol(X,u_out,M_out,n,tol,KKY,maxit,doPrint,isInitialized);
	}
	
	//Play safe and trivially allocate some variables
	gsl_vector* uold = gsl_vector_alloc(1);
	gsl_vector_int* smallAct = gsl_vector_int_alloc(1);
	gsl_vector* smallVar = gsl_vector_alloc(1);
	gsl_vector* u_e = gsl_vector_alloc(1);
	gsl_vector_int* acttemp = gsl_vector_int_alloc(1);
	gsl_vector* vartemp = gsl_vector_alloc(1);
	gsl_vector* utemp = gsl_vector_alloc(1);
	gsl_vector* uoldtemp = gsl_vector_alloc(1);
	
	//Set default inputs
 	gsl_vector* u = gsl_vector_alloc(m);
	gsl_vector* xj = gsl_vector_alloc(n);
	gsl_vector* Rxj = gsl_vector_alloc(n);
	gsl_vector* z = gsl_vector_alloc(n); //useful in the Cholesky update
	if (!isInitialized) {
		if (KKY>0) {
			gsl_vector_set_all(u, 1.0/m);
		}
		else {
			initwt6(X, u);
		}
	}
	else {
		gsl_vector_memcpy(u, u_out);
	}

	//print(u);
	
	//Create the Cholesky factor R of the matrix M = X*U*XT
	//inv(M) = factor * inv(R)*inv(R)^T
	gsl_matrix* M = gsl_matrix_calloc(n,n); getXUXT(X,u,M,0); //print(M);
	gsl_matrix* R = gsl_matrix_calloc(n,n); getChol(M,R); //print(R);//now R holds the upper triangular part
	double factor = 1.0;
	
	//act lists the mm non-eliminated columns of X, XX is the corresponding submatrix. 
	int mm=m; int oldmm = m;
	gsl_vector_int* act = gsl_vector_int_alloc(m); initAct(act);
	gsl_matrix* XX = gsl_matrix_alloc(n,m); gsl_matrix_memcpy (XX, X);
	
	
	// Initialize var: var(j) = x_j^T M^{-1} x_j,
	// and the number of drop, decrease, and add/increase iterations.
	// Statistics are printed every n or 100 iterations, and refactorization is
	// considered every n or 50000.
	gsl_matrix* RX = gsl_matrix_alloc(XX->size1,XX->size2); gsl_matrix_memcpy(RX, XX);
	My_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, R, RX); //print(RX);
	gsl_vector* var = gsl_vector_alloc(m); getVar(var, RX); //print(var);
	int iter=0;
	int n100 = GSL_MAX(n, 100);
	int n50000 = GSL_MAX(n,50000);
	gsl_vector* upos = gsl_vector_calloc(m); findPositive(upos,u); //print(u); print(upos);
	
	// maxvar is the maximum, and minvar the minimum over indices with u_i positive, variance.
	double maxvar; int maxj; getMax(var, &maxvar, &maxj);
	double minvar; int minj; minOnSubset(var,upos,&minvar,&minj);  double mvup = minvar; 	
	
	
	// Use the Harman-Pronzato test to see if columns of X can be eliminated.
	double ept = (maxvar-(double)n)/(double)n;
	double tresh=n*(1.0+0.5*ept-0.5*sqrt(ept*(4.0+ept-4.0/n)));
 	gsl_vector_int* e=gsl_vector_int_calloc(m); findThresh1or2(var,u,e, tresh,SMALLE);
	int length_e = sum(e); gsl_vector_int_free(smallAct); smallAct = gsl_vector_int_calloc(length_e); getSubVector(act,e,smallAct);
	gsl_vector_int_free(act); act = gsl_vector_int_alloc(smallAct->size); gsl_vector_int_memcpy(act, smallAct);
	gsl_matrix* newXX = gsl_matrix_alloc(n,act->size); getSubMatrixFromColumns(XX,act,newXX);
	gsl_matrix_free(XX); XX = gsl_matrix_alloc(n,act->size); gsl_matrix_memcpy(XX,newXX);
	mm = act->size;
	
	// If only n columns remain, recompute u, M, and R.
	if (mm==n) {
		gsl_vector_free(u); u=gsl_vector_alloc(n); gsl_vector_set_all(u,1.0/n);
		
		gsl_matrix_set_zero(M); getXUXT(XX,u,M,0);
		gsl_matrix_free(R); R = gsl_matrix_calloc(M->size1,M->size2); getChol(M,R); 
		factor = 1.0;
		gsl_matrix_free(RX); RX = gsl_matrix_alloc(XX->size1,XX->size2); gsl_matrix_memcpy(RX, XX);
		My_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, R, RX); //XX is RX in the MATLAB code
		gsl_vector_free(var); var = gsl_vector_alloc(n); getVar(var,RX);
	}
	else {
		gsl_vector_free(smallVar); smallVar = gsl_vector_calloc(mm); getSubVector(var, e, smallVar);
		gsl_vector_free(var); var = gsl_vector_alloc(mm); gsl_vector_memcpy (var, smallVar);
		
		gsl_vector_free(u_e); u_e = gsl_vector_calloc(mm); getSubVector(u, e, u_e); 
		gsl_vector_free(u); u = gsl_vector_alloc(mm); gsl_vector_memcpy (u, u_e);
		My_dscal(u, 1.0/blas_sum(u));
		
		gsl_matrix_free(RX); RX = gsl_matrix_alloc(XX->size1,XX->size2); gsl_matrix_memcpy(RX, XX);
	}	 
	
	oldmm = mm;
	getMax(var,&maxvar,&maxj); 
	minOnSubset(var,u,&minvar,&minj); mvup = minvar;
	
	//Start of loop. Check for termination.
	while (((maxvar > (1.0+tol)*n) || (mvup < (1.0-tol)*n)) && (iter < maxit)) 
	{
		//Find "furthest" and "closest" points using updated var
		iter = iter+1;
		int j=0; double mvar; int flag_decrease;
		if (maxvar + mvup > 2*n) {
			j = maxj;
			mvar = maxvar ;
			flag_decrease = 0;
		}
		else {
			j = minj;
			mvar = mvup;
			flag_decrease = 1;
		}
		//printf("j: %d\n",j);
		//printf("flag_decrease: %d\n",flag_decrease);
		//Compute Mxj = inv(M)*xj and recompute varj
		int flag_recompute = 0;
		gsl_vector_view column_j = gsl_matrix_column(XX,j);
		gsl_vector_memcpy(xj, &column_j.vector);
		gsl_vector_memcpy(Rxj, xj);

		//print(R);
		//print(Rxj);
		My_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, R, Rxj);
		gsl_vector_memcpy(z, Rxj);
		double mvarn = My_ddot(Rxj, Rxj);
		mvarn *= factor;
		My_dtrsv(CblasLower, CblasTrans, CblasNonUnit, R, Rxj);
		My_dscal(Rxj, factor);
		
		double mvarerror = mvarn - mvar; mvarerror = GSL_MAX(mvarerror,-mvarerror)/GSL_MAX(mvar,1.0);
		if (mvarerror>SMALLE) {
			//printf("%e %e\n",mvarn,mvar);
			flag_recompute = 1;
		}
		mvar = mvarn;
		
        // Compute stepsize tau (may be negative), epsilon, and improvement in logdet.
		double tau = 0.0;
		if (mvar<=1.0) {
			double uj = gsl_vector_get(u,j);
			tau = -uj/(1.0-uj);
		}
		else {
			tau = (mvar/n - 1.0)/(mvar - 1.0);
		}
		if (tau>(1.0 - 0.000001)) {
			tau = 1.0 - 0.000001;
		}
		
		//Update u and make sure it stays nonnegative.
		if (flag_decrease) {
			double uj = gsl_vector_get(u,j);
			if (tau< ( -uj/(1.0-uj))) {
				tau = -uj/(1.0-uj);
			}
		} 
		
		if ((iter == (int)(iter/n50000)*n50000)||(iter == (int)(iter/n100)*n100)) {
			gsl_vector_free(uold); uold = gsl_vector_alloc(u->size); gsl_vector_memcpy(uold,u);
		}
		My_dscal(u,1.0-tau);
		double * ujTemp = gsl_vector_ptr(u, j);
		*ujTemp += tau;
		if (flag_decrease) {
			if (gsl_vector_get(u,j)<0.0) {
				gsl_vector_set(u,j,0.0);
			}
		}
		
		//Update (or recompute) Cholesky factor and var
		if (iter == (int)(iter/n50000)*n50000) {
			
			gsl_matrix_set_zero(M); getXUXT(XX,uold,M,1);
			double normdiff = calcNormDiff(M,R,factor);
			if (normdiff> SMALLE)
			{	
				flag_recompute = 1;
			}
		}  
		if (flag_recompute) {
			numRecomputes++;
			//printf("%d\n",numRecomputes);
			gsl_matrix_set_zero(M); getXUXT(XX,u,M,0);
			gsl_matrix_set_zero(R); getChol(M,R); //now R holds the upper and lower triangular parts
			factor = 1.0;
			if (iter == (int)((double)(iter)/n100)*n100+1||iter==1||iter == (int)(iter/n50000)*n50000) {
				gsl_matrix_free(RX); RX = gsl_matrix_alloc(XX->size1, XX->size2); gsl_matrix_memcpy(RX, XX);
			}
			
			My_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, R, RX);
			
			gsl_vector_free(var); var = gsl_vector_alloc(XX->size2); getVar(var,RX);
			gsl_matrix_memcpy(RX, XX);
			
		}
		else {
			updateR(R, &factor, z,xj, &tau);
			double mult = tau / (1.0 - tau + tau * mvar);
			updateVar(var, &tau, &mult,Rxj, XX);
		}
		
		//Update maxvar.
		getAbsMax(var, &maxvar, &maxj);
		
		//Use the Harman-Pronzato test to see if further columns can be eliminated.
		if (iter == (int)((double)(iter)/n100)*n100) {
			ept = maxvar - n;
			tresh=n * (1 + 0.5*ept-0.5*sqrt(ept*(4.0+ept-4.0/n)));
			gsl_vector_int_free(e); e = gsl_vector_int_calloc(var->size);
			findThresh1or2(var, u,e,tresh, SMALLE);
			if (sum(e)<mm) {
				mm = sum(e);
				gsl_vector_int_free(acttemp); acttemp = gsl_vector_int_calloc(sum(e)); getSubVector(act,e,acttemp); 
				gsl_vector_int_free(act); act = gsl_vector_int_alloc(acttemp->size); gsl_vector_int_memcpy(act,acttemp);
				gsl_matrix_free(newXX); newXX = gsl_matrix_alloc(XX->size1,act->size); getSubMatrixFromColumns(X,act,newXX);
				gsl_matrix_free(XX); XX = gsl_matrix_alloc(newXX->size1,newXX->size2); gsl_matrix_memcpy(XX,newXX);
				gsl_matrix_free(RX); RX = gsl_matrix_alloc(XX->size1, XX->size2); gsl_matrix_memcpy(RX, XX);
				
				if (mm==n) {
					gsl_vector_free(u); u = gsl_vector_alloc(n); gsl_vector_set_all(u,1.0/n);
					gsl_vector_free(uold); uold = gsl_vector_alloc(u->size); gsl_vector_memcpy(uold,u);
					gsl_matrix_set_zero(M); getXUXT(XX,u,M,0);
					gsl_matrix_set_zero(R); getChol(M,R);
					factor = 1.0;
					My_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, R, XX); //Now XX is RX
					gsl_vector_free(var); var = gsl_vector_alloc(XX->size2); getVar(var,XX);
					getAbsMax(var,&maxvar,&maxj);
				}
				else {
					gsl_vector_free(vartemp); vartemp = gsl_vector_calloc(sum(e)); getSubVector(var,e,vartemp); 
					gsl_vector_free(var); var=gsl_vector_alloc(vartemp->size); gsl_vector_memcpy(var,vartemp);
					gsl_vector_free(utemp); utemp = gsl_vector_calloc(sum(e)); getSubVector(u,e,utemp); 
					gsl_vector_free(u); u=gsl_vector_alloc(utemp->size); gsl_vector_memcpy(u,utemp);
					double sumu = blas_sum(u); My_dscal(u,1.0/sumu);
					gsl_vector_free(uoldtemp); uoldtemp = gsl_vector_calloc(sum(e)); getSubVector(uold,e,uoldtemp); 
					gsl_vector_free(uold); uold=gsl_vector_alloc(uoldtemp->size); gsl_vector_memcpy(uold,uoldtemp);
					double sumuold = blas_sum(uold); My_dscal(uold,1.0/sumuold);
					getAbsMax(var,&maxvar,&maxj);
				}
				oldmm = mm;
			}
		}
		minOnSubset(var,u,&minvar,&minj); mvup = minvar;
		if (KKY == 1) {
			mvup = n;
		}
		
	}
	gsl_matrix_set_zero(M); 
	getXUXT(XX,u,M,1); 
	gsl_matrix_memcpy(M_out, M);
	if (doPrint)
		print(M);
	gsl_vector_set_all(u_out, 0.0);
	for (int ii=0; ii<act->size; ii++) {
		gsl_vector_set(u_out, gsl_vector_int_get(act, ii), gsl_vector_get(u, ii));
	}
	if (false) {
		printf("%d\n",iter);
		printf("%d\n",numRecomputes);
	}
	
	
	//do some cleaning....
	gsl_matrix_free(XX);
	gsl_matrix_free(newXX); 
	gsl_vector_free(xj); 
	gsl_vector_free(uold);
	gsl_matrix_free(M);
	gsl_vector_free(u);
	gsl_matrix_free(R);
	gsl_matrix_free(RX);
	gsl_vector_free(var);
	gsl_vector_int_free(act);
	gsl_vector_int_free(e);
	gsl_vector_int_free(smallAct);
	gsl_vector_free(smallVar);
	gsl_vector_free(u_e);
	gsl_vector_free(Rxj);
	gsl_vector_int_free(acttemp);
	gsl_vector_free(vartemp);
	gsl_vector_free(utemp);
	gsl_vector_free(uoldtemp);
	gsl_vector_free(upos);
	gsl_vector_free(z);
	if (doPrint) {
		printf("Number of iterations %d \n",iter);
	}
	return GSL_SUCCESS;
}

bool OptimalityTest(const gsl_matrix* X, const gsl_matrix* M,double tol, const gsl_vector* u)
{
	int n = X->size1;
	int m = X->size2;
	gsl_matrix* invM = gsl_matrix_calloc(n,n);
	gsl_matrix_memcpy(invM,M);	
	int info=0;
	char lower = 'U';
	int lda = invM->tda;
	dpotrf_(&lower, &n, invM->data, &lda, &info);
	dpotri_(&lower, &n, invM->data, &lda, &info);
	for (int i=0; i<n; i++) {
		for (int j=i+1 ; j<n; j++) {
			gsl_matrix_set(invM,i,j,gsl_matrix_get(invM,j,i));
		}
	}
	//My_dgemm(CblasNoTrans, CblasNoTrans, 1.0,M,invM, 1.0, eye);
	gsl_matrix* xi = gsl_matrix_calloc(n,1);
	gsl_matrix* xiTinvM = gsl_matrix_calloc(1,n);
	gsl_matrix* xiTinvMxi = gsl_matrix_calloc(1,1);
	for (int i=0; i<m; i++) {
		for (int j=0; j<n; j++) {
			gsl_matrix_set(xi,j,0,gsl_matrix_get(X,j,i));
		}
		My_dgemm(CblasTrans, CblasNoTrans, 1.0,xi,invM, 0.0,xiTinvM);
		My_dgemm(CblasNoTrans, CblasNoTrans, 1.0,xiTinvM,xi,0.0,xiTinvMxi);
		if (gsl_matrix_get(xiTinvMxi,0,0)>=n*(1.0-tol)) {		
			if (gsl_matrix_get(xiTinvMxi,0,0)<=(n*(1.0+tol))) {
				if (gsl_vector_get(u,i)==0.0) {
					return false;
				}
			}
			else {
				return false;
			}
		}
	}
	return true;
}

void initwt(const gsl_matrix* X, gsl_vector* u)
{
	int n = X->size1;
	int m = X->size2;
	gsl_vector_set_all(u,0.0);
	gsl_matrix* Q = gsl_matrix_alloc(n,n); gsl_matrix_set_identity(Q);
	gsl_vector_view d = gsl_matrix_column(Q,0);
	
	gsl_vector* dX = gsl_vector_alloc(m);
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_vector* y22 = gsl_vector_alloc(n);
	
	for (int j=0; j<n; j++) {
		My_dgemv(CblasTrans, 1.0, X, &d.vector, 0.0, dX);
		double maxdX; int indmax;
		double mindX; int indmin;
		getMax(dX, &maxdX, &indmax);
		getMin(dX, &mindX, &indmin);
		gsl_vector_set(u,indmax,1.0);
		gsl_vector_set(u,indmin,1.0);
		if (j==n-1) {
			break;
		}
		gsl_vector_const_view y1 = gsl_matrix_const_column(X, indmax);
		gsl_vector_memcpy(y, &y1.vector);
		gsl_vector_const_view y2 = gsl_matrix_const_column(X, indmin);
		gsl_vector_memcpy(y22, &y2.vector);
		gsl_vector_sub(y, y22);
		gsl_vector* z = gsl_vector_calloc(n-j);
		gsl_matrix_view Qjn = gsl_matrix_submatrix(Q, 0, j, n, n-j);
		My_dgemv(CblasTrans, 1.0, &Qjn.matrix, y, 0.0, z);
		double zeta = norm2(z); 
		if (gsl_vector_get(z, 0)<0.0) {
			gsl_vector_scale(z,-1.0);
		}
		gsl_vector* v = gsl_vector_alloc(z->size); gsl_vector_memcpy(v, z);
		gsl_vector_set(v,0,gsl_vector_get(v, 0)+zeta);
		double w = 1.0/gsl_vector_get(v,0)/zeta;
		gsl_vector* Qjnvw = gsl_vector_alloc(n);
		My_dgemv(CblasNoTrans, 1.0, &Qjn.matrix, v, 0.0, Qjnvw);
		My_dger(-w, Qjnvw, v, &Qjn.matrix);
		d = gsl_matrix_column(Q, j+1);
		gsl_vector_free(z);
		gsl_vector_free(v);
		gsl_vector_free(Qjnvw);
	}
	gsl_vector_scale(u,1.0/sum(u));
	gsl_matrix_free(Q);
	gsl_vector_free(dX);
	gsl_vector_free(y);
	gsl_vector_free(y22);
}

void initwt6(const gsl_matrix* X, gsl_vector* u)
{
	int n = X->size1;
	int m = X->size2;
	
	gsl_vector* ind = gsl_vector_alloc(m);
	for (int j=0; j<m; j++) {
		gsl_vector_set(ind, j, j);
	}
	gsl_vector* rho = gsl_vector_alloc(m);
	gsl_matrix* XX = gsl_matrix_alloc(n, m); gsl_matrix_memcpy(XX, X);
	gsl_matrix_mul_elements(XX, XX);

	gsl_vector* ones = gsl_vector_alloc(n);
	gsl_vector_set_all(ones, 1.0);
	My_dgemv(CblasTrans,1.0,XX,ones,0.0,rho);
	double meanrho = sqrt(sum(rho)/(double)m);

	// XX = [meanrho*ones(1,m); XX];
	gsl_matrix* XXtemp = gsl_matrix_alloc(n+1, m);
	for (int i=0; i<m; i++) {
		gsl_matrix_set(XXtemp,0,i,meanrho);
	}
	for (int i=0; i<n; i++) {
		for (int j=0; j<m; j++) {
			gsl_matrix_set(XXtemp,i+1,j,gsl_matrix_get(XX,i,j));
		}
	}
	gsl_matrix_free(XX); XX = gsl_matrix_alloc(n+1, m);
	gsl_matrix_memcpy(XX, XXtemp); gsl_matrix_free(XXtemp);
	
	for (int j=0; j<n+1; j++) {
		double maxrho=0.0;
		int maxind=0;
		getMax(rho, &maxrho, &maxind);
		gsl_vector* v = gsl_vector_alloc(n+1);
		gsl_matrix_get_col(v, XX, maxind);
		gsl_matrix_swap_columns(XX, maxind, j);

		int jj = gsl_vector_get(ind, maxind);
		gsl_vector_set(ind, maxind, gsl_vector_get(ind, j));
		gsl_vector_set(ind, j, jj);
		
		gsl_vector_set(rho, maxind, gsl_vector_get(rho, j));
		gsl_vector_set(rho, j, maxrho);
		
		gsl_vector* z = gsl_vector_alloc(n+1-j);
		for (int kk=j; kk<n+1; kk++) {
			gsl_vector_set(z, kk-j, gsl_vector_get(v, kk));
		}
		double zeta = norm2(z);
		if (gsl_vector_get(z, 0)<0.0) {
			gsl_vector_scale(z, -1.0);
		}
		gsl_vector_free(v); v = gsl_vector_alloc(z->size); gsl_vector_memcpy(v, z);		
		gsl_vector_set(v, 0, gsl_vector_get(v, 0)+zeta);
		double w = 1.0/zeta/gsl_vector_get(v, 0);
		gsl_matrix_view XXjnjm = gsl_matrix_submatrix(XX, j, j, n+1-j, m-j);
		gsl_vector* vXX = gsl_vector_alloc(m);
		My_dgemv(CblasTrans, 1.0, &XXjnjm.matrix, v, 0.0, vXX);
		My_dger(-w, v, vXX, &XXjnjm.matrix);
		gsl_vector_view rhojm = gsl_vector_subvector(rho, j, m-j);
		gsl_vector* XXjjm = gsl_vector_alloc(m-j);
		for (int kk=0; kk<m-j; kk++) {
			gsl_vector_set(XXjjm, kk, gsl_matrix_get(XX, j, j+kk)*gsl_matrix_get(XX, j, j+kk));
		}
		gsl_vector_sub(&rhojm.vector,XXjjm);
		gsl_vector_free(v);
		gsl_vector_free(z);
		gsl_vector_free(vXX);
		gsl_vector_free(XXjjm);
	}
	gsl_vector_set_all(u, 0.0);
	for (int i=0; i<n+1; i++) {
		gsl_vector_set(u, gsl_vector_get(ind, i), 1.0);
	}
	gsl_vector_scale(u, 1.0/(n+1.0));
	
	
	gsl_vector_free(rho);
	gsl_matrix_free(XX);
	gsl_vector_free(ones);
	gsl_vector_free(ind);
}

// Calculate what happens id you add point k to the set
// On output:
// 0 means that the point to be added is an interior point
// 1 otherwise
int calculateBoundAdd(const gsl_matrix* M, int k, const gsl_matrix* X, const double& optimal, double& bound)
{
	int n = X->size1;
	gsl_matrix* R = gsl_matrix_alloc(n, n);
	getChol(M,R);
	gsl_vector* RX = gsl_vector_alloc(n); gsl_matrix_get_col(RX, X, k);
	My_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, R, RX);
	double ksi_k = norm2(RX); ksi_k = ksi_k*ksi_k;
	gsl_vector_free(RX);
	gsl_matrix_free(R);
	if (ksi_k<=n) {
		bound = optimal;
		return 0;
	}
	else {
		double tau = ksi_k/(double)n - 1.0;
		tau =tau/(ksi_k - 1.0);
		bound = optimal + (n-1.0)*log(1.0-tau)+log(1+tau*ksi_k-tau);
		return 1;
	}	
}

//Calculate what happens if you remove point k from the set
void calculateBoundRemove(const gsl_matrix* M, int k, const gsl_vector* u, const double& optimal, double& bound)
{
	int n = M->size1;
	double u_k = gsl_vector_get(u, k);
	if (u_k==0.0) {
		bound = optimal;
		return;
	}
	else {
		bound = optimal - n*gsl_log1p(-u_k)+gsl_log1p(-n*u_k);
		return;
	}	
}

//Calculate what happens if you add point l and remove point k
void calculateBoundAddRemove(const gsl_matrix* M, int k , int l, const gsl_matrix* X, const gsl_vector* u, const double& optimal, double& bound)
{
	int n = X->size1;
	gsl_matrix* R = gsl_matrix_alloc(n, n);
	getChol(M,R);
	gsl_vector* RXk = gsl_vector_alloc(n); gsl_matrix_get_col(RXk, X, k);
	My_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, R, RXk);
	double ksi_k = norm2(RXk); ksi_k = ksi_k*ksi_k;
	gsl_vector* RXl = gsl_vector_alloc(n); gsl_matrix_get_col(RXl, X, l);
	My_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, R, RXl);
	double ksi_l = norm2(RXl); ksi_l = ksi_l*ksi_l;
	double ksi_kl = cblas_ddot(n, RXk->data, RXk->stride, RXl->data, RXl->stride);
	double u_k = gsl_vector_get(u, k);
	
	double delta = -1.0 + u_k*ksi_k + (1.0-u_k)*ksi_l - u_k*(1.0 - u_k)*(ksi_k * ksi_l - ksi_kl*ksi_kl);
	double tau2 = GSL_MAX(0.0,(delta - n + n*ksi_k-u_k*ksi_k)/n/delta);
	bound = optimal +(n-1.0)*gsl_log1p(-tau2)-n*gsl_log1p(-u_k)+gsl_log1p(-u_k*ksi_k+tau2*delta);
	gsl_vector_free(RXk);
	gsl_vector_free(RXl);
	gsl_matrix_free(R);
}

void calculateBoundAddRemove(const gsl_matrix* M, int k , int l,int i, const gsl_matrix* X, const gsl_vector* u, const double& optimal, double& bound)
{
	int n = X->size1;
	gsl_matrix* R = gsl_matrix_alloc(n, n);
	getChol(M,R);
	gsl_vector* RXk = gsl_vector_alloc(n); gsl_matrix_get_col(RXk, X, k);
	My_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, R, RXk);
	double ksi_k = norm2(RXk); ksi_k = ksi_k*ksi_k;
	gsl_vector* RXl = gsl_vector_alloc(n); gsl_matrix_get_col(RXl, X, l);
	My_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, R, RXl);
	double ksi_l = norm2(RXl); ksi_l = ksi_l*ksi_l;
	double ksi_kl = cblas_ddot(n, RXk->data, RXk->stride, RXl->data, RXl->stride);
	double u_k = gsl_vector_get(u, i); //The most important point!
	
	double delta = -1.0 + u_k*ksi_k + (1.0-u_k)*ksi_l - u_k*(1.0 - u_k)*(ksi_k * ksi_l - ksi_kl*ksi_kl);
	double tau2 = GSL_MAX(0.0,(delta - n + n*u_k*ksi_k+1.0-u_k*ksi_k)/n/delta);
	bound = optimal +(n-1.0)*log(1.0-tau2)-n*log(1.0-u_k)+gsl_log1p(-u_k*ksi_k+tau2*delta);
	gsl_vector_free(RXk);
	gsl_vector_free(RXl);
	gsl_matrix_free(R);
}



void getRandomMatrixMVEE(gsl_matrix* X, gsl_rng* r, double& optimalValue, int h)
{
	int n =X->size1;
	int m = X->size2;
	gsl_matrix* Y = gsl_matrix_alloc(n, m);
	getRandomMatrix(Y,r);
	gsl_matrix* Ytemp = gsl_matrix_alloc(n,m);
	orderMatrix(Y, Ytemp); gsl_matrix_free(Y);
	Y = gsl_matrix_alloc(n, h);
	for (int i=0; i<n; i++) {
		for (int j=0; j<h; j++) {
			gsl_matrix_set(Y, i, j, gsl_matrix_get(Ytemp, i, j));
		}
	} 
	gsl_matrix_free(Ytemp);
	gsl_vector_const_view Y0 = gsl_matrix_const_column(Y, h-1);
	double minNorm = norm2(&Y0.vector); printf("%e\n",minNorm);
	gsl_vector_const_view Y1 = gsl_matrix_const_column(Y, 0);
	double minNorm2 = norm2(&Y1.vector); printf("%e\n",minNorm2);
	//gsl_matrix_free(Ytemp);
	gsl_vector* uu = gsl_vector_alloc(h);
	gsl_matrix* M_out = gsl_matrix_alloc(n, n);
	minVol(Y, uu, M_out,0.0000001,0,10000,false);
	gsl_matrix* invM = gsl_matrix_calloc(n,n);
	gsl_matrix_memcpy(invM,M_out);	
	int info=0;
	char lower = 'U';
	int lda = invM->tda;
	dpotrf_(&lower, &n, invM->data, &lda, &info);
	dpotri_(&lower, &n, invM->data, &lda, &info);
	for (int i=0; i<n; i++) {
		for (int j=i+1 ; j<n; j++) {
			gsl_matrix_set(invM,i,j,gsl_matrix_get(invM,j,i));
		}
	}
	
	gsl_matrix* xi = gsl_matrix_calloc(n,1);
	gsl_matrix* xiTinvM = gsl_matrix_calloc(1,n);
	gsl_matrix* xiTinvMxi = gsl_matrix_calloc(1,1);
	int count = h;
	/*
	while (count!=m) {
		getRandomMatrix(xi,r);
		My_dgemm(CblasTrans, CblasNoTrans, 1.0,xi,invM, 0.0,xiTinvM);
		My_dgemm(CblasNoTrans, CblasNoTrans, 1.0,xiTinvM,xi,0.0,xiTinvMxi);
		if (gsl_matrix_get(xiTinvMxi,0,0)>=n) {
			for (int i=0; i<n; i++) {
				gsl_matrix_set(X, i, count, gsl_matrix_get(xi, i, 0));
			}
			count++;
		}
	}*/
	while (count!=m) {
		getRandomMatrix(xi,r);
		gsl_vector_const_view xivec = gsl_matrix_const_column(xi, 0);
	    double normxi = norm2(&xivec.vector);
		if (normxi>=minNorm) {
			for (int i=0; i<n; i++) {
				gsl_matrix_set(X, i, count, gsl_matrix_get(xi, i, 0));
			}
			count++;
		}
		else {
			//printf("failure\n");
		}

	}
	for (int i=0; i<n; i++) {
		for (int j=0; j<h; j++) {
			gsl_matrix_set(X, i, j, gsl_matrix_get(Y, i, j));
		}
	}
	gsl_matrix* LU=gsl_matrix_alloc(n,n);
	gsl_matrix_memcpy(LU, M_out);
	gsl_permutation* p = gsl_permutation_alloc(n);
	int signum = 0;
	gsl_linalg_LU_decomp (LU, p, &signum);
	optimalValue = log(gsl_linalg_LU_det(LU,signum));
	gsl_matrix_free(LU);
	gsl_permutation_free(p);
	gsl_matrix_free(Y);
	gsl_matrix_free(xi);
	gsl_matrix_free(xiTinvM);
	gsl_matrix_free(xiTinvMxi);
	gsl_matrix_free(invM);
	gsl_matrix_free(M_out);
}