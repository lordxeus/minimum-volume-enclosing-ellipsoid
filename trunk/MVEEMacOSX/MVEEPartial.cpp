#include "MVEEPartial.h"
#include "gsl/gsl_sort.h"
#include "gsl/gsl_sort_vector.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
#include "Cholesky.h"


//Basic constructor
MVEEPartial::MVEEPartial(int h, gsl_matrix* XX, const std::string& _solve_type): MVEE(XX), h(h), solve_type(_solve_type)
{
	heuristicSolutionIndices = gsl_vector_int_alloc(h);
	numTrials = 0;
	numBBTrials = 0;
}

//Detailed constructor
MVEEPartial::MVEEPartial(gsl_matrix* XX,int h, const double& tol, const int& KKY, const int& maxit, const bool& doPrint)
: MVEE(XX,tol,KKY,maxit,doPrint), h(h)
{
	heuristicSolutionIndices = gsl_vector_int_alloc(h);
	numBBTrials = 0;
}

//destructor - cleans up the only piece of memory on the heap
MVEEPartial::~MVEEPartial() 
{
	gsl_vector_int_free(heuristicSolutionIndices);
}

//Basic solution method
void MVEEPartial::solve()
{	
	if (solve_type=="COMPLETE_ENUMERATION")
		solveCompleteEnumeration();
	if (solve_type=="HEURISTIC_EID")
		solveHeuristic();
	if (solve_type=="HEURISTIC_PEELING")
		solveHeuristic();
	if (solve_type=="HEURISTIC_RANDOM_PEELING")
		solveHeuristic();
	if (solve_type=="BRANCH_AND_BOUND")
		solveBandB();
}

//solves the MVEE defined by X(:,c)
double MVEEPartial::solveSubset(MyCombination* c, gsl_matrix* M_out)
{
	int subsetSize = c->k;
	gsl_matrix* tempX = gsl_matrix_alloc(n, subsetSize);
	gsl_vector_int* indices = gsl_vector_int_alloc(subsetSize);
	for (int i=0; i<subsetSize; i++) {
		gsl_vector_int_set(indices, i, gsl_vector_int_get(c->data, i));
	} 
	//print(indices);
	getSubMatrixFromColumns(X,indices, tempX);
	MVEE tempMVEE(tempX,tol,KKY,maxit,doPrint);
	tempMVEE.solve();
	tempMVEE.calculateOptimalValue();
	gsl_matrix_memcpy(M_out, tempMVEE.M);
	//printf("Opt val: %e\n",tempMVEE.optimalValue);
	gsl_matrix_free(tempX);
	gsl_vector_int_free(indices);
	return tempMVEE.optimalValue;
	
}

//solves the MVEE defined by X(:,c)
double MVEEPartial::solveSubset(const gsl_vector_int* c, gsl_matrix* M_out, gsl_vector* u_out)
{
	int subsetSize = c->size;
	gsl_matrix* tempX = gsl_matrix_alloc(n, subsetSize);
	getSubMatrixFromColumns(X,c, tempX);
	MVEE tempMVEE(tempX,tol,KKY,maxit,doPrint);
	tempMVEE.solve();
	tempMVEE.calculateOptimalValue();
	gsl_matrix_memcpy(M_out, tempMVEE.M);
	tempMVEE.get_u(u_out);
	gsl_matrix_free(tempX);
	return tempMVEE.optimalValue;	
}

// solves the Partial MVEE using complete enumeration
// Only works for small instances
int MVEEPartial::solveCompleteEnumeration()
{
	double runningBest = GSL_POSINF; //set the running best to infinity
	MyCombination c(m,h,h);
	gsl_matrix* tempX = gsl_matrix_alloc(n, h);
	do {
		gsl_vector_int* indices = gsl_vector_int_alloc(h);
		for (int i=0; i<h; i++) {
			gsl_vector_int_set(indices, i, gsl_vector_int_get(c.data, i));
		} 
		getSubMatrixFromColumns(X,indices, tempX);
		MVEE tempMVEE(tempX);
		tempMVEE.solve();
		tempMVEE.calculateOptimalValue();
		runningBest = GSL_MIN(runningBest,tempMVEE.optimalValue);
		//printf("OptimalValue is:  %.5e, current was %.5e \n", runningBest, tempMVEE.optimalValue);
		gsl_vector_int_free(indices);
	} while (c.advance()!=0);
	
	this->optimalValue = runningBest;
	return 1;
}

void MVEEPartial::solveHeuristic()
{
	int opt = 0;
	gsl_vector* u_out = gsl_vector_alloc(h); gsl_vector* u_outtemp = gsl_vector_alloc(h);
	gsl_matrix* M_out = gsl_matrix_alloc(n, n);
	gsl_matrix* M_outtemp = gsl_matrix_alloc(n, n);
	gsl_vector_int* A = gsl_vector_int_alloc(h);
	double optimal_out = 0.0;
	double bound = 0.0;
	if (solve_type=="HEURISTIC_EID") {
		EID(A);
		optimal_out = solveSubset(A, M_out, u_out);
		initialOptimalValue = optimal_out;
		gsl_matrix_memcpy(M, M_out);
	}
	if (solve_type=="HEURISTIC_PEELING") {
		ellipsoidalPeeling(u_out,M_out,optimal_out,A);
		initialOptimalValue = optimal_out;
	}
	if (solve_type=="HEURISTIC_RANDOM_PEELING") {
		randomEllipsoidalPeeling(A);
		optimal_out = solveSubset(A, M_out, u_out);
	}
	gsl_vector_int* notA = gsl_vector_int_alloc(m-h);
	gsl_vector_int* Ahat = gsl_vector_int_alloc(h);
	getNotI(A, notA);
	gsl_vector_int* cA = gsl_vector_int_alloc(h); gsl_vector_int_memcpy(cA, A);
	while (opt==0) {
		int counti=0;
		int countj=0;
		for (int i=0; i<h; i++) {
			counti=i;
			for (int j=0; j<m-h; j++) {
				countj=j;
				int curri = gsl_vector_int_get(cA, i);
				int currj = gsl_vector_int_get(notA, j);
				if (gsl_vector_get(u_out, i)>0.00000000000001) {
					calculateBoundAddRemove(M_out, curri , currj,i, X, u_out, optimal_out, bound);
					if (bound<optimal_out) {
						gsl_vector_int_memcpy(Ahat, cA);
						gsl_vector_int_set(Ahat, i, currj);
						double cOptimal = solveSubset(Ahat,M_outtemp, u_outtemp);
						if (cOptimal<optimal_out) {
							optimal_out = cOptimal;
							gsl_permutation* p = gsl_permutation_alloc(h);
							gsl_sort_vector_int_index (p, Ahat);
							gsl_matrix_memcpy(M_out, M_outtemp);
							for (int kk=0; kk<h; kk++) {
								gsl_vector_int_set(cA, kk, gsl_vector_int_get(Ahat,gsl_permutation_get(p, kk)));
								gsl_vector_set(u_out, kk, gsl_vector_get(u_outtemp,gsl_permutation_get(p, kk)));
							}
							opt = 0;
							gsl_vector_int_memcpy(A, cA);
							getNotI(A, notA);
							gsl_permutation_free(p);
							break;
						}
						
					}
				}
			}
			if (counti==h-1 && countj==m-h-1) 
			{opt = 1;
				break; }
		}
		 
	}
	gsl_vector_int_memcpy(heuristicSolutionIndices, A);
	this->optimalValue = optimal_out;
	solveSubset(heuristicSolutionIndices,M_out, u_out);
	//gsl_matrix_memcpy(this->M, M_out);
	gsl_vector_free(u_out);
	gsl_vector_free(u_outtemp);
	gsl_matrix_free(M_out);
	gsl_matrix_free(M_outtemp);
	gsl_vector_int_free(A);
	gsl_vector_int_free(notA);
	gsl_vector_int_free(cA);
	gsl_vector_int_free(Ahat);
}

// "driver" method for solving using Branch and Bound
void MVEEPartial::solveBandB()
{	
	gsl_vector* u_out = gsl_vector_alloc(h);
	gsl_matrix* M_out = gsl_matrix_alloc(n, n);
	gsl_vector_int* A_out = gsl_vector_int_alloc(h);
	EID(A_out);
	print(A_out);
	upperBound = solveSubset(A_out, M_out, u_out);
	optimalValue = upperBound;
	//upperBound = -1.54;
	MyCombination* d = new MyCombination(m,n,h);
	this->solveRecursiveCore(d);
	printf("numtrials: %d\n",numTrials);
	delete d;
	gsl_vector_free(u_out);
	gsl_matrix_free(M_out);
	gsl_vector_int_free(A_out);
	
}

// main function of the Branch and Bound - part 1
// this is only for sequences of size n
void MVEEPartial::solveRecursiveCore(MyCombination* c)
{	
	//if (gsl_vector_int_get(c->data, 1)==9) {
	//	int asd = 2;
	//	printf("%d\n",asd);
	//}
	if (gsl_vector_int_get(c->data, c->k-1)>c->n-1+c->k-h) {
		if (c->advance()==0) {
			return;
		}
		else {
			solveRecursiveCore(c);
		}
	}
	numTrials++;
	//printf("%d\n",numTrials);
	//print(c->data);
	gsl_matrix* M_out = gsl_matrix_alloc(n, n);
	double currentValue = solveSubset(c,M_out);
	//print(c->data); printf("%e %e\n",currentValue,upperBound);
	
	//if the currentValue is bigger than the upper bound ignore that branch
	if (currentValue>=upperBound) {
		//if there is nowhere to go next, end it
		if	(c->advance()==0)
		{
			optimalValue = upperBound;
			gsl_matrix_free(M_out);
			return;
		}
		//otherwise move to the next node
		else {
			solveRecursiveCore(c);
		}
	}
	//if the currentValue is smaller than the upper bound we still have hope
	else {
		//if the current node is of the correct size
		if (c->k==h) {
			//update the upperBound
			upperBound = currentValue;
			optimalValue = upperBound;
			//print(c->data);
			//printf("%e  ",optimalValue); printf("%d\n",numTrials);
			//try to move
			if (c->advance()==0) {
				//if it fails go back
				gsl_matrix_free(M_out);
				return;
			}
			//otherwise go to the next node
			else {
				solveRecursiveCore(c);
			}
		}
		else {
			if (gsl_vector_int_get(c->data, c->k-1)>c->n-1+c->k-h) {
				//In this case we cannot move to higher dimension
				if (c->advance()==0) {
					gsl_matrix_free(M_out);
					return;
				}
				else {
					//printf("are we ever here? \n");
					solveRecursiveCore(c);
				}
			}
			else {
				MyCombination d = c->getHigherDim();
				double bound = 0.0;
				calculateBoundAdd(M_out,gsl_vector_int_get(d.data, d.k-1),X, currentValue, bound);
				//print(d.data); printf("bound--  %e %e\n",currentValue,bound);
				if (bound>=upperBound) {
					while (gsl_vector_int_get(d.data, d.k-1)<m-1+d.k-h) {
						gsl_vector_int_set(d.data, d.k-1, gsl_vector_int_get(d.data, d.k-1)+1);
						calculateBoundAdd(M_out,gsl_vector_int_get(d.data, d.k-1),X, currentValue, bound);
						//print(d.data); printf("bound--  %e %e\n",currentValue,bound);
						if (bound<upperBound) {
							solveRecursiveCore2(&d,currentValue);
							break;
						}
					}
					//the bounds is already too big, no need to bother entering
					if (c->advance()==0) {
						gsl_matrix_free(M_out);
						return;
					}
					else {
						solveRecursiveCore(c);
					}

				}
				else {
					solveRecursiveCore2(&d);
					if (c->advance()==0) {
						gsl_matrix_free(M_out);
						return;
					}
					else {
						solveRecursiveCore(c);
					}
				}
			}
		}
	}
	gsl_matrix_free(M_out);
}

// main function of the Branch and Bound - part 2
//The reason for this, is that when we are in a child, we don't advance at will,
// but we only go until we reach the last element
// and we should stop as soon as we reach the same value as the father
void MVEEPartial::solveRecursiveCore2(MyCombination* c, double fathersValue,bool isInteriorPoint)
{
	numTrials++;
	gsl_matrix* M_out = gsl_matrix_alloc(n,n);
	double currentValue = solveSubset(c,M_out);
	//print(c->data); printf("%e %e\n",currentValue,upperBound);
	//if the currentValue is bigger than the upper bound chop that branch
	if (currentValue>=upperBound) {
		//if there is nowhere to go next, end it
		if	(gsl_vector_int_get(c->data, c->k-1)>c->n-1+c->k-h)
		{
			optimalValue = upperBound;
			gsl_matrix_free(M_out);
			return;
		}
		//otherwise move to the next node
		else {
			if (c->advance()==1) {
				solveRecursiveCore2(c,fathersValue);
			}
			else {
				gsl_matrix_free(M_out);
				return;
			}

		}
	}
	//if the currentValue is smaller than the upper bound we still have hope
	else {
		//if the current node is of the correct size
		if (c->k==h) {
			//update the upperBound
			upperBound = currentValue;
			optimalValue = upperBound;
			//print(c->data);
			//printf("%e  ",optimalValue); printf("%d\n",numTrials);
			//try to move
			if (gsl_vector_int_get(c->data, c->k-1)>c->n-1+c->k-h||tolZero(currentValue, fathersValue, tol)) {
				//if it fails go back
				gsl_matrix_free(M_out);
				return;
			}
			//otherwise go to the next node
			else {
				if	(c->advance()==1) {
					solveRecursiveCore2(c,fathersValue);
				}
				else {
					gsl_matrix_free(M_out);
					return;
				}

			}
		}
		else {
			if (gsl_vector_int_get(c->data, c->k-1)>c->n-1+c->k-h) {
				//In this case we cannot move to higher dimension
				gsl_matrix_free(M_out);
				return;
			}
			else {
				MyCombination d = c->getHigherDim();
				double bound = 0.0;
				calculateBoundAdd(M_out,gsl_vector_int_get(d.data, d.k-1),X, currentValue, bound);
				//print(d.data); printf("bound--  %e %e\n",currentValue,bound);
				if (bound>upperBound) {
					// no need to solve it
					// just move ahead until we find something to solve
					while (gsl_vector_int_get(d.data, d.k-1)<m-1+d.k-h) {
						gsl_vector_int_set(d.data, d.k-1, gsl_vector_int_get(d.data, d.k-1)+1);
						calculateBoundAdd(M_out,gsl_vector_int_get(d.data, d.k-1),X, currentValue, bound);
						//print(d.data); printf("bound--  %e %e\n",currentValue,bound);
						if (bound<upperBound) {
							solveRecursiveCore2(&d,currentValue);
							break;
						}
					}
					if (gsl_vector_int_get(c->data, c->k-1)>c->n-1+c->k-h) {
						gsl_matrix_free(M_out);
						return;
					}
					else {
						if (c->advance()==1) {
							solveRecursiveCore2(c,fathersValue);
						}
						else {
							gsl_matrix_free(M_out);
							return;
						}

					}
				}
				else {
					// go ahead and solve.. it may be useful
					solveRecursiveCore2(&d,currentValue);
					if (gsl_vector_int_get(c->data, c->k-1)>c->n-1+c->k-h) {
						gsl_matrix_free(M_out);
						return;
					}
					else {
						if(c->advance()==1) {
							solveRecursiveCore2(c,fathersValue);
						}
						else {
							gsl_matrix_free(M_out);
							return;
						}

					}
				}
			}
		}
	}
	gsl_matrix_free(M_out);
}

//main driver for ellipsoidal peeling
void MVEEPartial::ellipsoidalPeeling(gsl_vector* u_out,gsl_matrix* M_out, double& optimal_out, gsl_vector_int* A_out)
{
	gsl_vector_int* indices= gsl_vector_int_alloc(m);
	for (int i=0; i<m; i++) {
		gsl_vector_int_set(indices, i, i);
	}
	for (int i=m; i>h; i--) {
		gsl_matrix* tempX = gsl_matrix_alloc(n, i);
		getSubMatrixFromColumns(X,indices, tempX);
		MVEE tempMVEE(tempX,tol,KKY,maxit,doPrint);
		tempMVEE.solve();
		gsl_vector* uu = gsl_vector_alloc(i);
		tempMVEE.get_u(uu);	
		//print(uu);
		double maxValue = 0.0;
		int maxPosition = 0;
		getMax(uu, &maxValue, &maxPosition);
		//gsl_vector_free(indices);
		gsl_vector_int* tempIndices = gsl_vector_int_alloc(i-1);
		for (int j=0; j<i-1; j++) {
			if (gsl_vector_int_get(indices, j)<maxPosition) {
				gsl_vector_int_set(tempIndices, j, gsl_vector_int_get(indices, j));
			}
			else {
				gsl_vector_int_set(tempIndices, j, gsl_vector_int_get(indices, j+1));
			}
		}
		gsl_vector_int_free(indices); indices = gsl_vector_int_alloc(i-1);
		gsl_vector_int_memcpy(indices, tempIndices); gsl_vector_int_free(tempIndices);
		gsl_vector_free(uu);
		gsl_matrix_free(tempX);
	}
	gsl_vector_int_memcpy(A_out, indices);
	gsl_matrix* tempX = gsl_matrix_alloc(n, h);
	getSubMatrixFromColumns(X,indices, tempX);
	MVEE tempMVEE(tempX,tol,KKY,maxit,doPrint);
	tempMVEE.solve();
	tempMVEE.calculateOptimalValue();
	this->upperBound = tempMVEE.optimalValue;
	tempMVEE.get_u(u_out);
	gsl_matrix_memcpy(M_out, tempMVEE.M);
	optimal_out = tempMVEE.optimalValue;
}

//main driver for ellipsoidal peeling
void MVEEPartial::ellipsoidalPeeling()
{
	gsl_vector_int* indices= gsl_vector_int_alloc(m);
	for (int i=0; i<m; i++) {
		gsl_vector_int_set(indices, i, i);
	}
	for (int i=m; i>h; i--) {
		gsl_matrix* tempX = gsl_matrix_alloc(n, i);
		getSubMatrixFromColumns(X,indices, tempX);
		MVEE tempMVEE(tempX,tol,KKY,maxit,doPrint);
		tempMVEE.solve();
		gsl_vector* uu = gsl_vector_alloc(i);
		tempMVEE.get_u(uu);	
		double maxValue = 0.0;
		int maxPosition = 0;
		getMax(uu, &maxValue, &maxPosition);
		//gsl_vector_free(indices);
		gsl_vector_int* tempIndices = gsl_vector_int_alloc(i-1);
		for (int j=0; j<i-1; j++) {
			if (gsl_vector_int_get(indices, j)<maxPosition) {
				gsl_vector_int_set(tempIndices, j, gsl_vector_int_get(indices, j));
			}
			else {
				gsl_vector_int_set(tempIndices, j, gsl_vector_int_get(indices, j+1));
			}
		}
		gsl_vector_int_free(indices); indices = gsl_vector_int_alloc(i-1);
		gsl_vector_int_memcpy(indices, tempIndices); gsl_vector_int_free(tempIndices);
		gsl_vector_free(uu);
		gsl_matrix_free(tempX);
	}
	gsl_matrix* tempX = gsl_matrix_alloc(n, h);
	getSubMatrixFromColumns(X,indices, tempX);
	MVEE tempMVEE(tempX,tol,KKY,maxit,doPrint);
	tempMVEE.solve();
	tempMVEE.calculateOptimalValue();
	this->upperBound = tempMVEE.optimalValue;
}

//main driver for EID
void MVEEPartial::EID(gsl_vector_int* A_out)
{
	gsl_vector_int* A = gsl_vector_int_alloc(m);
	for (int i=0; i<m; i++) {
		gsl_vector_int_set(A, i, i);
	}
	while (A->size!=h) {
		gsl_vector* x = gsl_vector_calloc(n);
		for (int i=0; i<A->size; i++) {
			gsl_vector_const_view Xi = gsl_matrix_const_column (X, gsl_vector_int_get(A, i));
			gsl_vector_add(x,&Xi.vector);
		}
		gsl_vector_scale(x, 1.0/A->size);
		gsl_matrix* XBar = gsl_matrix_alloc(n, A->size); //DONE
		getSubMatrixFromColumns(X, A, XBar);
		for (int i=0; i<A->size; i++) {
			gsl_vector_view Xi = gsl_matrix_column (XBar, i);
			gsl_vector_sub(&Xi.vector,x);
		}
		gsl_matrix* EIDMat = gsl_matrix_alloc(n, n); //DONE
		My_dgemm(CblasNoTrans,CblasTrans,1.0,XBar,XBar,0.0,EIDMat);
		gsl_matrix* R = gsl_matrix_calloc(n, n); //DONE
		getChol(EIDMat, R);
		gsl_matrix* Xtemp = gsl_matrix_calloc(n, m); gsl_matrix_memcpy(Xtemp, X); //DONE
		My_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, R, Xtemp);
		gsl_matrix_free(EIDMat); EIDMat = gsl_matrix_alloc(m, m);
		My_dgemm(CblasTrans, CblasNoTrans, 1.0, Xtemp, Xtemp, 0.0, EIDMat);
		gsl_vector* EIDVec = gsl_vector_alloc(m); //DONE
		for (int i=0; i<m; i++) {
			gsl_vector_set(EIDVec, i, gsl_matrix_get(EIDMat, i, i));
		}
		gsl_matrix_free(EIDMat);
		gsl_matrix_free(Xtemp);
		gsl_matrix_free(XBar);
		gsl_vector_free(x);
		double maxValue; int maxPosition;
		gsl_vector* EIDVecTemp= gsl_vector_alloc(A->size); //DONE
		for (int i=0; i<A->size; i++) {
			gsl_vector_set(EIDVecTemp,i,gsl_vector_get(EIDVec,gsl_vector_int_get(A, i)));
		}
		gsl_vector_free(EIDVec);
		getMax(EIDVecTemp, &maxValue, &maxPosition);
		gsl_vector_free(EIDVecTemp);
		gsl_vector_int* Atemp =gsl_vector_int_alloc(A->size-1);
		
		for (int i=0; i<maxPosition; i++) {
			gsl_vector_int_set(Atemp, i, gsl_vector_int_get(A, i));
		}
		for (int i=maxPosition+1; i<A->size; i++) {
			gsl_vector_int_set(Atemp, i-1, gsl_vector_int_get(A, i));
		}
		gsl_vector_int_free(A); A= gsl_vector_int_alloc(Atemp->size);
		gsl_vector_int_memcpy(A, Atemp); gsl_vector_int_free(Atemp);
		gsl_matrix_free(R);
	}
	//print(A);
	gsl_vector_int_memcpy(A_out, A);
	gsl_vector_int_free(A);
}

//main driver for random ellipsoidal peeling
void MVEEPartial::randomEllipsoidalPeeling(gsl_vector_int* A_out)
{
	const gsl_rng_type* T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	int* a=new int[m];
	int* b=new int[n+1];
	for (int i = 0; i < m; i++)
	{
		a[i] = i;
	}
	gsl_ran_choose (r,b,h, a,m,sizeof(int));
	gsl_vector_int* A = gsl_vector_int_alloc(n+1);
	for (int i = 0; i < n+1; i++)
	{
		gsl_vector_int_set(A, i, b[i]);
	}
	delete a;
	delete b;
	gsl_matrix* M_out = gsl_matrix_alloc(n, n);
	gsl_vector* u_out = gsl_vector_alloc(n+1);
	solveSubset(A, M_out, u_out);
	gsl_vector_int_free(A);
	std::vector<int> Atemp;
	std::vector<double> us;
	gsl_matrix* H = gsl_matrix_calloc(n,n);
	gsl_matrix_memcpy(H,M_out);	
	int info=0;
	char lower = 'U';
	int lda = H->tda;
	dpotrf_(&lower, &n, H->data, &lda, &info);
	dpotri_(&lower, &n, H->data, &lda, &info);
	for (int k=0; k<n; k++) {
		for (int j=k+1 ; j<n; j++) {
			gsl_matrix_set(H,k,j,gsl_matrix_get(H,j,k));
		}
	}
	for (int i=0; i<m; i++) {
		gsl_matrix* xi = gsl_matrix_calloc(n,1);
		gsl_matrix* xiTinvM = gsl_matrix_calloc(1,n);
		gsl_matrix* xiTinvMxi = gsl_matrix_calloc(1,1);
		for (int j=0; j<n; j++) {
			gsl_matrix_set(xi,j,0,gsl_matrix_get(X,j,i));
		}
		My_dgemm(CblasTrans, CblasNoTrans, 1.0,xi,H, 0.0,xiTinvM);
		My_dgemm(CblasNoTrans, CblasNoTrans, 1.0,xiTinvM,xi,0.0,xiTinvMxi);
		us.push_back(gsl_matrix_get(xiTinvMxi,0,0));
		if (gsl_matrix_get(xiTinvMxi,0,0)<=n*(1.0+tol)) {		
			Atemp.push_back(i);
		}
		gsl_matrix_free(xi);
		gsl_matrix_free(xiTinvM);
		gsl_matrix_free(xiTinvMxi);
	}
	gsl_matrix_free(H);
	int i=0;
	while (Atemp.size()!=h) {
		if (us[i]>n*(1+tol)) {
			Atemp.push_back(i);
		}
		i++;
	} 
	std::sort(Atemp.begin(), Atemp.end());
	for (int ii=0; ii<h; ii++) {
		gsl_vector_int_set(A_out, ii, Atemp[ii]);
	}
}


// "driver" method for solving using Branch and Bound
void MVEEPartial::solveBandB2()
{	
	gsl_vector* u_out = gsl_vector_alloc(h);
	gsl_matrix* M_out = gsl_matrix_alloc(n, n);
	gsl_vector_int* A_out = gsl_vector_int_alloc(h);
	EID(A_out);
	upperBound = solveSubset(A_out, M_out, u_out);
	gsl_vector_free(u_out);
	gsl_matrix_free(M_out);
	gsl_vector_int_free(A_out);
	
	optimalValue = upperBound;
	//upperBound = 1;
	MyCombination* d = new MyCombination(m,n,h);
	this->solveRecursiveCore(d);
	printf("numtrials: %d\n",numTrials);
	delete d;
	
}