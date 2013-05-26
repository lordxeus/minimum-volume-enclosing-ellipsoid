#ifndef MVEE_H
#define MVEE_H
/*
 *  MVEE.h
 *  LinAlg
 *
 *  Created by Spyros Schismenos on 27/06/2012.
 *
 */

/*! The MVEE class represents a minimum volume ellipsoid problem
 and when it is solved, its solution as well
 */
#include "minVolume.h"
#include "gsl/gsl_combination.h"

class MVEE 
{

public:	
	
	double tol; //tolerance to solve for
	int KKY; //type of solver
	int maxit; // maximum number of iterations
	bool doPrint; //whether to print output during the solve
	
	int n; //!< the dimension of the problem
	int m; //!< the number of points
	gsl_matrix* X; //!< the points in a [numOfDims][numOfPoints] format
	
	int isSolved; //!< is true when the problem is solved succesfully
	double optimalValue; //!< the optimal value of the problem
	
	gsl_matrix* M; //!< the optimal matrix M
	gsl_vector* u; //!< the optimal vector u
	gsl_vector_int* activeIndices; //!< the active points of the solution
	
	virtual void solve();
	virtual bool isOptimal(double tol);
	virtual int get_n() const;
	virtual int get_m() const;
	void get_u(gsl_vector* uu) const;
	virtual void set_doPrint(bool _doPrint);
	void printResults() const;
	virtual void calculateOptimalValue();
	
	MVEE() {};
	MVEE(gsl_matrix* XX);
	MVEE(gsl_matrix* XX,const double& tol, const int& KKY, const int& maxit, const bool& doPrint);
	~MVEE();
};

class MVEEGeneric: public MVEE
{
public:
	
	MVEE* mveeUpDimension;
	MVEEGeneric(gsl_matrix* XX);
	MVEEGeneric(gsl_matrix* XX,const double& tol, const int& KKY, const int& maxit, const bool& doPrint);
	void moveUpDimension(gsl_matrix* XXnew, const gsl_matrix* XX);
	virtual void solve();
	//virtual void moveDownDimension();
};

#endif