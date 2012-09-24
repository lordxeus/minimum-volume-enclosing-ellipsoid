#ifndef MVEE_PARTIAL_H
#define MVEE_PARTIAL_H
/*
 *  MVEEPartial.h
 *  LinAlg
 *
 *  Created by Spyros Schismenos on 15/07/2012.
 *
 */

#include "MVEE.h"

// This class represents a minimum volume partially enclosing ellipsoid problem
// It derived from the MVEE problem

class MVEEPartial : public MVEE {
public:
	int h; // the number of points needed
	std::string solve_type; // how to solve it
	gsl_vector_int* heuristicSolutionIndices;
	double upperBound;
	int numTrials;
	int numBBTrials;
	double initialOptimalValue;
public:
	MVEEPartial(int h, gsl_matrix* XX, const std::string& _solve_type); //default constructor
	MVEEPartial(gsl_matrix* XX,int h, const double& tol, const int& KKY, const int& maxit, const bool& doPrint);
	virtual void solve();
	virtual void calculateOptimalValue() {};
	virtual bool isOptimal(double tol) { printf("do not have an optimality test for partial inclusion\n"); return true;};
	~MVEEPartial();
	int solveCompleteEnumeration();
	void EID(gsl_vector_int* A_out);
	void get_M(gsl_matrix* M_out) {gsl_matrix_memcpy(M_out, M);};
private:
	void solveRecursiveCore(MyCombination* c);
	void solveRecursiveCore2(MyCombination* c, double fathersBound=0.0, bool isInteriorPoint=false);
	double solveSubset(MyCombination* c,gsl_matrix* M_out);
	double solveSubset(const gsl_vector_int* c,gsl_matrix* M_out, gsl_vector* u_out);
	void solveHeuristic();
	void solveBandB();
	void solveBandB2();
	void ellipsoidalPeeling(gsl_vector* u_out, gsl_matrix* M_out, double& optimal_out, gsl_vector_int* A_out);
	void ellipsoidalPeeling();
	void randomEllipsoidalPeeling(gsl_vector_int* A_out);

};

#endif