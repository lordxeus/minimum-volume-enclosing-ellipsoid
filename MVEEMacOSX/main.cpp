#include <iostream>
#include <fstream>
#include <stdio.h>
#include "MVEE.h"
#include "MVEEPartial.h"

#include <vector>
#include <numeric>

int testMVEEPartial()
{
	const gsl_rng_type* T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	int n = 2;
	int m = 50;//m=5;
	int h = 28;
	int N = 10; //number of solutions
	double tol = 0.0000001;
	tol = 0.01;
	int res; 
	bool doPrint = false;
	gsl_matrix* X = gsl_matrix_calloc(n,m);
	gsl_matrix* Y = gsl_matrix_calloc(n, m);
	//FILE* f = fopen("problematic.dat", "rb");
	//gsl_matrix_fread(f, X);
	//fclose(f);
	//print(X);
	//getRandomMatrix(X,r);
	//getRandomMatrix(X, r);
	//getRandomMatrix(X, r);
	//getRandomMatrix(X, r);
	
	//printMATLAB(X);
	//FILE* f = fopen("20-10K.txt", "wb");
	//gsl_matrix_fprintf(f, X, "%.10e");
	//fclose(f);
	
	//MVEE instance1(X); 
	//instance1.set_doPrint(doPrint);
	//gsl_vector* u = gsl_vector_calloc(m);
	std::vector<double> timeEID;
	std::vector<double> timePEELING;
	std::vector<double> timeBB;
	std::vector<double> optimalEID;
	std::vector<double> optimalPEELING;
	std::vector<double> initialOptimalEID;
	std::vector<double> initialOptimalPEELING;
	std::vector<double> optimalBB;
	std::vector<double> optimalConstruction;
	std::vector<double> bestHeuristicPeeling;
	int sucEID=0;
	int sucPEELING=0;
	int noSuc=0;
	for (int i=0; i<N; i++) {
		printf("------iteration %d------\n",i);
		double ov = 0;
		//getRandomMatrixMVEE(X, r,ov,h);
		getRandomMatrix(X, r);
		optimalConstruction.push_back(ov);
		orderMatrix(X, Y);
		gsl_matrix_memcpy(X, Y);
		//gsl_matrix_free(Y);
		clock_t begin=clock();
		MVEEPartial instance1(h,X, "HEURISTIC_EID"); 	
		instance1.tol=tol;
		instance1.solve();
		clock_t end=clock();
		timeEID.push_back((double)diffclock(end,begin));
		optimalEID.push_back(instance1.optimalValue);
		initialOptimalEID.push_back(instance1.initialOptimalValue);
		printf("optimal value EID: %e\n",optimalEID[i]);
		printf("time for EID: %e ms\n",timeEID[i]);
		begin=clock();
		gsl_matrix* M_out = gsl_matrix_alloc(n, n);
		instance1.get_M(M_out);
		orderMatrix(X, Y,M_out);
		gsl_matrix_free(M_out);
		gsl_matrix_memcpy(X, Y);
		MVEEPartial instance2(h,X, "BRANCH_AND_BOUND");
		//gsl_matrix* M_out = gsl_matrix_alloc(n, n);
		//instance1.get_M(M_out);
		//orderMatrix(X, Y,M_out);
		//gsl_matrix_memcpy(X, Y);
		instance2.tol=tol;
		begin=clock();
		instance2.solve();
	    end=clock(); 		
		timePEELING.push_back((double)diffclock(end,begin));
		optimalPEELING.push_back(instance2.optimalValue);
		initialOptimalPEELING.push_back(instance2.initialOptimalValue);
		printf("optimal value BB: %e\n",optimalPEELING[i]);
		printf("time for BB: %e ms\n",timePEELING[i]);
		
		if (optimalEID[i]<optimalPEELING[i])
		{
			bestHeuristicPeeling.push_back(0);
		}
		else {
			bestHeuristicPeeling.push_back(1);
		}
		if (optimalEID[i]<=ov) {
			sucEID++;
		}
		if (optimalPEELING[i]<=ov) {
			sucPEELING++;
		}
		if (optimalEID[i]>ov&&optimalPEELING[i]>ov) {
			noSuc++;
		}
		//begin=clock();

		//MVEEPartial instance3(h,X, "HEURISTIC_RANDOM_PEELING");
		//res = instance3.solve();

		//end=clock(); std::cout << "Time elapsed: " << double(diffclock(end,begin)) << " ms"<< std::endl;

		//begin=clock();
		//MVEEPartial instance4(h,X,"BRANCH_AND_BOUND");
		//res = instance4.solve();

		//end=clock();
		//timeBB.push_back(diffclock(end, begin));
		//optimalBB.push_back(instance4.optimalValue);
		//printf("ov BB: %e\n",optimalBB[i]);
	}
	double  avgTimeEID = Myavg(timeEID);
	double avgOptimalEID =Myavg(optimalEID);
	double avgTimePEELING =Myavg(timePEELING);
	double avgOptimalPEELING =Myavg(optimalPEELING);
	double avgInitialOptimalEID = Myavg(initialOptimalEID); 
	double avgInitialOptimalPEELING = Myavg(initialOptimalPEELING);
	double avgOptimalConstruction = Myavg(optimalConstruction);
	double freqBestHeuristicPeeling = Myavg(bestHeuristicPeeling);
	std::cout<<"----------------------------------"<<std::endl;
	std::cout<< "n, m, h " << n << " " << m << " " <<h <<std::endl;
	printf("avg time EID %e\n",avgTimeEID);
	printf("avg time BB %e\n",avgTimePEELING);
	printf("avg optimal value EID %e\n",avgOptimalEID);
	printf("avg initial value EID %e\n",avgInitialOptimalEID);
	printf("avg optimal value BB %e\n",avgOptimalPEELING);
	//printf("avg volume peeling %e\n",avgInitialOptimalPEELING);
	//printf("avg volume construction %e\n",avgOptimalConstruction);
	//printf("how often is best Heuristic the peeling? %e\n",freqBestHeuristicPeeling);
	//printf("how often is peeling better than the default ellipsoid? %e\n",sucPEELING/(double)N);
	//printf("how often is EID better than the default ellipsoid? %e\n",sucEID/(double)N);
	//printf("how often none of them is better than the default ellipsoid? %e\n",noSuc/(double)N);
	//instance1.calculateOptimalValue();
	//std::cout << "optimal Value is : " << instance1.optimalValue << std::endl;;
	//gsl_matrix* M = gsl_matrix_alloc(n, n);
	//gsl_vector* u = gsl_vector_alloc(m);
	//instance1.get_u(u);
	//printf("%e\n",sum(u));

	//gsl_matrix_memcpy(M, instance1.M);
	
	//if (instance1.isOptimal(tol))
	//{
	//	printf("success\n");
	//}
	//else {
	//	printf("failure\n");
	//}

	gsl_matrix_free(X); 
	gsl_matrix_free(Y);
	gsl_rng_free (r);    
	//gsl_vector_free(u);
	//std::ofstream myfile;
	//myfile.open ("example.txt");
	//myfile << "Writing this to a file.\n";
	//myfile.close();
	return GSL_SUCCESS;
}

int testMVEEPartialEID()
{
	const gsl_rng_type* T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	int n = 5;
	int m = 500;//m=5;
	int h = int((n+m+1)/2)+1;
	int N = 10; //number of solutions
	double tol = 0.01;
	int res; 
	bool doPrint = false;
	gsl_matrix* X = gsl_matrix_calloc(n,m);
	gsl_matrix* Y = gsl_matrix_calloc(n, m);
	//FILE* f = fopen("problematic.dat", "rb");
	//gsl_matrix_fread(f, X);
	//fclose(f);
	//print(X);
	//getRandomMatrix(X,r);
	//getRandomMatrix(X, r);
	//getRandomMatrix(X, r);
	//getRandomMatrix(X, r);
	
	//printMATLAB(X);
	//FILE* f = fopen("20-10K.txt", "wb");
	//gsl_matrix_fprintf(f, X, "%.10e");
	//fclose(f);
	
	//MVEE instance1(X); 
	//instance1.set_doPrint(doPrint);
	//gsl_vector* u = gsl_vector_calloc(m);
	std::vector<double> timeEID;
	std::vector<double> optimalEID;
	std::vector<double> initialOptimalEID;
	
	for (int i=0; i<N; i++) {
		printf("------iteration %d------\n",i);
		//getRandomMatrixMVEE(X, r,ov,h);
		getRandomMatrix(X, r);
		orderMatrix(X, Y);
		gsl_matrix_memcpy(X, Y);
		clock_t begin=clock();
		MVEEPartial instance1(h,X, "HEURISTIC_EID"); 	
		instance1.tol=tol;
		instance1.solve();
		clock_t end=clock();
		timeEID.push_back((double)diffclock(end,begin));
		optimalEID.push_back(instance1.optimalValue);
		initialOptimalEID.push_back(instance1.initialOptimalValue);
		printf("optimal value EID: %e\n",optimalEID[i]);
		printf("time for BB: %e ms\n",timeEID[i]);
	}
	double  avgTimeEID = Myavg(timeEID);
	double avgOptimalEID =Myavg(optimalEID);
	double avgInitialOptimalEID = Myavg(initialOptimalEID); 
	
	std::cout<<"----------------------------------"<<std::endl;
	std::cout<< "n, m, h " << n << " " << m << " " <<h <<std::endl;
	printf("avg time EID %e\n",avgTimeEID);
	printf("avg optimal value EID %e\n",avgOptimalEID);
	printf("avg initial value EID %e\n",avgInitialOptimalEID);
	//printf("avg volume peeling %e\n",avgInitialOptimalPEELING);
	//printf("avg volume construction %e\n",avgOptimalConstruction);
	//printf("how often is best Heuristic the peeling? %e\n",freqBestHeuristicPeeling);
	//printf("how often is peeling better than the default ellipsoid? %e\n",sucPEELING/(double)N);
	//printf("how often is EID better than the default ellipsoid? %e\n",sucEID/(double)N);
	//printf("how often none of them is better than the default ellipsoid? %e\n",noSuc/(double)N);
	//instance1.calculateOptimalValue();
	//std::cout << "optimal Value is : " << instance1.optimalValue << std::endl;;
	//gsl_matrix* M = gsl_matrix_alloc(n, n);
	//gsl_vector* u = gsl_vector_alloc(m);
	//instance1.get_u(u);
	//printf("%e\n",sum(u));
	
	//gsl_matrix_memcpy(M, instance1.M);
	
	//if (instance1.isOptimal(tol))
	//{
	//	printf("success\n");
	//}
	//else {
	//	printf("failure\n");
	//}
	
	gsl_matrix_free(X); 
	gsl_matrix_free(Y);
	gsl_rng_free (r);    
	//gsl_vector_free(u);
	//std::ofstream myfile;
	//myfile.open ("example.txt");
	//myfile << "Writing this to a file.\n";
	//myfile.close();
	return GSL_SUCCESS;
}

int testMVEE()
{
	const gsl_rng_type* T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	int n = 50;
	int m = 100000;
	double tol = 0.0000001;
	bool doPrint = false;
	gsl_matrix* X = gsl_matrix_calloc(n,m);
	//gsl_matrix* Y = gsl_matrix_calloc(n, m);
	//FILE* f = fopen("problematic.dat", "rb");
	//gsl_matrix_fread(f, X);
	//fclose(f);
	//print(X);
	getRandomMatrix(X,r);
	//getRandomMatrix(X, r);
	//getRandomMatrix(X, r);
	//getRandomMatrix(X, r);
	
	//printMATLAB(X);
	//FILE* f = fopen("20-10K.txt", "wb");
	//gsl_matrix_fprintf(f, X, "%.10e");
	//fclose(f);
	
	MVEE instance1(X);
	clock_t begin=clock();
	instance1.solve();
	clock_t end = clock(); std::cout << "Time elapsed: " << double(diffclock(end,begin)) << " ms"<< std::endl;
	//if (instance1.isOptimal(tol))
	//{
	//	printf("success\n");
	//}
	
	gsl_rng_free(r);
	gsl_matrix_free(X);
	return GSL_SUCCESS;
	
	 
}

int testBadUposCase()
{
	int n = 2;
	int m = 4;
	double tol = 0.0000001;
	bool doPrint = false;
	gsl_matrix* X = gsl_matrix_calloc(n,m);
	gsl_matrix_set(X, 0, 0, -0.226119);
	gsl_matrix_set(X, 0, 1,-0.440694);
	gsl_matrix_set(X, 0, 2,0.948226 );
	gsl_matrix_set(X, 0, 3,-1.3212 );
	gsl_matrix_set(X, 1, 0,-0.0804288 );
	gsl_matrix_set(X, 1, 1,-0.211832 );
	gsl_matrix_set(X, 1, 2,0.520374 );
	gsl_matrix_set(X, 1, 3, -0.776972);
	MVEE instance1(X);
	clock_t begin=clock();
	instance1.solve();
	clock_t end = clock(); std::cout << "Time elapsed: " << double(diffclock(end,begin)) << " ms"<< std::endl;
	if (instance1.isOptimal(tol))
	{
		printf("success\n");
	}
	gsl_matrix_free(X);
	return GSL_SUCCESS;	
}

int main(int argc , char* argv[])
{
	//testMVEEPartial();
	//testBadUposCase();
	testMVEE();
	//testMVEEPartialEID();
	
}