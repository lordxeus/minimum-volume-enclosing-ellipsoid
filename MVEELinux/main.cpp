//============================================================================
// Name        : MVEE.cpp
// Author      : Spyros
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "MVEE.h"
#include <cblas.h>

using namespace std;

int testMVEE()
{
        const gsl_rng_type* T;
        gsl_rng * r;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        int n = 10;
        int m = 1000;
        double tol = 0.0000001;
        gsl_matrix* X = gsl_matrix_calloc(n,m);
        getRandomMatrix(X,r);

        MVEE instance1(X);
        clock_t begin=clock();
        instance1.solve();
        clock_t end = clock(); std::cout << "Time elapsed: " << double(diffclock(end,begin)) << " ms"<< std::endl;
        if (instance1.isOptimal(tol))
        {
              printf("success\n");
        }

        gsl_rng_free(r);
        gsl_matrix_free(X);
        return GSL_SUCCESS;


}

int main(int argc , char* argv[])
{

        testMVEE();

}
