#include <iostream>
#include<armadillo>
#include <stdlib.h> /* atof function */
#include <math.h>   /* sine function */
#include <stdio.h>  /* printf function */
#include <fstream>
#include <iomanip>
#include <time.h>


using namespace std;
using namespace arma;

double f(double x){             //The function f(x)
    return 100*exp(-10*x);
}

int main()
{
    int N;

    cout<< "Enter a value for N: ";
    cin >> N;
    double h = 1./(N+1);

    mat A = zeros<mat>(N,N);

    for (int i=0; i<N; i++){  //Makes the b-vector. I call it x.
        for (int j=0; j<N; j++){
            if (i ==j){A(i,j)=2;}
            if (i == j+1 || i == j-1){A(i,j)=-1;}
        }
    }

    vec b = vec(N);
    for (int i=0; i<N; i++){  //Makes the b-vector. I call it x.
        b(i) = pow(h, 2) * f(i*h);
    }

    clock_t start, finish; // declare start and final time
    start = clock();

    mat L, U;
    vec x = vec(N);
    vec y = vec(N);
    lu(L, U, A);
    solve(y,L,b);
    solve(x, U, y);

    finish = clock();
    cout<< "Time:"<< double(finish-start)/(double)CLOCKS_PER_SEC <<endl;


    return 0;
}

