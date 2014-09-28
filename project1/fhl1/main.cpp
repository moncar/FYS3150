#include <iostream>
#include<armadillo>
#include <stdlib.h> /* atof function */
#include <math.h>   /* sine function */
#include <stdio.h>  /* printf function */
#include <fstream>
#include <iomanip>

using namespace std;
using namespace arma;

double f(double x){             //The function f(x)
    return 100*exp(-10*x);
}

double eps(double u, double v){ //Function for calculating relative error
    return log10(abs((v-u)/u));
}

int main()
{
    int N;
    double h;
    cout<< "Enter a value for N: ";
    cin >> N;
    h = 1/double(N+1);
    vec v = vec(N); // Making bunch of vectors...
    vec x = v;
    vec a = v;
    vec b = v;
    vec c = v;
    vec g = v;
    vec p = v;
    vec e = v;
    vec u = v;

    for (int i=0; i<N; i++){  //Makes the b-vector. I call it x.
        x(i) = pow(h, 2) * f(i*h);
        a(i) = -1;
        b(i) = 2;
        c(i) = -1;
        }

    a(0) = 0;   //These values are not in the matrix.
    c(N-1) = 0;

    clock_t start, finish; // declare start and final time
    start = clock();       // starts the timer.

    g(0) = c(1)/b(1);
    p(0) = x(1)/b(1);
    //int FLOPS = 0;
    for (int i=1; i<N; i++){    //Forward substitution.
        g(i) = c(i)/(b(i)-(a(i)*g(i-1)));
        p(i) = (x(i)-a(i)*p(i-1))/(b(i)-a(i)*g(i-1));
    //FLOPS += 8;
    }

    v(0)=0; //Initial condition
    v(N-1) = p(N-1);
    for (int i = N-2; i>=1; i--){   //Backward substitution.
        v(i) = p(i)-g(i)*v(i+1);
    //    FLOPS += 2;
    }

    finish = clock();
    cout<< "Time:"<< double(finish-start)/(double)CLOCKS_PER_SEC <<endl;

    //cout<<"#FLOPS: "<<FLOPS<<endl;

    //Writes to file in order to plot using python.
    ofstream ofile;
    ofile.open("/home/filiphl/Desktop/v-er.txt");
    for (int n=0;n<N;n++){
        ofile<<v(n)<<endl;
    }
    ofile<<N<<endl;
    ofile.close();
    system("python /home/filiphl/Desktop/fhl_p2.py");    //Starts the python program.

    for (int i=0; i<N; i++){    //calculates the relative error.
        e(i) = eps(1-(1-exp(-10))*(i*h)-exp(-10*(i*h)),v(i));
    }

    cout<<"Maximum relativ error: "<<max(e)<<endl;


    return 0;
}

