#include <iostream>
#include<armadillo>
#include <stdlib.h> /* atof function */
#include <math.h>   /* sine function */
#include <stdio.h>  /* printf function */
#include <fstream>
#include <iomanip>
#include <plot.h>
#include <schemes.h>
using namespace std;
using namespace arma;


int main()
{

    int N = 10;
    double dx = 1./N;
    double dt = 0.001;
    int L = 10;  //Can be chosen freely.
    vec x = linspace(0,L,N);
    vec v = x/L -1;     //Initial condition/state

    Schemes One(dx, dt);
    int c = 0;
    for (double t=0; t<1;t+=dt)
    {
        v = One.Crank_Nicolson(v);
        c+=1;
        if (c==100)
        {
            Plot t(x,v,"r");
            t.Labels("x", "u(x)");
            t.Axis(0,L,-1.1,0.01);
            t.Show();
            c=0;
        }
    }



/*
    vec v = zeros<vec>(N);
    for (int n=1; n<10000; n++)
    {
        v+= (-2/(n*M_PI))*sin(M_PI*n*x/L);
    }

    Plot plt(x,v, "r");
    plt.Legend("N = 100");
    plt.Labels("x", "u(x)");
    plt.Title("Fourier-tilnaerming");
    plt.Show();
*/











    return 0;
}

