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

vec Analytical(double t, double L , vec x)
{
    vec A = zeros(x.n_elem);
    for (int n=1; n<50; n++)
    {
        A +=  (2/(n*M_PI)) * (cos(n*M_PI)*(1-L) - 1) * sin(M_PI*n*x/L) * exp(-pow((M_PI*n/L),2)*t);
    }
    return A;
}

int main()
{

    int N = 10;
    double dx = 1./N;
    double dt = 0.5*dx*dx; //Limit of the explicit schemes stability.
    double L = 1;  //Can be chosen freely.
    vec x = linspace(0,L,N);
    vec v = x/L -1;     //Initial condition/state
    vec Anal;   //haha...
    vec CN, E, I, CNu, Eu, Iu;
    Schemes One(dx, dt);    //Just some object. Random name.
    int c = 0;
    E = One.Explicit(v);
    I = One.Implicit(v);
    CN = One.Crank_Nicolson(v);
    for (double t=0; t<0.5;t+=dt)
    {
        E = One.Explicit(E);
        I = One.Implicit(I);
        CN = One.Crank_Nicolson(CN);
        //Transforming back v→u
        Eu = E+(1-x/L);
        Iu = I+(1-x/L);
        CNu = CN +(1-x/L);
        fstream plotdata4c;
        char str[100];
        c+=1;
        if (c == 10 || c== 20|| c==30)  //Writes some cases to files. Used to make plots.
        {
            Anal = Analytical(t,L,x) + (1-x/L);   // Calculates the analytical solution.
            sprintf(str, "/home/filiphl/Desktop/FYS3150/project4/project4/plotdata/plotdata4c%d.txt", c);
            plotdata4c.open(str, ios::out);
            plotdata4c << setw(20) <<setprecision(10)<<"x"
                       << setw(20) <<setprecision(10)<< "Explicit"
                       << setw(20) <<setprecision(10)<< "Implicit"
                       << setw(20) <<setprecision(10)<< "Crank Nicolson"
                       << setw(20) <<setprecision(10)<< "Analytical"
                       << setw(20) <<setprecision(10)<< "Diff. Explicit"
                       << setw(20) <<setprecision(10)<< "Diff. Implicit"
                       << setw(20) <<setprecision(10)<< "Diff. Crank Nicolson"<<endl;
            for (int i=0; i<N; i++)
            {
                plotdata4c << setw(20) <<setprecision(10)<<x(i)
                           << setw(20) <<setprecision(10)<<Eu(i)
                           << setw(20) <<setprecision(10)<<Iu(i)
                           << setw(20) <<setprecision(10)<<CNu(i)
                           << setw(20) <<setprecision(10)<<Anal(i)
                           << setw(20) <<setprecision(10)<<abs( (Anal(i)-Eu(i) ) / Anal(i))
                           << setw(20) <<setprecision(10)<<abs( (Anal(i)-Iu(i) ) / Anal(i))
                           << setw(20) <<setprecision(10)<<abs( (Anal(i)-CNu(i) ) / Anal(i)) <<endl;
            }
            plotdata4c.close();
        }
    }

    return 0;
}

