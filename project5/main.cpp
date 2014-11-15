#include <iostream>
#include<armadillo>
#include <stdlib.h> /* atof function */
#include <math.h>   /* sine function */
#include <stdio.h>  /* printf function */
#include <fstream>
#include <iomanip>
#include <plot.h>
//#include <schemes.h>
#include <schemes2d.h>

using namespace std;
using namespace arma;

mat Analytical(double t , vec x, vec y)
{
    int N = x.n_elem;
    double Fnm;
    mat Z = zeros<mat>(N,N);
            for (int i=0;i<N;i++)
            {
                for (int j=0; j<N;j++)
                {
                    for (int n=1; n<50; n++)
                    {
                        for (int m=1; m<50; m++)
                        {
                             Fnm = (4/(M_PI*M_PI*n*m)) * (cos(M_PI*m) - 1);
                             Z(j,i) +=  Fnm * sin(n*M_PI*x(i)) * sin(m*M_PI*y(j)) * exp(-(n*n*M_PI*M_PI + m*m*M_PI*M_PI )*t) ;
                        }
                    }
                    Z(j,i) += 1-x(i);
                 }
            }
    return Z;
}

int main()
{
    vec r;
    int nr_of_bins = 100;
    int nr_of_particles = 10000;
    vec particles = zeros<vec>(nr_of_bins);
    vec new_dist;
    vec a = particles;
    particles(0) = nr_of_particles;
    fstream hist1;
    hist1.open("/home/filiphl/Desktop/figs/hist1.txt", ios::out);

    for (int t = 0; t<1000;t++)
    {
        for (int k = 0; k<nr_of_bins; k++)
        {
            hist1 << setw(10) <<  particles(k)/nr_of_particles;
        }
        hist1<<endl;
        new_dist = particles;
        for (int i = 0; i<nr_of_bins-1; i++)
        {
            r = randu<vec>(particles(i));

            for (int j=0; j<particles(i); j++)
            {
                if (r(j) > 0.5)
                {
                    new_dist(i+1)+= 1;
                }
                else
                {
                    if (i != 0)
                    {
                        new_dist(i-1)+= 1;
                    }
                }
                new_dist(i)-= 1;
            }
        }
        particles = new_dist;
        particles(0) = nr_of_particles;
        particles(nr_of_bins - 1) = 0;
    }
    hist1.close();

    system("python /home/filiphl/Desktop/figs/hist.py 999");
    system("rm /home/filiphl/Desktop/figs/hist1.txt");
    //Plot p(linspace(0,1,nr_of_bins), particles);
    //p.Show();



    /*
    int N = 21;
    double dt = 0.004;
    double dx = 0.2;
    vec x = linspace(0,1,N);
    vec y = linspace(0,1,N);
    mat A;
    mat V0=zeros<mat>(N,N);
    for (int i = 1; i<N-1; i++)
    {
        V0.col(i) = x-1;
    }
    V0.row(0) = zeros<vec>(N).t();
    V0.row(N-1) = zeros<vec>(N).t();


    Schemes2d Two(dx,dt);
    mat E, EU = ones<mat>(N,N);
    E = Two.Explicit2d(V0);

    double t = 0;
    double timelimit = 2;
    int c = 0;
    int inst = 0;

    while ( t<timelimit)
    {
        fstream plot1;
        char str[100];
        E = Two.Explicit2d(E);

        for (int i=0; i<N; i++)
        {
            EU.col(i) = E.col(i) + 1-x;
        }

        if (inst == 10)
        {
            sprintf(str, "/home/filiphl/Desktop/figs/plotE%d.txt", c);
            plot1.open(str, ios::out);
            for (int i=0; i<N; i++)
            {
                for (int j=0; j<N; j++)
                {
                    plot1 << setw(20) <<setprecision(10)<<EU(i,j);
                }
                plot1 << endl;
            }
            inst = 0;
        }


        inst +=1;
        c+=1;
        t+=dt;

        /*A = Analytical(t,x,y);   // Calculates the analytical solution.

        sprintf(str, "/home/filiphl/Desktop/figs/plot%d.txt", c);
        plot1.open(str, ios::out);
        for (int i=0; i<N; i++)
        {
            for (int j=0; j<N; j++)
            {
                plot1 << setw(20) <<setprecision(10)<<A(i,j);
            }
            plot1 << endl;
        }
        */


    return 0;
}


