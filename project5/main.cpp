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
#include </home/filiphl/Desktop/FYS3150/project5/cppLibrary/lib.h>




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


/*
    // Monte Carlo a)

    vec r;
    int nbins = 100;
    int npart = 10000;
    vec part = zeros<vec>(nbins);
    vec new_M;
    vec a = part;
    part(0) = npart;
    fstream hist1;
    hist1.open("/home/filiphl/Desktop/figs/hist1.txt", ios::out);

    for (int t = 0; t<2000;t++)
    {
        for (int k = 0; k<nbins; k++)
        {
            hist1 << setw(10) <<  part(k)/npart;
        }
        hist1<<endl;
        new_M = part;
        for (int i = 0; i<nbins-1; i++)
        {
            r = randu<vec>(part(i));

            for (int j=0; j<part(i); j++)
            {
                if (r(j) > 0.5)
                {
                    new_M(i+1)+= 1;
                }
                else
                {
                    if (i != 0)
                    {
                        new_M(i-1)+= 1;
                    }
                }
                new_M(i)-= 1;
            }
        }
        part = new_M;
        part(0) = npart;
        part(nbins - 1) = 0;
    }
    hist1.close();

    system("python /home/filiphl/Desktop/figs/hist.py 1999");
    system("rm /home/filiphl/Desktop/figs/hist1.txt");
    //Plot p(linspace(0,1,nbins), part);
    //p.Show();

    */



    //********Monte Carlo b)*********//
    /*Since there is no longer an integer times l_0 steplength,
    we can no longer use the configuration above. We can no longer
    use a vector telling how many part are in the one x-position,
    because most likely none of the part are at the same place
    at all! Except off course at the boundary x=0. Instead, we let
    each vector element represent the number of part in a
    region of the space. This does demand more if-tests, but this
    is much more efficient than having a vector containing the
    position of each particle and expand and erase elements as
    part enter and exit...
    */



    int nbins = 100;
    int npart = 10000;
    mat new_M = zeros<mat>(nbins, nbins);
    mat M(nbins, nbins);
    vec rx, ry;

    // Set up initial state
    for (int i = 0; i<nbins; i++)
    {
        M(i,0) = npart/nbins;
    }


    //mat1.open("/home/filiphl/Desktop/figs/mat1.txt", ios::out);

    int m;
    int c=0;
    double prosess;
    double timelimit = 5000;
    // Start iterations
    for (int t = 0; t<timelimit; t++)
    {
        new_M = zeros<mat>(nbins, nbins);

        for (int i = 0; i<nbins; i++)
        {
            for (int j=0; j<nbins; j++)
            {
                m = M(i, j);
                if (m!=0)
                {
                    rx = randu<vec>(m);
                    ry = randu<vec>(m);
                    for ( int n=0; n<m; n++ )   // Loop over all particles in the position (i,j)
                    {
                        if ( (rx(n) > 0.5) && (ry(n) > 0.5) && (j!=nbins-1) )
                        {
                            if (i==nbins-1) { new_M(0, j+1) += 1;}
                            else {new_M(i+1,j+1) += 1;}
                        }

                        if ( (rx(n) > 0.5) && (ry(n) < 0.5) && (j!=nbins-1) )
                        {
                            if (i==0) { new_M(nbins-1, j+1) += 1;}
                            else {new_M(i-1, j+1) += 1;}
                        }

                        if ( (rx(n) < 0.5) && (ry(n) > 0.5) && (j!=0) )
                        {
                            if (i==nbins-1) { new_M(0, j-1) += 1;}
                            else {new_M(i+1, j-1) += 1;}
                        }

                        if ( (rx(n) < 0.5) && (ry(n) < 0.5) && (j!=0) )
                        {
                            if (i==0) { new_M(nbins-1, j-1) += 1;}
                            else {new_M(i-1, j-1) += 1;}
                        }
                    }
                }
            }
        }
        M = new_M;
        for (int i=0; i<nbins; i++)
        {
            M(i, 0) = npart/nbins;
            M(i, nbins-1) = 0;
        }
        if (c==100)
        {
            fstream mat1;
            char str[100];
            sprintf(str, "/home/filiphl/Desktop/figs/mat%d.txt", t);
            mat1.open(str, ios::out);


            for (int i=0; i<nbins; i++)
            {
                for (int j=0; j<nbins; j++)
                {
                    mat1 <<  setw(20) << M(i,j);   //writes the matrix.
                }
                mat1 << endl;
            }
            c=0;
            prosess += 100;
            cout << 100*prosess/timelimit<<"% finished"<<endl;
        }
        c+=1;


    }





/*


            {
                rx = randu<vec>(partx(i));
                for (int j=0; j<partx(i); j++)
                {
                    if ( (rx(j) > 0.5) && (i!=nbins-1) )
                    {
                        new_Mx(i+1)+= 1;
                    }
                    else
                    {
                        if ( (i != 0) && (i!=nbins-1) )
                        {
                            new_Mx(i-1)+= 1;
                        }
                    }
                }
            }

            if (party(i) > 0)
            {
                ry = randu<vec>(party(i));
                for (int k = 0; k<party(i); k++)
                {
                    if (ry(k) > 0.5)
                    {

                        if (i == nbins-1)
                       {
                           new_My(0)+= 1;
                       }

                        else
                        {
                            new_My(i+1)+= 1;
                        }
                    }
                    else
                    {
                        if (i == 0)
                        {
                            new_My(nbins-1)+= 1;
                        }
                        else
                        {
                            new_My(i-1)+= 1;
                        }
                    }
                }
            }
        }

        partx = new_Mx;
        partx(0) = npart;
        partx(nbins - 1) = 0;
        party = new_My;

        if (c == 50)
        {
            fstream mat1;
            char str[100];
            sprintf(str, "/home/filiphl/Desktop/figs/mat%d.txt", t);
            mat1.open(str, ios::out);


            for (int i =0; i<nbins; i++)
            {
                for (int j=0; j<nbins; j++)
                {
                    mat1 <<  setw(20) << partx(i)*party(j)/(npart*npart);   //writes a normalized matrix.
                }
                mat1 << endl;
            }
            c=0;
        }
        c+=1;
    }


*/







/*
    int N = 5;
    double dt = 0.004;
    double dx = 0.2;
    vec x = linspace(0,1,N);
    vec y = linspace(0,1,N);
    mat U=zeros<mat>(N,N);
    U.col(0) = 1-x;
    U.col(N-1)=1-x;
    U.row(0) = ones<mat>(1,N);
    mat V0(N,N);
    for (int i = 0; i<N; i++)
    {
        V0.col(i)=U.col(i)+x-1;
    }

    Schemes2d Two(dx,dt);
    mat E, EU = ones<mat>(N,N);
    E = Two.Explicit2d(V0);

    Two.Implicit2d(V0);
*/

/*
    double t = 0;
    double timelimit = 2;
    int c = 0;
    int inst = 10;

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
    }
*/







/*
        A = Analytical(t,x,y);   // Calculates the analytical solution.
        char str[50];
        sprintf(str, "/home/filiphl/Desktop/figs/plot%d.txt", c);
        fstream plot1;
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


