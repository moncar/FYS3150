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
#include "cppLibrary/lib.h"
#include "random.h"
#include <algorithm>



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



    // Monte Carlo a)

    /*
    int L = 1;
    double dt = 0.01;
    double l0= sqrt(2*dt);
    int npart = 10000;
    long int idum = -1;
    double r;
    vector<double>dpos(npart);
    vector<double>new_dpos;
    double pastdpos;
    double nextdpos;

    fstream hista;
    hista.open("/home/filiphl/Desktop/figs/hist_a.txt", ios::out);

    for (double t = 0; t<1; t+=dt)
    {
        cout << t<<endl;
        //Write to file
        for (int k = 0; k<dpos.size(); k++)   {hista << setw(15) <<  dpos[k];}
        hista<<endl;

        for (int i=0; i<dpos.size(); i++)
        {
                pastdpos = dpos[i];
                r = ran0(&idum);
                if (r > 0.5){nextdpos = pastdpos + l0;}
                else {nextdpos = pastdpos - l0;}
                if ((nextdpos>0)&&(nextdpos<L)) {new_dpos.push_back(nextdpos);}
                if (pastdpos==0) {new_dpos.push_back(0);
            }
        }
        dpos = new_dpos;
        new_dpos.clear();
    }
    hista.close();
    //Run python script for histogram.
    system("python /home/filiphl/Desktop/figs/hist.py hist_a.txt 90");

    */






    //********Monte Carlo b)*********//
    /* Since there is no longer an integer times l_0 steplength,
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
    /*
    int L = 1;
    double dt = 0.00004;
    double l0= sqrt(2*dt);
    double dx = 0.01;
    Random r(-1);
    int npart = 10000;

    vector<double> pos;
    vector<double> updatedpos(npart); // FIlled with zeroes as default.
    double sd =1/sqrt(2); //standard deviation
    double pastpos;
    double nextpos;

    fstream hist1;
    hist1.open("/home/filiphl/Desktop/figs/hist1.txt", ios::out);

    int c = 0;
    int cc = 0;
    for (double t=0; t<1; t+=dt)
    {

        //cout<<t<<endl;

        if (c==1000)
        {
            cout << 100*t <<"% done"<<endl;
            for (int k = 0; k<pos.size(); k++)  //Writes data to file.
            {
                hist1 << setw(25) << setprecision(8) <<  pos[k];
            }
            hist1 << endl;
            c = 0;
        }
        for (int i=0; i<pos.size(); i++)
        {
            pastpos = pos[i];
            pos[i] += r.nextGauss(0.0, sd)*l0;
            nextpos = pos[i];

            //Includes relevant values.
            if ( (nextpos>=0)&&(nextpos<L) ) { updatedpos.push_back(nextpos); }

            // Partical moving from inside the dx-interval to outside
            if ( (0<=pastpos) && (pastpos<dx) )
            {
                if ( (nextpos<0) || (nextpos>dx) ) {cc--;}
            }

            //Partical moving from outside the dx-interval to inside
            if ( (pastpos>dx) && (nextpos<dx) && (nextpos>0) ) {cc++;}
        }

        while (cc < 0)
        {
            updatedpos.push_back(dx/2); // Notice the position!
            cc++;
        }
        pos = updatedpos;
        updatedpos.clear();
        c++;
    }


    system("python /home/filiphl/Desktop/figs/hist.py hist1.txt 24");
*/


    /*
                // *********monte carlo 2D************

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

                */



    // Monte Carlo 2D

    int L = 1;
    double dt = 0.001;
    double l0= sqrt(2*dt);
    double dx = 0.05;   // dy=dx
    Random r(-1);
    int npart = 100000;

    vector<float> xpos;  // FIlled with zeroes as default.
    vector<float> ypos;
    for (float j=0; j<npart; j++) { xpos.push_back(dx/2); ypos.push_back(j/npart); }
    vector<float> new_xpos;
    vector<float>new_ypos;
    float sd =1/sqrt(2); //standard deviation
    float pastxpos, pastypos, nextxpos, nextypos;
    long int idum = -1;
    int nbins = L/dx-1;
    mat M;
    int c = 1;
    int a = 0;
    int cc = 0;
    float y;
    for (float t=0; t<0.25; t+=dt)
    {
        if (c==1)
        {
            cout << 100*t/0.25 <<"% done"<<endl;

            M = zeros<mat>(nbins, nbins);
            for (int i=0; i<xpos.size(); i++)
            {
                for (int row=0; row<nbins; row++)
                {
                    for (int col=0; col<nbins; col++)
                    {
                        if ( (xpos[i]>=row*dx) && (xpos[i]<(row+1)*dx) )
                        {
                            if ( (ypos[i]>=col*dx) && (ypos[i]<(col+1)*dx) )  {M(col,row)+=1;}
                        }
                    }
                }
            }
            fstream mc2;
            char str[100];
            sprintf(str, "/home/filiphl/Desktop/figs/mc2d%d.txt", a);
            mc2.open(str, ios::out);

            for (int row=0; row<nbins; row++)
            {
                for (int col=0; col<nbins; col++)
                {
                    mc2<< setw(15)<<M(row,col);
                }
                mc2<<endl;
            }
            a+=c;
            c = 0;
        }

        for (int i=0; i<xpos.size(); i++)   //xpos and ypos should have equal size.
        {
            pastxpos = xpos[i];
            pastypos = ypos[i];
            xpos[i] += r.nextGauss(0.0, sd)*l0;
            ypos[i] += r.nextGauss(0.0, sd)*l0;
            nextxpos = xpos[i];
            nextypos = ypos[i];


            //Includes relevant values.
            if ( (nextxpos>0)&&(nextxpos<L) )
            {
                new_xpos.push_back(nextxpos);
                if ( (nextypos>=0)&&(nextypos<=L) ) { new_ypos.push_back(nextypos); }
                // Particles moving for inside y=[0,1] to outside: periodical boundary conditions.
                else if ( (pastypos<=L) && (nextypos>L) ) { new_ypos.push_back(nextypos-L); }
                else if ( (pastypos>=0) && (nextypos<0) ) { new_ypos.push_back(nextypos+L); }
            }

            // Partical moving from inside the dx-interval to outside
            if ( (0<=pastxpos) && (pastxpos<dx) )
            {
                if ( (nextxpos<0) || (nextxpos>dx) ) {cc--;}
            }

            //Partical moving from outside the dx-interval to inside
            if ( (pastxpos>dx) && (nextxpos<dx) && (nextxpos>0) ) {cc++; }

        }



        // Add new elements at x = dx/2 and y = [0,1]
        y = cc;
        while (cc < 0)
        {
            new_xpos.push_back(dx/2.0);
            new_ypos.push_back(cc/y);
            cc++;
        }

        xpos = new_xpos;
        ypos = new_ypos;
        new_xpos.clear();
        new_ypos.clear();

        c++;
    }




































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


