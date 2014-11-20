#include "schemes2d.h"
using namespace arma;

Schemes2d::Schemes2d(double dx, double dt)
{
    alpha = dt/(dx*dx);
}

mat Schemes2d::Explicit2d(mat U)
{
    int n = U.n_cols;   //Assuming nxn-matrix
    mat V = zeros<mat>(n,n);
    for (int i=1; i<n-1; i++)
    {
        for (int j=1; j<n-1; j++)
        {
            V(i,j) = U(i,j) + alpha*(  U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) - 4*U(i,j)  );
        }
    }
    //The boundary conditions are conserved as they are not changed at all.
    return V;
}

mat Schemes2d::Implicit2d(mat U)
{
    int n = U.n_cols;   //Assuming nxn-matrix
    int a = n-2;
    int N = a*a;
    mat V = zeros<mat>(N,N);
    vec u = zeros<vec>(N);
    vec x = u;
    for (int i=2; i<=a; i++)
    {
        for (int j=2; j<=a; j++)
        {

            V(i, (i-1)*(n-2)+j) = 1;
        /*
            if (i == 1 ||  i== n-1)
            {

            }
            */
        }
    }
    cout <<V<<endl;
    
    
    
    return V;
}
