#include "schemes2d.h"
using namespace arma;

Schemes2d::Schemes2d(double dx, double dt)
{
    alpha = dt/(dx*dx);
    beta = 1/(1+4*alpha);
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
    int n = U.n_cols;    //Assuming nxn-matrix
    mat V1 = U;          // Past iteration
    mat V2;                 // Precent iteration
    for (int k=0; k<30; k++)
    {
        V2 = zeros<mat>(n,n);
        for (int i=1; i<n-1; i++)
        {
            for (int j=1; j<n-1; j++)
            {
                V2(i,j) = beta*( alpha*( V1(i+1,j) + V1(i-1,j) + V1(i,j+1) + V1(i,j-1) ) + U(i,j)  );
            }
        }
        V1 = V2;
    }
    return V2;
}
