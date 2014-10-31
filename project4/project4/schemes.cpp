#include "schemes.h"
using namespace arma;

Schemes::Schemes(double dx, double dt)
{
    alpha = dt/(dx*dx);
}

vec Schemes::Explicit(vec u)
{
    int n = u.n_elem;
    double a = alpha;
    double b = 1-2*alpha;
    vec y(n);
    for (int i=1; i<n-1; i++)
    {
        y(i) = u(i-1)*a + u(i)*b + u(i+1)*a;
    }
    y(0) = 0;     //Boundary condition
    y(n-1) = 0;  //Boundary condition

    return y;
}

vec Schemes::Implicit(vec y)
{
    int n = y.n_elem;
    // Au=y   , A is a matrix and u and y are vectors.
    // Elements of A are all constants so there is no need of a matrix.
    double a = -alpha;
    double b = 1+2*alpha;
    //Make one time iteration:
    vec u = Tridiag(y, a, b, a);    //Solves the tridiagonal matrix equation. Returning the u in Au=y.
    u(0) = 0;      //Boundary condition
    u(n-1) = 0;   //Boundary condition

    return u;
}

vec Schemes::Crank_Nicolson(vec u)
{
    // Setting up the output vector
    int n = u.n_elem;
    vec y(n);

    // Calculate v0' in v0'=(2I-aB)v0
    double a = alpha;
    double b = 2-2*alpha;
    for (int i=1; i<n-1; i++)
    {
        y(i) = u(i-1)*a + u(i)*b + u(i+1)*a;
    }
    y(0) = 0;     //Boundary conditions
    y(n-1) = 0;

    //Calculate v1 in (2I+aB)v1=v0'
    y = Tridiag(y, -alpha, 2+2*alpha, -alpha);
    y(0)=0;      //Boundary conditions
    y(n-1) =0;

    return y;
}


vec Schemes::Tridiag(vec x, double a, double b, double c)
{
    int N = x.n_elem;
    vec v = vec(N); // Making bunch of vectors...
    vec g = v;
    vec p = v;

    g(0) = c/b; //a(0) = 0
    p(0) = x(0)/b;

    for (int i=1; i<N-1; i++){    //Forward substitution.
        g(i) = c/(b-(a*g(i-1)));
        p(i) = (x(i)-a*p(i-1))/(b-a*g(i-1));
    }
    g(N-1) = 0; //c(N-1)=0
    p(N-1) = (x(N-1)-a*p(N-2))/(b-a*g(N-2));


    //v(0)=0; //Initial condition
    v(N-1) = p(N-1);
    for (int i = N-2; i>=0; i--){   //Backward substitution.
        v(i) = p(i)-g(i)*v(i+1);
    }
    return v;
}

double u(double x, double t)
{
    double L =10;
    double v = 0;
    for (int n=1; n<100; n++)
    {
        v+= (-2/(n*M_PI))*sin(M_PI*n*x/L)*exp(-pow((M_PI*n/L),2)*t);
    }
    return v;
}
