#ifndef SCHEMES_H
#define SCHEMES_H
#include <armadillo>

//using namespace arma;
class Schemes
{
public:
    Schemes(){}
    Schemes(double dx, double dt);
    arma::vec Explicit(arma::vec u);
    arma::vec Implicit(arma::vec y);
    arma::vec Crank_Nicolson(arma::vec u);
    arma::vec Tridiag(arma::vec x, double a, double b, double c);
    double alpha;
};

#endif // SCHEMES_H
