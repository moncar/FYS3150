#ifndef SCHEMES2D_H
#define SCHEMES2D_H

#include <armadillo>

//using namespace arma;
class Schemes2d
{
public:
    Schemes2d(){}
    Schemes2d(double dx, double dt);
    arma::mat Explicit2d(arma::mat U);
    arma::mat Implicit2d(arma::mat U);
    arma::vec Tridiag(arma::vec x, double a, double b, double c);
    double alpha;
    double beta;
};
#endif // SCHEMES2D_H
