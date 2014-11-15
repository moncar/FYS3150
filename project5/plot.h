#ifndef PLOT_H
#define PLOT_H

#include <iostream>
#include<armadillo>
using namespace std;
using namespace arma;

class Plot
{
public:
    Plot();
    Plot(vec a, vec b);
    Plot(vec a, vec b, string x);
    void basicstructure(vec a,vec b);
    void Title(string T);
    void Labels(string A, string B);
    void Legend(string L);
    void Axis(string E);
    void Axis(double xmin, double xmax, double ymin, double ymax);
    void Figure(int f);
    void Show();

    double n;
};

#endif // PLOT_H
