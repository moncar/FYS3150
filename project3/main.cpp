#include <iostream>
#include <iostream>
#include<armadillo>
#include <stdlib.h> /* atof function */
#include <math.h>   /* sine function */
#include <stdio.h>  /* printf function */
#include <fstream>
#include <iomanip>

using namespace std;
using namespace arma;


class System
{
private:
    int dt;       //Timestep length.
    vec A;      // vector containing all positions and velocities.
public:
    System(){}
    vec forces(Planet a, Planet b);
    int n;        //Number of objects.
    Planet planets[10];     //List of planets.

};
vec System::forces(Planet a, Planet b)
{
    vec F(4);
    double r = sqrt(pow((a.position(0) - b.position(0)),2) + pow((a.position(1)-b.position(1)),2));

}


class Planet
{
private:
    double m_sun = 2e30;
public:
    Planet() {}
    Planet(double x0, double y0, double vx0, double vy0, double m0);
    vec position;
    vec velocity;
    double mass;
    System::planets[n] = planet;
    System::n += 1;


    friend ostream& operator<< (ostream &out, Planet &planet);
};

Planet::Planet(double x0, double y0, double vx0, double vy0, double m0)    //Setter mulighet for syntax planet(p,v,m)
{
    position = vec(2);
    velocity = vec(2);
    position(0) = x0;
    position(1) = y0;
    velocity(0) = vx0;
    velocity(1) = vy0;
    mass = m0/m_sun;
}

//How cout<< will work with this class.
ostream& operator<< (ostream &out, Planet &planet)
{
    // Since operator<< is a friend of the Point class, we can access
    // Point's members directly.
    out << "(" << planet.position(0) << ", " << planet.position(1) << ", " <<
        planet.velocity(0) << ", " << planet.velocity(1)  << ", " <<
        planet.mass << ")"<<endl;
    return out;
}


int main()
{
//    vec A(2);
//    A(0) = 1;
//    A(1) = 2;
//    cout << A<<endl;
    Planet earth(1.0, 0.0, 0.0, 4.0, 6e24);
    Planet mars(1.0,  2.0, 3.0, 4.0 ,5.0);
    Planet* planets = new Planet[2];

    planets[0] = earth;

    cout << planets[0] <<  endl;
    cout <<earth<<endl;
    //Planet mars(double 1, double 2, double 3);
    //mars.bigger(20);
    //cout << mars<<endl;
/*
    cout<<"pos: "<<mars.pos<<endl;
    cout<<"vel: "<<mars.vel<<endl;
    cout<<"m: "<<mars.m<<endl;
*/
    return 0;
}

