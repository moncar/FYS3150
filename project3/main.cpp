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

/*------------------------------------------------------------------------------------*/

class Planet
{
private:
    double m_sun = 2e30;
public:
    Planet() {}
    Planet(string nm, double x0, double y0, double vx0, double vy0, double m);
    string name;
    vec position;
    vec velocity;
    double mass;



    friend ostream& operator<< (ostream &out, Planet &planet);
};

Planet::Planet(string nm, double x0, double y0, double vx0, double vy0, double m)    //Setter mulighet for syntax planet(p,v,m)
{
    position = vec(2);
    velocity = vec(2);
    position(0) = x0;
    position(1) = y0;
    velocity(0) = vx0;
    velocity(1) = vy0;
    mass = m/m_sun;
    name = nm;
}


//How cout<< will work with Planet class.
ostream& operator<< (ostream &out, Planet &planet)
{
    // Since operator<< is a friend of the Point class, we can access
    // Point's members directly.
    out << setw(8) <<"name: " << setw(10) <<  planet.name <<endl;
    out << setw(8) <<"Position: " << setw(4) <<planet.position(0) << ", " << planet.position(1) <<endl;
    out << setw(8) <<"Velocity: " << setw(4) << planet.velocity(0) << ", " << planet.velocity(1)  << endl;
    out << setw(8) <<"mass: " << setw(10) <<planet.mass<<endl;
    return out;
}


/*-----------------------------------------------------------------------------------*/


class System
{
private:
    double G =4*M_PI*M_PI; //39.478; (4*pi^2 AU^3 * M_sun^-1 * yr^-2)
    double dt;       //Time-step in years.
    vec P;      // vector containing all positions.
    vec V;      // vector containing all velocities.
public:
    System(){}
    int n = 0;        //Number of objects.
    Planet planets[10];     //List of objects.
    vec Acceleration(vec A);    //Calculates acceleration of  every object.
    void Addplanet(string nm, double x0, double y0, double vx0, double vy0, double m);
    void Setup();
    void Solve(double);
    vec Centerofmass(vec C);
    void Euler();
    void RK4();
    void Verlet();
    vec P2;
    vec nextP;
    vec nextV;

};

void System::Addplanet(string nm, double x0, double y0, double vx0, double vy0, double m)
{
    Planet x(nm, x0,  y0,  vx0,  vy0,  m);
    planets[n]  = x;
    n += 1;
}

void System::Setup()
{
    P = vec(n*2);   // Saves space for position of all objects.
    V = vec(n*2);   // ------------//--------- velocities.
    for (int i = 0; i<n; i++)
    {   // This algo stores position and velocity data of every planet, one at a time.
        P(2*i) = planets[i].position(0);
        V(2*i) = planets[i].velocity(0);
        P(2*i + 1) = planets[i].position(1);
        V(2*i + 1) = planets[i].velocity(1);
    }
}


vec System::Acceleration(vec A)
{
    vec B = zeros<vec>(2*n);
    for (int i=0; i<n; i++)
    {
        for (int j=i+1; j<n; j++)
        {

            double dx = A(2*j) - A(2*i);
            double dy = A(2*j+1) - A(2*i+1);
            double dr = sqrt(dx*dx + dy*dy);
            double F_tot = G*planets[i].mass*planets[j].mass/(dr*dr);
            double Fx = F_tot*(dx/dr);  //F_tot*cos(angle)
            double Fy = F_tot*(dy/dr);  //F_tot*sin(angle)
            B(2*i) +=  Fx/planets[i].mass;
            B(2*i+1) += Fy/planets[i].mass;
            B(2*j) += -Fx/planets[j].mass;
            B(2*j+1) += -Fy/planets[j].mass;
        }
    }
    return B;
}

vec System::Centerofmass(vec C)
{
    vec CM = zeros<vec>(2);
    double M=0;
    for (int i=0;i<n;i++)
    {
        M+= planets[i].mass;
        CM(1) += planets[i].mass*C(2*i);
        CM(2) += planets[i].mass*C(2*i+1);
    }
    CM = CM/M;
    return CM;
}

void System::Euler()    //Advancing one time-step.
{

    V += Acceleration(P)*dt;
    P += V*dt;
}

void System::RK4()  //Advancing one time-step.
{
    vec k1p(2*n), k2p(2*n), k3p(2*n), k4p(2*n), k1v(2*n), k2v(2*n), k3v(2*n), k4v(2*n); //Declaring vectors
    k1p = V * dt;
    k1v = Acceleration(P) * dt;
    k2p = (V + 0.5 * k1v) * dt;
    k2v = Acceleration(P + 0.5 * k1p) * dt;
    k3p = (V+ 0.5 * k2v) * dt;
    k3v = Acceleration(P+ 0.5 * k2p) * dt;
    k4p = (V + k3v ) * dt;
    k4v = Acceleration(P + k3p) * dt;

    P += (1/6.)*(k1p + 2 * (k2p+k3p) + k4p);    //Updating position-vector
    V += (1/6.)*(k1v + 2 * (k2v+k3v) + k4v);    //Updating velocity-vector
}

void System::Verlet()
{
    P2 = nextP;
    V = nextV;
    nextP = 2*nextP - P + Acceleration(nextP)*dt*dt;
    nextV = (nextP - P)/(2*dt);
    P = P2;
}


void System::Solve(double time_period)
{
    dt = 1/(365.0*24.); //Time-step one Hour.
    double t = 0;          //Set initial time to zero years.
    int Hour = 0;          //Counts itterations.

    nextP = P+V*dt + 0.5*Acceleration(P)*dt*dt; // Needed for verlet.
    nextV = V + (Acceleration(P) + Acceleration(nextP))*dt / 2.; // Needed for verlet.

    fstream planetpos, planetvel;
    planetpos.open("/home/filiphl/Desktop/project3/PlanetData.txt", ios::out);
    planetvel.open("/home/filiphl/Desktop/project3/PlanetVelocites.txt", ios::out);
    planetpos <<"Number_of_objects:"<<"       " << n<<endl; // I use this in python...
    planetpos <<"Hours";
    planetvel << "Hours";
    for (int i=0; i<n; i++)
    {
        planetpos << "        "<< setw(8)<<planets[i].name + "_x"<< "        "<<setw(8)<< planets[i].name + "_y";
        planetvel << "        "<< setw(8)<<planets[i].name + "_vx"<< "        "<<setw(8)<< planets[i].name + "_vy";
    }
    planetpos << endl;
    planetvel  << endl;
    while (t<time_period)
    {
        planetpos << Hour;
        planetvel  << Hour;
        for (int i=0; i<2*n; i++)
        {
            planetpos <<"        "<< setw(8)<<setprecision(5)<<P(i);
            planetvel <<"        "<< setw(8)<<setprecision(5)<<V(i);
        }
        planetpos << endl;
        planetvel  << endl;

        // Optional methods:
        //Euler();
        RK4();
        //Verlet();

        Hour+=1;
        t+=dt;
    }
    planetpos.close();
    planetvel.close();
    cout<<P(2) <<endl;
    cout<<P(3) <<endl;
    cout<<V(2) <<endl;
    cout<<V(3) <<endl;

}

/*-----------------------------------------------------------------------------------*/


int main()
{
    double g = 4*M_PI*M_PI;
    System Solarsystem;
    Solarsystem.Addplanet("Sun", 0.0, 0.0, 0.0, 0.0, 2e30);
    Solarsystem.Addplanet("Earth", 1.0, 0.0, 0.0, sqrt(g), 6e24);
    Solarsystem.Addplanet("Mars", 1.52, 0.0, 0.0, sqrt(g/1.52), 6.6e23);
    Solarsystem.Addplanet("Jupiter",5.20,0.0,0.0, sqrt(g/5.20), 1.9e27);


    //for (int i=0; i<Solarsystem.n; i++)
    //{
    //    cout << Solarsystem.planets[i]<<endl;
    //}

    Solarsystem.Setup();
    Solarsystem.Solve(1.);

    system("python /home/filiphl/Desktop/project3/p3plot.py");  //Runs the python script for plotting.
    return 0;
}
