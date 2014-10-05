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
    double G = 39.478;  //4*pi^2 AU^3 * M_sun^-1 * yr^-2
    double dt;       //Time-step in years.
    vec P;      // vector containing all positions.
    vec V;      // vector containing all velocities.
public:
    System(){}
    void Addplanet(string nm, double x0, double y0, double vx0, double vy0, double m);
    void Setup();
    vec Acceleration(vec A);
    void Solve(double);
    void Euler();
    void RK4();
    int n = 0;        //Number of objects.
    Planet planets[10];     //List of planets.
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

void System::Euler()    //Advancing one time-step.
{
    P += V*dt;
    V += Acceleration(P)*dt;
}

void System::RK4()  //Advancing one time-step.
{
    vec k1(2*n), k2(2*n), k3(2*n), k4(2*n); //Declaring vectors
    k1 = Acceleration(P) * dt;
    k2 = Acceleration(P + 0.5 * k1) * dt ;
    k3 = Acceleration(P+ 0.5 * k2) * dt ;
    k4 = Acceleration(P + k3 ) * dt ;
    V+=  (1/6.)*(k1 + 2 * (k2+k3) + k4);    //Updating velocity-vector
    P +=V*dt;
}

void System::Solve(double time_period)
{
    dt = 1/(365.0*24.); //Time-step one Hour.
    double t = 0;          //Set initial time to zero years.
    int Hour = 0;          //Counts itterations.

    fstream outfile;
    outfile.open("/home/filiphl/Desktop/project3/PlanetData.txt", ios::out);
    outfile <<"Number_of_objects:"<<"       " << n<<endl; // I use this in python...
    outfile <<"Hours";
    for (int i=0; i<n; i++)
    {
        outfile << "        "<< setw(8)<<planets[i].name + "_x"<< "        "<<setw(8)<< planets[i].name + "_y";
    }
    outfile << endl;
    while (t<time_period)
    {
        outfile <<Hour;
        for (int i=0; i<2*n; i++)
        {
            outfile <<"        "<< setw(8)<<setprecision(5)<<P(i);
        }
        outfile<<endl;

        // Optional methods:
        //Euler();
        RK4();

        Hour+=1;
        t+=dt;
    }
    outfile.close();
    cout<<P(2) <<endl;
}

/*-----------------------------------------------------------------------------------*/


int main()
{
//    vec Data(2);
//    Data(0) = 1;
//    Data(1) = 2;
//    cout << Data<<endl;
//    cout << Data+3.0<<endl;
//    vec B(2);
//    B(0) = 4;
//    B(1) = 10;

//    Planet Earth("earth", 1.0, 0.0, 0.0, 4.0, 6e24);
//    Planet Mars("mars", 1.0,  2.0, 3.0, 4.0 ,5.0);
System Solarsystem;
Solarsystem.Addplanet("Sun", 0.0, 0.0, 0.0, 0.0, 2e30);
Solarsystem.Addplanet("Earth", 1.0, 0.0, 0.0, 2*3.14, 6e24);
//Solarsystem.Addplanet("Moon", 0.99744, 0.0,0.0, -0.21498+ 2*3.14, 7.3477e22);
//Solarsystem.Addplanet("Mars", 1.52, 0.0, 0.0, 5.0963, 6.6e23);
//Solarsystem.Addplanet("Jupiter",5.20,0.0,0.0,-2.7554,1.9e27);
//for (int i=0; i<Solarsystem.n; i++)
//{
//    cout << Solarsystem.planets[i]<<endl;
//}

Solarsystem.Setup();
Solarsystem.Solve(10.0);

system("python /home/filiphl/Desktop/project3/p3plot.py");  //Runs the python script for plotting.
    return 0;
}
