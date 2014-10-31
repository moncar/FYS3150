#include "plot.h"

Plot::Plot()
{
    cout<< "Please give two vectors as arguments. \n Exampe: Plot(vec x, vec y)"<<endl;
}

void Plot::basicstructure(vec a, vec b)
{
    if (a.n_elem != b.n_elem)
    {
        cout << "Vectors must be of the same size."<<endl;
        exit(0);
    }

    n = a.n_elem;
    fstream outFile;
    outFile.open("plotdata.txt", ios::out );
    for (int i=0; i<n; i++)
    {
        outFile << a(i)<<"      "<<b(i)<<endl;
    }
    outFile.close();

    fstream pyscript;
    pyscript.open("pythonplot.py", ios::out);
    pyscript << "from numpy import *"<<endl;
    pyscript << "from matplotlib import pyplot as plt"<<endl;
    pyscript << "a=[]"<<endl;
    pyscript << "b=[]"<<endl;
    pyscript << "infile = open('plotdata.txt', 'r')"<<endl;
    pyscript << "for line in infile:"<<endl;
    pyscript << "   words = line.split()"<<endl;
    pyscript << "   a.append(float(words[0]))"<<endl;
    pyscript << "   b.append(float(words[1]))"<<endl;
}

Plot::Plot(vec a, vec b)
{
    basicstructure(a,b);
    fstream pyscript;
    pyscript.open("pythonplot.py", ofstream::out | ofstream::app);
    pyscript << "plt.plot(a,b)"<<endl;
    pyscript.close();
}

Plot::Plot(vec a, vec b, string x)
{
    basicstructure(a,b);
    fstream pyscript;
    pyscript.open("pythonplot.py", ofstream::out | ofstream::app);
    pyscript << "plt.plot(a,b,'"<<x<<"')"<<endl;
    pyscript.close();
}

void Plot::Title(string T)
{
    fstream pyscript;
    pyscript.open("pythonplot.py", ofstream::out | ofstream::app);
    pyscript << "plt.title(' "<<T<<" ')"<<endl;
    pyscript.close();
}

void Plot::Labels(string A, string B)
{
    fstream pyscript;
    pyscript.open("pythonplot.py", ofstream::out | ofstream::app);
    pyscript << "plt.xlabel(' "<<A<<" ')"<<endl;
    pyscript << "plt.ylabel(' "<<B<<" ')"<<endl;
    pyscript.close();
}

void Plot::Legend(string L)
{
    fstream pyscript;
    pyscript.open("pythonplot.py", ofstream::out | ofstream::app);
    pyscript << "plt.legend([' "<<L<<" '])"<<endl;
    pyscript.close();
}

void Plot::Axis(double xmin, double xmax, double ymin, double ymax)
{
    fstream pyscript;
    pyscript.open("pythonplot.py", ofstream::out | ofstream::app);
    pyscript << "plt.axis(["<<xmin<<","<<xmax<<","<<ymin<<","<<ymax<<"])"<<endl;
    pyscript.close();
}

void Plot::Axis(string E)
{
    fstream pyscript;
    pyscript.open("pythonplot.py", ofstream::out | ofstream::app);
    pyscript << "plt.axis(' "<<E<<" ')"<<endl;
    pyscript.close();
}


void Plot::Figure(int f)
{
    fstream pyscript;
    pyscript.open("pythonplot.py", ofstream::out | ofstream::app);
    pyscript << "plt.figure("<<f<<")"<<endl;
    pyscript.close();
}


void Plot::Show()
{
    fstream pyscript;
    pyscript.open("pythonplot.py", ofstream::out | ofstream::app);
    pyscript << "plt.show()"<<endl;
    pyscript.close();
    system("python pythonplot.py");
    system("rm plotdata.txt pythonplot.py");
}

