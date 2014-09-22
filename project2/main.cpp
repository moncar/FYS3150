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


//Function for finding the indexes of the largest offdiagonal element, and return the value.
double max_offdiagonal(int n, mat A,int& k, int& l)
{
    double max = 0.0;
    for(int i=0; i<n;i++){
        for(int j=i+1; j<n; j++){
            if (fabs(A(i,j))>max){
                max = fabs(A(i,j));
                l = j;
                k = i;
            }
        }
    }
    return A(k,l);
}

//Funcion impelementing the rotation algorithm.
void rotate(mat& A, mat& R, int k, int l, int n)
{
    double c ,s, t, tau;
    if (A(k,l) != 0){
        tau = (A(l,l)-A(k,k))/(2*A(k,l));
        if (tau>0){
            t = 1./(tau + sqrt(1+tau*tau));}
        else{ t = 1./(tau - sqrt(1+tau*tau));}
        c = 1/sqrt(1+t*t);
        s = c*t;}
    else{c=1; s=0;}

    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);

    A(k,k) = c*c*a_kk - 2.*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.*c*s*A(k,l) + c*c*a_ll;
    A(l,k) = 0;
    A(k,l) = 0;

    for (int i=0; i<n; i++){
        if (i!=k && i!=l){
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
    return;
}

mat Jacobi(mat A, int n)
{
    //Setting up eigenvector matrix
    mat R = zeros<mat>(n,n);
    for (int i=0; i<n; i++){
        R(i,i) = 1;
    }

    int k, l;
    int iterations = 0;
    double epsilon = 1.0e-10; //Conditon for presission.
    double max_iterations = (double)n*(double)n*(double)n;
    double max_offdiag = max_offdiagonal(n, A, k, l); //Stores the largest offdiagonal elementvalue.

    while (fabs(max_offdiag)>epsilon && (double) iterations < max_iterations){
        max_offdiag = max_offdiagonal(n,A,k,l); //Changes k,l and max_offdiag.
        rotate(A,R,k,l,n);  //Does the rotation turning A(k,l) to zero.
        iterations++;
    }

    //Plot






/*
    Normering av egenvektorer.
    vec sR;
    cout<<"før: \n"<<R.col(0)<<endl;
    for (int i=0; i<n; i++){
        sR = sqrt(R.col(i).t()*R.col(i));
        double lR = sR(0);
        cout<<"lR = "<<lR<<endl;
        R.col(i) = R.col(i)/lR;
    }
    cout<<"etter \n"<<R.col(0)<<endl;

    ALLEREDE NORMERT!!!! <:o
*/

    //cout<<"iterations: "<<iterations<<endl;
    return A;
}


//--------------------------------------------------------------

int main()
{
    int n_step = 150;   //280 in b)
    int p_max = 5;
    double h = (double)p_max/n_step;

      /*-------------------------------*/
//    // Part B
//    //Setting up the tridiagonal matrix
//    mat A = zeros<mat>(n_step,n_step);

//    for (int i=1; i<n_step; i++)
//    {
//        A(i-1,i-1) = (2/(h*h))+pow((i*h),2);
//        A(i-1,i) = -1./(h*h);
//    }
//    //Adding the last diagonal element.
//    A(n_step-1,n_step-1)= (2./(h*h))+pow(((n_step)*h),2);
//    A = symmatu(A); //EXPLAIN THIS IN THE REPORT!

//    //Finding the eigenvalues
//    mat Eigval = Jacobi(A,n_step);
//    vec D = sort(Eigval.diag()); //Vector containing the eigenvalues sorted by size.

//    cout<<"Part b) \n\n";
//    cout<<"lambda0 = " << setw(5) << setprecision(4) << D(0)<<endl;
//    cout<<"lambda1 = "<< setw(5) << setprecision(4) << D(1)<<endl;
//    cout<<"lambda2 = "<< setw(5) << setprecision(4) << D(2)<<endl;


//    /*-----------------------------------------------------------*/
//    //oppgave c)
//    cout<<"\n\n\nPart c) \n \n";
//    double w[3];
//    w[0] = 0.01;
//    w[1] = 0.5;
//    w[2] = 1;
//    n_step = 150;

//    for (int j=0;j<3;j++){
//        //Setting up the new tridiagonal matrix (with the new potentials)
//        mat B = zeros<mat>(n_step,n_step);

//        double Vc;
//        for (int i=1; i<n_step; i++)
//        {
//            Vc = pow(w[j],2)*pow((i*h),2) + 1./(i*h);
//            B(i-1,i-1) = (2/(h*h))+Vc;
//            B(i-1,i) = -1./(h*h);
//        }
//        //Adding the last diagonal element.
//        B(n_step-1,n_step-1)= (2/(h*h))+pow(w[j],2)*pow((n_step*h),2) + 1./(n_step*h);
//        B = symmatu(B); //EXPLAIN THIS IN THE REPORT!

//        //Finding the eigenvalues
//        mat Eigval2 = Jacobi(B,n_step);
//        vec D2 = sort(Eigval2.diag()); //Vector containing the eigenvalues sorted by size.

//        cout<< "omega = "<< setw(4) << setprecision(5)<<w[j];
//        cout<<setw(15)<<"lambda0 = "<< setw(5) << setprecision(5) << D2(0)<<endl;
//    }


    /*---------------------------------------------------*/
    //Part d)

    //Setting up the tridiagonal matrix
    double w = 0.01;
    mat A = zeros<mat>(n_step,n_step);
    double Vd;
    for (int i=1; i<n_step; i++)
    {
        Vd = pow(w,2)*pow((i*h),2); // + 1./(i*h) if with repultion
        A(i-1,i-1) = (2/(h*h))+Vd;
        A(i-1,i) = -1./(h*h);
    }
    //Adding the last diagonal element.
    A(n_step-1,n_step-1)= (2./(h*h))+pow(w,2)*pow((n_step*h),2);
    A = symmatu(A); //EXPLAIN THIS IN THE REPORT!

    //Finding the eigenvalues
    mat Eigval = Jacobi(A,n_step);
    vec D = sort(Eigval.diag()); //Vector containing the eigenvalues sorted by size.

    cout<<"lambda0 ="<<setw(5)<<setprecision(4)<<D(0)<<endl;
    cout<<"lambda1 ="<<setw(5)<<setprecision(4)<<D(1)<<endl;
    cout<<"lambda2 ="<<setw(5)<<setprecision(4)<<D(2)<<endl;




    return 0;
}
