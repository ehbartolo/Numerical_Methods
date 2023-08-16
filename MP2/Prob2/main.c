#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define maxerr 0.00001

/*
    ES 204 MP2 Problem #2 - Adaptive Quadrature with Simpson's 1/3 Rule
*/

double fn(double t){
    double temp;

    temp = 20+10*(1-exp(-2*t))*sin(exp(t)-1);

    return temp;
}
double compFull(double h,double xi, double xii){

    return (h/6)*(fn(xi)+4*fn(xi+(h/2)) + fn(xii));
}
double compHalf(double h,double xi, double xii){

    double hlf;

    hlf = fn(xi)+4*fn(xi+(h/4)) + 2*fn(xi+(h/2));
    hlf = hlf + 4*fn(xi+((3*h)/4)) + fn(xii);
    hlf = (h/12)*hlf;

    return hlf;
}
double adaptivequad(double a0, double b0, double a, double b,double maxsubint, double *subint){
    double h,xi,xii;
    double Ifull,Ihalf,err,tol,mean,Iapp;

    Iapp = 0;
    h = b-a;
    xi = a;
    xii = b;

    Ifull = compFull(h,xi,xii);
    Ihalf = compHalf(h,xi,xii);

    err = fabs(Ifull - Ihalf)/15;
    tol = ((h)/(b0-a0))*(maxerr);

    if(err<=tol){
        Iapp = Iapp + Ihalf;
        //printf("Interval [%0.5lf,%0.5lf], error: %0.7lf\n",xi,xii,err);
    }else if(*subint >=maxsubint){
        printf("The maximum subinterval is exceeded\n");
        Iapp = Iapp + Ihalf;
        return Iapp;
    }else{
        mean = (a+b)/2;
        //Call again the recursive adaptive quadrature
        (*subint) = (*subint) + 1;
        Iapp = Iapp + adaptivequad(a0,b0,a,mean,maxsubint,subint); // Left Sub Interval
        Iapp = Iapp + adaptivequad(a0,b0,mean,b,maxsubint,subint); // Right Sub interval
    }

    return Iapp;

}

int main()
{
    double a,ai,b,bi,h,xi,xii;
    double Ifull,Ihalf,err,tol,mean,Iapp;
    int converge;
    double maxsub, subinter;

    subinter = 1;
    maxsub = pow(2,10);
    Iapp = adaptivequad(0,4,0,4,maxsub,&subinter);

    printf("\n\n Answer = %0.9lf\n",Iapp);
    printf(" Number of sub intervals = %0.0lf\n",subinter);
    return 0;
}
