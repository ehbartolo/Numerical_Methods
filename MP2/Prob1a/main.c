#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
    ES 204 MP2 Problem #1a - Newton's Interpolating Polynomial order 3
*/

double evaluate(double a[],double xi[],double x ){
    int i,j;
    double fx = 0;
    double prod;

    for(i=0;i<4;i++){
        prod = a[i];
        if(i>0){
            for(j=i-1;j>=0;j--){
                prod = prod*(x-xi[j]);
            }
        }
        fx = fx + prod;
    }
    return fx;
}

int main()
{
    double xi[4] = {0.4,0.6,0.8,1};
    double fxi[4] = {20,50,105,180};
    double a[4];
    double derivative1[3], derivative2[2], derivative3;
    double res,num;
    int i,j,k ;

    for(i=0;i<3;i++){
        derivative1[i] = (fxi[i+1] - fxi[i])/(xi[i+1] - xi[i]);
    }
    for(i=0;i<2;i++){
        derivative2[i] = (derivative1[i+1] - derivative1[i])/(xi[i+2] - xi[i]);
    }
    derivative3 = (derivative2[1] - derivative2[0])/(xi[3] - xi[0]);

    a[0] = fxi[0];
    a[1] = derivative1[0];
    a[2] = derivative2[0];
    a[3] = derivative3;

    printf("\n fn(x) = \n");
    for(i=0;i<4;i++){
        printf(" + ");
        printf("%0.6lf",a[i]);
        if(i>0){

            for(j=i-1;j>=0;j--){
                printf("(");
                printf("x-%0.1lf",xi[j]);
                printf(")");
            }
        }
        printf("\n");
    }
    num = 0.9;
    res = evaluate(a,xi,num);
    printf("\n Value at %0.1lf m = %0.9lf MW\n",num,res);

    return 0;
}
