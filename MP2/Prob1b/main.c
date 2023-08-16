#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
    ES 204 MP2 Problem #1b - Cubic Splines with Natural End Conditions
*/

//Computes the value of System of linear equations through LU Decomposition
double *LUDecomp(double a[4][4],double b[4]){
    double U[4][4], L[4][4];
    static double y[4], x[4];
    int total, col, row, k, i, j;

    col = sizeof(a[0])/sizeof(a[0][0]);
    row = col;

    //Initialize the L and U Matrix to zero
    for (j=0; j<col; j++){
        for (i=0; i <row; i++){
            U[i][j] = 0;
            L[i][j] = 0;
        }
    }

    for (j=0; j<col; j++){
        U[j][j] = 1;
        for (i=0; i <row; i++){
            // Compute for elements of L
            if (i>=j){
                L[i][j] = a[i][j];
                for(k = 0; k < j; k ++){
                    L[i][j] = L[i][j] - L[i][k]*U[k][j];
                }
            }
        }
        for (i=0; i <row; i++){
            // Compute for elements of U
            if (i>j) {
                U[j][i] = a[j][i];
                for(k = 0; k<j; k++){
                    U[j][i] = U[j][i] - L[j][k]*U[k][i];
                }
                U[j][i] = U[j][i]/L[j][j];
            }
        }
    }

    for(i=0;i<row;i++){
        //Forward Substitution
        y[i] = b[i];
        for(k = 0; k < i; k ++) {
            y[i] = y[i] - L[i][k]*y[k];
        }
        y[i] = y[i]/L[i][i];
    }
    for(i=(row-1); i >= 0;i--){
        //Backward Substitution
        x[i] = y[i];
        for(k=(row-1); k > i; k --) {
            x[i] = x[i] - U[i][k]*x[k];
        }
    }

    return x;

}

int main()
{
    double x[4] = {0.4,0.6,0.8,1};
    double f[4] = {20,50,105,180};
    double z[4][4];
    double c[4], u[4];
    double h[3],a[3],b[3],d[3];
    double divdiff[2];
    double *p;
    double res, num;
    long i, j;


    //Compute for Values of hi
    for(i=0;i<3;i++){
        h[i] = x[i+1] - x[i];
    }
    //Compute for Values of ui
    for(i=0;i<4;i++){
        if (i==0 || i==3){
            u[i] = 0; // condition for Natural Ends
        }else{
            divdiff[1] = (f[i+1]-f[i])/h[i];
            divdiff[0] = (f[i]-f[i-1])/h[i-1];
            u[i] = 3*(divdiff[1] - divdiff[0]);
        }
    }

    // Initialize to zero the z matrix
    for(i=0;i<4;i++){
        for(j=0;j<4;j++){
            z[i][j] = 0;
        }
    }

    //Compute values for z matrix
    for(i=0;i<4;i++){
        for(j=0;j<4;j++){
            if(i==j){
                //Middle Diagonal Elements
                if(i==0 || i==3){
                    z[i][i] = 1;
                }else{
                    z[i][i] = 2*(h[i-1]+h[i]);
                }
            }else if(i>0 && i <3){
                if(i == j-1){
                //Upper diagonal elements
                    z[i][j] = h[i];
                }else if(i==j+1){
                //Lower diagonal elements
                    z[i][j] = h[i-1];
                }
            }
        }
    }

    // Get the value of ci's through LU Decomposition
    p = LUDecomp(z,u);
    for(i=0;i<4;i++){
        c[i] = *(p+i);
    }

    printf("\n The cubic splines equations are: \n\n");
    for(i=0;i<3;i++){
        // Compute for ai,bi,di and print results
        a[i] = f[i];
        b[i] = ( (f[i+1]-f[i])/h[i] ) - (h[i]/3)*(2*c[i]+c[i+1]);
        d[i] = (c[i+1]-c[i])/(3*h[i]);

        printf(" S%d(x) = ",i);
        printf("%0.3lf",a[i]);
        printf(" + %0.3lf(x-%0.1lf)",b[i],x[i]);
        printf(" + %0.3lf(x-%0.1lf)^2",c[i],x[i]);
        printf(" + %0.3lf(x-%0.1lf)^3",d[i],x[i]);
        printf("    where %0.1lf <= x <= %0.1lf",x[i],x[i+1]);
        printf("\n\n");
    }

    //Evaluate D = 0.9
    num = 0.9;
    res = a[2]+b[2]*(num-x[2]) + c[2]*pow(num-x[2],2)+d[2]*pow(num-x[2],3);
    printf(" Value at %0.1lf m = %0.9lf MW\n",num,res);

    free(p);
    return 0;
}
