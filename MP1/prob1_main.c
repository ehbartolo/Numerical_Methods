#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAXERR 0.0000001

/*  Problem # 1:
    Newton Raphson Method
*/

double therm_eff(double s)
{
    double n;
    double p1, p2, p3;
    double lambda, alpha;
    double k, l, h;
    double A, P;

    k = 240;
    h = 9;
    l = 0.1;

    A = s*s;
    P = s*4;
    lambda = sqrt((k*A)/(h*P));
    alpha = sqrt((h*A)/(k*P));

    p1 = sinh(l/lambda) + alpha*cosh(l/lambda);
    p2 = cosh(l/lambda) + alpha*sinh(l/lambda);
    p3 = (lambda/(l+(A/P)));

    n = (p1/p2)*p3;

    return n;

}
double derivative(double s)
{
    double dfx;
    double z, dz;
    double Alpha, dAlpha, Lambda, dLambda;
    double k, l, h;
    double A, B, C, D, E, F;
    double dA,dB,dC,dD,dE,dF;

    k = 240;
    h = 9;
    l = 0.1;

    z = 0.25*s;      /** Area over Perimeter **/
    dz = 0.25;       // Derivative of Area over Perimeter

    Alpha = sqrt(h/k) * sqrt(z);
    Lambda = sqrt(k/h) * sqrt(z);


    dAlpha = 0.5*sqrt(h/k)*(1/sqrt(z))*dz;
    dLambda = 0.5*sqrt(k/h)*(1/sqrt(z))*dz;

    F = l + z;
    dF = dz;

    E = Lambda;
    dE = dLambda;

    D = cosh(l/Lambda) + Alpha*sinh(l/Lambda);
    dD = (-0.1)*(l/(Lambda*Lambda))*(sinh(l/Lambda)+ \
    Alpha*cosh(l/Lambda))*dLambda + sinh(l/Lambda)*dAlpha;

    C = sinh(l/Lambda) + Alpha*cosh(l/Lambda);
    dC = (-0.1)*(l/(Lambda*Lambda))*(cosh(l/Lambda)+ \
    Alpha*sinh(l/Lambda))*dLambda + cosh(l/Lambda)*dAlpha;

    B = E/F;
    dB = (dE*F-E*dF)/(F*F);

    A = C/D;
    dA = (dC*D-C*dD)/(D*D);

    return dfx = dA*B + A*dB;
}
int main()
{
    double xk, xk_old, diff, abserr,root, c;
    int k = 0;

    xk = 0.1;                             //Initial Guess
    xk_old = xk;
    c = 0.95;                             //Desired Thermal Efficiency

    abserr = 100;

    printf("Iteration        Root        \
           Thermal_Efficiency     Absolute Error   \n");
    while (abserr > MAXERR) {
        if ((therm_eff(xk)-c) == 0) {     // xk is the root
            root = xk;
            break;
        }else if ( derivative(xk) == 0 ) { // provide another initial guess
            printf("Provide another initial Guess\n");
            return 0;
        }else {                            // Update value of xk
            xk = xk - ((therm_eff(xk)-c)/derivative(xk));
            root = xk;
            k++;                           // Update iteration number
        }

        diff = xk - xk_old;
        abserr = sqrt(diff*diff);          // Determine the absolute error
        //abserr = sqrt((therm_eff(xk)-c)*(therm_eff(xk)-c));
        if (k<10) {
            printf("  ");
        }else if(k<100) {
            printf(" ");
        }
        printf("   %d         %0.9lf       %0.9lf         \
               %0.9lf\n",k,root,therm_eff(root),abserr);

        xk_old = xk;                      // Update the previous iteration value
    }

    printf("\nConverged at %0.9lf meters, after %d iterations\n",root,k);
    printf("\nThus the square cross section has\n");
    printf("Side Length of:  %0.9lf meters\n",root);
    printf("Area of       :  %0.9lf square meter\n",root*root);
    printf("Perimeter of  :  %0.9lf meters\n",root*4);

    return 0;
}