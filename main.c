
#include <stdio.h>
#include <math.h>
#include "solve.h"

int main(){
    int n;
    double delta_x;/**< Step length in x direction*/
    double delta_y;/**< Step length in y direction*/
    double a, b;
    double E;

    n = (2*N-1)*(N-1);
    delta_x = 2*sqrt(3)/(3*N);
    delta_y = sqrt(3)*delta_x;

    a = 1/(delta_x*delta_x);
    b = 1/(delta_y*delta_y); 

    E = 0.001;

    jacobi(n, E, a, b);
}


   
