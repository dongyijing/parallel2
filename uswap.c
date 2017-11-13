
#include <stdio.h>
#include <stdlib.h>

void swapU(double *a, double *temp, int n){
    int i;
    double tmp;

    for (i = 0; i < n; i++){
        tmp = temp[i];
        temp[i] = a[i];
        a[i] = tmp;
    }
   /* for (i = 0; i < n; i++){
        printf("%f\n",a[i]);
    }*/
    return;
}
        

