/**
* @file jacobi.c
* @Brief  Using the Jacobi iterection to calculate the solutions
* @author Dong Yijing,1701110028@pku.edu.cn
* @version 1
* @date 2017-11-04
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solve.h"


/**
* @Param n 
*/
void jacobi( int n, double E, double b, double c){
		int i;/**< for circle*/
		int k = 0;/**< record the number of iterations*/
		double maxE;/**< the error between two iterations*/;
    double *a;
    double *temp;
    FILE *fp;
		
    printf("b=%f,c=%f\n",b,c);
		
		/**< initial a */ 
		a = (double *)malloc(sizeof(double)*n);
		for (i=0;i<n;i++){
				a[i] = 0;
		}
   /* printf("\n");*/

		/**< store the new value after a iteration*/
    printf("n= %d\n",n);
		temp =(double *)malloc(sizeof(double)*n);
		do{
				calculateJ(a,temp,n,b,c);
        for (i = 0; i < n; i++){
            printf("temp[%d] = %f ",i,temp[i]);
        } 
        printf("\n");
				k++;
				swapU(temp,a,n);
				maxE=fabs(temp[0]-a[0]);
				for (i=0;i<n;i++){
						if(maxE<fabs(temp[i]-a[i])){
							maxE=fabs(temp[i]-a[i]);
              }
				}
		}while(maxE >=E);
    fp = fopen("/home/egg/pictures/a.txt","w");
    for (i = 0; i < (2*N-1)*(N-1); i++){
        if (fp != NULL){
            fprintf(fp,"%f ",a[i]);
            if ((i+1)%(2*N-1) == 0){
                fprintf(fp,"\n");
            }
        }
            printf("a[%d]=%f\n", i, a[i]);
    }
    fclose(fp);
		return;
}
	


