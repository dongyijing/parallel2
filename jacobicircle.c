/**
* @file jacobicircle.c
* @Brief  Do a circle by using Jacobi method
* @author Dong Yijing,1701110028@pku.edu.cn
* @version 1
* @date 2017-11-05
*/

#include <stdio.h>
#include <stdlib.h>
#include "solve.h"

/**
* @Brief To do a circle by using Jacobi method.
*
* @Param a The initial value of the circle. 
* @Param temp The new value after a circle.
*/
void calculateJ(double *a,double *temp, int n,double b ,double c){
		int i,j;

    /*********************************************************/
    /**< To calculate f$ AU=f f$,fisrt,we should get A and f.*/
    /*********************************************************/

		/**< declar a \f$(2N-1)\times(5)\f& matrix d,where d contains all none zero elements in the matrix */
		double **d;
		d = malloc(sizeof(double*)*n);
		for(i = 0; i < n; i++){
			d[i] = malloc(sizeof(double)*n);
		}
		/**< use the function calculatea to fulfill the matrix A，only save none zero elements by diagonal*/
    matrixa(d, b, c);
    /**< declar the right hand side vector f with (2N-1)dimention*/
    double *f;
		f = (double*)malloc(sizeof(double)*n);
		vectorf(f, b, c);
    /**<Start the isolation and put the result into vector temp*/
    /**< for \f[ eu_{i,j}-au_{i-1,j}-bu_{i+1,j}-cu_{i,j-1}-du_{i,j+1} \]
     * where we let \f$e = d_{i,2},a = d_{i,1},b=d_{i,3},c=d_{i,0},d=d_{i,4}\f$
     * Then we have the iterative scheme 
     * \f[ u_{i,j}=\frac{1}{d((2*N-1)*(j-1)+i-1,2)(d((2*N-1)*(j-1)+i-1,1)u_{i-1,j}+d((2*N-1)*(j-1)+i-1,3)u_{i+1,j}+d((2*N-1)*(j-1)+i-1,4)u_{i,j+1}+d((2*N-1)*(j-1)+i-1,0)u_{i,j+1}+f_{i,j})\f]
     * where \f$ u_{i,j} = a[(2N-1)(j-1)+i-1], f_{i,j} = f[(2N-1)(j-1)+i-1]/f$
     */

    /*Jacobi isolate*/
    temp[0] = 0;
    temp[1] = 0;
    for (i = 2; i < 2*N - 1; i++){
        if (d[i][2] != 0){
            temp[i] = - (d[i][1]*a[i - 1] + d[i][3]*a[i + 1] + d[i][4]*a[i + 2*N - 1] - f[i])*(double)(1/d[i][2]);
        }
        else{
            temp[i] = 0;
        }
    }

    for (j = 1; j < N - 2; j++){
        for (i = 0; i < 2*N - 1; i++){
            if(d[(2*N-1)*j + i][2] != 0){
                temp[(2*N - 1)*j + i] = -(d[(2*N - 1)*j + i][1]*a[(2*N - 1)*j + i-1]
                    + d[(2*N - 1)*j + i][3] * a[(2*N - 1)*j + i + 1]
                    + d[(2*N - 1)*j + i][4] * a[(2*N - 1)*(j + 1) + i]
                    + d[(2*N - 1)*j + i][0] * a[(2*N - 1)*(j - 1) + i]
                    - f[(2*N - 1)*j + i])*(double)(1/d[(2*N - 1)*j + i][2]);
            }
            else{
                temp[(2*N - 1)*j + i] = 0;
            }
        }
    }
    /**< On the last block matrix*/
    for (i = 0; i < 2*N - 3; i++){
        if (d[(2*N - 1)*(N - 2) + i][2] != 0){
            temp[(2*N - 1)*(N - 2) + i] =  -(d[(2*N - 1)*(N - 2) + i][1]*a[(2*N - 1)*(N - 2) + i - 1] 
                                   + d[(2*N - 1)*(N - 2) + i][3]*a[(2*N - 1)*(N - 2) + i + 1] 
                                   + d[(2*N - 1)*(N - 2) + i][0]*a[(2*N - 1)*(N - 3) + i] 
                                   - f[(2*N - 1)*(N - 2) + i])*(double)(1/d[(2*N - 1)*(N - 2) + i][2]);
        }
        else{
            temp[i] = 0;
        }
    }
    temp[2*N - 2] = 0;
    temp[2*N - 3] = 0;
    /******************/
    /**< SOR isolate***/
    /******************/
   /** temp[0] = 0;
    temp[1] = 0;
    for (i = 2; i < 2*N-1; i++){
        if (d[i][2] != 0){
            temp[i] = (1 - w)*a[i]-w*(d[i][1]*a[i - 1] + d[i][3]*a[i + 1] + d[i][4]*a[i + 2*N - 1] - f[i])*(double)(1/d[i][2]);
        }
        else{
            temp[i] = 0;
        }
    }

    for (j = 1; j < N - 2; j++){
        for (i = 0; i < 2*N - 1; i++){
            if(d[(2*N-1)*j + i][2] != 0){
                temp[(2*N - 1)*j + i] = (1 - w)*a[(2*N - 1)*j + i] - w*(d[(2*N - 1)*j + i][1]*a[(2*N - 1)*j + i-1]
                    + d[(2*N-1)*j + i][3] * a[(2*N-1)*j + i + 1]
                    + d[(2*N-1)*j + i][4] * a[(2*N-1)*(j + 1) + i]
                    + d[(2*N-1)*j + i][0] * a[(2*N-1)*(j - 1) + i]
                    - f[(2*N-1)*j + i])*(double)(1/d[(2*N-1)*j + i][2]);
            }
            else{
                temp[(2*N-1)*j + i] = 0;
            }
        }
    }
    for (i = 0; i < 2*N-3; i++){
        if (d[(2*N - 1)*(N-2) + i][2] != 0){
            temp[(2*N-1)*(N-2) + i] = (1 - w)*a[(2*N - 1)*(N - 2) + i] -w*(d[(2*N - 1)*(N - 2) + i][1]*a[i - 1] 
                                   + d[(2*N - 1)*(N - 2) + i][3]*a[i + 1] 
                                   + d[(2*N - 1)*(N - 2) + i][0]*a[(2*N - 1)*(N - 3) + i] 
                                   - f[(2*N - 1)*(N - 2) + i])*(double)(1/d[(2*N - 1)*(N-2) + i][2]);
        }
        else{
            temp[i] = 0;
        }
    }
    temp[2*N - 2] = 0;
    temp[2*N - 3] = 0;
   */ 

   return;
}

