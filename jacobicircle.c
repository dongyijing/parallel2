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
void calculateJ(double *a,double *temp, int n, double b ,double c, double **d, double *f, int size, int myid, int *count, double *buff){
		int i,j;
		int size1;
		for (i = 0; i < 4; i++){
			count[i] = 0;
		}

    /**<Start the isolation and put the result into vector temp*/
    /**< for \f[ eu_{i,j}-au_{i-1,j}-bu_{i+1,j}-cu_{i,j-1}-du_{i,j+1} \]
     * where we let \f$e = d_{i,2},a = d_{i,1},b=d_{i,3},c=d_{i,0},d=d_{i,4}\f$
     * Then we have the iterative scheme 
     * \f[ u_{i,j}=\frac{1}{d((2*N-1)*(j-1)+i-1,2)(d((2*N-1)*(j-1)+i-1,1)u_{i-1,j}+d((2*N-1)*(j-1)+i-1,3)u_{i+1,j}+d((2*N-1)*(j-1)+i-1,4)u_{i,j+1}+d((2*N-1)*(j-1)+i-1,0)u_{i,j+1}+f_{i,j})\f]
     * where \f$ u_{i,j} = a[(2N-1)(j-1)+i-1], f_{i,j} = f[(2N-1)(j-1)+i-1]/f$
     */


	/**< parallel iteration*/
	size1 = size - 1;
	if (myid == 0){ /**< put the value on the first row and the last row into process 0*/
		temp[0] = 0;
		count[0] = count[0] + 1;
		buff[count[0] - 1] = temp[0];
		temp[1] = 0;
		count[0] = count[0] + 1;
		buff[count[0] - 1] = temp[1];
      for (i = 2; i < 2*N - 1; i++){
          if (d[i][2] != 0){
              temp[i] = - (d[i][1]*a[i - 1] + d[i][3]*a[i + 1] + d[i][4]*a[i + 2*N - 1] - f[i])*(double)(1/d[i][2]);
          }
          else {
              temp[i] = 0;
          }
		  count[0] = count[0] + 1;
		  buff[count[0] - 1] = temp[i];
      }

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
		  count[0] = count[0] + 1;
		  buff[count[0] - 1] = temp[(2*N - 1)*(N - 2) + i];
      }
		temp[2*N - 2] = 0;
		count[0] = count[0] + 1;
		buff[count[0] - 1] = temp[2*N - 2];
		temp[2*N - 3] = 0;
		count[0] = count[0] + 1;
		buff[count[0] - 1] = temp[2*N - 3];
	}        
	if (myid == 1){ /**< put values to process 1*/ 
		for (j = myid; j < N - 2; j += size1){
          for (i = 0; i < 2*N - 1; i++){
			/**< Fill the value of temp*/
            if(d[(2*N-1)*j + i][2] != 0){
                temp[(2*N - 1)*j + i] = -(d[(2*N - 1)*j + i][1]*a[(2*N - 1)*j + i - 1]
                    + d[(2*N - 1)*j + i][3] * a[(2*N - 1)*j + i + 1]
                    + d[(2*N - 1)*j + i][4] * a[(2*N - 1)*(j + 1) + i]
                    + d[(2*N - 1)*j + i][0] * a[(2*N - 1)*(j - 1) + i]
                    - f[(2*N - 1)*j + i])*(double)(1/d[(2*N - 1)*j + i][2]);
            }
            else{
                temp[(2*N - 1)*j + i] = 0;
            }

			  count[1] = count[1] +1;
			  buff[count[1] - 1] = temp[(2*N - 1)*j + i];
		  }
		}
		/*for( i = 0; i < count[1]; i++){
		printf("In process 1,buff[%d] = %f\n",i,buff[i]);
		}*/
	}
	if (myid == 2){/**< Put values to process 2*/
		for (j = myid; j < N - 2; j += size1){
          for (i = 0; i < 2*N - 1; i++){
            if(d[(2*N-1)*j + i][2] != 0){
                temp[(2*N - 1)*j + i] = -(d[(2*N - 1)*j + i][1]*a[(2*N - 1)*j + i - 1]
                    + d[(2*N - 1)*j + i][3] * a[(2*N - 1)*j + i + 1]
                    + d[(2*N - 1)*j + i][4] * a[(2*N - 1)*(j + 1) + i]
                    + d[(2*N - 1)*j + i][0] * a[(2*N - 1)*(j - 1) + i]
                    - f[(2*N - 1)*j + i])*(double)(1/d[(2*N - 1)*j + i][2]);
            }
            else{
                temp[(2*N - 1)*j + i] = 0;
            }
			  count[2] = count[2] +1;
			  buff[count[2] - 1] = temp[(2*N - 1)*j + i];
			}
		}		
		/*for (i = 0; i < count[2]; i++){
		  printf("In process 2,a[%d] = %f\n",i,a[i]);
	    }*/
 	}
	if (myid == 3){/**<Put values to process 3*/ 
		for (j = myid; j < N - 3; j += size1){
          for (i = 0; i < 2*N - 1; i++){
            if(d[(2*N-1)*j + i][2] != 0){
                temp[(2*N - 1)*j + i] = -(d[(2*N - 1)*j + i][1]*a[(2*N - 1)*j + i - 1]
                    + d[(2*N - 1)*j + i][3] * a[(2*N - 1)*j + i + 1]
                    + d[(2*N - 1)*j + i][4] * a[(2*N - 1)*(j + 1) + i]
                    + d[(2*N - 1)*j + i][0] * a[(2*N - 1)*(j - 1) + i]
                    - f[(2*N - 1)*j + i])*(double)(1/d[(2*N - 1)*j + i][2]);
            }
            else{
                temp[(2*N - 1)*j + i] = 0;
            }
			  count[3] = count[3] + 1;
			  buff[count[3] - 1] = temp[(2*N - 1)*j + i];
          }
        }
		/*for (i = 0; i < count[3]; i++){
		printf("In process 1,buff[%d] = %f\n",i,buff[i]);
		}*/
    }
	
   return; 
}

