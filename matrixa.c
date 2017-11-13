/**
* @file matrixa.c
* @Brief  Return the coefficient matrix of the linear equations 
* @author Dong Yijing,1701110028@pku.edu.cn
* @version 1
* @date 2017-11-04
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solve.h"

void matrixa(double **d, double b,double c){
		int i, j, k;
    FILE *fp;
    /*******************************************************/
		/**< Initial the value without using the boundary value*/
    /*******************************************************/

		/**< Input the inital value among the diagonal*/
		for(i = 0; i < (2*N - 1)*(N - 1); i++){
				d[i][2] = 2*(b+c);
		}
		/**< Input the second diagonal value on the right of the main dioagonal*/
		for(i = 0; i < (2*N - 1)*(N - 1)-1; i++){
				d[i][3] = -b;
		}
		for(j = 2*N-2; j < (2*N-1)*(N-1); j += 2*N - 1){/**< change some elements into 0 because of the block matrix*/
				d[j][3] = 0;
		}
		/**< Input the second diagonal value on the left of the main dioagonal*/
		for(i = 1;i < (2*N-1)*(N-1);i++){
				d[i][1] = -b;
		}
		for(j = 2*N-1;j < (2*N-1)*(N-1);j += 2*N-1){
				d[j][1]=0;
		}
		/**< Input the value of the third diagonal on the right*/
		for(i = 0;i < (2*N-1)*(N-2); i++){
				d[i][4] = -c;
		}
		/**< Input the value of the third diagonal on the left*/
		for(i = 2*N-1; i < (2*N-1)*(N-1); i++){
				d[i][0] = -c;
		}
    
    /************************************************/
		/** using boundary conditions                   */
		/***********************************************/

		/**< using boundary conditions on f=1 and f=4*/

		/* On f=4*/
		for(i = N/2;i <= 3*N/2;i++){
				d[i - 1][2]=d[i - 1][2] - c;
		}
		/* On f=1 */
		for(i = N/2 + (2*N - 1)*(N - 2);i <= 3*N/2 + (2*N - 1)*(N - 2); i++){
				d[i - 1][2]=d[i - 1][2] - c;
		}

		for(i = 1; i < N/2; i++){
				/**< using boundary condition on f=0*/
				d[(N/2 + i - 1)*(2*N - 1) + 2*N - i - 1][2] = b + c;
				d[(N/2 + i - 1)*(2*N - 1) + 2*N - i - 1][3] = 0;
				d[(N/2 + i - 1)*(2*N - 1) + 2*N - i - 1][4] = 0;
				/**< using boundary condition on f=2*/
				d[(N/2 + i - 1)*(2*N - 1) + i - 1][2] = b + c;
				d[(N/2 + i - 1)*(2*N - 1) + i - 1][1] = 0;
				d[(N/2 + i - 1)*(2*N - 1) + i - 1][4] = 0;
			  /**< using boundary condition on f=3*/
				d[(i - 1)*(2*N - 1) + N/2 - i - 1][2] = b + c;
				d[(i - 1)*(2*N - 1) + N/2 - i - 1][1] = 0;
				d[(i - 1)*(2*N - 1) + N/2 - i - 1][0] = 0;
        /**< using boundary condition on f=5*/
				d[(i - 1)*(2*N - 1) + 3*N/2 + i - 1][2] = b + c;
				d[(i - 1)*(2*N - 1) + 3*N/2 + i - 1][3] = 0;
				d[(i - 1)*(2*N - 1) + 3*N/2 + i - 1][0] = 0;
		}
    
    /*************************************************************/
    /**< Change  coefficients of points out of the area into zero*/
    /**************************************************************/

    /**< Left down corner*/
    for (k = 0; k < 5; k++){
        for (j = 1; j < N/2-1; j++){
            for (i = 1; i < N/2 - j; i++){
                d[(2*N - 1)*(j - 1) + i - 1][k] = 0;
              /* if ((2*N - 1)*(j - 1) + i - 1 + 2*N - 1 < (2*N - 1)*(N - 1)){
                d[(2*N - 1)*(j - 1) + i - 1 + 2*N - 1][4] = 0;
                }
                if ((2*N - 1)*(j - 1) + i - 1 - (2*N - 1) >= 0 ){
                d[(2*N - 1)*(j - 1) + i - 1 - (2*N - 1)][0] = 0;
                }
                if ((2*N - 1)*(j - 1) + i  < (2*N - 1)*(N - 1)){
                d[(2*N - 1)*(j - 1) + i ][3] = 0;
                }
                if ((2*N - 1)*(j - 1) + i - 2 >= 0){
                d[(2*N - 1)*(j - 1) + i - 2][1] = 0;
                }*/
            }
        }
    }
    /**< Right down corner*/
    for (k = 0; k < 5; k++){
        for (j = 1; j < N/2 -1; j++){
            for ( i = 3*N/2 + j + 1; i < 2*N; i++){
                d[(2*N-1)*(j-1) + i - 1][k] = 0;
               /* if ((2*N - 1)*(j - 1) + i -1 + 2*N - 1 < (2*N-1)*(N-1)){
                d[(2*N - 1)*(j - 1) + i -1 + 2*N - 1][4] = 0;
                }
                if ((2*N - 1)*(j - 1) + i - 1 - (2*N - 1) >= 0 ){
                d[(2*N - 1)*(j - 1) + i - 1 - (2*N - 1)][0] = 0;
                }
                if ((2*N - 1)*(j - 1) + i  < (2*N-1)*(N-1)){
                d[(2*N - 1)*(j - 1) + i ][3] = 0;
                }
                if ((2*N - 1)*(j - 1) + i - 2 >= 0){
                d[(2*N - 1)*(j - 1) + i - 2][1] = 0;
                }*/
            }
        }
    }
    /**< Left up corner*/
    for (k = 0; k < 5; k++){
        for (j = N/2 + 2; j < N; j++){
            for (i = 1; i < j - N/2; i++){
                d[(2*N-1)*(j-1) + i - 1][k] = 0;
               /* if ((2*N - 1)*(j - 1) + i -1 + 2*N - 1 < (2*N-1)*(N-1)){
                d[(2*N - 1)*(j - 1) + i - 1 + 2*N - 1][4] = 0;
                }
                if ((2*N - 1)*(j - 1) + i - 1 - (2*N - 1) >= 0 ){
                d[(2*N - 1)*(j - 1) + i - 1 - (2*N - 1)][0] = 0;
                }
                if ((2*N - 1)*(j - 1) + i  < (2*N - 1)*(N - 1)){
                d[(2*N - 1)*(j - 1) + i ][3] = 0;
                }
                if ((2*N - 1)*(j - 1) + i - 2 >= 0){
                d[(2*N - 1)*(j - 1) + i - 2][1] = 0;
                }*/
            }
        }
    }
    /**< Right up corner*/
    for (k = 0; k < 5; k++){
        for ( j = N/2 + 2; j < N; j++){
            for ( i = 5*N/2 - j + 1; i < 2*N; i++){
                d[(2*N - 1)*(j - 1) + i - 1][k] = 0;
                /*if ((2*N - 1)*(j - 1) + i -1 + 2*N - 1 < (2*N-1)*(N-1)){
                d[(2*N - 1)*(j - 1) + i -1 + 2*N - 1][4] = 0;
                }
                if ((2*N - 1)*(j - 1) + i -1 - (2*N - 1) >= 0 ){
                d[(2*N - 1)*(j - 1) + i -1 - (2*N - 1)][0] = 0;
                }
                if ((2*N - 1)*(j - 1) + i  < (2*N-1)*(N-1)){
                d[(2*N - 1)*(j - 1) + i ][3] = 0;
                }
                if ((2*N - 1)*(j - 1) + i -2 >= 0){
                d[(2*N - 1)*(j - 1) + i - 2][1] = 0;
                }*/
            }
        }
    }

    fp = fopen("/home/egg/pictures/matrix.txt","w");
    for (j = 0; j < (2*N - 1)*(N - 1); j++){
        for (i = 0; i < 5; i++){
           if( fp != 0){
               fprintf(fp,"%f ",d[j][i]);
             if (i == 4){
               fprintf (fp,"\n");
             }
           }
        }
    }
    fclose(fp);
  	return;
}
