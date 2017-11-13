/**
* @file vectorf.c
* @Brief   
* @author Dong Yijing,1701110028@pku.edu.cn
* @version 1
* @date 2017-11-04
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solve.h"
void vectorf(double *f, double a, double b){
		int i,j;
    double c = -1.25;
    FILE *fp;
    /***********************************************************************************/
		/**<Initial the value of f according to the right hand side of the Possion equation*/
    /***********************************************************************************/
		for (i = 0; i < (2*N-1)*(N-1); i++){
				f[i]=0;
		}

		for (j = N/2 + 1;j < N; j++){
				/**< Initial the f=1 area*/
				for(i=3*N/2-j + 1; i <= N/2 + j; i++){
						f[(j - 1)*(2*N - 1) + i - 1] = 1;
				}
				/**< Initial the f=2 area*/
				for(i = j - N/2; i <= 3*N/2 - j; i++){
						f[(j - 1)*(2*N - 1) + i - 1] = 2;
				}
		}
		/**< Initial the f=3 area*/
		for (j = 1; j <= N/2; j++){
				for(i = N/2-j; i < N/2 + j; i++){
						f[(j - 1)*(2*N - 1) + i - 1] = 3;
				}
		}
		/**< Initial the f=4 area*/
		for (j = 1; j < N/2; j++){
				for(i = N/2 + j; i <= 3*N/2 - j - 1; i++){
						f[(j - 1)*(2*N - 1) + i - 1] = 4;
				}
		}
		/**< Initial the f=5 area*/
		for (j = 1; j < N/2; j++){
				for(i = 3*N/2 - j; i <= 3*N/2 + j; i++){
						f[(j - 1)*(2*N-1) + i - 1] = 5;
				}
		}
    fp = fopen("/home/egg/pictures/f1.txt","w");
    for (i = 0; i < (2*N - 1)*(N - 1); i++){
        if (fp != 0){
            fprintf (fp,"%f\n",f[i]);
        }
    }
    fclose(fp);
    /**************************************************************/
    /**< Change the value of f according to the boundary condition*/
    /**************************************************************/
    /**< Using the boundary conditon on f=1 and f=4*/
    for (j =  N/2 - 1 ; j < 3*N/2; j++){
        f[j] = f[j] + sqrt(b)*c;
        f[(2*N - 1)*(N - 2) + j] = f[(2*N - 1)*(N - 2) + j] + sqrt(b)*c;
    }

    fp = fopen("/home/egg/pictures/f2.txt","w");
    for (i = 0; i < (2*N - 1)*(N - 1); i++){
        if (fp != 0){
            fprintf (fp,"%f\n",f[i]);
        }
    }
    fclose(fp);
    for (i = 1; i < N/2; i++){
        f[(2*N - 1)*(N/2 + i - 1) + 2*N - i - 1] = sqrt(a)*2.0*(sqrt(3.0)*1.0)/3.0*c;/**< f=0*/
        f[(2*N - 1)*(N/2 + i - 1) + i - 1] = sqrt(a)*2.0*(sqrt(3.0)*1.0)/3.0*c;/**< on  f=2*/
        f[(i - 1)*(2*N - 1) + N/2 - i - 1] =  sqrt(a)*2.0*(sqrt(3.0)*1.0)/3.0*c;/**< on f=3*/
        f[(i - 1)*(2*N - 1) + 3*N/2 + i - 1] =  sqrt(a)*2.0*(sqrt(3.0)*1.0)/3.0*c;/**< on f=5*/
    }
    f[(2*N - 1)*(N/2 - 1)] = f[(2*N - 1)*(N/2 - 1)] - sqrt(a)*c;
    f[(2*N - 1)*(N/2 - 1) + 2*N - 1 -1] = f[(2*N - 1)*(N/2 - 1) + 2*N - 1 -1] + sqrt(a)*c;

    fp = fopen("/home/egg/pictures/f.txt","w");
    for (i = 0; i < (2*N - 1)*(N - 1); i++){
        if (fp != 0){
            fprintf (fp,"%f\n",f[i]);
        }
    }
    fclose(fp);

    return;
}
					  
