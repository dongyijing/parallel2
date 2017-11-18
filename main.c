
#include <stdio.h>
#include <math.h>
#include "solve.h"
#include "mpi.h"

int main(int argc, char* argv[]){
    int n;
	int m1;
	int m2;
	int m3;
	int i,j;
	int k = 1;
    double delta_x;/**< Step length in x direction*/
    double delta_y;/**< Step length in y direction*/
    double b, c;
    double E;
	double *a;
	double *temp;
	double **d;
	double *f;
	double maxE;
	FILE *fp;
	int num ,num1,num2,num3;
	int *count;
	double *buff;
	double *a1;
	double *a2;
	double *a3;
	double wall_time;

    int size, myid;
    /**< MPI Initialization*/
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);


    n = (2*N-1)*(N-1);
    delta_x = 2*sqrt(3)/(3*N);
    delta_y = sqrt(3)*delta_x;

    b = 1/(delta_x*delta_x);
    c = 1/(delta_y*delta_y); 

    E = 0.000001;

		
    /**< initial a */ 
    a = (double *)malloc(sizeof(double)*n);
	for (i = 0; i < n; i++){
			a[i] = 0;
	}
    /* printf("\n");*/

    /**< store the new value after a iteration*/

    temp = (double *)malloc(sizeof(double)*n);
	count = (int *)malloc(sizeof(int)*4);/**< Record the number of every process*/
	for (i = 0; i < 4; i++){
		count[i] = 0;
    }
    /**< declar a \f$(2N-1)\times(5)\f& matrix d,where d contains all none zero elements in the matrix */
	d = malloc(sizeof(double*)*n);
	for(i = 0; i < n; i++){
		d[i] = malloc(sizeof(double)*n);
	}
	/**< use the function calculatea to fulfill the matrix A，only save none zero elements by diagonal*/
    matrixa(d, b, c);
    /**< declar the right hand side vector f with (2N-1)dimention*/
	f = (double*)malloc(sizeof(double)*n);
    vectorf(f, b, c);
	num = (N - 3)/(size - 1);

	buff = (double *)malloc(sizeof(double)*n);

	/**< Define vectors in process 0 to receive values from different process*/
	a1 = (double*)malloc(sizeof(double)*n);
	a2 = (double*)malloc(sizeof(double)*n);
	a3 = (double*)malloc(sizeof(double)*n);
	

    wall_time =MPI_Wtime();
    do{
		calculateJ(a, temp, n, b, c, d, f, size, myid, count, buff);
		k++;
		printf("k = %d",k);
        /*for (i = 0; i < n; i++){
            printf("temp[%d] = %f ",i,temp[i]);
        } 
        printf("\n");*/
		swapU(temp,a,n);
        /**< Exchange the value on each process*/

		/**< Send the value on other process to 0 process*/
		if (myid != 0){
			    MPI_Send(buff, n , MPI_DOUBLE, 0, 10 + myid, MPI_COMM_WORLD);
		}
		else{
		    	MPI_Recv(a1, n, MPI_DOUBLE, 1, 11, MPI_COMM_WORLD,&status);
		    	MPI_Recv(a2, n, MPI_DOUBLE, 2, 12, MPI_COMM_WORLD,&status);
		      	MPI_Recv(a3, n, MPI_DOUBLE, 3, 13, MPI_COMM_WORLD,&status);
        }
		m1 = 0;
		m2 = 0;
		m3 = 0;

		if (myid == 0){
           for (i = 0; i < 2*N - 1; i++){
				a[i] = buff[i];
				a[(2*N - 1)*(N - 2) + i] = buff[i + 2*N - 1];
		   }
		   for (j = 1; j < N - 2; j += (size - 1)){
			  if (j%(size - 1) == 1){
				 m1 = m1 + 1;
				 for ( i = 0; i < 2*N - 1; i++){
				  a[(2*N - 1)*j + i] = a1[(m1 - 1)*(2*N - 1) + i];
			     }
			  }
		   }
		   for (j = 2; j < N - 2; j += (size - 1)){
			  if (j%(size - 1) == 2){
			     m2 = m2 + 1;
			     for (i = 0; i < 2*N - 1; i++){
			      a[(2*N - 1)*j + i] = a2[(m2 - 1)*(2*N - 1) + i];
			     }
		      }
		   }
		   for (j = 3; j < N - 2; j += (size - 1)){
			   if (j%(size - 1) == 0){
			     m3 = m3 + 1;
			     for (i = 0; i < 2*N - 1; i++){
				        a[(2*N - 1)*j + i] = a3[(m3 - 1)*(2*N - 1) + i];
				 }
			  }
		   }
		}
		  /* for (i = 0; i < n; i++){
			   printf("a[%d] = %f\n",i,a[i]);
		   }*/
		/**<Send the value on 0 process to other process*/
        MPI_Bcast(a,n,MPI_DOUBLE,0,MPI_COMM_WORLD);

			/**< Calculate the max error*/
		maxE=fabs(temp[0]-a[0]);
		  for (i = 0;i < n; i++){
			if(maxE < fabs(temp[i]-a[i])){
				maxE = fabs(temp[i]-a[i]);
            }
	      }
	}while(maxE >= E && k < 200000);

	wall_time = MPI_Wtime() - wall_time;
	if( myid == 0){
		printf("Wall clock time = %f secs\n", wall_time);
	}

	MPI_Finalize();

	/**< Print solutions to the function*/
    fp = fopen("/home/egg/pictures/u256.txt","w");
    for (i = 0; i < (2*N-1)*(N-1); i++){
         if (fp != NULL){
            fprintf(fp,"%f ",a[i]);
            if ((i + 1)%(2*N - 1) == 0){
               fprintf(fp,"\n");
            }
         }
            /*printf("a[%d]=%f\n", i, a[i]);*/
    }
    fclose(fp);
	return 0;
}


   
