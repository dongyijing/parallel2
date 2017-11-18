
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"
#define N 256
#ifndef __SOLVE_H
#define __SOLVE_H

void calculateJ(double *a, double *temp, int n, double b, double c, double **d,double *f, int size, int myid, int *count, double *buff);
void matrixa(double **d, double b, double c);
void vectorf(double *f, double a, double b);
void swapU(double *a, double *temp, int n);
#endif
