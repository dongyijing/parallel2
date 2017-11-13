
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define N 24
#ifndef __SOLVE_H
#define __SOLVE_H

void jacobi(int n, double E, double b, double c);
void calculateJ(double *a, double *temp, int n, double b, double c);
void matrixa(double **d, double b, double c);
void vectorf(double *f, double a, double b);
void swapU(double *a, double *temp, int n);
#endif
