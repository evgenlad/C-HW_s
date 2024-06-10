#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <time.h>
using namespace std;

const int f = 0;
void Error();
int fill(double* matrix, int n, int k, FILE* fin);
double nev(double *A, int n, int m);
int cr(double *, int, int);
double form (int, int, int, int);
int tr12 (double *A, int n, double eps);
int printM(double *A, int n, int s, int m);
int tr1(double *A, int n);
int tr2(double *A, int n, double eps);
int LR(double *A, int n, int m);
int printP(double *A, int n, int s, int m);
int IterLR(double *A, int n, int m, double epss);
int ev(double *A, int n, double epss);
int EigenPrint(double *A, int n, int m, int k);
double nev1(double *A, int n);
double nev2(double *A, int n);
double nev3(double *A, int n, int fl);
