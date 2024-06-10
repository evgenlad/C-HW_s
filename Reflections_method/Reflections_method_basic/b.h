#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <time.h>
using namespace std;

//const int f = 0;
void Error();
int mult_m_v(double *A, double *b, double *x, int n);
int fill(double* matrix, int n, FILE* fin);
int copydata(double *A, double *B, int n);
int makeE(double *A, int n);
double form(int k, int n, int i, int j);
double nev(double *, double *, int);
double nev4(double* A, int n);
double nev5(double* A, int n);
int cr(double *, int, int, FILE* fin);
int divM(double *A, int n, double d, int flag);
int printM(double *A, int n, int s, int m);
int tr12(double *A, double *A1, int n);
int tr1(double *A, double *A1, int n);
int tr2(double *A, double *A1, int n);
int normalize(double* A, double n10, int n, int flag);
double nev9(double *A, int n);
double nev10(double* A, int n);
double nev11(double *A, int n);
double nev0(double *A, double *B, int n);
