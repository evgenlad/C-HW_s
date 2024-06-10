#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include "mpi.h"
using namespace std;


void Error();
double ab(double a);
void printT (int k);
void printPM(double *A, int n, int m);
int fill(double* matrix, int n, FILE* fin, int p, int pn);
int nev(double *, double *, int);
int cr(double *A, int n, int k, int p, int pn, FILE* fin);
int makeE(double *B, int n, int p, int pn);
int printM(double *A, int n, int m, int p, int pn);
void printMes (int p, int pn, int code, double m);
int tr1(double *A, double *B, int n, int p, int pn);
int tr2(double *A, double *B, int n, int p, int pn);
int tr12(double *A, double *B, int n, int p, int pn);
int fill0(double* matrix, int n, int p, int pn);
double nev5 (double *A, int n, int p, int pn);
double nev10 (double *A, int n, int p, int pn);
double nev(double *A, double *B, int n, int p, int pn);

int normalize (double *A, int n, double n11, int p, int pn);
//void *tr0(void *ae);
//void synchronize(int threads_count);
