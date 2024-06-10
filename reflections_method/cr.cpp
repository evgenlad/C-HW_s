#include "b.h"

void Error() {
	printf("Incorrect data in the terminal line\n");
}

int mult_m_v(double *A, double *b, double *x, int n)
{
    for (int i = 0; i < n; i++)
        x[i] = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            x[i] += A[i*n + j] * b[j];
        }
    }
    return 0;
}

int makeE(double *A, int n) //try-catch сюда!
{
   A[0]=1.0;
    for (int i = 1; i < n * n; i++)
        A[i] = ((i%(n+1))==0);
    return 0;
}

int copydata(double *A, double *B, int n)
{
    for (int k1 = 0; k1 < n*n; k1++)
        A[k1] = B[k1];
    return 0;
}

double form(int k, int n, int i, int j)
{
    if (k == 1)
        return (n - max(i, j) + 1);
    if (k == 2)
        return (max(i, j));
    if (k == 3)
        return (max(i - j, j - i));
    if (k == 4)
        return (1.0/(max(i + j - 1, - j - i + 1)));
    return -1;
}

int fill(double* matrix, int n, FILE* fin) 
{    
	int s = 0;
	double x;
	while (fscanf(fin, "%lf", &x) == 1) 
    {
        matrix[s] = x;
        s++;
    }
    if (s != n * n) 
    {
        return -1;
    }
	return 0;
}

int cr(double *A, int n, int k, FILE* fin)
{
    if ((k == 0)&&(fin != nullptr))
    {
        if (fill(A, n, fin)!=0) 
		    return -3;
    }
    else if (k == 0)
        return -1;
    if ((k != 0)&&(fin != nullptr))
        return -1;
    if ((k == 1)||(k == 2)||(k == 3)||(k == 4))
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                A[i * n + j] = form (k, n, i + 1, j + 1);
    return 0;
}

int printM(double *A, int n, int s, int m)
{
    for (int i=0; i < min(n, m); i++)
    {
        for (int j=0; j<min(n, m); j++)
            printf("%e\t", A[i*n + j + s]);
        cout<<endl;
    }
    cout<<endl<<endl<<endl;
    return 0;
}

int divM(double *A, int n, double d, int flag)
{
    if (flag == 1)
    {
        d = 1/d;
        for(int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                A[i*n+j] = A[i*n+j] * d;
            }
        }
    return 0;
    }
    if (flag == 0)
    {
        for(int i= 0; i < n; i++)
        {
            for (int j= 0; j < n; j++)
            {
                A[i*n+j] *= d;
            }
        }
    return 0;
    }   
    return 0;
}

int normalize(double* A, double n10, int n, int flag) // начиная с A[0], нормализовать m строк длиной n
{
    double m6;
    if (flag == 1)
    {
        for(int i = 0; i < n; i++)
        {
            m6 = A[n*n + n*n + i];
            for(int j = 0; j < n; j++)
            {
                A[j*n + i + n*n] *= m6;
            }
        }
    }
    if (flag == 2)
    {
        m6 = 1/n10;
        for (int i = 0; i < n*n; i++)
        {
            A[i] *= m6;
        }
    }
    return 0;
}
