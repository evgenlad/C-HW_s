#include "b.h"

double nev(double *A, int n, int m) //строчная норма для алгоритма(максимальная сумма модулей элементов строки)
{
    double m6, m7;
    m6 = 0;
    m7 = 0;
    for(int i = 0; i < m; i++)
    {
        m6 = 0;
        for (int j = 0; j < m; j++)
        {
            if (A[i*n + j] > 0)
                m6 += A[i*n + j];
            else
                m6 += -A[i*n + j];
        }
        if (m7 < m6)
            m7 = m6;
    }
    return m7;
}

double nev1(double *A, int n) //след
{
    double res = 0;
    for (int i=0; i<n; i++)
    {
        res+=A[i*n + i];
    }
    return res;
}

double nev2(double *A, int n) //сумма квадратов
{
    double res=0;
    for (int i = 0; i < n*n; i++)
    {
        res+=A[i]*A[i];
    }
    return sqrt(res);
}

double nev3(double *A, int n, int fl)
{
    double res = 0;
    if (fl == 1)
    {    
        for(int i = 0; i < n; i++)
        {
            res += A[i*n+i]*A[i*n+i];
        }
    }
    else if (fl == -1)
    {
        for(int i = 0; i < n; i++)
        {
            res += A[i * n + i] * A[i * n + i];
        }
        res += 2 * A[1] * A[1];
    }
    //res = sqrt(res);
    return sqrt(res);
}
