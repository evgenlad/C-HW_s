#include "b.h"

double nev(double *A, double *B, int n)
{
    //Невязочка
    double m6, m7;
    m6 = 0;
    m7 = 0;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            m6 = 0;
            for(int k = 0; k < n; k++)
            {
                m6 += B[i * n + k] * A[k * n + j]; // считаем (ij) элемент невязки
            } 
            if (i == j)
                m6 = m6 - 1;
            m7 += m6*m6; //сразу норму
        }
    }
    m7 = sqrt(m7);

    cout<<"Невязка= "<<m7<<endl;

    return m7;
}

double nev0(double *A, double *B, int n)
{
    //Невязочка
    double m6, m7, m8;
    m6 = 0;
    m7 = 0;
    m8 = 0;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            m6 = 0;
            for(int k = 0; k < n; k++)
            {
                m6 += B[i * n + k] * A[k * n + j]; // считаем (ij) элемент невязки
            } 
            if (i==j)
                m6--;
            m7 += max(m6, -m6); //сразу норму
        }
        m8 = max(m7, m8);
    }
    return m8;
}

double nev4(double* A, int n) // верхняя оценка определителя
{
    double m6, m7;
    m6 = 0;
    m7 = 1;
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            m6 = max(m6, max(-A[i*n + j], A[i*n + j]));
        }
        m7 *= m6;
    }
    m7 *= n;
    return m7;
}

double nev5(double* A, int n)
{
    double m6, m7;
    m6 = 0;
    m7 = 1;
    for (int i = 0; i < n; i++)
    {
        m6 = 0;
        for (int j = 0; j < n; j++)
        {
            m6 += max(-A[i*n + j], A[i*n + j]);
        }
        m7 = min(m7, m6);
    }
    return m7;
}
double nev10(double* A, int n)
{
    double m6, m7;
    m6 = 0;
    m7 = 1;
    for (int j = 0; j < n; j++)
    {
        m6 = 0;
        for (int i = 0; i < n; i++)
        {
            m6 += max(-A[i*n + j], A[i*n + j]);
        }
        m7 = min(m7, m6);
    }
    return m7;
}

double nev9(double *A, int n)
{
    double m6;
    m6 = 0;
    for (int i = 0; i < n; i++)
    {
        m6 = min(-A[i], A[i]);
    }
    return m6; 
}

double nev11(double *A, int n)
{
    double m6, m7;
    m6 = 0;
    m7 = 0;
    for(int i = 0; i < n; i++)
    {
        m6 = 0;
        for (int j = 0; j < n; j++)
        {
            m6 = max (m6, max(A[j*n+i], -A[j*n+i]));
            //m6 += max(A[j*n+i], -A[j*n+i]);
            cout << "m6 = " << m6 << endl;
        }
        //m6 = m6/n;
        if (i == 0)
            m7 = m6;
        m7 = min(m7, m6);
    }
    return m7;
}
