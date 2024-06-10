#include "b.h"

double nev(double *A, double *B, int n, int p, int pn)
{
    int an, bn, apn, bpn, dp1, d, k;
    double *buf;
    double *a, *b;
    double m6, m7, m8, m77;
    
    an = (n/p)*pn + min(pn, n%p);
    bn = (n/p)*(pn + 1) + min((pn + 1), n%p) - 1;
    buf = new double[n];

    m77 = 0;
    m8 = 0;
    m7 = 0;
    for (d = 0; d < n; d++)
    {
        for (int i = an; i <= bn; i++)
        {
            buf[i] = A[(i - an)*n + d]; //d-тый столбец
        }
        for (int j = 0; j < p; j++)
        {
            apn = (n/p)*j + min(j, n%p);
            bpn = (n/p)*(j + 1) + min((j + 1), n%p) - 1;
            dp1 = bpn - apn + 1; 
            MPI_Bcast(buf + apn, dp1, MPI_DOUBLE, j, MPI_COMM_WORLD); // каждый делится своим кусочком d-того столбца 
        }
        for (int i = an; i <= bn; i++)
        {            
            m6 = 0;
            for (k = i*n, b = B + i*n - an*n, a = buf; k < (i+1)*n; k++)
            {
                m6 += (*(b++)) * (*(a++)); //элемент B*A[i, d]
            }
            if (d == i)
            {
               m6 -= 1;
            }
            m7 += max(m6, -m6);
        }
        //printMes(p, pn, 7, m7);
        MPI_Allreduce (&m7, &m77, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        m7 = m77;
        //printMes(p, pn, 77, m77);
        m8 = max (m8, m7);
        //printMes(p, pn, 8, m8);
        m7 = 0;
        m77 = 0;
    }
    delete[] buf;
    return m8;
}



double nev5 (double *A, int n, int p, int pn)
{
    int an, bn, d1;
    double m6, m7, m66;
    m6 = 0;
    m7 = 0;
    m66 = 0;
    
    an = (n/p)*pn + min(pn, n%p);
    bn = (n/p)*(pn + 1) + min((pn + 1), n%p) - 1;
    d1 = bn - an + 1;
    //минимóм максимóмов в каждом столбце
    //было б хорошо конечно сделать грóппó компов с d1 > 0 на слóчай, когда их больше n.. Предположим, что их не больше n, тогда все ок
    
    for (int j = 0; j < n; j++)
    {
        m6 = 0;
        for (int i = 0; i < d1; i++)
        {
            m6 = max (m6, ab(A[i*n + j]));
        }
        //нашли локальный максмóм в этом столбце
        MPI_Allreduce (&m6, &m66, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        m6 = m66;
        if (j == 0)
            m7 = m6;
        else
            m7 = max (m7, m6);
    }
    return m7;
}

double nev10 (double *A, int n, int p, int pn)
{
    int an, bn, d1;
    double m6, m7, m77;
    m6 = 0;
    m7 = 0;
    m77 = 0;
    
    an = (n/p)*pn + min(pn, n%p);
    bn = (n/p)*(pn + 1) + min((pn + 1), n%p) - 1;
    d1 = bn - an + 1;
    //минимóм максимóмов в каждой строке
    
    for (int i = 0; i < d1; i++)
    {
        m6 = 0;
        for (int j = 0; j < n; j++)
        {
            m6 = max (m6, ab(A[i*n + j]));
        }
        //нашли локальный максмóм в этой строке
        if (i == 0)
            m7 = m6;
        else
            m7 = max (m7, m6);
    }
    MPI_Allreduce (&m7, &m77, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    m7 = m77;
    return m7;
}
     
    
