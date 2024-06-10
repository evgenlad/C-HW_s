#include "b.h"

int tr12(double *A, double *B, int n, int p, int pn)
{
    int res;
    double n11;
    makeE(B, n, p, pn);
    n11 = min (nev5(A, n, p, pn), nev10(A, n, p, pn));
    res = normalize(A, n, n11, p, pn);
    res = tr1 (A, B, n, p, pn);
    cout << "OK 2.1" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    if (res < 0 )
    {
        return -1;
    }
    tr2 (A, B, n, p, pn);
    cout << "OK 2.2" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    normalize(B, n, n11, p, pn);  //здесь присоединенная матрица
    return 0;
}

double ab(double a)
{
    if (a > 0)
        return a;
    if (a < 0)
        return (-a);
    return 0;
}

int tr1(double *A, double *B, int n, int p, int pn)
{
    int sn, s1n, fi, flag, flag1, buf, d1; //номер первого элемента данной полосы в матрице шага индукции (размера m*m) 
    double eps, m1, m2, m3, m11, m33;
    double *a;
    double *b;
    MPI_Status status;
    int an;
    int bn;
    eps = 1;
    
    an = (n/p)*pn + min(pn, n%p);
    bn = (n/p)*(pn + 1) + min((pn + 1), n%p) - 1;
    
    while (1.0 + eps/2 > 1.0)
    {
        eps = eps/2;
    }
    eps = eps * 100;
    
    d1 = bn - an + 1;
    fi = 0;
    buf = 0;
    m1 = 0;
    m2 = 0;
    m3 = 0;
    m11 = 0;
    m33 = 0;
    //flag -- индикатор вырожденности
    //flag1 -- индикатор использования этой полосы на данном шаге
    // pn -> sn // sn - первый элементы матрицы шага индóкции в пересечении с кóском компа.
    for (int m = 0; m < n; m++)
    {
    	//if (pn == p-1) cout << (m+1) << endl;
        //cout << "OK 2.1.1" << endl;
        m11 = 0;
        m33 = 0;
        m1 = 0;
        m2 = 0;
        m3 = 0;
        flag = 0;
        flag1 = 0;
        s1n = 0;
        sn = 0;
        if (( m <= bn )&&( d1 > 0)) 
        {
            flag1 = 1; //полоса компа задевает матрицó шага индóкции
            if (m >= an)
            {
                flag1 = 2; //эта полоса - первая в матрице шага индукции
                sn = (m - an)*n + m;
        
            }
            if (an > m)
            {
                sn = m;
            }
        }
        // sn - локальный номер первого элемента данной полосы в матрице шага индóкции
        
        if ((ab(A[sn]) > eps)&&(flag1 == 2))
        {
            flag = 1;
        }
        s1n = sn;
        if (flag1 == 2)
        {
            s1n += n;
        }
        if (flag1 > 0)
        {
            while (s1n < d1 * n)
            {
                if (ab(A[s1n]) > eps)
                {
                    flag = 2;
                }
                s1n += n;
            }
            if (flag1 == 2)
            {
                fi = pn;
                if (pn != 0)
                {
                    MPI_Send (&pn, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);        
                }   
            }
        }
        if (flag1 != 2)
        {
            if (pn == 0)
            {
                MPI_Recv(&fi, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                fi = status.MPI_SOURCE;
            }

        }
        MPI_Bcast(&fi, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (pn > fi)
        {
            MPI_Send(&flag, 1, MPI_INT, fi, 1, MPI_COMM_WORLD);
        }
        if (pn == fi)
        {
            for (int k = fi + 1; k < p; k++)
            {
                MPI_Recv(&buf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                flag += buf;
            }
        }
        MPI_Bcast(&flag, 1, MPI_INT, fi, MPI_COMM_WORLD);
        buf = 0;
        if (flag == 0)
        {
            return -1;
        }
        if (flag == 1)
        {
            //можно скипнуть это m
        }
        
    if (flag >= 2)
    {
        if (flag1 > 0)
        {
            for (a = A + sn, s1n = sn; s1n < d1*n; s1n += n)
            {
                m1 += (*(a)) * (*(a));
                a += n;
            }
        }
        MPI_Allreduce (&m1, &m11, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        m1 = m11;
        if (flag1 == 2)
        {
            m2 = 2*m1 - 2*A[sn]*sqrt(m1);
            m2 = sqrt(m2);
            m2  = 1/m2;
        }
        m1 = sqrt(m11);
        if (ab(m1) < eps)
        {
            return -1;
        }
        MPI_Bcast(&m2, 1, MPI_DOUBLE, fi, MPI_COMM_WORLD);
        for (int d = 1; d < n - m; d++)
        {
            m3 = 0;
            m33 = 0;
            if (flag1 > 0)
            {
                s1n = sn;
                m3 = 0;
                for (a = A + sn, b = A + sn + d, s1n = sn; s1n < d1 * n; s1n += n)
                {
                    if ((flag1 == 2)&&(s1n == sn))
                    {
                        m3 -= m1 * m2 * (*(b));
                    }   
                    m3 += (*(a)) * m2 * (*(b));
                    a += n;
                    b += n;
                }
            }
            MPI_Allreduce (&m3, &m33, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            m3 = m33;
            if (flag1 > 0)
            {
                for (a = A + sn, b = A + sn + d, s1n = sn; s1n < d1 * n; s1n += n)
                {
                    if ((flag1 == 2)&&(s1n == sn))
                    {
                        *(b) += 2 * m1 * m2 * m3;
                    }
                    *(b) -= 2 * m2 * m3 * (*(a));
                    b += n;
                    a += n;
                }
            }        
        }
   
        for (int d = 0; d < n; d++)
        {
            m3 = 0;
            m33 = 0;
            if (flag1 > 0)
            { //b = A + sn - m + d1 * n + d
                for (a = A + sn, b = B + sn - m + d, s1n = sn; s1n < d1 * n; s1n += n) //кажется, здесь b - это присоединенная матрица(фрагмент)
                {
                    if ((flag1 == 2)&&(s1n == sn))
                    {
                        m3 -= m1 * m2 * (*b);
                    }
                    m3 += (*(a)) * m2 * (*(b));
                    a += n;
                    b += n;
                }
            }
            MPI_Allreduce (&m3, &m33, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            m3 = m33;
        
            s1n = sn; 
            if (flag1 > 0)
            { //b = A + sn - m + d1 * n + d
                for (a = A + sn, b = B + sn - m + d, s1n = sn; s1n < d1 * n; s1n += n) //кажется, b - присоединенная.
                {
                    if ((flag1 == 2)&&(s1n == sn))
                    {
                        *b += 2 * m1 * m2 * m3;
                    }
                    *(b) -= 2 * m2 * m3 * (*(a));
                    b += n;
                    a += n;
                    
                }
                
            }
        }
        
        m3 = 0;
        m33 = 0;
        s1n = sn;

        if (flag1 > 0)
        {
            for (a = A + sn, s1n  = sn; s1n < d1 * n; s1n += n)
            {
                *(a) = 0;
                a += n;
            }
            if (flag1 == 2)
            {
                A[sn] = m1;
            }
        }

    }
    }
    return 0;
}

int tr2(double *A, double *B, int n, int p, int pn)
{
    double *buf;
    double *a;
    double *b;
    double *c;
    double *e;
    double r;
    double eps;
    double buf0;
    int an;
    int bn;
    int d1;
    int i;
    int j;
    
    an = (n/p)*pn + min(pn, n%p);
    bn = (n/p)*(pn + 1) + min((pn + 1), n%p) - 1;
    d1 = bn - an + 1;
    eps = 1.e-15;
    
    MPI_Request request, request1;
    MPI_Status status, status1;

    r = 0;
    buf = new double[n];
    
   
    for (int d = n - 1; d >= 0; d--)
    {
    	//if (pn == 0) cout << n - d << endl;
        if (d > bn) //ща скоро пришлют строку ГАУССА, ЛОВИТЕ ЕЁ!
        {

            MPI_Recv (buf, n, MPI_DOUBLE, MPI_ANY_SOURCE, d, MPI_COMM_WORLD, &status);
            MPI_Recv (&buf0, 1, MPI_DOUBLE, MPI_ANY_SOURCE, d + n, MPI_COMM_WORLD, &status1);//Гениально! Представим, что отправилось 

            //куча сообщений в MPI-буфер с тегами = номер шага индукции. Тогда эта невероятно умная MPI_хрень заставит принимать сообщения в нужной последовательности!  
            // Ща как Гаусну тебя, строка!
            if (ab(buf0) < eps ) cout << "ALARM " << endl;
            for (a = A + d, j = 0; j < d1; j++)
            {
                if (ab(*a) > eps)
                {
                    r = (-1) * (*(a))/buf0; //здесь должна быть зависимость от d// j*n + d
                    //b = a - d + d1 * n
                    // b = B
                    for (b = B + j*n, c = buf, i = 0; i < n; i++)
                    {
                        *(b++) += r * (*(c++));
                    }
                    a += n;
                }
            }
        }
        else if (d >= an)
        {
            for (i = pn - 1; i >= 0; i--)
            {
                //A + (d1 + d - an)*n
                MPI_Isend (B + (d - an)*n, n, MPI_DOUBLE, i, d, MPI_COMM_WORLD, &request);
                MPI_Isend (A + (d - an)*n + d, 1, MPI_DOUBLE, i, d + n, MPI_COMM_WORLD, &request1);// ОТПРАВИЛ И ДАЛЬШЕ ПОШЕЛ
            }
            //b = A + d + (d - an) * n
            for (a = A + d, b = A + d + (d - an) * n, j = 0; j < (d - an); j++)
            {
                r = (-1) * (*(a)) / *(b);
                //c = A + d1 * n + j*n, e = A + d1 * n + (d - an) * n
                for (c = B + j*n, e = B + (d - an) * n, i = 0; i < n; i++)
                {
                    *(c++) += r * (*(e++));
                }
                a += n;
            }

        }
    
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (a = A - an*n, i = 0; i < n; i++)
    {
        if ((i*(n+1) >= an*n)&&(i*(n+1) <= bn*n + n - 1))
        {
            r = *a;
            //b = A + d1*n + (i - an)*n
            for (b = B + (i - an)*n, j = 0; j < n; j++)
            {
                *(b++) /= r;
            }
        }
        a += (n + 1);
    }

    delete[] buf;
    
    // добавить проверку того, что сообщения дошли через MPI_REQUEST/ хотя может забить? ладно, пофиг, предположим, что дойдóт
    return 0;
}

//надо добавить оптимизации, идя вдоль по строке
