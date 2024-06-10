#include "b.h"

int tr12(double *A, double *A1, int n)
{
    double n11;
    int res;
    makeE(A1, n);  
    n11 = min (nev5(A, n), nev10(A, n));
    normalize(A, n11, n, 2);
    res = tr1(A, A1, n);
    printM(A, n, 0, n);
    if (res == 0)
    {
        tr2(A, A1, n);
        normalize(A1, n11, n, 2);
    }
    if (res == 1)
    {
        return 1;
    }
    return 0;
}


int tr1(double *A, double *A1, int n)
{
    double m1, m2, m3;
    double eps;
    int d, flag, k;
    double *a;
    double *b;
    
    m1 = 0;
    m2 = 0;
    d = 0;
    m3 = 0;
    eps = 1;
    
    while (1.0 + eps/2 > 1.0)
    {
        eps = eps/2;
    }
    eps = eps * 100;

    for (int m = 0; m < n; m++) //A_m = (a_{k*n+m}, m <= k <= n-1)
    {
        flag = 0;
        m1 = 0; //norm(y)
        m2 = 0; // 1/norm(y-m1*e1)
        d = 0; //
        m3 = 0; //scalar product
        
        if (max(A[m * (n + 1)], (-1) * A[m * (n + 1)]) > eps)
        { flag = 1; }
        for (int k = 1; k <= (n - 1 - m); k++)
        {
            if (max(-A[k * n + m * (n + 1)], A[k * n + m * (n + 1)]) > eps)
                flag = 2;
        }

        if (flag == 1)
        {
            if (max(A[m * (n + 1)], (-1) * A[m * (n + 1)]) < eps)
            {
                return 1;
            }
        }
        if (flag == 0)
        {
            return 1;
        } 
        if (flag == 2)
        {
             
        for (a = A + m * (n + 1), k = 0; k <= (n - 1 - m); k++)
        {
            m1 += (*a) * (*a);
            a += n;
        }

        m2 = m1;
        m2 = 2 * m2 - 2 * A[m * n + m] * sqrt(m1); // ||A_m - norm(A_m)*e1||^2 = ||A_m||^2 - 2*norm(A_m)*(A_m^0) + norm(A_m)^2 = 2*norm(A_m)^2 - 2*norm(A_m)*(A_m^0) = 2*norm(A_m)*(norm(A_m) - A_m^0)   
        m1 = sqrt(m1);
        
        if (max(m1, -m1) < eps)
        {
            //cout<<"no A^{-1} or too low accuracy max m1 < eps";
            return 1;
        }
        
        m2 = sqrt(m2);
        m2 = 1/m2;
        
        for (d = 1; d < n-m; d++) //подействовать на все столбцы кроме первого
        {
            
            m3 = 0;
            for (k = 0, a = A + m * (n + 1) - n, b = a + d; k <= (n - 1 - m); k++)
            {
                m3 += (*(a += n)) * m2 * (*(b += n));
            }
            m3 -= m1 * m2 * A[m*(n+1) + d];
            for (k = 0, a = A + m*(n + 1) - n, b = a + d; k <= (n - 1 - m); k++)
            {
                if (k == 0)
                    A[m * (n + 1) + d] = A[m * (n + 1) + d] + 2 * m1 * m2 * m3;
                *(b += n) -= 2 * m2 * m3 * (*(a += n)); 
        
            } 
            //умножить d-ый столбец m-обрезанной матрицы на U(первый столбец)
        }
        for (d = 0; d < n; d++) //подействовать на все столбцы кроме первого
        {
            
            m3 = 0;
            for (k = m, a = A + m * n + m - n, b = A1 + m * n + d - n; k <= (n-1); k++)
            {
                m3 += (*(a += n)) * m2 * (*(b += n)); 
            }
            m3 -= m1 * m2 * A1[m*n + d];
            for (k = m, a = A + m*n + m - n, b = A1 + m*n + d - n; k <= (n-1); k++)
            {
                *(b += n) -= 2 * m2 * m3 * (*(a += n));
                if (k == m)
                    A1[m*n + d] += 2*m1*m2*m3;
                
            }
        }
       
        m3 = 0;
        d = 0;
     
        for (k = 0, a = A + m*(n + 1) - n; k <= (n - 1 - m); k++)
        {
            (*(a += n)) = 0;
        }
        
            A[m*(n+1)] = m1;
        
        }
    }
    return 0;
}

int tr2(double *A, double *A1, int n)
{
    //Обратный Гаусс
    
    double m4, m5; 
    double *a;
    double *b;
    double *c;
    int i;
    int j;
    int k;
    
    for (i = 1; i <= n; i++)
    {
        m5 = A[(n - i) * (n + 1)]; //Диагональный элемент
        for(j = 1, a = A - i; j <= n - i; j++) //Элемент нужного столбца = A[j*n - i]
        {
            m4 = (-1)*(*(a += n))/m5; //коэффициент метода Гаусса для строчки
            *a = 0;
            for (k = 0, b = A1 + (j - 1) * n, c = A1 + (n - i) * n; k < n; k++)
            {
                *(b++) += m4 * (*(c++));
            } 

        }
             
        for (k = 0, a = A1 + (n - i)*n; k < n; k++)
        {
            *(a++) /= m5;
        }
     
        A[(n-i)*(n+1)] = 1; 
    }
    
 
    return 0;
}

