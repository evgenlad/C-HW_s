#include "b.h"

int tr12 (double *A, int n, double epss)
{
    int res;
    res = tr1(A, n);
    //printM(A, n, 0, m);
    res += ev(A, n, epss);
    return res;
}

int tr1(double *A, int n) //Почти треугольная
{
    double cos;
    double sin;
    double l;
    double temp1, temp2;
    int flag;
    flag = 1;
    double eps;
    eps = 1;
    while (1.0 + eps/2 > 1.0)
    {
        eps = eps/2;
    }

    for (int k = 0; k < n-2; k++) //количество столбцов/строк для изменения
    {
        for (int j = k + 2; j < n; j++) //поворотов в столбце/строке
        {
            flag = 1;
            l = sqrt(A[(k+1)*n+k]*A[(k+1)*n+k] + A[j*n+k]*A[j*n+k]); //длина 
            if (l < eps)
            {    flag = 0; }
            cos = A[(k+1)*n+k]/l; //Угол поворота
            sin = (-1)*A[j*n+k]/l;
            
            if (flag)
            {
                for(int m = 0; m < n-k; m++)
                {
                    
                    temp1 = A[(k+1)*n+k+m]; //Значения двух элементов столбца, которые поворачиваем
                    temp2 = A[j*n+k+m];
       
                    A[(k+1)*n+k+m] = cos*temp1 - sin*temp2; // Его величество поворот
                    A[j*n+k+m] = sin*temp1 + cos*temp2;
                    
                }   
                for(int m = 0; m < n-k; m++)
                {
                    //cout<<"------------"<<endl;
                    //cout<< " Поворот в строке : ПОИХАЛЫ" << endl;
                    //printM(A, n, 0, n);
                    
                    temp1 = A[(k+m)*n + k + 1]; //Значения двух элементов строки, которые поворачиваем
                    temp2 = A[(k+m)*n + j];
                    //cout<< " Поворот в строке : ГОТОВО" << endl;
                    //printM(A, n, 0, n);
                    A[(k+m)*n + k + 1] = temp1*cos - temp2*sin;
                    A[(k+m)*n + j] = temp1*sin + temp2*cos;
                }
            }
        } 
    }
    return 0;
}

int LR(double *A, int n, int m)
{
    double eps;
    eps = 1;
    
    while (1.0 + eps/2 > 1.0)
    {
        eps = eps/2;
    }
    for (int k = 0; k < m; k++) // Заполняется первая строка R
    {
        //cout<< "Заполняется нулевая строка R" <<endl;
        //r_{ik} = A[k*n + i + 2n]
        //r_{0k} = A[k*n + 2n]
        A[k*n + 2*n] = A[k];
    }
    //printP(A, n+3, 0, n);
    for (int s = 1; s < m; s++) //s -- Номер строчки R и диаг. элемента L, которые заполняем
    {
        //l_{s, s-1} = A[(n+2)*n + s]
        //l_{s, s-1} := a_{s, s-1}/r_{s-1, s-1} 
        if (max(A[(s-1)*n + s-1 + 2*n], -A[(s-1)*n + s-1 + 2*n]) < eps)
        {
            //cout<<"dividing by zero in LR, modify A - s" <<endl;
            return 1;
        }   
        A[(n+2)*n + s] = A[s*n + s-1]/A[(s-1)*n + s-1 + 2*n]; //добавили элемент в почти-диагональ L      
        for (int k = s; k < m; k++) //Заполнение s-той R строчки в нумерации s <= k <= m-1
        {
            //cout<< "Заполняется "<<s << " строка R" <<endl;
            //r_{sk} := a_{sk} - l_{s, s-1}*r_{s-1, k}, k>=s
            A[k*n + s + 2*n] = A[s*n + k] - A[(n+2)*n + s]*A[k*n + s - 1 + 2*n];
        }
    }   
    
    for (int i = 0; i < m; i++) //Вычисляем произведение RL. Треугольная часть без последнего столбца.
    {
        //a_{ik} = r_{ik} + r_{i, k+1}*l_{k+1, k}
        for (int k = i; k < m - 1; k++)
        {
            A[i*n + k] = A[k*n + i + 2*n] + A[(k+1)*n + i + 2*n]*A[n*(n+2) + k+1];
        }
    }
    for (int i = 0; i < m; i++) //Вычисляем последний столбец k=n-1
    {
        A[i*n + m - 1] = A[(m-1)*n + i + 2*n]; 
    }
    for (int i = 1; i < m; i++) //Вычисляем почти-диагональ
    {
        A[i*n + i - 1] = A[i*n + i + 2*n]*A[n*(n+2) + i];
    }
    //cout<< "LR-> RL done" <<endl;
    //printP(A, n+3, 0, n);
    return 0;
}

int IterLR(double *A, int n, int m, double epss)
{
    int res, res1;
    double s; 
    double norm;
    double eps;
    eps = 1;
    
    norm = nev(A, n, m);
    while (1.0 + eps/2 > 1.0)
    {
        eps = eps/2;
    }

    res = 0;
    if (norm < eps)
        return 1;
    while (max(A[n*(m-1) + m-2], -A[n*(m-1) + m-2]) > epss*norm)
    {
        if (res == 0)
        {
            s = A[(m-1)*n + m-1];
    
            for (int i = 0; i < m; i++)
            {
                A[i*n + i] -= s; 
            }
        }
        if (res >= 1)
        {
            s = A[(m-1)*n + m-1]*res + eps*100*res;
    
            for (int i = 0; i < m; i++)
            {
                A[i*n + i] -= s; 
            }
        }

        for (int i = 0; i < m; i++)
        {
            A[i*n + i] += s; 
        }   
        res1 = LR(A, n, m);
        if (res1 == 0)
            res = 0;
        else res++;
        if(res == 16500) // << ..хватит уже терпеть! >>
            return 1;

    }
 
    return 0;
}

int ev(double *A, int n, double epss)
{
    double trace;
    double det;
    double D;
    int res;
    for(int k = n; k > 2; k--)
    {
        res = IterLR(A, n, k, epss);
        //cout<< "value found" <<endl;
        //printM(A, n, 0, n);
        if (res == 1)
        {
            cout<<"LR error"<<endl;
            return 1;
        }
    }
    
    trace = A[0] + A[n + 1];
    det = A[0]*A[n+1] - A[1]*A[n];
    if (trace*trace - 4*det < 0)
    {
        D = sqrt(-trace*trace + 4*det);
        A[0] = (trace)/2;
        A[n+1] = (trace)/2;
        A[1] = -D/2;
        A[n] = D/2;
        cout<<"disc is negative"<<endl;
        /*for (int i = 1; i < n; i++)
        {
            for(int k = i-1; k < i; k++)
            {
                A[n*i + k] = 0;
            }
        }*/
        return -1;
    }
    D = sqrt(trace*trace - 4*det);
    A[0] = (trace + D)/2;
    A[n+1] = (trace - D)/2;
    A[1] = 0;
    A[n] = 0;
    for (int i = 2; i < n; i++)
    {
        for(int k = 0; k < i-1; k++)
        {
            A[n*i + k] = 0;
        }
    }

    //printM(A, n, 0 , m);
    return 0;
}
