#include "b.h"

void Error() {
	printf("Incorrect приколы in the terminal line\n");
}

void printMes (int p, int pn, int code, double m)
{
    MPI_Status status;
    if ((pn > 0)&&(pn < p))
    {
        MPI_Send (&pn, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send (&code, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send (&m, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);        
    }
    if (pn == 0)
    {
        printT(1);
        for (int i = 0; i < p; i++)
        {
            if (i > 0)
            {
                MPI_Recv(&pn, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
                MPI_Recv(&code, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(&m, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);
            }
            printf("pn:\t%d\tcode:\t%d\tmessage:\t%f\n", pn, code, m);
        }
        printT(1);
        cout << endl;
    }
}


int fill0(double* matrix, int n, int p, int pn)
{
    int an, bn, tag, bcast;
    double x;
    bcast = 0;

    an = (n/p)*pn + min(pn, n%p);
    bn = (n/p)*(pn + 1) + min((pn + 1), n%p) - 1;
    MPI_Status status;
    for (int i = 0; i < (bn - an + 1)*n; i++)
    {
        MPI_Recv (&x, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        
        tag = status.MPI_TAG;
        if (tag == 0){
            break; }
        matrix[tag-1] = x;
    }
    MPI_Bcast(&bcast, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return bcast;
}
int fill(double* matrix, int n, FILE* fin, int p, int pn) 
{
    int an, bn, bcast;
    int s;
    double x;
    int amn, bmn, pmn;
    s = 0;
    bcast = 0;
    an = (n/p)*pn + min(pn, n%p);
    bn = (n/p)*(pn + 1) + min((pn + 1), n%p) - 1;
    x = 0;
    
    if (pn == 0)
    {
        amn = an;
        bmn = bn;
        pmn = pn;
	while (fscanf(fin, "%lf", &x) == 1) 
        {
        	cout << "in s= " << s << endl;
            if (s >= n * n)
            {
                bcast = -1;
                break;
            }
            if (pmn == 0)
            {
                matrix[s] = x;
            }   
            else 
            {
                MPI_Send (&x, 1, MPI_DOUBLE, pmn, s - amn*n + 1, MPI_COMM_WORLD);
            }
            s++;
            if (!((s >= amn*n)&&(s <= bmn*n + n - 1)))
            {
                pmn ++;
                amn = (n/p)*pmn + min(pmn, n%p);
                bmn = (n/p)*(pmn + 1) + min((pmn + 1), n%p) - 1;
            }
        }
    if ((fgetc(fin) != EOF )) //не только цифры
    {
        bcast = -1; 
    }
	if (s != n * n) //много или мало цифр
    {
		bcast = -1;
        if (fgetc(fin) == EOF)
            bcast = -2;
	}
    if ((bcast != 0)&&(s < n * n))
    {
        for (int i = max(pmn, 1); i < p; i ++)
            MPI_Send (&x, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(&bcast, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    return bcast;
}

int cr(double *A, int n, int k, int p, int pn, FILE* fin)
{
    int an;
    int bn;
    int res;
    MPI_Request request;
    MPI_Status status; 
    res = 0;
    an = (n/p)*pn + min(pn, n%p);
    bn = (n/p)*(pn + 1) + min((pn + 1), n%p) - 1;
    if (pn == 0){
    if ((fin != nullptr)&&(k != 0))
        res = -1;
    if ((k == 0)&&(fin == nullptr))
        res = -1; 
    for (int i = 1; i < p; i++)
    {
        MPI_Isend (&res, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
    } }
    if (pn != 0)
    {
        MPI_Recv (&res, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }
    if (res != 0) return res;

    if (k ==0)
    {
        if (pn != 0)
        {
            res = fill0 (A, n, p, pn);
        }
        if (pn == 0)
        {   
            res = fill (A, n, fin, p, pn);
        }

        if (res != 0)
            return res;
    }
    if (k == 1)
    {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                {
                    if ((i*n + j >= an*n)&&(i*n + j <= bn*n + n - 1))
                    {    A[i*n + j - an*n] = n - max(i+1, j+1) + 1; 
                    }
                }
    }
    if (k == 2)
    {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                {
                    if ((i*n + j >= an*n)&&(i*n + j <= bn*n + n - 1))
                        A[i*n + j - an*n] = max(i+1, j+1);
                }
    }
    if (k == 3)
    {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                {
                    if ((i*n + j >= an*n)&&(i*n + j <= bn*n + n - 1))
                        {
                            A[i*n + j - an*n] = max(i-j, j-i);
                        }
                }
    }
    if (k == 4)
    {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                {
                    if ((i*n + j >= an*n)&&(i*n + j <= bn*n + n - 1))
                        A[i*n + j - an*n] = 1.0/(i + j + 1);
                }
    }
    return 0;
}
int makeE(double* B, int n, int p, int pn)
{
    int an, bn;

    an = (n/p)*pn + min(pn, n%p);
    bn = (n/p)*(pn + 1) + min((pn + 1), n%p) - 1;
    
    if (pn == 0)
        B[0]=1.0;
    for (int i=1; i < n*n; i++)
    {
        if ((i >= an*n)&&(i <= bn*n + n - 1))
            B[i - an*n] = ((i%(n+1))==0); 
    }     

    return 0;
}

void printT(int k)
{
    for (int i = 0; i <= k; i++)
        cout << "---------------------" << endl;
}

void printPM(double *A, int n, int m)
{
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++)
            printf("%e\t", A[i*n + j]);
        printf ("\n"); }
}

int printM(double *A, int n, int m, int p, int pn)
{
    int an, bn;
    MPI_Status status;
    double x;

    an = (n/p)*pn + min(pn, n%p);
    bn = (n/p)*(pn + 1) + min((pn + 1), n%p) - 1;
    //d1 = bn - an + 1;
    ///cout << "a" << pn << " = " << an << endl;
    //cout << "b" << pn << " = " << bn << endl;
    if (bn >= an)
    {//y - глобальные координаты, x - локальные, тогда:
        // y = an*n + x + Ind(x из левой части)*(n*n - d1*n) -> x = y - an*n - Ind(y из левой части)*(n*n - d1*n)
        for (int i = 0; i < min(n, m); i++)
        {
            for (int j = 0; j < min(n, m); j++)
            {
                if ((i*n + j >= an*n ) && (i*n + j <= bn*n + n - 1))
                {
                    if (pn == 0)
                    {
                        printf("%e\t", A[i*n + j - an*n]);
                    }
                    else
                    {
                        MPI_Send(A + i*n + j - an*n, 1, MPI_DOUBLE, 0, i*n + j, MPI_COMM_WORLD);
                    }
                }
                else if (pn == 0)
                {
                    MPI_Recv(&x, 1, MPI_DOUBLE, MPI_ANY_SOURCE, i*n + j, MPI_COMM_WORLD, &status);
                    printf("%e\t", x);
                }
            }
            if (pn == 0)
                cout << endl;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (pn == 0)
        cout << endl << endl << endl;
    return 0;
}


int normalize (double *A, int n, double n11, int p, int pn)
{
    int an, bn, i, j, d1;
    
    an = (n/p)*pn + min(pn, n%p);
    bn = (n/p)*(pn + 1) + min((pn + 1), n%p) - 1;
    d1 = bn - an + 1;
    
    for (i = 0; i < d1; i++)
        for (j = 0; j < n; j++)
            A[i*n + j] /= n11;
        
    return 0;
}
