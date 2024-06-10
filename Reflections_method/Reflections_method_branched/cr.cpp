#include "b.h"

void Error() {
	printf("Incorrect data in the terminal line\n");
}

double form(int n, int k, int i, int j)
{
    if (k == 1)
        return n - max(i, j) + 1;
    if (k == 2)
        return (i == j) ? 2 : (((i - j == 1)||(j - i == 1)) ? (-1) : 0);
    if (k == 3)
        return (i == j)&&(i < n) ? 1 : ((i == n)||(j == n) ? min(i, j) : 0);
    if (k == 4)
        return (1.0/(i + j - 1));
    return 0;
}

int fill(double* matrix, int n, int k, FILE* fin) {
 if (k == 0){
		int s = 0;
		double x;
		while (fscanf(fin, "%lf", &x) == 1) {
			matrix[s] = x;
			s++;
		}
		if (s != n * n) {
			return -1;
		}
	}
	return 0;
}

int cr(double *A, int n, int k)
{
    if ((k == 1)||(k == 2)||(k == 3)||(k == 4))
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                A[i * n + j] = form(n, k, i + 1, j + 1);
    for (int i=n*n; i < n*(n+3); i++)
        A[i] = 0;   
    return 0;
}

int printM(double *A, int n, int s, int m)
{
    for (int i=0; i<min(n, m); i++)
    {
        for (int j=0; j<min(n, m); j++)
            printf("%e\t", A[i*n + j + s]);
        cout<<endl;
    }
    cout<<endl<<endl<<endl;
    return 0;
}
int printP(double *A, int n, int s, int m) // P - типо Прямоугольная матрица
{
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            printf("%e\t", A[i*m + j + s]);
        }
        cout<<endl;
    }
    cout<<endl<<endl<<endl;
    return 0;
}

int EigenPrint(double *A,  int n, int m, int k)
{
    for (int i = 0; i < min(n, m); i++)
    {
        cout<<"Value № "<< i + k << " = "<< A[(i + k)*n + i + k] << endl;
    }
    return 0;
}
