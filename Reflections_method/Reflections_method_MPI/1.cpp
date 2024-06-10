#include "b.h"



int main(int argc, char* argv[])
{
    
    int res;
    int n, r, s, pn, p, an, bn;
    double residual, elapsed;
    double *A;
    double *B;
    double t;
    int la;
    FILE* fin;
    MPI_Status status;
    MPI_Request request;

    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &(pn));
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    if (argc < 4 || argc > 5) {
		if (pn == 0) Error();
		return -1;
	}

	

	if (sscanf(argv[1], "%d", &n) != 1 || n <= 0) {
        if (pn == 0) Error();
		return -1;
	}

	if (sscanf(argv[2], "%d", &r) != 1 || r <= 0 || r>n) {
		if (pn == 0) Error();
		return -1;
	}

	if (sscanf(argv[3], "%d", &s) != 1 || s < 0 || s>4) {
		if (pn == 0) Error();
		return -1;
	}

    if (p > n)
    {
        if ( pn == 0)
        cout << "Nonsense" << endl;
        MPI_Finalize();
        return 0;
    }
    
    an = 0;
    bn = 0;

    an = (n/p)*pn + min(pn, n%p);
    bn = (n/p)*(pn + 1) + min((pn + 1), n%p) - 1;
    
    la = bn - an + 1;
    
    B = new double[n*la]; 
    A = new double[n*la];

    fin = nullptr;

    res = 0;

    if (s == 0)
    {
        if (pn == 0)
        {
            fin = ((s == 0) && (argv[4] != NULL)) ? fopen(argv[4], "r") : nullptr;            
            res = 0;
            if (argv[4] == NULL)
            {
                printf ("Empty file name\n");
                res = 1;
            }
            if (fin == nullptr)
            {
                printf ("The matrix cannot be loaded\n");
                res = 1;
            }
            for (int i = 1; i < p; i++)
            {
                MPI_Send (&res, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Isend (&res, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
            }
        }
        if(pn != 0)
        {
            res = 0;  
            MPI_Recv (&res, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        }
    }
    if (res != 0)
    {
        if (fin != 0)
            fclose(fin);
        delete[] A;
        delete[] B;
        MPI_Finalize();
        return 0;
        cout << "res = " << res << endl;
        cout << "OH SHIT" << endl;
        
    }
    res = 0;
    res = cr(A, n, s, p, pn, fin); //заполнить: в матрицу А данную и единичную, в B скопировать данную.
    cout << "OK1" << endl;
    if (fin != 0)
        fclose(fin);
    if (res != 0)
    {
        if (pn == 0) {
        if (res == -1) cout << "Take this matrix away from my algorhythm!" << endl;
        if (res == -2) cout << "Empty matrix = nice matrix" << endl; }
        delete[] A;
        delete[] B;
        MPI_Finalize();
        return 0;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime(); 
    printM(A, n, r, p, pn);
    MPI_Barrier(MPI_COMM_WORLD);
    res = tr12(A, B, n, p, pn);
    if (res == -1)
    {
        if (pn == 0)
            cout << "det A = 0" << endl;
        t = MPI_Wtime() - t;
        if (pn == 0) printf("Working time is %lf\n", t);
        delete[] A;
        delete[] B;
        MPI_Finalize();
        return 0;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime() - t;
    elapsed = t;
    if ((s == 0)&&(pn == 0))
        fin = ((s == 0) && (argv[4] != NULL)) ? fopen(argv[4], "r") : nullptr;
    res = cr(A, n, s, p, pn, fin);
    residual = nev(B, A, n, p, pn);
    if (pn == 0) {printf("%s : residual = %e elapsed = %.2f s = %d n = %d r = %d p = %d\n",
argv[0], residual, elapsed, s, n, r, p); }
    MPI_Barrier(MPI_COMM_WORLD);
    delete[] A;
    delete[] B;
    if (pn == 0)
    {
        if (fin != nullptr) 
        {
			fclose(fin);
		}
    }
    MPI_Finalize();    
    return 0;
}

