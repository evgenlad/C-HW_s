#include "b.h"

int main(int argc, char* argv[])
{
    double *A, *A1, *b, *x, *x1;
    int res;
	double y;
    int s = 0;

    if (argc < 4) {
		Error();
		return -1;
	}
    A = 0;
    x1 = 0;
    x = 0;
    A1 = 0;
    b = 0;

	int n, m, k;

	if (sscanf(argv[1], "%d", &n) != 1 || n <= 0) {
		Error();
		return -1;
	}

	if (sscanf(argv[2], "%d", &m) != 1 || m <= 0 || m > n) {
		Error();
		return -1;
	}

	if (sscanf(argv[3], "%d", &k) != 1 || k < 0 || k>4) {
		Error();
		return -1;
	}
    FILE* fin = ((k == 0) && (argv[4] != NULL)) ? fopen(argv[4], "r") : nullptr;
	if ((k == 0) && (argv[4] == NULL)) 
    {
		cout << "stupid file name" << endl;
		delete[] A;
		return -1;
	}

    if ((fin==nullptr)&&(k==0))
    {
        cout << "disgusting file data, can't open" << endl;
		delete[] A;
		return -1;
    }
    try{
        A = new double[n*n]; 
    }
    catch(...){
        cout << "buy some memory pls" << endl;
        delete[] A;
        return -5;
    }
    res = cr(A, n, k, fin); 
    if (res != 0)
    {
        if (res == -3)
        {
            cout << "cannot fill your stupid data" << endl;
		    delete[] A;
		    return -3;
        }
        if (res == -1)
        {

            Error();
            return -1;
        }
    }
    if (fin != nullptr) 
    {
		fclose(fin);
	}
    //A1 - адрес, в который погрузится обратная матрица
    try{
        A1 = new double[n*n]; 
    }
    catch(...){
        cout << "buy some memory pls" << endl;
        delete[] A;
        delete[] A1;
        return -5;
    }
    printM(A, n, 0, m);

    clock_t chasi1 = clock();
	double chasi1_ = (double)chasi1;
    
    res = tr12(A, A1, n);
    
    clock_t chasi2 = clock();
    double chasi2_ = (double)chasi2;
    
    if (res != 0)
    {
        cout << "Take this bad matrix away from my algorithm!" << endl;
    }
    
    if (res == 0)
    {
        printM(A1, n, 0, m);
        
        FILE* fin = ((k == 0) && (argv[4] != NULL)) ? fopen(argv[4], "r") : NULL;

	    if ((k == 0) && (argv[4] == NULL)) 
        {
		    cout << "stupid file name" << endl;
		    delete[] A;
		    return -1;
	    }

        if ((fin==NULL)&&(k==0))
        {
            cout << "disgusting matrix data" << endl;
		    delete[] A;
		    return -1;
        }

        cr(A, n, k, fin);
        if (fin != NULL) 
        {
			fclose(fin);
		}
        cout << "Невязочка" << nev0(A1, A, n) << endl;
        
        //столбец правой части
        FILE* fin1 = ((k == 0) && (argv[4] != NULL)) ? fopen("b.txt", "r") : nullptr;
        b = new double[n];
        x = new double[n];
	    while (fscanf(fin1, "%lf", &y) == 1) 
        {
            b[s] = y;
            s++;
        }

        mult_m_v(A1, b, x, n);
        
        x1 = new double[n];
        mult_m_v(A, x, x1, n);
        printf("Working time is %lf\n", (chasi2_ - chasi1_)/(CLOCKS_PER_SEC));
        y = 0;
        for (int i = 0; i < n; i++)
        {
            //y = x1[i] - b[i];
            //cout << y << endl;
            y += max(x1[i] - b[i], b[i] - x1[i]);
        }
        delete[] x1;
        cout << y << endl;
        

        FILE* fout = fopen("/home/evgen/Documents/progr/p/num_m_spring/x.txt", "w");
        for (int i = 0; i < n; i++)
            fprintf(fout, "%lf\n", x[i]);
            
        FILE* fout1 = fopen("/home/evgen/Documents/progr/p/num_m_spring/a-.txt", "w");
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                fprintf(fout1, "%lf ", A1[i*n + j]);
            fprintf(fout1, "\n");
        }     
        
        fclose(fout1);
        fclose(fin1);
        fclose(fout);

    }

    printf("Working time is %lf\n", (chasi2_ - chasi1_)/(CLOCKS_PER_SEC));
    delete[] A;
    delete[] b;
    delete[] x;
    delete[] A1;
    
    return 0;
}
