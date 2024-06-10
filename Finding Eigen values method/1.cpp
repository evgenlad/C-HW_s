#include "b.h"



int main(int argc, char* argv[])
{

    double *A;
    int res;
    int n, m, k;
    double epss; 
    double trace1, trace2, sq1, sq2;
    clock_t chasi1;
	double chasi1_;
    clock_t chasi2;
	double chasi2_;

    if (argc < 5) {
		Error();
		return -1;
	}

	

	if (sscanf(argv[1], "%d", &n) != 1 || n <= 0) {
		Error();
		return -1;
	}

	if (sscanf(argv[2], "%d", &m) != 1 || m <= 0 || m > n) {
		Error();
		return -1;
	}
    if (sscanf(argv[3], "%lf", &epss) != 1) {
		Error();
		return -1;
	}
	if (sscanf(argv[4], "%d", &k) != 1 || k < 0 || k>4) {
		Error();
		return -1;
	}
    
    //epss = 1.e-10;
    //cout << "eps = "<<epss<<endl;


    

    A = new double[n*(n+3)];

    FILE* fin = ((k == 0) && (argv[5] != NULL)) ? fopen(argv[5], "r") : NULL;

	if ((k == 0) && (argv[5] == NULL)) {
		printf("The matrix cannot be loaded\n");
		delete[] A;
		return -1;
	}

        if ((fin==NULL)&&(k==0)){
             printf("The matrix cannot be loaded\n");
		delete[] A;
		return -1;
        }

    if (k == 0)
	{
        if (fill(A, n, k, fin)!=0) 
        {
		    printf("The matrix cannot be loaded\n");
		    if (fin != NULL) 
            {
			    fclose(fin);
		    }
		    delete[] A;
		    return -3;
	    }  
    }





    cr(A, n, k);

    trace1 = nev1(A, n);
    sq1 = nev2(A, n);
    cout << "sq1 = " << sq1 << endl;
    res = printM(A, n, 0, m);

    chasi1 = clock();
	chasi1_ = (double)chasi1;

    res = tr12(A, n, epss); 

    chasi2 = clock();
	chasi2_ = (double)chasi2;

    if ((res == 0)||(res == -1))
    {
    printM(A, n, 0, m);
    //printP(A, n+3, 0, n);
    cout<<endl;
    
    if (res == -1)
    {
        cout << A[0] << " + I* " << A[1] << endl;
        cout << A[n+1] << " + I* " << A[n] << endl;
        
        EigenPrint(A, n, m - 2, 2);
        trace2 = nev1(A, n);
        sq2 = nev3(A, n, -1);
        cout << "sq2 = " << sq2 << endl;
    }
    else
    {
        trace2 = nev1(A, n);
        sq2 = nev3(A, n, 1);
        EigenPrint(A, n, m, 0);
    }
    
    cout << "NEV1 = "<< max((trace1 - trace2), (trace2-trace1))<<endl;
    cout<< "NEV2 = " << max((sq1-sq2), (sq2-sq1))<<endl;
    printf("Working time is %lf\n", (chasi2_ - chasi1_)/(CLOCKS_PER_SEC));
    delete[] A; 
    return 0;
    }
    cout << "matrix is degenerate or tricky for algorithm" << endl;
    printM(A, n, 0, m);
    delete[] A;
    printf("Working time is %lf\n", (chasi2_ - chasi1_)/(CLOCKS_PER_SEC)); 
    return 0;
}
