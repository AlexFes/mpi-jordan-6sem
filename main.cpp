#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "methods.h"

int main (int argc, char* argv[]) {
    MPI_Init (&argc, &argv);

    int r;
    int np;
    int n = atoi(argv[1]);

    MPI_Comm_rank (MPI_COMM_WORLD, &r);
    MPI_Comm_size (MPI_COMM_WORLD, &np);

    // Before algorithm
    double* a = new double[(n/np + 1)*n];
    double* b = new double[n/np + 1];

    // After algorithm
    double* a1 = new double[(n/np + 1)*n];
    double* b1 = new double[n/np + 1];

    double res;
    double time;
    char* input = NULL;
    if (argc == 3) {
        input = argv[2];
    }

    for (int i = 0; i < n/np + 1; ++i) {
        b[i] = 0;
        b1[i] = 0;
    }

    if (argc == 3) {
        readMatrix(n, r, np, a, b, a1, b1, input);
    } else {
        createMatrix(n, r, np, a, b, a1, b1);
    }
    //writeMatrix(n, r, np, a, b);

    time = MPI_Wtime();
    solveMatrix(n, r, np, a, b);
    time = MPI_Wtime() - time;

    if (n <= 10) {
        writeMatrix(n, r, np, a, b);      
    }

    res = discrepancy(n, r, np, a1, b1, b);

    if (r == 0) {
        printf("\n res=%e \n", res);
        printf("\n time=%.2f \n", time);
    }

    delete[]a;
    delete[]b;
    delete[]a1;
    delete[]b1;

    MPI_Finalize();
    return 0;
}
