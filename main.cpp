#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "methods.cpp"

int main (int argc, char* argv[]) {
    MPI_Init (&argc, &argv);

    int r;
    int np;
    int n = atoi(argv[1]);

    MPI_Comm_rank (MPI_COMM_WORLD, &r);
    MPI_Comm_size (MPI_COMM_WORLD, &np);

    double* a = new double[(n/np + 1)*n];
    double* b = new double[n/np + 1];
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
        readMatrix(n, r, np, a, b, input);
    } else {
        createMatrix(n, r, np, a, b);
    }
    writeMatrix(n, r, np, a, b);

    time = MPI_Wtime();
    solveMatrix(n, r, np, a, b);
    time = MPI_Wtime() - time;

    if (argc == 3) {
        readMatrix(n, r, np, a, b1, input);
    } else {
        createMatrix(n, r, np, a, b1);
    }
    res = discrepancy(n, r, np, a, b1, b);

    if (r == 0) {
        printf("\n res=%e \n", res);
        printf("\n time=%.2f \n", time);
    }

    delete[]a;
    delete[]b;
    delete[]b1;

    MPI_Finalize();
    return 0;
}
