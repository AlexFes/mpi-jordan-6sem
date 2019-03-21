#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "methods.cpp"

int main (int argc, char* argv[]) {
    MPI_Init (&argc, &argv);

    int r, np, n = atoi(argv[1]);

    MPI_Comm_rank (MPI_COMM_WORLD, &r);
    MPI_Comm_size (MPI_COMM_WORLD, &np);

    double* a = new double[(n/np + 1)*n];
    double* b = new double[n/np + 1];

    createMatrix(n, r, np, a, b);
    writeMatrix(n, r, np, a, b);

    MPI_Finalize();
    return 0;
}
