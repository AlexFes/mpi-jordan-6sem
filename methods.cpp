#include <mpi.h>
#include <stdio.h>
#include <algorithm>

using namespace std;

void createMatrix (int n, int r, int np, double *a, double *b);
void writeMatrix (int n, int r, int np, double *a, double *b);

// Auto generates matrix shape n x n for np processes
void createMatrix (int n, int r, int np, double *a, double *b) {
    for (int i = 0; i < n/np + (r + 1 <= n%np ? 1 : 0); ++i) {
        b[i] = 0;

        for (int j = 0; j < n; ++j) {
            a[i*n + j] = max(i*np + r, j);
            // a[i*n + j] = 1.0/(i*np + r + j);
            // a[i*n + j] = abs(i*np + r - j);
            b[i] += j%2 ? 0 : a[i*n + j];
        }
    }
}

// Outputs matrix to terminal
void writeMatrix (int n, int r, int np, double *a, double *b) {
    MPI_Status st;
    int N = min(n, 10);
    double* buf = new double[n + 1];

    if (r == 0) {
        printf("\n\n");

        for (int i = 0; i < N; ++i) {
            if (i%np == r) {
                for (int j = 0; j < N; ++j) {
                    printf(" %f", a[(i/np)*n + j]);
                }

                printf("    %f\n", b[(i/np)]);
            } else {
                MPI_Recv(buf, n + 1, MPI_DOUBLE, i%np, MPI_ANY_TAG, MPI_COMM_WORLD, &st);

                for (int j = 0; j < n; ++j) {
                    printf(" %f", buf[j]);
                }

                printf("    %f\n", buf[n]);
            }
        }

        printf("\n");
    } else {
        for (int i = 0; i < N/np + (r + 1 <= N%np ? 1 : 0); ++i) {
            memcpy(buf, a + i*n, n*sizeof(double));
            buf[n] = b[i];
            MPI_Send(buf, n + 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

    delete[]buf;
    return;
}

// Solve system of linear equations (Gauss–Jordan elimination); array b holds the result vector
void solveMatrix (int n, int r, int np, double* a, double* b) {
    // Numer of rows; each process has matrix shape N x n
    int nLocal = n/np + (r + 1 <= n%np ? 1 : 0);

    for (int iGlobal = 0; iGlobal < n; ++iGlobal) {
        // Get main element
        for (int ii = i; ii < N; ++ii) {

        }
    }
}
