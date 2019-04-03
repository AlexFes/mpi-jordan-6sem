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
            // a[i*n + j] = max(i*np + r, j);
            // a[i*n + j] = 1.0/(i*np + r + j);
            a[i*n + j] = abs(i*np + r - j);
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
    // Number of rows; each process has matrix shape nLocal x n
    int nLocal = n/np + (r + 1 <= n%np ? 1 : 0);
    double* mainElement = new double[2*np];

    // Global iteration
    for (int iGlobal = 0; iGlobal < n; ++iGlobal) {
        int firstRow = iGlobal%np <= r ? iGlobal/np : iGlobal/np + 1;
        double tmp = fabs(firstRow >= nLocal ? 0 : a[firstRow*n + iGlobal]);
        int tmpRow = firstRow >= nLocal ? -1 : firstRow;
        int tmpRank = rank;

        // Get main element for process
        for (int iLocal = firstRow + 1; iLocal < nLocal; ++iLocal) {
            if (fabs(a[iLocal*n + iGlobal]) > tmp) {
                tmp = fabs(a[iLocal*n + iGlobal]);
                tmpRow = iLocal;
            }
        }

        double mainElementTmp[2] = {tmp, (double) tmpRow};
        MPI_Gather(mainElementTmp, 2, MPI_DOUBLE, mainElement, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Get main element for iteration
        if (r == 0) {
            for (int i = 0; i < 2*np; i+=2) {
                if (fabs(mainElement[i] > tmp)) {
                    tmp = fabs(mainElement[i]);
                    tmpRow = mainElement[i + 1];
                    tmpRank = i;
                }
            }

            if (tmpRank > 0) {

            }

            printf("\n%f  %f\n", mainElement[0], mainElement[0]);
        }
    }
}
