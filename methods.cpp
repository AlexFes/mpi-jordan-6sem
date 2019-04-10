#include <mpi.h>
#include <stdio.h>
#include <algorithm>

using namespace std;

void createMatrix (int n, int r, int np, double *a, double *b);
void writeMatrix (int n, int r, int np, double *a, double *b);

// Generate matrix shape n x n for np processes
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

// Output matrix to stdout
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

    MPI_Status st;
    double* mainElement = new double[2*np];
    double* buf = new double[n + 1];
    double* buf1 = new double[n + 1];

    // Global iteration
    for (int iGlobal = 0; iGlobal < n; ++iGlobal) {
        int firstRow = iGlobal%np <= r ? iGlobal/np : iGlobal/np + 1;
        double tmp = fabs(firstRow >= nLocal ? 0 : a[firstRow*n + iGlobal]);
        int tmpRow = firstRow >= nLocal ? -1 : firstRow;
        int tmpRank = r;

        // Get main element for process
        for (int iLocal = firstRow + 1; iLocal < nLocal; ++iLocal) {
            if (fabs(a[iLocal*n + iGlobal]) > tmp) {
                tmp = fabs(a[iLocal*n + iGlobal]);
                tmpRow = iLocal;
            }
        }

        double mainElementTmp[2] = {tmp, (double) tmpRow};
        MPI_Allgather(mainElementTmp, 2, MPI_DOUBLE, mainElement, 2, MPI_DOUBLE, MPI_COMM_WORLD);

        // Get main element for iteration
        for (int i = 0; i < 2*np; i+=2) {
            if (fabs(mainElement[i]) > tmp || (fabs(mainElement[i]) >= tmp && fabs(mainElement[i]) <= tmp && tmpRank > i/2)) {
                tmp = fabs(mainElement[i]);
                tmpRow = mainElement[i + 1];
                tmpRank = i/2;
            }
        }

        // Write the main row to buffer (root)
        if (r == tmpRank) {
            memcpy(buf, a + tmpRow*n + iGlobal, (n - iGlobal)*sizeof(double));
            buf[n - iGlobal] = b[tmpRow];
        }

        // Broadcast the main row to all processes
        MPI_Bcast(buf, n - iGlobal + 1, MPI_DOUBLE, tmpRank, MPI_COMM_WORLD);

        // Get replacement for the main row (root)
        if (r == tmpRank && r != iGlobal%np) {
            MPI_Recv(buf1, n - iGlobal + 1, MPI_DOUBLE, iGlobal%np, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
            memcpy(a + tmpRow*n + iGlobal, buf1, (n - iGlobal)*sizeof(double));
            b[tmpRow] = buf1[n - iGlobal];
        } else if (r == tmpRank && r == iGlobal%np) {
            memcpy(buf1, a + (iGlobal/np)*n + iGlobal, (n - iGlobal)*sizeof(double));
            buf1[n - iGlobal] = b[iGlobal/np];
            memcpy(a + (iGlobal/np)*n + iGlobal, buf, (n - iGlobal)*sizeof(double));
            b[iGlobal/np] = buf[n - iGlobal];
            memcpy(a + tmpRow*n + iGlobal, buf1, (n - iGlobal)*sizeof(double));
            b[tmpRow] = buf1[n - iGlobal];
        }

        // Send replacement for the main row to the root
        if (r != tmpRank && r == iGlobal%np) {
            memcpy(buf1, a + (iGlobal/np)*n + iGlobal, (n - iGlobal)*sizeof(double));
            buf1[n - iGlobal] = b[iGlobal/np];
            MPI_Send(buf1, n - iGlobal + 1, MPI_DOUBLE, tmpRank, 0, MPI_COMM_WORLD);
            memcpy(a + (iGlobal/np)*n + iGlobal, buf, (n - iGlobal)*sizeof(double));
            b[iGlobal/np] = buf[n - iGlobal];
        }
    }
}
