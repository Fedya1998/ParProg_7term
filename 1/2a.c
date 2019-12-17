#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define IJ i*jsize+j
double *a = NULL;

int main(int argc, char **argv)
{
    if (argc != 3) {
        printf("too few arguments!\n");
        return 0;
    }

    MPI_Init(NULL, NULL);

    int isize = atoi(argv[1]);
    int jsize = atoi(argv[2]);

    int i = 0, j = 0;
    int n_workers = 0;
    int rank = 0;
    double *a = NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_workers);

    int cell_size = isize / n_workers;
    if (cell_size <= 4) {
        printf("Too many workers or too few strings!\n");
        return -1;
    }
    int starti = cell_size * rank;

    if (rank == n_workers - 1) {
        cell_size += isize % n_workers;
    }

    a = (double*) calloc(isize*jsize, sizeof(double));

#define IJ i*jsize+j
    double start_time = MPI_Wtime();
    for (i = starti; i < cell_size + starti + 1 & i < isize; i++) {
        for (j = 0; j < jsize; j++){
            a[IJ] = 10*i + j;
        }
    }

    for (i = starti; i < cell_size + starti & i < isize - 1; i++) {
        for (j = 1; j < jsize; j++) {
            a[IJ] = sin(0.00001*a[(i+1)*jsize + j - 1]);
        }
    }

    int *recvcounts = (int*)calloc(n_workers, sizeof(int));
    int *displ = (int*)calloc(n_workers, sizeof(int));

    if (!rank) {
        for (i = 0; i < n_workers - 1; i++) {
            recvcounts[i] = jsize * cell_size;
            displ[i] = i * cell_size * jsize;
        }
        recvcounts[n_workers-1] = jsize * (cell_size + isize % n_workers);
        displ[n_workers-1] = (n_workers - 1) * cell_size * jsize;
    }
    double end_time = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(&a[starti*jsize], jsize*cell_size, MPI_DOUBLE, a, recvcounts, displ, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        FILE *ff;
        ff = fopen("result.txt", "w");

        for (i = 0; i < isize; i++) {
            for (j = 0; j < jsize; j++) {
                fprintf(ff, "%3.5f ", a[IJ]);
            }
            fprintf(ff, "\n");
        }
        fclose(ff);
    }

    MPI_Finalize();

    if (rank == 0) {
        printf("time: %lg\n", end_time - start_time);
    }

    return 0;
}
