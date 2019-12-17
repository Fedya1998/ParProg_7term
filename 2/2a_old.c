#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
//#define ISIZE 40
//#define JSIZE 40

int main(int argc, char **argv)
{
	int ISIZE = atoi(argv[1]);
	int JSIZE = atoi(argv[2]);
	
	double * a = (double *) calloc(ISIZE * JSIZE, sizeof(double));
    int i, j;
    FILE *ff;
    for (i=0; i<ISIZE; i++){
        for (j=0; j<JSIZE; j++){
            a[i * JSIZE + j] = 10*i +j;
        }
    }
    for (i=1; i<ISIZE; i++){
        for (j = 3; j < JSIZE - 1; j++){
            a[i * JSIZE + j] = sin(0.00001*a[(i - 1) * JSIZE + j - 3]);
        }
    }
    ff = fopen("result_old.txt","w");
    for(i=0; i < ISIZE; i++){
        for (j=0; j < JSIZE; j++){
            fprintf(ff,"%3.5f ",a[i * JSIZE + j]);
        }
        fprintf(ff,"\n");
    }
    fclose(ff);
    return 0;
}
