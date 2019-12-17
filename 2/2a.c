#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

double *a = NULL;

int diag_index_from_norm(int i, int j, int isize, int jsize, int diag0);
int need_sin_on_element(int index, int isize, int jsize, int diag0);

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
	double *b = NULL;
    double *a = NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_workers);
	
	int n_in_snake = 4 * jsize; // because 3 * jsize + levels
	int snakes = jsize * isize / n_in_snake;
	int left_without_snake = jsize * isize % n_in_snake;
	int my_snakes = snakes / n_workers;
	int snakes_left = snakes % n_workers;
	int my_elements = my_snakes * n_in_snake;
	if (!rank)
		my_elements += left_without_snake + snakes_left * n_in_snake;

	if (snakes_left)
		if (!rank)
			my_snakes += snakes_left;
	

    a = (double*) calloc(2*isize*jsize, sizeof(double));
	

	int pos = 0;
	int level = 0;
	if (rank) 
		level = jsize - my_snakes * (n_workers - rank) * 4;
	
	int levels_to_jump = 0;
	int starti = 0;
	if (rank)
		starti = left_without_snake + snakes_left * n_in_snake
				+ (rank) * my_snakes * n_in_snake;
/*
	int init_pos;
	if (jsize % 3 == 0)
		init_pos = jsize - 3;
	else
		init_pos = jsize - (jsize % 3);

	pos = init_pos;
	if (init_pos % 3 == 0)
		level -= 3;
	else
		level -= init_pos % 3;*/
	int level0 = level;
	int n_in_diags[3];
	n_in_diags[0] = jsize / 3;
	if (jsize % 3)
		n_in_diags[0]++;
	n_in_diags[1] = jsize / 3;
	if (jsize % 3 == 2)
		n_in_diags[1]++;
	n_in_diags[2] = jsize / 3;
		

	//printf("My rank %i, my snakes %i, my elements %i starti %i\n-------------------\n", rank, my_snakes, my_elements, starti);
    double start_time = MPI_Wtime();
	for (int i_in_arr = starti; i_in_arr < starti + my_elements; i_in_arr++) {
		int index = i_in_arr;
		/*printf("first %i %i %i\n", i, level % jsize, pos);
		if (level >= isize & rank == 0)
			level -= isize;
		if (level < 0 & rank == 0)
			level += isize;

		a[i] = 10 * level + (jsize - pos - 1);
	
		if (pos - 3 >= 0) {
			level++;
			levels_to_jump++;
			pos -= 3;
		}
		else {
			try_jump_back:
			if (pos + 3 * levels_to_jump + 1 < jsize) {
				level -= levels_to_jump;
				pos += 3 * levels_to_jump + 1;
				levels_to_jump = 0;
			}
			else {
				levels_to_jump -= 1;
				goto try_jump_back;
			}
		}*/

		int n_in_line = jsize;
		int line = index / n_in_line;
		int index_in_line = index % n_in_line;
	
		// a diag may contain jsize / 3 or jsize / 3 + 1 elements
		int left = jsize % 3;
		int diag_in_line;
		int in_diag;
	
		int k = 0;
		for (k = 0; k < 3; k++) {
			if (index_in_line - n_in_diags[k] < 0) {
				in_diag = index_in_line;
				diag_in_line = k;
				/*
				if (jsize % 3 == 0)
					diag_in_line = k;
				if (jsize % 3 == 1) {
					if (k == 0)
						diag_in_line = 0;
					if (k == 1)
				}
					diag_in_line = ;
				
				if (jsize % 3 == 2)
					diag_in_line = (k + 1) % 3;
				*/
				break;
			}
			else {
				index_in_line -= n_in_diags[k];
			}
		}
	
		
		int pos = (n_in_diags[k] - in_diag - 1) * 3 + diag_in_line;
		int j = jsize - pos - 1;
		int i = line - (n_in_diags[k] - in_diag - 1);	
		
		if (i < 0)
			i += isize;
		if (i >= isize)
			i -= isize;
		//printf("first %i %i %i\n", i_in_arr, i % jsize, pos);
		a[i_in_arr] = 10 * i + (jsize - pos - 1);
	
    }
	//sleep(1000);
    for (i = starti; i < starti + my_elements; i++) {
        if (need_sin_on_element(i, isize, jsize, level0)) {
			if (i < 1 || i > isize * jsize)
				printf("Pizdec in sinus\n\n-----------\n\n");
			a[i] = sin(0.00001 * a[i - 1]);
		}
    }

    int *recvcounts = (int*)calloc(n_workers, sizeof(int));
    int *displ = (int*)calloc(n_workers, sizeof(int));

    if (1 || !rank) {
        for (i = 1; i < n_workers; i++) {
            recvcounts[i] = snakes / n_workers * n_in_snake;
			//recvcounts[i] = left_without_snake + (my_snakes + snakes_left) * n_in_snake;
            displ[i] = left_without_snake + snakes_left * n_in_snake
			                + i * snakes / n_workers * n_in_snake;
        }
    	recvcounts[0] = left_without_snake + (snakes / n_workers + snakes_left) * n_in_snake;
    	displ[0] = 0;
	}
	double end_time = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(&a[starti], recvcounts[rank], MPI_DOUBLE, a, recvcounts, displ, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        FILE *ff;
        ff = fopen("result.txt", "w");

        for (i = 0; i < isize; i++) {
            for (j = 0; j < jsize; j++) {
				int index = diag_index_from_norm(i, j, isize, jsize, level0);
				if (index > isize * jsize - 1 || index < 0)
					printf("PIZDEC!!!!!!\n----------------------------\n");

                fprintf(ff, "%3.5f ", a[index]);
            }
            fprintf(ff, "\n");
        }
		//getchar();
        fclose(ff);
		//getchar();
    }
    MPI_Finalize();

    if (rank == 0) {
        printf("time: %lg\n", end_time - start_time);
    }

    return 0;
}

int diag_index_from_norm(int i, int j, int isize, int jsize, int diag0) {
	int pos = jsize - j - 1;
	int level_adj = pos / 3;
	int line = (i + level_adj) % isize;
	int diag_in_line;
	diag_in_line = pos % 3;
		
	int n_in_line = jsize;
	
	int prev_in_diags = 0;



	while(diag_in_line--) {
		int n_in_diag = jsize / 3;
		if ((jsize % 3) > diag_in_line) {
			n_in_diag++;
		}
		prev_in_diags += n_in_diag;
	}
	int i_in_diag = j / 3;
	int index = line * n_in_line + prev_in_diags + i_in_diag;
	//printf("index %i from i%i pos%i    line %i p_i_d %i i_i_d %i\n", index, i, pos, line, prev_in_diags, i_in_diag);
	//printf("index %i from i%i pos%i j %i\n", index, i, pos, j);
	return index;
}


int need_sin_on_element(int index, int isize, int jsize, int diag0) {
	int n_in_line = jsize;
	int line = index / n_in_line;
	int index_in_line = index % n_in_line;

	// a diag may contain jsize / 3 or jsize / 3 + 1 elements
	int left = jsize % 3;
	int diag_in_line;
	int in_diag;

	int n_in_diags[3];
	n_in_diags[0] = jsize / 3;
	if (jsize % 3)
		n_in_diags[0]++;
	n_in_diags[1] = jsize / 3;
	if (jsize % 3 == 2)
		n_in_diags[1]++;
	n_in_diags[2] = jsize / 3;
	
	int k;
	for (k = 0; k < 3; k++) {
		if (index_in_line - n_in_diags[k] < 0) {
			in_diag = index_in_line;
			diag_in_line = k;
			/*
			if (jsize % 3 == 0)
				diag_in_line = k;
			if (jsize % 3 == 1) {
				if (k == 0)
					diag_in_line = 0;
				if (k == 1)
			}
				diag_in_line = ;
			
			if (jsize % 3 == 2)
				diag_in_line = (k + 1) % 3;
			*/
			break;
		}
		else {
			index_in_line -= n_in_diags[k];
		}
	}
	

	
	int pos = (n_in_diags[k] - in_diag - 1) * 3 + diag_in_line;
	int j = jsize - pos - 1;
	int i = line - (n_in_diags[k] - in_diag - 1);	
	
	if (i < 0)
		i += isize;
	if (i >= isize)
		i -= isize;
	// printf("index %i %i %i diag0 %i line %i k %i n_i_d %i i_i_d %i\n", index, i, pos, diag0, line, k, n_in_diags[k], in_diag);

	//printf("index %i %i %i\n", index, i, pos);
	/*
	if (j > 2 & i > 0 & j < jsize - 1)
		printf(" sin ");
	puts("");*/
	return j > 2 & i > 0 & j < jsize - 1;
}





