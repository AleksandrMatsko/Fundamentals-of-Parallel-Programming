#include <iostream>
#include <malloc.h>
#include <stdlib.h>
#include <mpi.h>

int N = 200;
double TAU = 0.000001;
double EPSILON = 0.00001;
const unsigned int MAX_ITERATIONS = 1000000;
const int BUFF_SIZE = 5;

void Fill_A(double *A, FILE *in) {
    for (int i = 0; i < N * N; i++) {
        fscanf(in, "%lf", &A[i]);
    }
}

void Fill_b(double *b, FILE *in) {
    for (int i = 0; i < N; i++) {
        fscanf(in, "%lf", &b[i]);
    }
}

double ScalarProduct(double *v, double *u, int len) {
    double res = 0;
    for (int i = 0; i < len; i++) {
        res += v[i] * u[i];
    }
    return res;
}

void MulMatrixVector(double *matrix, double *vector, double *res) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            res[i] += matrix[i * N + j] * vector[j];
        }
    }
}

void SubVectors(double *vector_1, double *vector_2, int len, double *res) {
    for (int i = 0; i < len; i++) {
        res[i] = vector_1[i] - vector_2[i];
    }
}

void MulScalarAndVector(double scalar, double *vector, int len,  double *res) {
    for (int i = 0; i < len; i++) {
        res[i] = vector[i] * scalar;
    }
}

void CalcTMP(double *A, int num_of_strings, double *x, double *b, int shift, double *tmp) {
    for (int i = 0; i < num_of_strings; i++) {
        tmp[i] = ScalarProduct(&A[i * N], x, N);
        tmp[i] -= b[shift + i];
    }
}

double SumCoordinates(double *v, int len) {
    double res = 0;
    for (int i = 0; i < len; i++) {
        res += v[i];
    }
    return res;
}

void PrintMatrix(double *matrix, size_t size) {
    if (size % N != 0) {
        std::cout << "Problems with size of matrix" << std::endl;
	return;
    }
    int lines = size / N;
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << matrix[i * N + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int SumBuff(int *buff) {
    int res = 0;
    for (int i = 0; i < BUFF_SIZE; i++) {
        res += buff[i];
    }
    return res;
}

void ZeroBuff(int *buff) {
    for (int i = 0; i < BUFF_SIZE; i++) {
        buff[i] = 0;
    }
}

void ZeroVector(double *vector, int len) {
    for (int i = 0; i < len; i++) {
        vector[i] = 0;
    }
}

int main(int argc, char **argv) {
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    FILE *in = NULL;
    if (rank == 0) {
	in = fopen("data.txt", "r");
	fscanf(in, "%d %lf %lf", &N, &TAU, &EPSILON);
    }

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&TAU, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&EPSILON, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

    int *strings_in_tasks = (int *) calloc(size, sizeof(int));
    int *shifts = (int *) calloc(size, sizeof(int));

    int num_of_strings = N / size;
    int remainder = N % size;
    for (int i = 0; i < size; i++) {
        strings_in_tasks[i] = num_of_strings;
	if (i < remainder) {
	    strings_in_tasks[i] += 1;
	}
    }
    
    int start_index = 0;
    for (int i = 0; i < size; i++) {
        shifts[i] = start_index;
	start_index += strings_in_tasks[i] * N;
    }

    int *shifts_in_x = (int *)calloc(size, sizeof(int));
    for (int i = 0; i < size; i++) {
	shifts_in_x[i] = shifts[i] / N;
    } 

    double *A = (double *)calloc(N * strings_in_tasks[rank], sizeof(double));
    double *b = (double *)calloc(N, sizeof(double));
    double module_b = 0;

    double *full_matrix_A = (double *)calloc(N * N, sizeof(double));
    
    if (rank == 0) {
	Fill_A(full_matrix_A, in);
	Fill_b(b, in);
	fclose(in);
    }
    
    int *to_send = (int *)calloc(size, sizeof(int));
    for (int i = 0; i < size; i++) {
        to_send[i] = N * strings_in_tasks[i]; 
    }
    MPI_Scatterv(full_matrix_A, to_send, shifts, MPI_DOUBLE, A, N * strings_in_tasks[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    free(full_matrix_A);
    free(to_send);
    
    double *x = (double *)calloc(N, sizeof(double));
    int *buffer = (int *)calloc(BUFF_SIZE, sizeof(int));
    int buff_counter = 0;
    double *tmp = (double *)calloc(strings_in_tasks[rank], sizeof(double));
    unsigned int iteration = 0;
    double g_x = 0;

    double start_time = MPI_Wtime();
    double part_module_b = ScalarProduct(&b[shifts_in_x[rank]], &b[shifts_in_x[rank]], strings_in_tasks[rank]);
    MPI_Allreduce(&part_module_b, &module_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double exit_val = EPSILON * EPSILON * module_b; 

    while (iteration < MAX_ITERATIONS) {
	double module_tmp = 0;
        CalcTMP(A, strings_in_tasks[rank], x, b, shifts_in_x[rank], tmp);
	double tmp_val = ScalarProduct(tmp, tmp, strings_in_tasks[rank]);
        MPI_Allreduce(&tmp_val, &module_tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	g_x = module_tmp / module_b;
        if (module_tmp < exit_val) {
            buffer[buff_counter] = 1;
            buff_counter = (buff_counter + 1) % BUFF_SIZE;
        }
        else {
            ZeroBuff(buffer);
        }
        if (SumBuff(buffer) == BUFF_SIZE) {
            break;
        }

        MulScalarAndVector(TAU, tmp, strings_in_tasks[rank], tmp);
        SubVectors(&x[shifts_in_x[rank]], tmp, strings_in_tasks[rank], tmp);
        MPI_Allgatherv(tmp, strings_in_tasks[rank], MPI_DOUBLE, x, strings_in_tasks, shifts_in_x, MPI_DOUBLE, MPI_COMM_WORLD);
	ZeroVector(tmp, strings_in_tasks[rank]);
        iteration += 1;
    }
    
    double end_time = MPI_Wtime();
    if (rank == 0) {
	std::cout << "Parallel 1" << std::endl;
	std::cout << "size = " << size << std::endl;
        if (iteration == MAX_ITERATIONS) {
            std::cout << "solution wasn't find after " << MAX_ITERATIONS << " iterations" << std::endl;
	    std::cout << "g(x)^2 = " << g_x << std::endl;
        }
        else {
            std::cout << "!!!solution found after " << iteration << " iteration!!!" << std::endl;
        }
        std::cout << "Time taken " << end_time - start_time << std::endl << std::endl;
    }

    free(shifts);
    free(shifts_in_x);
    free(strings_in_tasks);
    free(buffer);
    free(tmp);
    free(A);
    free(x);
    free(b);
    MPI_Finalize();
    return 0;
}
