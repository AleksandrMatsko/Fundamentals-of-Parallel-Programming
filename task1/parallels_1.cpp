#include <iostream>
#include <malloc.h>
#include <stdlib.h>
#include <mpi.h>

const int N = 200;
const double TAU = 0.00001;
const double EPSILON = 0.001;
const unsigned int MAX_ITERATIONS = 1000000;

const int size_buff = 5;

void Fill_A(double *A) {
    srand(1);
    for (int i = 0; i < N; i++) {
        A[i * N + i] = rand() / (double)RAND_MAX * 50 + 400;
        for (int j = i + 1; j < N; j++) {
            A[i * N + j] = rand() / (double)RAND_MAX * 50;
            A[j * N + i] = A[i * N + j];
        }
    }
}

void Fill_b(double *b) {
    for (int i = 0; i < N; i++) {
        b[i] = N + 1;//rand() / (double)RAND_MAX * 200 - 100;
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
        std::cerr << "Problems with size of matrix" << std::endl;
	return;
    }
    int lines = size / N;
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << matrix[i * N + j] << " ";
        }
        std::cerr << std::endl;
    }
    std::cerr << std::endl;
}

int SumBuff(int *buff) {
    int res = 0;
    for (int i = 0; i < size_buff; i++) {
        res += buff[i];
    }
    return res;
}

void ZeroBuff(int *buff) {
    for (int i = 0; i < size_buff; i++) {
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

    int *strings_in_tasks = (int *) calloc(size, sizeof(int));
    int *shifts = (int *) calloc(size, sizeof(int));

    int num_of_strings = N / size;
    int remainder = N % size;
    for (int i = 0; i < remainder; i++) {
        if (rank == i) {
            num_of_strings += 1;
        }
    }
    MPI_Allgather(&num_of_strings, 1, MPI_INT, strings_in_tasks, 1, MPI_INT, MPI_COMM_WORLD);

    int start_index = 0;
    for (int i = 0; i < size; i++) {
        if (rank == i) {
            start_index = i * num_of_strings * N;
            if (rank >= remainder) {
                start_index += N * remainder;
            }
        }
    }
    MPI_Allgather(&start_index, 1, MPI_INT, shifts, 1, MPI_INT, MPI_COMM_WORLD);

    int *shifts_in_x = (int *)calloc(size, sizeof(int));
    for (int i = 0; i < size; i++) {
	    shifts_in_x[i] = shifts[i] / N;
    } 

    double *A = (double *)calloc(N * num_of_strings, sizeof(double));
    double *b = (double *)calloc(N, sizeof(double));
    double module_b = 0;

    double *full_matrix_A = (double *)calloc(N * N, sizeof(double));
    Fill_A(full_matrix_A);
    
    int *to_send = (int *)calloc(size, sizeof(int));
    for (int i = 0; i < size; i++) {
        to_send[i] = N * strings_in_tasks[i]; 
    }
    MPI_Scatterv(full_matrix_A, to_send, shifts, MPI_DOUBLE, A, N * num_of_strings, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    Fill_b(b);
    module_b = ScalarProduct(b, b, N);
    MPI_Bcast(b, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&module_b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(full_matrix_A);
    free(to_send);
    
    double *x = (double *)calloc(N, sizeof(double));
    int *buffer = (int *)calloc(size_buff, sizeof(int));
    int buff_counter = 0;
    double *tmp = (double *)calloc(num_of_strings, sizeof(double));
    unsigned int iteration = 0;
    double g_x = 0;

    double start_time = MPI_Wtime();
    while (SumBuff(buffer) != size_buff && iteration < MAX_ITERATIONS) {
	    double module_tmp = 0;
        CalcTMP(A, num_of_strings, x, b, shifts_in_x[rank], tmp);
	    double tmp_val = ScalarProduct(tmp, tmp, num_of_strings);
        MPI_Allreduce(&tmp_val, &module_tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	
	    g_x = module_tmp / module_b;
        if (g_x < EPSILON * EPSILON) {
            buffer[buff_counter] = 1;
            buff_counter = (buff_counter + 1) % size_buff;
        }
        else {
            ZeroBuff(buffer);
        }
        if (SumBuff(buffer) == size_buff) {
            break;
        }

        MulScalarAndVector(TAU, tmp, num_of_strings, tmp);
        SubVectors(&x[shifts_in_x[rank]], tmp, num_of_strings, tmp);
        MPI_Allgatherv(tmp, num_of_strings, MPI_DOUBLE, x, strings_in_tasks, shifts_in_x, MPI_DOUBLE, MPI_COMM_WORLD);
	    ZeroVector(tmp, num_of_strings);
        iteration += 1;
    }

    if (rank == 0) {
	std::cerr << "size = " << size << std::endl;
        if (iteration == MAX_ITERATIONS) {
            std::cerr << "solution wasn't find after " << MAX_ITERATIONS << " iterations" << std::endl;
	    std::cerr << "g(x)^2 = " << g_x << std::endl;
        }
        else {
            std::cerr << "solution found after " << iteration << " iteration" << std::endl;
        }
        std::cerr << "Time taken " << MPI_Wtime() - start_time << std::endl;
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
