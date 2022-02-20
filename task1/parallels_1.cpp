#include <iostream>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <mpi.h>

const int N = 4;
const double TAU = 0.001;
const double EPSILON = 0.00001;
const unsigned int MAX_ITERATIONS = 100000;

const int size_buff = 5;

void Fill_A(double *A) {
    srand(1);
    for (int i = 0; i < N; i++) {
        A[i * N + i] = rand() / (double)RAND_MAX * 200 + 200;
        for (int j = i + 1; j < N; j++) {
            A[i * N + j] = rand() / (double)RAND_MAX * 200 - 100;
            A[j * N + i] = A[i * N + j];
        }
    }
}

void FillVector(double *x) {
    for (int i = 0; i < N; i++) {
        x[i] = rand() / (double)RAND_MAX * 200 - 100;
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

double CalcModuleOfVector(double *vector) {
    double res = 0;
    for (int i = 0; i < N; i++) {
        res = vector[i] * vector[i];
    }
    res = sqrt(res);
    return res;
}

void SubVectors(const double *vector_1, const double *vector_2, int len, double *res) {
    for (int i = 0; i < len; i++) {
        res[i] = vector_1[i] - vector_2[i];
    }
}

void MulScalarAndVector(double scalar, const double *vector, int len,  double *res) {
    for (int i = 0; i < len; i++) {
        res[i] = vector[i] * scalar;
    }
}

void CalcTMP(double *A, int num_of_strings, double *x, double *b, double *tmp) {
    for (int i = 0; i < num_of_strings; i++) {
        tmp[i] = ScalarProduct(&A[i * N], x, N);
        tmp[i] -= b[i];
    }
}

double SumCoordinates(double *v, int len) {
    double res = 0;
    for (int i = 0; i < len; i++) {
        res += v[i];
    }
    return res;
}

/*void PrintMatrix(double *matrix, size_t size) {
    if (size % N != 0) {
        std::cerr << "Problems with size of matrix" << std::endl;
    }
    int lines = size / N;
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << matrix[i * N + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}*/

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
    strings_in_tasks[rank] = num_of_strings;
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
    shifts[rank] = start_index;
    MPI_Allgather(&start_index, 1, MPI_INT, shifts, 1, MPI_INT, MPI_COMM_WORLD);

    double *A = (double *)calloc(N * num_of_strings, sizeof(double));
    double *b = (double *)calloc(N, sizeof(double));
    double module_b = 0;

    if (rank == 0) {
        double *full_matrix_A = (double *)calloc(N * (N + 1), sizeof(double));
        Fill_A(full_matrix_A);
        MPI_Scatterv(full_matrix_A, strings_in_tasks, shifts, MPI_DOUBLE, A, sizeof(A), MPI_DOUBLE, rank, MPI_COMM_WORLD);
        //memcpy(A, full_matrix_A, sizeof(double) * N * num_of_strings);
        FillVector(&full_matrix_A[N * N]);
        double module_B = CalcModuleOfVector(&full_matrix_A[N * N]);
        MPI_Scatter(&full_matrix_A[N * N], N, MPI_DOUBLE, b, N, MPI_DOUBLE, rank, MPI_COMM_WORLD);
        MPI_Scatter(&module_B, 1, MPI_DOUBLE, &module_b, 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
        //module_b = module_B;
        //memcpy(b, B, N * sizeof(double));
        free(full_matrix_A);
    }

    double *x = (double *)calloc(N, sizeof(double));

    //std::cout << "vector x_0:" << std::endl;
    //PrintMatrix(x, N);


    int *buffer = (int *)calloc(size_buff, sizeof(int));
    int buff_counter = 0;

    //double *to_calc_exit_condition = (double *) calloc(size, sizeof(double));
    double *tmp = (double *)calloc(num_of_strings, sizeof(double));

    unsigned int iteration = 0;

    double start_time = MPI_Wtime();
    while (SumBuff(buffer) != size_buff && iteration < MAX_ITERATIONS) {
        double module_tmp = 0;
        CalcTMP(A, num_of_strings, x, &b[start_index / N], tmp);
        double tmp_val = ScalarProduct(tmp, tmp, num_of_strings);
        MPI_Allreduce(&tmp_val, &module_tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        module_tmp = sqrt(module_tmp);

        if (module_tmp / module_b < EPSILON) {
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
        SubVectors(&x[start_index / N], tmp, num_of_strings, tmp);
        MPI_Allgatherv(tmp, num_of_strings, MPI_DOUBLE, x, strings_in_tasks, shifts, MPI_DOUBLE, MPI_COMM_WORLD);
        ZeroVector(tmp, num_of_strings);
        iteration += 1;
    }

    if (rank == 0) {
        if (iteration == MAX_ITERATIONS) {
            std::cout << "solution wasn't find after " << MAX_ITERATIONS << "iterations" << std::endl;
        }
        else {
            std::cout << "solution found after " << iteration << "iteration" << std::endl;
        }
        std::cout << "Time taken" << MPI_Wtime() - start_time << std::endl;
    }


    /*std::cout << std::endl << "result:" << std::endl;
    for (int i = 0; i < N; i++) {
        std::cout << x[i] << " ";
    }
    std::cout << std::endl;*/

    free(shifts);
    free(strings_in_tasks);
    free(buffer);
    free(tmp);
    //(to_calc_exit_condition);
    free(A);
    free(x);
    free(b);

    MPI_Finalize();
    return 0;
}

