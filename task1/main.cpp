#include <iostream>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>

const int N = 200;
const double TAU = 0.000000002;
const double EPSILON = 0.000001;
const unsigned int MAX_ITERATIONS = 1000000;

const int size_buff = 5;

void Fill_A(double *A) {
    srand(2);
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
        b[i] = N + 1;//rand() / (double)RAND_MAX * 20 - 10;
    }
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

void SubVectors(const double *vector_1, const double *vector_2, double *res) {
    for (int i = 0; i < N; i++) {
        res[i] = vector_1[i] - vector_2[i];
    }
}

void MulScalarAndVector(double scalar, const double *vector, double *res) {
    for (int i = 0; i < N; i++) {
        res[i] = vector[i] * scalar;
    }
}

void CalcTMP(double *A, double *x, double *b, double *tmp) {
    MulMatrixVector(A, x, tmp);
    SubVectors(tmp, b, tmp);
}

void PrintMatrix(double *matrix, size_t size) {
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

void ZeroVector(double *vector) {
    for (int i = 0; i < N; i++) {
        vector[i] = 0;
    }
}

int main() {
    struct timespec start, end;

    double *A = (double *)calloc(N * N, sizeof(double));
    double *x = (double *)calloc(N, sizeof(double));
    double *b = (double *)calloc(N, sizeof(double));
    double *tmp = (double *)calloc(N, sizeof(double));

    Fill_A(A);
    //std::cout << "matrix A:" << std::endl;
    //PrintMatrix(A, N * N);

    Fill_b(b);
    //std::cout << "vector b:" << std::endl;
    //PrintMatrix(b, N);
    double module_b = CalcModuleOfVector(b);

    //std::cout << "vector x_0:" << std::endl;
    //PrintMatrix(x, N);

    int *buffer = (int *)calloc(size_buff, sizeof(int));
    int buff_counter = 0;
    unsigned int iteration = 0;
    double g_x = 0;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    while (SumBuff(buffer) != size_buff && iteration < MAX_ITERATIONS) {
        CalcTMP(A, x, b, tmp);
        //std::cout << std::endl << std::endl << "vector tmp:" << std::endl;
        //PrintMatrix(tmp, N);
        double module_tmp = CalcModuleOfVector(tmp);
        g_x = module_tmp / module_b;
        if (g_x < EPSILON) {
            buffer[buff_counter] = 1;
            buff_counter = (buff_counter + 1) % size_buff;
        }
        else {
            ZeroBuff(buffer);
        }

        if (SumBuff(buffer) == size_buff) {
            break;
        }
        MulScalarAndVector(TAU, tmp, tmp);
        SubVectors(x, tmp, x);
        //std::cout << "vector x:" << std::endl;
        //PrintMatrix(x, N);
        ZeroVector(tmp);
	iteration += 1;
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    std::cerr << "Without MPI:" << std::endl;
    if (iteration == MAX_ITERATIONS) {
	std::cerr << "Solution wasn't found after " << iteration << " iterations" << std::endl;
	std::cerr << "g(x) = " << g_x << std::endl;
    }
    else {
	std::cerr << "!!!Solution was found after " << iteration << " iteratons!!!" << std::endl;
    }
    std::cerr << "Time taken: " << end.tv_sec-start.tv_sec + 0.000000001 * (end.tv_nsec-start.tv_nsec) << " sec" << std::endl;;
    /*std::cout << std::endl << "result:" << std::endl;
    for (int i = 0; i < N; i++) {
        std::cout << x[i] << " ";
    }
    std::cout << std::endl;*/

    free(A);
    free(x);
    free(b);
    free(tmp);
    free(buffer);
    return 0;
}
