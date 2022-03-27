#include <iostream>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>

int N = 200;
double TAU = 0.0000001;
double EPSILON = 0.000001;
const unsigned int MAX_ITERATIONS = 1000000;
const int MAX_NUM_CHECKS= 5;

void Fill_A(double *A, int N, FILE *in) {
    for (int i = 0; i < N * N; i++) {
	fscanf(in, "%lf", &A[i]);
    }
}

void Fill_b(double *b, int N, FILE *in) {
    for (int i = 0; i < N; i++) {
        fscanf(in, "%lf", &b[i]);
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
        res += vector[i] * vector[i];
    }
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

void ZeroVector(double *vector) {
    for (int i = 0; i < N; i++) {
        vector[i] = 0;
    }
}

int main() {
    struct timespec start, end;

    FILE *in = fopen("data.txt", "r");
    fscanf(in, "%d %lf %lf", &N, &TAU, &EPSILON);

    double *A = (double *)calloc(N * N, sizeof(double));
    double *x = (double *)calloc(N, sizeof(double));
    double *b = (double *)calloc(N, sizeof(double));
    double *tmp = (double *)calloc(N, sizeof(double));
    
    Fill_A(A, N, in);
    Fill_b(b, N, in);

    int check_counter = 0;
    unsigned int iteration = 0;
    double g_x = 0;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    double module_b = CalcModuleOfVector(b);
    std::cerr << "module b = " << module_b << std::endl;
    double exit_val = EPSILON * EPSILON * module_b;
    while (iteration < MAX_ITERATIONS) {
        CalcTMP(A, x, b, tmp);
        double module_tmp = CalcModuleOfVector(tmp);
        g_x = module_tmp / module_b;
        if (module_tmp < exit_val) {
            check_counter += 1;
        }
        else {
            check_counter = 0;
        }
        if (check_counter == MAX_NUM_CHECKS) {
            break;
        }
        MulScalarAndVector(TAU, tmp, tmp);
        SubVectors(x, tmp, x);
        ZeroVector(tmp);
	    iteration += 1;
        if (iteration % 1000 == 0) {
            std::cerr << iteration << " g(x) = " << g_x << std::endl;
        }
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    std::cout << "Without OpenMP:" << std::endl;
    if (iteration == MAX_ITERATIONS) {
	std::cout << "Solution wasn't found after " << iteration << " iterations" << std::endl;
	std::cout << "g(x) = " << g_x << std::endl;
    }
    else {
	std::cout << "!!!Solution was found!!!" << std::endl;
    }
    std::cout << "Time taken: " << end.tv_sec-start.tv_sec + 0.000000001 * (end.tv_nsec-start.tv_nsec) << " sec" << std::endl << std::endl;
    
    free(A);
    free(x);
    free(b);
    free(tmp);
    return 0;
}
