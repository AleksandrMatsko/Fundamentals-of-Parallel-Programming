#include <iostream>
#include <malloc.h>
#include <stdlib.h>
#include <omp.h>

int N = 200;
double TAU = 0.0000001;
double EPSILON = 0.000001;
unsigned int MAX_ITERATIONS = 1000000;
const int MAX_NUM_CHECKS = 5;

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
#pragma omp for
    for (int i = 0; i < N; i++) {
        res[i] = 0;
        for (int j = 0; j < N; j++) {
            res[i] += matrix[i * N + j] * vector[j];
        }
    }
}

/*double CalcModuleOfVector(double *vector) {
    double res = 0;
#pragma omp shared(res) for reduction(+: res)
    for (int i = 0; i < N; i++) {
        res += vector[i] * vector[i];
    }
    return res;
}*/

void SubVectors(const double *vector_1, const double *vector_2, double *res) {
#pragma omp for
    for (int i = 0; i < N; i++) {
        res[i] = vector_1[i] - vector_2[i];
    }
}

void MulScalarAndVector(double scalar, const double *vector, double *res) {
#pragma omp for
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
    for (int i = 0; i < MAX_NUM_CHECKS; i++) {
        res += buff[i];
    }
    return res;
}

void ZeroBuff(int *buff) {
#pragma omp for
    for (int i = 0; i < MAX_NUM_CHECKS; i++) {
       buff[i] = 0;
    }
}

void ZeroVector(double *vector) {
#pragma omp for
    for (int i = 0; i < N; i++) {
        vector[i] = 0;
    }
}

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Not enough arguments" << std::endl;
        return 0;
    }
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
    unsigned int last_iteration = 0;
    double g_x = 0;

    int n = atoi(argv[1]);

    double start = omp_get_wtime();
    double module_b;
    double module_tmp;
    #pragma omp parallel num_threads(n) firstprivate(iteration, MAX_NUM_CHECKS, MAX_ITERATIONS, check_counter)
    {
        module_b = 0;
        #pragma omp for reduction(+: module_b)
        for (int i = 0; i < N; i++) {
            module_b += b[i] * b[i];
        }
        double exit_val = EPSILON * EPSILON * module_b;
        while (check_counter < MAX_NUM_CHECKS && iteration < MAX_ITERATIONS) {
            CalcTMP(A, x, b, tmp);
            module_tmp = 0;
            #pragma omp for reduction(+: module_tmp)
            for (int i = 0; i < N; i++) {
                module_tmp += tmp[i] * tmp[i];
            }
            g_x = module_tmp / module_b;
            if (module_tmp < exit_val) {
                check_counter += 1;
                #pragma omp master
                last_iteration = iteration;
            } else {
                check_counter = 0;
            }
            MulScalarAndVector(TAU, tmp, tmp);
            SubVectors(x, tmp, x);
            ZeroVector(tmp);
            //#pragma omp atomic
            iteration += 1;
            /*if (iteration % 1000 == 0) {
                std::cerr << iteration << " g(x) = " << g_x << std::endl;
            }*/
            #pragma omp barrier
        }
    };
    double end = omp_get_wtime();

    std::cout << "With OpenMP" << std::endl;
    if (last_iteration == MAX_ITERATIONS) {
	    std::cout << "Solution wasn't found after " << MAX_ITERATIONS << " iterations" << std::endl;
	    std::cout << "g(x) = " << g_x << std::endl;
    }
    else {
	    std::cout << "!!!Solution was found after " << last_iteration << " iteratons!!!" << std::endl;
    }
    std::cout << "Time taken: " << end - start << " sec" << std::endl << std::endl;
    
    free(A);
    free(x);
    free(b);
    free(tmp);
    //free(buffer);
    return 0;
}
