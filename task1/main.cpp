#include <iostream>
#include <math.h>

const int N = 4;
const double tau = 0.001;
const double epsilon = 0.00001;

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

void Fill_Vector(double *x) {
    for (int i = 0; i < N; i++) {
        x[i] = rand() / (double)RAND_MAX * 200 - 100;
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
    double *A = (double *)calloc(N * N, sizeof(double));
    double *x = (double *)calloc(N, sizeof(double));
    double *b = (double *)calloc(N, sizeof(double));
    double *tmp = (double *)calloc(N, sizeof(double));

    Fill_A(A);
    std::cout << "matrix A:" << std::endl;
    PrintMatrix(A, N * N);

    Fill_Vector(b);
    std::cout << "vector b:" << std::endl;
    PrintMatrix(b, N);
    double module_b = CalcModuleOfVector(b);

    //Fill_Vector(x);
    std::cout << "vector x_0:" << std::endl;
    PrintMatrix(x, N);

    int *buffer = (int *)calloc(size_buff, sizeof(int));
    int buff_counter = 0;
    while (SumBuff(buffer) != size_buff) {
        CalcTMP(A, x, b, tmp);
        std::cout << std::endl << std::endl << "vector tmp:" << std::endl;
        PrintMatrix(tmp, N);
        double module_tmp = CalcModuleOfVector(tmp);
        std::cout << "g(x):" << std::endl;
        std::cout << module_tmp / module_b << std::endl;
        if (module_tmp / module_b < epsilon) {
            buffer[buff_counter] = 1;
            buff_counter = (buff_counter + 1) % size_buff;
        }
        else {
            ZeroBuff(buffer);
        }

        if (SumBuff(buffer) == size_buff) {
            break;
        }
        MulScalarAndVector(tau, tmp, tmp);
        SubVectors(x, tmp, x);
        std::cout << "vector x:" << std::endl;
        PrintMatrix(x, N);
        ZeroVector(tmp);
    }

    std::cout << std::endl << "result:" << std::endl;
    for (int i = 0; i < N; i++) {
        std::cout << x[i] << " ";
    }
    std::cout << std::endl;

    free(A);
    free(x);
    free(b);
    free(tmp);
    free(buffer);
    return 0;
}
