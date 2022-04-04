#include <iostream>
#include <malloc.h>
#include <stdio.h>

int N;
int M;
int K;

void FillMatrix(double *matrix, int lines, int columns, FILE *in) {
    for (int i = 0; i < lines * columns; i++) {
        fscanf(in, "%lf", &matrix[i]);
    }
}

void MulMatrix(double *A, double *B, double *C) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < K; k++) {
                C[i * K + k] += A[i * M + j] * B[j * K + k];
            }
        }
    }
}

void PrintMatrix(double *matrix, int lines, int columns) {
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j < columns; j++) {
            std::cout << matrix[i * columns + j] << " ";
        }
        std::cout << std::endl;
    }
}

int main() {
    FILE *in = fopen("data.txt", "r");
    fscanf(in, "%d %d %d", &N, &M, &K);
    double *A = (double *)calloc(N * M, sizeof(double));
    double *B = (double *)calloc(M * K, sizeof(double));
    FillMatrix(A, M, N, in);
    FillMatrix(B, M, K, in);
    double *C = (double *)calloc(N * K, sizeof(double));
    MulMatrix(A, B, C);
    PrintMatrix(C, N, K);
    return 0;
}
