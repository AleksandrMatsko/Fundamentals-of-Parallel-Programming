#include <iostream>
#include <stdlib.h>
#include <stdio.h>

void FillMatrix(double *matrix, int lines, int columns, int key) {
    srand(key);
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j < columns; j++) {
            matrix[i * columns + j] = rand() / (double)RAND_MAX * 50 - 25;
        }
    }
}


void PrintMatrix(double *m, int lines, int columns, FILE *out) {
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j < columns; j++) {
            fprintf(out, "%f ", m[i * columns + j]);
        }
        fprintf(out, "\n");
    }
}

int main(int argc, char ** argv) {
    if (argc != 5) {
        std::cerr << "Not enough arguments" << std::endl;
        return -1;
    }
    int key = atoi(argv[1]);
    int N = atoi(argv[2]);
    int M = atoi(argv[3]);
    int K = atoi(argv[4]);


    double *A = (double *)calloc(N * M, sizeof(double));
    double *B = (double *)calloc(M * K, sizeof(double));

    FILE *out = fopen("data.txt", "w");
    fprintf(out, "%d %d %d\n", N, M, K);
    FillMatrix(A, N, M, key);
    FillMatrix(B, M, K, key + 1);

    PrintMatrix(A, N, M, out);
    PrintMatrix(B, M, K, out);

    fclose(out);
    free(A);
    free(B);
    return 0;
}