#include <iostream>
#include <stdlib.h>
#include <stdio.h>

void Fill_A(double *A, int N, int key) {
    srand(key);
    for (int i = 0; i < N; i++) {
	A[i * N + i] = rand() / (double)RAND_MAX * 50 + 500;
	for (int j = i + 1; j < N; j++) {
	    A[i * N + j] = rand() / (double)RAND_MAX * 50 - 25;
	    A[j * N + i] = A[i * N + j];
	}
    }
}

void Fill_b(double *b, int N) {
    for (int i = 0; i < N; i++) {
	b[i] = N + 1;
    }
}

void PrintMatrix(double *m, int lines, int N, FILE *out) {
    for (int i = 0; i < lines; i++) {
	for (int j = 0; j < N; j++) {
	    fprintf(out, "%f ", m[i * N + j]);
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
    double TAU = atof(argv[3]);
    double EPSILON = atof(argv[4]);
    
    double *A = (double *)calloc(N * N, sizeof(double));
    double *b = (double *)calloc(N * N, sizeof(double));
    
    FILE *out = fopen("data.txt", "w");
    fprintf(out, "%d %lf %lf\n", N, TAU, EPSILON);
    Fill_A(A, N, key);
    Fill_b(b, N);
    
    PrintMatrix(A, N, N, out);
    PrintMatrix(b, 1, N, out);
    
    fclose(out);
    free(A);
    free(b);
    return 0;
}