#include <iostream>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

int N;
int M;
int K;
int LINES_LATTICE = 3;
int COLUMNS_LATTICE = 8;

void FillMatrix(double *matrix, int lines, int columns, FILE *in) {
    for (int i = 0; i < lines * columns; i++) {
        fscanf(in, "%lf", &matrix[i]);
    }
}

void MulMatrix(double *A, double *B, double *C) {
    for (int i = 0; i < N / LINES_LATTICE; i++) {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < K / COLUMNS_LATTICE; k++) {
                C[i * K / COLUMNS_LATTICE + k] += A[i * M + j] * B[j * K / COLUMNS_LATTICE + k];
            }
        }
    }
}

void PrintMatrix(double *matrix, int lines, int columns) {
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j < columns; j++) {
            std::cerr << matrix[i * columns + j] << " ";
        }
        std::cerr << std::endl;
    }
}

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Wrong number of arguments" << std::endl;
        return 0;
    }

    LINES_LATTICE = atoi(argv[1]);
    COLUMNS_LATTICE = atoi(argv[2]);

    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double start = MPI_Wtime();
    MPI_Comm COMM_2D;
    int dims[2] = {LINES_LATTICE, COLUMNS_LATTICE};
    int periods[2] = {0, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, false, &COMM_2D);
    int coords[2];
    MPI_Cart_coords(COMM_2D, rank, 2, coords);

    MPI_Comm COMM_VERTICAL;
    MPI_Comm COMM_HORIZONTAL;
    int remain_dims_ver[2] = {true, false};
    int remain_dims_hor[2] = {false, true};
    MPI_Cart_sub(COMM_2D, remain_dims_hor, &COMM_HORIZONTAL);
    MPI_Cart_sub(COMM_2D, remain_dims_ver, &COMM_VERTICAL);

    FILE *in = fopen("data.txt", "r");
    if (in == NULL) {
        std::cerr << "File not open" << std::endl;
    }
    fscanf(in, "%d %d %d", &N, &M, &K);
    double *A = (double *)calloc(N * M, sizeof(double));
    double *B = (double *)calloc(M * K, sizeof(double));
    if (coords[0] == 0 && coords[1] == 0) {
        FillMatrix(A, M, N, in);
        FillMatrix(B, M, K, in);
    }
    fclose(in);
    double *A_part = (double *)calloc(N / LINES_LATTICE * M, sizeof(double));
    MPI_Scatter(A, M * N / LINES_LATTICE, MPI_DOUBLE, A_part, M * N / LINES_LATTICE, MPI_DOUBLE, 0, COMM_VERTICAL);
    MPI_Bcast(A_part, M * N / LINES_LATTICE, MPI_DOUBLE, 0, COMM_HORIZONTAL);

    MPI_Datatype B_TYPE;
    MPI_Type_vector(M, K / COLUMNS_LATTICE, K, MPI_DOUBLE, &B_TYPE);
    MPI_Type_commit(&B_TYPE);

    double *B_part = (double *) calloc(K / COLUMNS_LATTICE * M,sizeof(double));
    int coords_to_receive[2] = {0, 0};
    int root;
    MPI_Cart_rank(COMM_2D, coords_to_receive, &root);
    for (int i = 1; i < COLUMNS_LATTICE; i++) {
        int rank_receive;
        coords_to_receive[1] = i;
        MPI_Cart_rank(COMM_2D, coords_to_receive, &rank_receive);
        if (coords[0] == 0 && coords[1] == 0) {
            MPI_Send(&B[i * K / COLUMNS_LATTICE], 1, B_TYPE, rank_receive, i, COMM_2D);
        }
    }
    if (coords[0] == 0 && coords[1] != 0) {
        MPI_Recv(B_part, M * K / COLUMNS_LATTICE, MPI_DOUBLE, root, coords[1], COMM_2D, MPI_STATUS_IGNORE);
    }
    if (coords[0] == 0 && coords[1] == 0) {
        for (int i = 0; i < M; i++) {
            memcpy(&B_part[i * K / COLUMNS_LATTICE], &B[i * K], K / COLUMNS_LATTICE * sizeof(double));
        }
    }
    MPI_Bcast(B_part, M * K / COLUMNS_LATTICE, MPI_DOUBLE, 0, COMM_VERTICAL);

    double *C_part = (double *)calloc(N / LINES_LATTICE * K / COLUMNS_LATTICE, sizeof(double));
    MulMatrix(A_part, B_part, C_part);

    MPI_Datatype C_TYPE;
    MPI_Type_vector(N / LINES_LATTICE, K / COLUMNS_LATTICE, K, MPI_DOUBLE, &C_TYPE);
    MPI_Type_commit(&C_TYPE);

    double *C = (double *)calloc(N * K, sizeof(double));

    if (rank != root) {
        MPI_Send(C_part, N / LINES_LATTICE * K / COLUMNS_LATTICE, MPI_DOUBLE, root, rank, COMM_2D);
    }
    int coords_src[2] = {0, 0};
    for (int i = 0; i < LINES_LATTICE; i++) {
        for (int j = 0; j < COLUMNS_LATTICE; j++) {
            if ((i != 0 || j != 0) && rank == root) {
                int rank_src;
                coords_src[0] = i;
                coords_src[1] = j;
                MPI_Cart_rank(COMM_2D, coords_src, &rank_src);
                MPI_Recv(&C[i * N / LINES_LATTICE * K + j * K / COLUMNS_LATTICE], 1, C_TYPE, rank_src, rank_src, COMM_2D, MPI_STATUS_IGNORE);
            }
        }
    }
    if (rank == root) {
        for (int i = 0; i < N / LINES_LATTICE; i++) {
            memcpy(&C[i * N], &C_part[i], K / COLUMNS_LATTICE * sizeof(double));
        }
        //PrintMatrix(C, N, K);
    }

    free(C);
    free(B);
    free(A);
    free(B_part);
    free(A_part);
    free(C_part);
    MPI_Type_free(&C_TYPE);
    MPI_Type_free(&B_TYPE);
    MPI_Comm_free(&COMM_VERTICAL);
    MPI_Comm_free(&COMM_HORIZONTAL);
    MPI_Comm_free(&COMM_2D);

    double end = MPI_Wtime();

    if (rank == 0) {
        std::cout << "MPI " << LINES_LATTICE << " * " << COLUMNS_LATTICE << std::endl;
        std::cout << "Time taken: " << end - start << std::endl << std::endl;
    }

    MPI_Finalize();
    return 0;
}

