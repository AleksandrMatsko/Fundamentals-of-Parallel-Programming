#include <iostream>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include <malloc.h>
#include <string.h>
#include <limits.h>

int NUM_LINES;
int NUM_COLUMNS;
int NUM_UINTS;
const int MAX_ITERATION = 100000;

typedef unsigned int u_int;
const size_t UINT_BIT_SIZE = sizeof(u_int) * 8;

u_int getBit(u_int line, u_int index) {
    return (line >> index) & 1;
}

/*
 * 0 0 -> 0
 * 0 1 -> 1
 * 1 0 -> 0
 * 1 1 -> 1
 */

u_int setBit(u_int line, u_int index, u_int val) {
    return (line & (UINT_MAX - (1 << index))) + (val << index);
}

void printUINT(u_int num) {
    for (int i = 0; i < sizeof(u_int) * 8; i++) {
        std::cerr << ((num >> i) & 1);
    }
}

void init(u_int *stage) {
    stage[0] = setBit(stage[0], 1, 1);
    stage[NUM_UINTS] = setBit(stage[NUM_UINTS], 2, 1);
    stage[2 * NUM_UINTS] = setBit(stage[2 * NUM_UINTS], 0, 7);
}

void calcNextStageCenter(u_int *stage, u_int *next_stage, int num_lines) {
    for (int i = 1; i < num_lines - 1; i++) {
        for (int j = 0; j < NUM_COLUMNS; j++) {
            /*
             *      order of checking neighbours
             *
             *      3 4 5
             *      2 c 1
             *      6 7 8
             */
            u_int num_alive = 0;
            num_alive += getBit(stage[i * NUM_UINTS + (j + 1) * (j != NUM_COLUMNS - 1) / UINT_BIT_SIZE],
                                (j + 1) * (j != NUM_COLUMNS - 1) % UINT_BIT_SIZE);
            num_alive += getBit(stage[i * NUM_UINTS + (j - 1 + (!j) * NUM_COLUMNS) / UINT_BIT_SIZE],
                                (j - 1 + (!j) * NUM_COLUMNS) % UINT_BIT_SIZE);

            num_alive += getBit(stage[(i - 1) * NUM_UINTS + (j - 1 + (!j) * NUM_COLUMNS) / UINT_BIT_SIZE],
                                (j - 1 + (!j) * NUM_COLUMNS) % UINT_BIT_SIZE);
            num_alive += getBit(stage[(i - 1) * NUM_UINTS + j / UINT_BIT_SIZE],
                                j % UINT_BIT_SIZE);
            num_alive += getBit(stage[(i - 1) * NUM_UINTS + (j + 1) * (j != NUM_COLUMNS - 1) / UINT_BIT_SIZE],
                                (j + 1) * (j != NUM_COLUMNS - 1) % UINT_BIT_SIZE);

            num_alive += getBit(stage[(i + 1) * NUM_UINTS + (j - 1 + (!j) * NUM_COLUMNS) / UINT_BIT_SIZE],
                                (j - 1 + (!j) * NUM_COLUMNS) % UINT_BIT_SIZE);
            num_alive += getBit(stage[(i + 1) * NUM_UINTS + j / UINT_BIT_SIZE],
                                j % UINT_BIT_SIZE);
            num_alive += getBit(stage[(i + 1) * NUM_UINTS + (j + 1) * (j != NUM_COLUMNS - 1) / UINT_BIT_SIZE],
                                (j +  1) * (j != NUM_COLUMNS - 1) % UINT_BIT_SIZE);

            u_int bit = getBit(stage[i * NUM_UINTS + j / UINT_BIT_SIZE], j % UINT_BIT_SIZE);
            if ((num_alive < 2 || num_alive > 3) && bit == 1) {
                next_stage[i * NUM_UINTS + j / UINT_BIT_SIZE] = setBit(next_stage[i * NUM_UINTS + j / UINT_BIT_SIZE],
                                                                       j % UINT_BIT_SIZE, 0);
            }
            else if (num_alive == 3 && bit == 0) {
                next_stage[i * NUM_UINTS + j / UINT_BIT_SIZE] = setBit(next_stage[i * NUM_UINTS + j / UINT_BIT_SIZE],
                                                                       j % UINT_BIT_SIZE, 1);
            }
            else {
                next_stage[i * NUM_UINTS + j / UINT_BIT_SIZE] = setBit(next_stage[i * NUM_UINTS + j / UINT_BIT_SIZE],
                                                                       j % UINT_BIT_SIZE, bit);
            }
        }
    }
}

void calcNextStageUp(u_int *stage, u_int *next_stage, u_int *cells_up) {
    for (int j = 0; j < NUM_COLUMNS; j++) {
        u_int num_alive = 0;
        num_alive += getBit(stage[(j + 1) * (j != NUM_COLUMNS - 1) / UINT_BIT_SIZE],
                            (j + 1) * (j != NUM_COLUMNS - 1) % UINT_BIT_SIZE);
        num_alive += getBit(stage[(j - 1 + (!j) * NUM_COLUMNS) / UINT_BIT_SIZE],
                            (j - 1 + (!j) * NUM_COLUMNS) % UINT_BIT_SIZE);

        num_alive += getBit(cells_up[(j - 1 + (!j) * NUM_COLUMNS) / UINT_BIT_SIZE],
                            (j - 1 + (!j) * NUM_COLUMNS) % UINT_BIT_SIZE);
        num_alive += getBit(cells_up[j / UINT_BIT_SIZE],
                            j % UINT_BIT_SIZE);
        num_alive += getBit(cells_up[(j + 1) * (j != NUM_COLUMNS - 1) / UINT_BIT_SIZE],
                            (j + 1) * (j != NUM_COLUMNS - 1) % UINT_BIT_SIZE);

        num_alive += getBit(stage[NUM_UINTS + (j - 1 + (!j) * NUM_COLUMNS) / UINT_BIT_SIZE],
                            (j + UINT_BIT_SIZE - 1 + (!j) * NUM_COLUMNS) % UINT_BIT_SIZE);
        num_alive += getBit(stage[NUM_UINTS + j / UINT_BIT_SIZE],
                            j % UINT_BIT_SIZE);
        num_alive += getBit(stage[NUM_UINTS + (j + 1) * (j != NUM_COLUMNS - 1) / UINT_BIT_SIZE],
                            (j + 1) * (j != NUM_COLUMNS - 1) % UINT_BIT_SIZE);

        u_int bit = getBit(stage[j / UINT_BIT_SIZE], j % UINT_BIT_SIZE);
        if ((num_alive < 2 || num_alive > 3) && bit == 1) {
            next_stage[j / UINT_BIT_SIZE] = setBit(next_stage[j / UINT_BIT_SIZE],
                                                   j % UINT_BIT_SIZE, 0);
        }
        else if (num_alive == 3 && bit == 0) {
            next_stage[j / UINT_BIT_SIZE] = setBit(next_stage[j / UINT_BIT_SIZE],
                                                   j % UINT_BIT_SIZE, 1);
        }
        else {
            next_stage[j / UINT_BIT_SIZE] = setBit(next_stage[j / UINT_BIT_SIZE],
                                                   j % UINT_BIT_SIZE, bit);
        }
    }
}

void calcNextStageDown(u_int *stage, u_int *next_stage, int num_lines, u_int *cells_down) {
    for (int j = 0; j < NUM_COLUMNS; j++) {
        u_int num_alive = 0;
        num_alive += getBit(stage[(num_lines - 1) * NUM_UINTS + (j + 1) * (j != NUM_COLUMNS - 1) / UINT_BIT_SIZE],
                            (j + 1) * (j != NUM_COLUMNS - 1) % UINT_BIT_SIZE);
        num_alive += getBit(stage[(num_lines - 1) * NUM_UINTS + (j - 1 + (!j) * NUM_COLUMNS) / UINT_BIT_SIZE],
                            (j - 1 + (!j) * NUM_COLUMNS) % UINT_BIT_SIZE);

        num_alive += getBit(stage[(num_lines - 2) * NUM_UINTS + (j - 1 + (!j) * NUM_COLUMNS) / UINT_BIT_SIZE],
                            (j - 1 + (!j) * NUM_COLUMNS) % UINT_BIT_SIZE);
        num_alive += getBit(stage[(num_lines - 2) * NUM_UINTS + j / UINT_BIT_SIZE],
                            j % UINT_BIT_SIZE);
        num_alive += getBit(stage[(num_lines - 2) * NUM_UINTS + (j + 1) * (j != NUM_COLUMNS - 1) / UINT_BIT_SIZE],
                            (j + 1) * (j != NUM_COLUMNS - 1) % UINT_BIT_SIZE);

        num_alive += getBit(cells_down[(j - 1 + (!j) * NUM_COLUMNS) / UINT_BIT_SIZE],
                            (j - 1 + (!j) * NUM_COLUMNS) % UINT_BIT_SIZE);
        num_alive += getBit(cells_down[j / UINT_BIT_SIZE],
                            j % UINT_BIT_SIZE);
        num_alive += getBit(cells_down[(j + 1) * (j != NUM_COLUMNS - 1) / UINT_BIT_SIZE],
                            (j +  1) * (j != NUM_COLUMNS - 1) % UINT_BIT_SIZE);

        u_int bit = getBit(stage[(num_lines - 1) * NUM_UINTS + j / UINT_BIT_SIZE], j % UINT_BIT_SIZE);
        if ((num_alive < 2 || num_alive > 3) && bit == 1) {
            next_stage[(num_lines - 1) * NUM_UINTS + j / UINT_BIT_SIZE] = setBit(next_stage[(num_lines - 1) * NUM_UINTS + j / UINT_BIT_SIZE],
                                                                                 j % UINT_BIT_SIZE, 0);
        }
        else if (num_alive == 3 && bit == 0) {
            next_stage[(num_lines - 1) * NUM_UINTS + j / UINT_BIT_SIZE] = setBit(next_stage[(num_lines - 1) * NUM_UINTS + j / UINT_BIT_SIZE],
                                                                                 j % UINT_BIT_SIZE, 1);
        }
        else {
            next_stage[(num_lines - 1) * NUM_UINTS + j / UINT_BIT_SIZE] = setBit(next_stage[(num_lines - 1) * NUM_UINTS + j / UINT_BIT_SIZE],
                                                                                 j % UINT_BIT_SIZE, bit);
        }
    }
}


bool checkPrevStage(u_int *prev_stage, u_int *stage, int num_lines) {
    for (int i = 0; i < num_lines; i++) {
        for (int j = 0; j < NUM_UINTS; j++) {
            if (stage[i * NUM_UINTS + j] != prev_stage[i * NUM_UINTS + j]) {
                return false;
            }
        }
    }
    return true;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Wrong amount of arguments: expected 3" << std::endl;
        return 1;
    }

    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    NUM_LINES = atoi(argv[1]);
    NUM_COLUMNS = atoi(argv[2]);
    NUM_UINTS = NUM_COLUMNS / (UINT_BIT_SIZE) + (NUM_COLUMNS % (UINT_BIT_SIZE) != 0);

    int *strings_in_tasks = (int *)calloc(size, sizeof(int));
    int *shifts = (int *) calloc(size, sizeof(int));

    int num_of_strings = NUM_LINES / size;
    int remainder = NUM_LINES % size;
    for (int i = 0; i < size; i++) {
        strings_in_tasks[i] = num_of_strings;
        if (i < remainder) {
            strings_in_tasks[i] += 1;
        }
    }

    int start_index = 0;
    for (int i = 0; i < size; i++) {
        shifts[i] = start_index;
        start_index += strings_in_tasks[i] * NUM_UINTS;
    }

    int *to_send = (int *)calloc(size, sizeof(int));
    for (int i = 0; i < size; i++) {
        to_send[i] = NUM_UINTS * strings_in_tasks[i];
    }

    u_int *prev_stages[MAX_ITERATION];
    size_t full_field_size = NUM_LINES * NUM_UINTS;
    size_t field_size = strings_in_tasks[rank] * NUM_UINTS;
    u_int *stage = (u_int *)calloc(field_size, sizeof(u_int));
    u_int *next_stage = (u_int *)calloc(field_size, sizeof(u_int));

    double start_time = MPI_Wtime();
    u_int *full_field;
    if (rank == 0 ) {
        full_field = (u_int *)calloc(full_field_size, sizeof(u_int));
        init(full_field);
    }

    MPI_Scatterv(full_field, to_send, shifts, MPI_UNSIGNED, stage, NUM_UINTS * strings_in_tasks[rank], MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        free(full_field);
    }

    u_int *cells_up = (u_int *)calloc(NUM_UINTS, sizeof(u_int));
    u_int *cells_down = (u_int *)calloc(NUM_UINTS, sizeof(u_int));
    int period[MAX_ITERATION];
    int global_period[MAX_ITERATION];

    int current_iteration = 0;
    int is_global_period = 0;
    while (current_iteration < MAX_ITERATION) {
        MPI_Request req_recv_up, req_recv_down;
        MPI_Irecv(cells_up, NUM_UINTS, MPI_UNSIGNED, (rank - 1 + size) % size, 0, MPI_COMM_WORLD, &req_recv_up);
        MPI_Irecv(cells_down, NUM_UINTS, MPI_UNSIGNED, (rank + 1) % size, 1, MPI_COMM_WORLD, &req_recv_down);

        MPI_Request req_send[2];
        MPI_Isend(stage, NUM_UINTS, MPI_UNSIGNED, (rank - 1 + size) % size, 1, MPI_COMM_WORLD, req_send);
        MPI_Isend(&stage[(strings_in_tasks[rank] - 1) * NUM_UINTS], NUM_UINTS, MPI_UNSIGNED, (rank + 1) % size, 0, MPI_COMM_WORLD, &req_send[1]);

        prev_stages[current_iteration] = (u_int *) calloc(field_size, sizeof(u_int));
        memcpy(prev_stages[current_iteration], stage, field_size * sizeof(u_int));

        calcNextStageCenter(stage, next_stage, strings_in_tasks[rank]);
        bool up_done = false;
        bool down_done = false;
        int suc_up;
        int suc_down;
        while (!up_done || !down_done) {
            MPI_Test(&req_recv_up, &suc_up, MPI_STATUS_IGNORE);
            if (suc_up && !up_done) {
                calcNextStageUp(stage, next_stage, cells_up);
                up_done = true;
            }
            MPI_Test(&req_recv_down, &suc_down, MPI_STATUS_IGNORE);
            if (suc_down && !down_done) {
                calcNextStageDown(stage, next_stage, strings_in_tasks[rank], cells_down);
                down_done = true;
            }
        }
        memcpy(stage, next_stage, field_size * sizeof(u_int));
        MPI_Waitall(2, req_send, MPI_STATUS_IGNORE);

        for (int i = current_iteration; i > -1; i--) {
            if (checkPrevStage(prev_stages[i], stage, strings_in_tasks[rank])) {
                period[i] = 1;
            }
            else {
                period[i] = 0;
            }
        }
        is_global_period = 0;
        MPI_Allreduce(&period, &global_period, current_iteration + 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        for (int i = current_iteration; i >= 0; i--) {
            if (global_period[i] == size) {
                is_global_period = global_period[i];
                break;
            }
        }
        if (is_global_period == size) {
            break;
        }
        current_iteration += 1;
    }
    for (int i = 0; i < current_iteration; i++) {
        free(prev_stages[i]);
    }
    free(cells_up);
    free(cells_down);
    double end_time = MPI_Wtime();

    if (is_global_period && rank == 0) {
        std::cout << "size = " << size << std::endl;
        std::cout << "Periodic" << std::endl;
        std::cout << "Time taken: " << end_time - start_time << std::endl << std::endl;
    }

    free(stage);
    free(next_stage);
    free(strings_in_tasks);
    free(shifts);
    free(to_send);
    MPI_Finalize();
    return 0;
}
