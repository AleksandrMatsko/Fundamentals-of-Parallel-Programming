#include <iostream>
#include <stdlib.h>
#include <stdbool.h>

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
        std::cout << ((num >> i) & 1);
    }
    //printf("\n");
}

void init(u_int *stage) {
    stage[0] = setBit(stage[0], 1, 1);
    stage[NUM_UINTS] = setBit(stage[NUM_UINTS], 2, 1);
    stage[2 * NUM_UINTS] = setBit(stage[2 * NUM_UINTS], 0, 7);

}

void calcNextStage(u_int *stage, u_int *next_stage, int num_lines) {
    for (int i = 0; i < num_lines; i++) {
        for (int j = 0; j < NUM_COLUMNS; j++) {
            u_int num_alive = 0;
            num_alive += getBit(stage[i * NUM_UINTS + (j + 1) * (j != NUM_COLUMNS - 1) / UINT_BIT_SIZE],
                                (j + 1) * (j != NUM_COLUMNS - 1) % UINT_BIT_SIZE);
            num_alive += getBit(stage[i * NUM_UINTS + (j - 1 + (!j) * NUM_COLUMNS) / UINT_BIT_SIZE],
                                (j - 1 + (!j) * NUM_COLUMNS) % UINT_BIT_SIZE);

            num_alive += getBit(stage[(i - 1 + (!i) * num_lines) % num_lines * NUM_UINTS + (j - 1 + (!j) * NUM_COLUMNS) / UINT_BIT_SIZE],
                                (j - 1 + (!j) * NUM_COLUMNS) % UINT_BIT_SIZE);
            num_alive += getBit(stage[(i - 1 + (!i) * num_lines) % num_lines * NUM_UINTS + j / UINT_BIT_SIZE],
                                j % UINT_BIT_SIZE);
            num_alive += getBit(stage[(i - 1 + (!i) * num_lines) % num_lines * NUM_UINTS + (j + 1) * (j != NUM_COLUMNS - 1) / UINT_BIT_SIZE],
                                (j + 1) * (j != NUM_COLUMNS - 1) % UINT_BIT_SIZE);

            num_alive += getBit(stage[(i + 1) * (i != num_lines - 1) % num_lines * NUM_UINTS + (j - 1 + (!j) * NUM_COLUMNS) / UINT_BIT_SIZE],
                                (j - 1 + (!j) * NUM_COLUMNS) % UINT_BIT_SIZE);
            num_alive += getBit(stage[(i + 1) * (i != num_lines - 1) % num_lines * NUM_UINTS + j / UINT_BIT_SIZE],
                                j % UINT_BIT_SIZE);
            num_alive += getBit(stage[(i + 1) * (i != num_lines - 1) % num_lines * NUM_UINTS + (j + 1) * (j != NUM_COLUMNS - 1) / UINT_BIT_SIZE],
                                (j + 1) * (j != NUM_COLUMNS - 1) % UINT_BIT_SIZE);

            u_int bit = getBit(stage[i * NUM_UINTS + j / UINT_BIT_SIZE], j % UINT_BIT_SIZE);
            if ((num_alive < 2 || num_alive > 3) && bit == 1) {
                next_stage[i * NUM_UINTS + j / UINT_BIT_SIZE] = setBit(next_stage[i * NUM_UINTS + j / UINT_BIT_SIZE], j % UINT_BIT_SIZE, 0);
            }
            else if (num_alive == 3 && bit == 0) {
                next_stage[i * NUM_UINTS + j / UINT_BIT_SIZE] = setBit(next_stage[i * NUM_UINTS + j / UINT_BIT_SIZE], j % UINT_BIT_SIZE, 1);
            }
            else {
                next_stage[i * NUM_UINTS + j / UINT_BIT_SIZE] = setBit(next_stage[i * NUM_UINTS + j / UINT_BIT_SIZE], j % UINT_BIT_SIZE, bit);
            }
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
    NUM_LINES = atoi(argv[1]);
    NUM_COLUMNS = atoi(argv[2]);
    NUM_UINTS = NUM_COLUMNS / (UINT_BIT_SIZE) + (NUM_COLUMNS % (UINT_BIT_SIZE) != 0);
    u_int *prev_stages[MAX_ITERATION];
    int current_iteration = 0;
    size_t field_size = NUM_LINES * NUM_UINTS;
    u_int *stage = (u_int *)calloc(field_size, sizeof(u_int));
    u_int *next_stage = (u_int *)calloc(field_size, sizeof(u_int));

    init(stage);
    for (int i = 0; i < NUM_LINES; i++) {
        for (int j = 0; j < NUM_UINTS; j++) {
            printUINT(stage[i * NUM_UINTS + j]);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    bool period = false;
    while (current_iteration < MAX_ITERATION && !period) {
        prev_stages[current_iteration] = (u_int *) calloc(field_size, sizeof(u_int));
        memcpy(prev_stages[current_iteration], stage, field_size * sizeof(u_int));
        calcNextStage(stage, next_stage, NUM_LINES);
        memcpy(stage, next_stage, field_size * sizeof(u_int));
        /*for (int i = 0; i < NUM_LINES; i++) {
            for (int j = 0; j < NUM_UINTS; j++) {
                printUINT(stage[i * NUM_UINTS + j]);
            }
            std::cout << std::endl;
        }*/
        for (int i = current_iteration; i > -1; i--) {
            if (checkPrevStage(prev_stages[i], stage, NUM_LINES)) {
                period = true;
                break;
            }
        }
        current_iteration += 1;
        std::cout << current_iteration << std::endl;
    }

    if (period) {
        std::cout << "Periodic" << std::endl;
    }

    for (int i = 0; i < current_iteration; i++) {
        free(prev_stages[i]);
    }
    free(stage);
    free(next_stage);
    return 0;
}