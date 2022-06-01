#include <iostream>
#include <mpi.h>
#include <string.h>
#include <pthread.h>
#include <stdlib.h>
#include <math.h>

int MAX_ITERATIONS = 5;
int TOTAL_TASK_NUM = 300;
int TASK_MULTIPLIER = 1000;
int TASKS_IN_PROCESS;

const int STOP_CODE = -1;
const int PERCENT = 40;
const int MIN_TASK_SEND = 2;

int left_task = 0;
int right_task = 0;

pthread_mutex_t task_mutex;


struct arguments_t {
    int *tasks;
    int rank;
    int size;
};

void *answerTask(void *argument) {
    MPI_Status status;
    arguments_t parameters = *((arguments_t *)argument);
    int *tasks = parameters.tasks;

    int status_rank = 0;

    while (true) {
        MPI_Recv(&status_rank, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        if (status_rank == STOP_CODE) {
            break;
        }

        int send_tasks = 0;

        pthread_mutex_lock(&task_mutex);
        send_tasks = ((right_task - left_task) * PERCENT) / 100;
        if (send_tasks > MIN_TASK_SEND) {
            right_task -= send_tasks;
        }
        else {
            send_tasks = 0;
        }
        pthread_mutex_unlock(&task_mutex);

        MPI_Send(&send_tasks, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
        if (send_tasks > MIN_TASK_SEND) {
            if (right_task / TASKS_IN_PROCESS != (right_task + send_tasks) / TASKS_IN_PROCESS) {
                MPI_Send(&tasks[right_task % TASKS_IN_PROCESS], TASKS_IN_PROCESS - right_task % TASKS_IN_PROCESS,
                         MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
                send_tasks -= TASKS_IN_PROCESS - right_task % TASKS_IN_PROCESS;
                MPI_Send(tasks, send_tasks, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
            }
            else {
                MPI_Send(&tasks[right_task % TASKS_IN_PROCESS], send_tasks, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
            }
        }

    }


}

void calculationTask(arguments_t *argument) {
    int *tasks = argument->tasks;
    int rank = argument->rank;
    int size = argument->size;
    MPI_Status status;
    int *tasks_buffer = (int *) calloc(TASKS_IN_PROCESS * PERCENT / 100 + 1, sizeof(int));

    pthread_mutex_lock(&task_mutex);
    while (true) {
        int task_size = tasks[left_task % TASKS_IN_PROCESS];
        tasks[left_task % TASKS_IN_PROCESS] = 0;
        left_task += 1;
        int remaining_tasks = right_task - left_task;
        pthread_mutex_unlock(&task_mutex);

        double sum = 0;
        for (int i = 0; i < task_size; i++) {
            sum += sin(i);
        }
        if (sum >= task_size) {
            std::cerr << "wow" << std::endl;
        }

        pthread_mutex_lock(&task_mutex);
        if (left_task == right_task) {
            pthread_mutex_unlock(&task_mutex);

            int num_refused = 0;
            int i = (rank + 1) % size;


            while (i != rank) {
                int rank_to_ask = i % size;
                if (rank_to_ask == rank) {
                    continue;
                }
                MPI_Send(&rank_to_ask, 1, MPI_INT, rank_to_ask, 0, MPI_COMM_WORLD);
                int recv_tasks = -1;
                MPI_Recv(&recv_tasks, 1, MPI_INT, rank_to_ask, 1, MPI_COMM_WORLD, &status);
                if (recv_tasks > 0) {
                    MPI_Recv(tasks_buffer, recv_tasks, MPI_INT, rank_to_ask, 1, MPI_COMM_WORLD, &status);
                    int recv_after_first_send = 0;
                    MPI_Get_count(&status, MPI_INT, &recv_after_first_send);
                    if (recv_after_first_send < recv_tasks) {
                        MPI_Recv(&tasks_buffer[recv_after_first_send], recv_tasks - recv_after_first_send,
                                 MPI_INT, rank_to_ask, 1, MPI_COMM_WORLD, &status);
                    }

                    pthread_mutex_lock(&task_mutex);
                    if (right_task / TASKS_IN_PROCESS != (right_task + recv_tasks) / TASKS_IN_PROCESS) {
                        memcpy(&tasks[right_task % TASKS_IN_PROCESS], tasks_buffer,
                               sizeof(int) * (TASKS_IN_PROCESS - right_task % TASKS_IN_PROCESS));
                        memcpy(tasks, &tasks_buffer[TASKS_IN_PROCESS - right_task % TASKS_IN_PROCESS],
                               sizeof(int) * (recv_tasks - (TASKS_IN_PROCESS - right_task % TASKS_IN_PROCESS)));
                    }
                    else {
                        memcpy(&tasks[right_task % TASKS_IN_PROCESS], tasks_buffer, sizeof(int) * recv_tasks);
                    }
                    right_task += recv_tasks;
                    pthread_mutex_unlock(&task_mutex);
                    break;
                }
                else {
                    num_refused += 1;
                }
                i = (i + 1) % size;
            }


            pthread_mutex_lock(&task_mutex);
            if (num_refused == size - 1) {
                break;
            }
        }

    }
    free(tasks_buffer);
    pthread_mutex_unlock(&task_mutex);
}


int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Invalid amount of arguments" << std::endl;
        return 1;
    }
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *tasks_completed = (int *)calloc(size, sizeof(int));
    double *time_executed = (double *)calloc(size, sizeof(double));

    MAX_ITERATIONS = atoi(argv[1]);
    TOTAL_TASK_NUM = atoi(argv[2]);

    pthread_attr_t attributes;
    if (pthread_attr_init(&attributes) != 0) {
        perror("Cannot initialize attributes");
        exit(2);
    }

    pthread_mutex_init(&task_mutex, NULL);

    if (0 != pthread_attr_setdetachstate(&attributes, PTHREAD_CREATE_JOINABLE)) {
        perror("Error in setting attributes");
        exit(3);
    }

    TASKS_IN_PROCESS = TOTAL_TASK_NUM / size;
    if (TOTAL_TASK_NUM % size > rank) {
        TASKS_IN_PROCESS += 1;
    }
    right_task = TASKS_IN_PROCESS;

    int *tasks = (int *)calloc(TASKS_IN_PROCESS, sizeof(int));

    arguments_t argument = {tasks, rank, size};

    pthread_t answerer;
    if (pthread_create(&answerer, &attributes, answerTask, &argument) != 0) {
        perror("Error in creating thread answerer");
        exit(4);
    }

    pthread_attr_destroy(&attributes);

    for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
        for (int i = 0; i < TASKS_IN_PROCESS; i++) {
            tasks[i] = abs(rank - (iteration % size)) * TASK_MULTIPLIER * 2 * TASKS_IN_PROCESS * size;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        double start = MPI_Wtime();
        calculationTask(&argument);
        double end = MPI_Wtime() - start;

        double max_time = 0;
        double min_time = 0;
        MPI_Reduce(&end, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&end, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            std::cout << "Iteration = " << iteration << std::endl;
            double disbalance = max_time - min_time;
            std::cout << "Disbalance = " << disbalance << std::endl;
            std::cout << "Disbalance in % = " << disbalance / max_time * 100 << std::endl;
        }
        MPI_Gather(&right_task, 1, MPI_INT, tasks_completed, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&end, 1, MPI_DOUBLE, time_executed, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            std::cout << "Tasks completed" << std::endl;
            for (int i = 0; i < size; i++) {
                std::cout << tasks_completed[i] << " ";
            }
            std::cout << std::endl;
            std::cout << "Time executed" << std::endl;
            for (int i = 0; i < size; i++) {
                std::cout << time_executed[i] << " ";
            }
            std::cout << std::endl << std::endl;
        }
        left_task = 0;
        right_task = TASKS_IN_PROCESS;
    }
    MPI_Send(&STOP_CODE, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);

    if (pthread_join(answerer, NULL) != 0) {
        perror("Cannot join answerer");
        exit(5);
    }

    pthread_mutex_destroy(&task_mutex);
    free(tasks_completed);
    free(tasks);
    free(time_executed);
    MPI_Finalize();
    return 0;
}
