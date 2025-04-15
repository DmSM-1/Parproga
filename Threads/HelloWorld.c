#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#define NUM_THREADS 4


typedef struct _thread_data_t{
    int tid;
    int rank;
    double stuff;
}thread_data_t;


void* thr_func(void* arg){
    thread_data_t *thr_data = (thread_data_t*)arg;

    printf("Hello World| thr:%2d| total:%2d| rank:%2d|\n", thr_data->tid, NUM_THREADS, thr_data->rank);
    pthread_exit(NULL);
}


int main(int argc, char **argv){
    int commsize, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    pthread_t thr[NUM_THREADS];
    thread_data_t thr_data[NUM_THREADS];
    
    
    for(int i = 0; i < NUM_THREADS; i++){
        thr_data[i].rank = rank;
        thr_data[i].tid = i;
        pthread_create(&thr[i], NULL, thr_func, &thr_data[i]);
    }

    for(int i = 0; i < NUM_THREADS; i++){
        pthread_join(thr[i], NULL);
    }

    MPI_Finalize();
}